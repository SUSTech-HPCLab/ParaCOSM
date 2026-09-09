# Task A: 2-Phase Binary Search 优化

## 文件路径
- 本文档: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/TASK_A.md`
- 总览: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/README.md`
- 目标文件: `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/core/gpu/gpu_bfs_search.cu`
- 参考实现: `/home/v-haibinlai/haibin/GraphMiner/include/search.cuh` (binary_search_2phase)

## 问题描述

当前 `bfs_bsearch` 是纯 global memory binary search：
```cuda
__device__ __forceinline__ int bfs_bsearch(
    const uint32_t* arr, int size, uint32_t target
) {
    int lo = 0, hi = size - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        uint32_t val = arr[mid];          // ← 每次比较都是 global memory load (~100 cycles)
        if (val == target) return mid;
        if (val < target) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}
```

对 degree=1000 的邻居列表，需要 ~10 次 global memory 比较。ncu 显示这是 #1 延迟源。

## 优化方案

### 核心思想
在搜索前，将目标列表均匀采样 32 个 pivot 加载到 shared memory。搜索时先在 shared memory 中定位区间（~5 次比较，~1 cycle each），再在 global memory 的小区间中精确搜索（~5 次比较）。Global memory 访问从 ~10 次降至 ~5 次。

### 适用场景分析

在 `bfs_expand_count_lj_kernel` 的 INNER 循环中：
```cuda
// INNER: 32 lanes split the w candidate list
for (uint32_t k = lane; k < nc2; k += 32) {
    ...
    for (uint32_t j = q2s; j < q2e; j++) {
        // ← 所有 32 lanes 搜索 **同一条邻居列表** csr_neighbors[s..s+nc)
        //   但每个 lane 有不同的 key (w)
        int pos = bfs_bsearch(csr_neighbors + s, nc, w);
        ...
    }
}
```

**关键**：在 LJ kernel 中，joinability 内层循环的每次迭代 j，所有 32 个 lane 搜索的是**同一条列表** `csr_neighbors[s..s+nc)`（因为 `m[uo]` 对所有 lane 相同——来自 shared memory 的 warp-uniform 值），只是搜索的 key (w) 不同。这**完美匹配** GraphMiner 的 2-phase 模式。

### 实现步骤

#### Step 1: 添加 shared memory cache 空间

在 `bfs_expand_count_lj_kernel` 已有的 shared memory 声明旁边添加：
```cuda
// 已有:
__shared__ uint32_t s_m_block[8][BFS_MAX_Q];

// 新增: 2-phase search 的 pivot cache，每个 warp 独占 32 个 uint32_t
__shared__ uint32_t s_search_cache[8][32];  // 8 warps × 32 pivots = 1024 bytes
```

注意：`__launch_bounds__(256, 8)` 下 8 warps/block × 32 × 4B = 1KB，加上已有的 s_m_block (8×20×4=640B)，总 shared memory ~1.6KB，远低于 48KB 限制。

#### Step 2: 添加 2-phase binary search 函数

```cuda
// 2-phase binary search: 先在 shared memory cache 中粗定位，再在 global memory 小区间精确搜索
// cache 必须在调用前由 warp 协作加载: cache[i] = list[i * size / 32]
__device__ __forceinline__ int bfs_bsearch_2phase(
    const uint32_t* list, const uint32_t* cache, int size, uint32_t target
) {
    // Phase 1: shared memory 中 O(log 32) = 5 次比较
    int bottom = 0, top = 32;
    while (top > bottom + 1) {
        int mid = (top + bottom) / 2;
        uint32_t y = cache[mid];
        if (target == y) return mid * size / 32;  // 精确命中 pivot (返回近似 index)
        if (target < y) top = mid;
        else bottom = mid;
    }
    // Phase 2: global memory 中 O(log(size/32)) 次比较
    int lo = bottom * size / 32;
    int hi = top * size / 32 - 1;
    if (hi >= size) hi = size - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        uint32_t val = list[mid];
        if (val == target) return mid;
        if (val < target) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}
```

#### Step 3: 在 joinability 循环前加载 cache

对 `bfs_expand_count_lj_kernel` 的 INNER joinability 循环，在**每次进入新的搜索列表**前协作加载：

```cuda
// 在 joinability 循环内，当要搜索 csr_neighbors[s..s+nc) 时：
// 加载 32 个等间距 pivot 到 warp 的 cache
uint32_t* my_cache = s_search_cache[warp_in_block];

// 在 inner loop 的 joinability check 中：
for (uint32_t j = q2s; j < q2e; j++) {
    uint32_t uo = bfs_q_neighbors[j];
    if (uo == u2_min) continue;
    uint32_t muo = (uo == u) ? v : m[uo];
    if (muo == UINT32_MAX) continue;
    uint32_t s = csr_offsets[muo];
    uint32_t nc = ...;
    
    // ★ 关键改动：协作加载 pivot cache（只需要 lane < 32 即可，即全 warp）
    // 但注意！在 INNER 循环中，不同 lane 有不同的 w，但搜索的 LIST 相同
    // 所以 cache 加载是 warp-uniform 的（all lanes 加载同一条列表的 pivots）
    if (nc >= 32) {  // 只对较长列表启用 2-phase（短列表直接 bsearch 更快）
        my_cache[lane] = csr_neighbors[s + lane * nc / 32];
        __syncwarp();
        int pos = bfs_bsearch_2phase(csr_neighbors + s, my_cache, nc, w);
        ...
    } else {
        int pos = bfs_bsearch(csr_neighbors + s, nc, w);
        ...
    }
}
```

#### Step 4: 对 OUTER joinability 也应用（warp-uniform case）

在 LJ kernel 的 OUTER 循环中，joinability check 也调用 bfs_bsearch：
```cuda
// OUTER joinability: 所有 lanes 搜索同一 list 同一 key (v is warp-uniform)
// 这里 2-phase 仍有效（减少 global 访问），但由于只有一个 key，
// 收益主要来自将 Phase 1 缩小搜索范围
```

这里可以用相同的 cache 机制。

### 注意事项

1. **`__syncwarp()` 放置**：cache 加载后必须有 `__syncwarp()` 才能被其他 lane 读取
2. **短列表 fallback**：`nc < 32` 时 pivot 采样无意义，直接用原始 `bfs_bsearch`
3. **Phase 1 精确命中**：`target == cache[mid]` 时，返回的 index 是 `mid * size / 32`，这不一定是精确位置。需要返回后仍然验证 elabel。但当前 bfs_bsearch 的调用方已经检查 `csr_elabels[s + pos]`，所以可以在 Phase 1 命中时做一次精确定位。
4. **寄存器压力**：2-phase 函数本身很短，`__forceinline__` 确保零调用开销。
5. **Shared memory 预算**：1KB cache + 640B m[] = 1.6KB/block，极低。

### 改动范围

仅修改 `gpu_bfs_search.cu` 中：
1. 添加 `bfs_bsearch_2phase` 设备函数（在 `bfs_bsearch` 下方）
2. `bfs_expand_count_lj_kernel` 中：添加 `s_search_cache` 声明 + 修改 INNER joinability 循环
3. 可选：`bfs_count_kernel`、`bfs_expand_kernel` 中类似应用

### 预期收益
- Global memory 访问减少 ~50%（从 ~10 次/search 降至 ~5 次）
- 直接缓解 ncu 报告的 long_scoreboard stall
- 预期 LJ 9v heavy queries +15-25%

### 验证
- 所有 query 的 match count 必须与优化前 bit-identical
- 用 `run.ipynb` 中的标准 benchmark 验证
