# Task C: BlockReduce 优化 + Early Termination (实验性)

## 文件路径
- 本文档: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/TASK_C.md`
- 总览: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/README.md`
- 目标文件: `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/core/gpu/gpu_bfs_search.cu`
- 参考: `/home/v-haibinlai/haibin/GraphMiner/include/set_intersect.cuh` (early termination pattern)

## 任务 1: CUB BlockReduce 替代 per-warp atomicAdd

### 现状
Count kernel 的归约模式：
```cuda
// 每个 warp 做 warp_reduce，然后 lane 0 发一次 atomicAdd
local_count = warp_reduce_add_u64(local_count);
if (lane == 0 && local_count > 0) {
    atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
              static_cast<unsigned long long>(local_count));
}
```

256 threads/block = 8 warps → 每个 block 发 8 次 atomicAdd。

### 优化方案
用 CUB BlockReduce 先做 block 内归约，只发 1 次 atomicAdd/block：

```cuda
#include <cub/cub.cuh>  // 需要在文件头添加

// 在 kernel 内：
typedef cub::BlockReduce<uint64_t, 256> BlockReduce;
__shared__ typename BlockReduce::TempStorage temp_storage;

// 替换原来的 warp_reduce + atomicAdd：
uint64_t block_total = BlockReduce(temp_storage).Sum(local_count);
if (threadIdx.x == 0 && block_total > 0) {
    atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
              static_cast<unsigned long long>(block_total));
}
```

### 适用 kernel
- `bfs_count_kernel`
- `bfs_count_versioned_kernel`
- `bfs_expand_count_kernel`
- `bfs_expand_count_lj_kernel`
- `bfs_expand_count_lj_v2_kernel`

### 注意事项
1. **CUB 头文件**：需要 `#include <cub/cub.cuh>`。确认编译环境有 CUB（CUDA 11+ 自带）。
2. **Shared memory 开销**：`BlockReduce::TempStorage` 对 256 threads + uint64_t 约 2KB。需确认与 Agent B 的 shared memory 不冲突（Agent B 用 640B，CUB 用 ~2KB，总 ~2.6KB，仍远低于限制）。
3. **`__syncthreads()` 隐含**：CUB BlockReduce 内部会调用 `__syncthreads()`。需确保所有线程都到达此点（不能在之前有 divergent return）。
   - 当前 kernel 的 early return `if (warp_id >= in_count) return;` 会导致部分线程不参与 BlockReduce。**这是一个问题**。
   - **解决方案**：让 early-return 的线程不 return，而是将 `local_count = 0` 然后继续到 BlockReduce。或者保持现有的 warp_reduce + atomicAdd 模式（8 次 atomic 对性能影响极小）。

### 风险评估
- 收益较小：atomic 竞争不是当前瓶颈（已有 warp-level 聚合）
- 实现复杂度中等：需处理 early return 问题
- **建议**：如果 early return 改造复杂，可跳过此优化，专注 Task C-2

---

## 任务 2: Early Termination (排序邻居列表)

### 前提条件验证
首先需要确认 CSR 邻居列表是否已排序。当前 `bfs_bsearch` 用了 binary search，**隐含邻居列表已排序**。需验证 BuildCSR 是否保证排序（检查 `Graph::GetNeighbors()` 返回值是否排序）。

### 优化场景 1: expand_kernel 的 strided sweep 中 label filter early-exit

在 `bfs_expand_kernel` 的 strided sweep 中：
```cuda
for (uint32_t base = 0; base < nbr_count; base += 32) {
    uint32_t i = base + lane;
    ...
    if (i < nbr_count) {
        v = csr_neighbors[nbr_start + i];
        // 如果邻居按 ID 排序，且我们知道合法 vertex label 的 ID 范围，
        // 可以在整 warp 的 v 都不满足 label 时提前退出。
        // 但 vertex label 和 ID 无关，所以这里 early exit 不适用。
    }
}
```

**结论**：对 label-based filter，early termination 不适用（label 和 vertex ID 无相关性）。

### 优化场景 2: joinability bsearch 中的上界剪枝

在 joinability 检查中，如果邻居列表排序且当前候选 v 比列表最大值大，可直接判定不存在：
```cuda
// 在 bfs_bsearch 中添加 early-out:
__device__ __forceinline__ int bfs_bsearch(
    const uint32_t* arr, int size, uint32_t target
) {
    if (size == 0) return -1;
    // ★ 新增: 如果 target > 列表最大值，直接返回（利用排序性质）
    // 但这需要一次额外的 global load (arr[size-1])，不一定划算
    // 因为大部分 binary search 在前几次比较就会走到 hi < lo
    ...
}
```

**结论**：对单次 bsearch 添加 max-check 不一定划算（额外一次 load）。但如果配合 2-phase cache（Task A），Phase 1 的 `cache[31]` 就是列表的近似最大值，可以零开销判断。**此优化依赖 Task A 先完成**。

### 优化场景 3: bfs_expand_count_lj_kernel OUTER 循环 early-exit

在 LJ kernel 的 OUTER 循环中：
```cuda
for (uint32_t i = 0; i < nbr_count; i++) {
    uint32_t v = csr_neighbors[nbr_start + i];
    if (vlabels[v] != bfs_q_vlabels[u]) continue;
    if (csr_elabels[nbr_start + i] != u_min_label) continue;
    ...
}
```

如果 `csr_elabels` 也按某种顺序存储（不一定），可以做 early termination。但当前没有证据表明 elabels 排序。**不适用**。

### 优化场景 4: visited check 用 bitmask 替代线性扫描

当前 visited check：
```cuda
bool visited = false;
for (uint32_t d = 0; d < Q; d++) {
    if (m[d] == v) { visited = true; break; }
}
```

对 Q ≤ 20，这是 20 次比较。可以用 bitmask + warp vote 优化：
```cuda
// 如果 Q ≤ 32，可以用 warp 协作做 parallel compare
// 但这需要 m[] 已在 shared memory（Agent B 的工作）
// 且需要确保不增加寄存器压力

// 更简单的优化：预计算 m[] 的 min/max，快速排除大部分 v
uint32_t m_min = UINT32_MAX, m_max = 0;
for (uint32_t d = 0; d < Q; d++) {
    if (m[d] != UINT32_MAX) {
        m_min = min(m_min, m[d]);
        m_max = max(m_max, m[d]);
    }
}
// 然后:
if (v < m_min || v > m_max) { /* definitely not visited */ }
else { /* do the full scan */ }
```

**预期收益**：对 power-law 图，大部分候选 v 的 ID 范围远超 m[] 中的值，可跳过 ~80% 的 full scan。

### 推荐实现优先级

1. **场景 4 (visited bitmask/range check)** — 简单、独立、无依赖 ⭐⭐
2. ~~场景 1~~ — 不适用
3. **场景 2 (上界剪枝)** — 依赖 Task A ⭐
4. ~~场景 3~~ — 不适用

## 改动范围

### 场景 4 实现位置
在以下 kernel 的 visited check 前添加 range check：
- `bfs_expand_kernel` (line ~230 区域)
- `bfs_count_kernel` (line ~310 区域)
- `bfs_expand_count_kernel` (line ~400 区域)
- `bfs_expand_count_lj_kernel` (line ~600 区域)

模式：
```cuda
// 在 kernel 入口处（pivot 计算后、主循环前）预计算 m[] 的 min/max
uint32_t m_min = UINT32_MAX, m_max = 0;
for (uint32_t d = 0; d < Q; d++) {
    uint32_t mv = m[d];
    if (mv != UINT32_MAX) {
        if (mv < m_min) m_min = mv;
        if (mv > m_max) m_max = mv;
    }
}

// 然后在 visited check 处替换为:
bool visited = false;
if (v >= m_min && v <= m_max) {
    for (uint32_t d = 0; d < Q; d++) {
        if (m[d] == v) { visited = true; break; }
    }
}
```

## 预期收益
- BlockReduce: ~0-3%（atomic 不是瓶颈）
- Visited range check: ~3-8%（减少无效 m[] 扫描）
- 总计: 可能 +3-10%

## 验证
- 所有 query 的 match count 必须与优化前 bit-identical
- range check 绝对不能 false-negative（只能 false-positive → fall through to full scan）
