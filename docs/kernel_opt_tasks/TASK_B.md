# Task B: Shared Memory m[] 扩展 + __launch_bounds__ 全覆盖

## 文件路径
- 本文档: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/TASK_B.md`
- 总览: `/home/v-haibinlai/haibin/ParaCOSM/docs/kernel_opt_tasks/README.md`
- 目标文件: `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/core/gpu/gpu_bfs_search.cu`

## 问题描述

### 问题 1: m[] 从 global memory 反复读取
只有 `bfs_expand_count_lj_kernel` 将 partial match m[] 缓存到 shared memory：
```cuda
__shared__ uint32_t s_m_block[8][BFS_MAX_Q];
uint32_t* s_m = s_m_block[warp_in_block];
for (uint32_t t = lane; t < Q; t += 32) s_m[t] = pm[1 + t];
__syncwarp();
const uint32_t* m = s_m;
```

其他 kernel（`bfs_expand_kernel`, `bfs_count_kernel`, `bfs_expand_count_kernel`）仍从 global memory 读 m[]。每个候选 v/w 都要扫描 m[] 做 visited check（Q 次读取），joinability 循环也读 m[uo]。

### 问题 2: 缺少 __launch_bounds__
只有 LJ kernel 有 `__launch_bounds__(256, 8)`。其他 kernel 编译器自由分配寄存器，可能用 40+ regs 导致 occupancy 不足。

## 优化方案

### Part 1: 为 bfs_expand_kernel 添加 shared memory m[]

在 kernel 开头添加：
```cuda
__global__ void __launch_bounds__(256, 8) bfs_expand_kernel(  // 也加 launch_bounds
    ...
) {
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 1;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];

    // ★ 新增: cache m[] in shared memory
    __shared__ uint32_t s_m_block[8][BFS_MAX_Q];
    uint32_t warp_in_block = threadIdx.x >> 5;
    uint32_t* s_m = s_m_block[warp_in_block];
    for (uint32_t t = lane; t < Q; t += 32) s_m[t] = pm[1 + t];
    __syncwarp();
    const uint32_t* m = s_m;  // 替代原来的 pm + 1
    
    // 以下代码保持不变，m 已指向 shared memory
    ...
}
```

### Part 2: 为 bfs_count_kernel 添加 shared memory m[]

同样模式：
```cuda
__global__ void __launch_bounds__(256, 8) bfs_count_kernel(  // 加 launch_bounds
    ...
) {
    ...
    // ★ 新增
    __shared__ uint32_t s_m_block[8][BFS_MAX_Q];
    uint32_t warp_in_block = threadIdx.x >> 5;
    uint32_t* s_m = s_m_block[warp_in_block];
    for (uint32_t t = lane; t < Q; t += 32) s_m[t] = pm[1 + t];
    __syncwarp();
    const uint32_t* m = s_m;
    ...
}
```

### Part 3: 为 bfs_expand_count_kernel 添加 shared memory m[]

```cuda
__global__ void __launch_bounds__(256, 8) bfs_expand_count_kernel(  // 加 launch_bounds
    ...
) {
    ...
    // ★ 新增
    __shared__ uint32_t s_m_block[8][BFS_MAX_Q];
    uint32_t warp_in_block = threadIdx.x >> 5;
    uint32_t* s_m = s_m_block[warp_in_block];
    for (uint32_t t = lane; t < Q; t += 32) s_m[t] = pm[1 + t];
    __syncwarp();
    const uint32_t* m = s_m;
    ...
}
```

### Part 4: Versioned kernels

对 `bfs_expand_versioned_kernel` 和 `bfs_count_versioned_kernel` 同样添加：
```cuda
__global__ void __launch_bounds__(256, 8) bfs_expand_versioned_kernel(...)
__global__ void __launch_bounds__(256, 8) bfs_count_versioned_kernel(...)
```

注意 versioned kernel 的 m 指针偏移是 `pm + 2`（而非 `pm + 1`），缓存逻辑相应调整：
```cuda
for (uint32_t t = lane; t < Q; t += 32) s_m[t] = pm[2 + t];
```

## 改动清单

| Kernel | 添加 shared m[] | 添加 __launch_bounds__ |
|--------|----------------|----------------------|
| bfs_expand_kernel | ✅ | ✅ (256, 8) |
| bfs_count_kernel | ✅ | ✅ (256, 8) |
| bfs_expand_count_kernel | ✅ | ✅ (256, 8) |
| bfs_expand_versioned_kernel | ✅ | ✅ (256, 8) |
| bfs_count_versioned_kernel | ✅ | ✅ (256, 8) |
| bfs_expand_count_lj_kernel | 已有 | 已有 |
| bfs_expand_count_lj_v2_kernel | ✅ | ✅ (256, 8) |
| bfs_init_kernel | 不需要 | 不需要（thread-per-task，无 m[] 循环） |
| bfs_init_versioned_kernel | 不需要 | 不需要 |

## Shared Memory 预算验证

每个 block 用量: `8 warps × BFS_MAX_Q × 4B = 8 × 20 × 4 = 640 bytes`

A100 shared memory limit: 48KB/block（默认）或 164KB（可配置）。
640B 远低于限制。与 `__launch_bounds__(256, 8)` 兼容（8 blocks/SM × 640B = 5KB/SM）。

## 注意事项

1. **Early return 后不能访问 shared memory**：`if (warp_id >= in_count) return;` 必须在 shared memory 声明之前（已满足）
2. **`__syncwarp()` 而非 `__syncthreads()`**：因为 shared memory 只在 warp 内共享（每个 warp 有独立的 `s_m_block[warp_in_block]` 行）
3. **不影响其他 agent 的改动**：Agent A 修改 bsearch 函数体，Agent C 修改归约出口；都不冲突

## 预期收益
- m[] 每次 visited check 从 global (~100 cycles) 变为 shared (~1 cycle)
- 每个候选 v 节省 Q 次 global load（Q ≤ 20）
- `__launch_bounds__` 降低寄存器用量，提高 occupancy
- 预期 +5-15%（取决于 Q 和候选数量）

## 验证
- 所有 query 的 match count 必须与优化前 bit-identical
- 编译时检查 ptxas 输出的寄存器数（应 ≤ 32）
