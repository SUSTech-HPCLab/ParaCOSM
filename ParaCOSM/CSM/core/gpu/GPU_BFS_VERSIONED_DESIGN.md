# GPU Versioned BFS Search 设计文档

> 基于 CPU versioned search 的成功验证，将时间戳机制扩展到 GPU BFS pipeline

## 1. 动机

### 当前 GPU 模式的问题

| 模式 | 语义 | 性能 | 问题 |
|---|---|---|---|
| `gpu_bfs` (batch) | inter-update | 快 | 匹配数 ≠ inner-update（多算了） |
| `gpu_bfs_single` (流式) | inner-update | 极慢 | 每条边 1 次 kernel launch + CSR 增量更新 |

**目标**：`gpu_bfs_versioned` = inner-update 正确语义 + batch 级性能

### CPU versioned 验证结果

| 查询 | single | batch_all | versioned | versioned 加速比 |
|---|---:|---:|---:|---:|
| 8v Q_1 | 4,246ms | 780ms | 1,158ms | **3.7x** vs single |
| 8v Q_80 | 5,834ms | 968ms | 1,309ms | **4.5x** vs single |
| 9v Q_10 | 23,531ms | 4,145ms | 4,114ms | **5.7x** vs single |

正确性：4/4 查询匹配数与 `single` 完全一致 ✅

---

## 2. 设计

### 核心思想

1. 一次性把**所有边**（safe + unsafe）加入 GPU CSR，每条边带 timestamp
2. 每个 unsafe 边的搜索任务带自己的 `max_timestamp`
3. BFS kernel 中遍历邻居时，跳过 `timestamp > max_timestamp` 的边
4. 所有 unsafe 边的搜索**完全并行**（只读 CSR，无写冲突）

### GPU CSR 扩展

```
现有:
  d_csr_offsets_[V+1]      // 顶点偏移
  d_csr_neighbors_[2E]     // 邻居列表
  d_csr_elabels_[2E]       // 边标签

新增:
  d_csr_timestamps_[2E]    // 每条边的时间戳 (uint32_t)
```

内存开销：Amazon 2.2M 边 × 2 × 4B = ~17MB，LiveJournal 43M 边 × 2 × 4B = ~344MB（A100 80GB 充裕）

### BFS Kernel 修改

#### Init Kernel

```cuda
// 现有: 检查 data edge 的 label 是否匹配 query edge
// 新增: 每个 init task 带 max_timestamp

__global__ void bfs_init_versioned_kernel(
    ...,
    const uint32_t* unsafe_edge_timestamps,  // 每个 unsafe 边的 max_ts
    uint32_t num_unsafe_edges,
    ...
) {
    // task_id = unsafe_edge_idx * num_query_edges * 2 + ...
    uint32_t unsafe_idx = task_id / (num_query_edges * 2);
    uint32_t max_ts = unsafe_edge_timestamps[unsafe_idx];
    
    // partial match 中额外存储 max_ts
    // partial_match = [order_idx, max_ts, m[0], ..., m[Q-1]]
    //                  stride = Q + 2
}
```

#### Expand Kernel

```cuda
__global__ void bfs_expand_versioned_kernel(
    ...,
    const uint32_t* d_csr_timestamps,  // 边时间戳
    ...
) {
    // 从 partial match 取出 max_ts
    uint32_t max_ts = cur_buf[pm_idx * stride + 1];
    
    // 遍历 u_min 的邻居时:
    for (uint32_t i = 0; i < degree; i++) {
        uint32_t neighbor = d_csr_neighbors[offset + i];
        uint32_t edge_ts = d_csr_timestamps[offset + i];
        
        if (edge_ts > max_ts) continue;  // 版本过滤 ← 唯一新增的检查
        
        // ... 正常的 label check, joinability, visited ...
        
        // joinability 中的 binary search 也要检查 timestamp:
        // if (other_ts[dis] > max_ts) → not joinable
    }
}
```

#### Count Kernel

同 Expand，加 timestamp 过滤。

### Partial Match 格式变化

```
现有: [order_idx, m[0], ..., m[Q-1]]     stride = Q+1
新增: [order_idx, max_ts, m[0], ..., m[Q-1]]  stride = Q+2
```

每个 partial match 多 4 bytes。对于 9v 查询，stride 从 10 → 11 (10% 增长)。

### Inter-Executor 流程

```
BatchUpdates_GPU_BFS_Versioned:

1. CPU: 顶点更新 + OMP 并行分类 (safe/unsafe)
2. CPU: 为所有边分配时间戳 (ts=1,2,3,...)
3. CPU→GPU: BuildCSR_Versioned(data)  // 含 timestamps 数组
4. GPU: SearchBatchEdgesBFS_Versioned(
         unsafe_edges[],      // 只搜 unsafe 边
         timestamps[],        // 每个 unsafe 边的 max_ts
         num_query_edges, Q)
5. GPU→CPU: 返回 per-unsafe-edge match counts
6. CPU: 汇总结果
```

---

## 3. 实现计划

### 文件修改清单

| 文件 | 修改 |
|---|---|
| `gpu_bfs_search.h` | 新增 `BuildCSR_Versioned()`, `SearchBatchEdgesBFS_Versioned()` |
| `gpu_bfs_search.cu` | 新增 `d_csr_timestamps_`, 3 个 versioned kernel, 2 个 host 函数 |
| `inter_executor.h/cpp` | 新增 `BatchUpdates_GPU_BFS_Versioned()` |
| `main.cpp` | 新增 `-m gpu_bfs_versioned` 入口 |

### 预期开销

| 操作 | gpu_bfs (现有) | gpu_bfs_versioned (新) |
|---|---|---|
| CSR 构建 | 35ms | ~40ms (+timestamps 传输) |
| Kernel 额外比较 | 0 | 1 次 uint32 比较/邻居 (~0%) |
| Partial match 大小 | Q+1 | Q+2 (+10%) |
| Buffer 容量 | 200M | ~180M (stride 增大) |

**预期性能：与 gpu_bfs 几乎相同，但匹配数与 single 一致。**

---

## 4. 正确性论证

对于更新流 `[e1, e2, e3, ...]`，timestamp 分配为 `ts=1, 2, 3, ...`。

搜索 unsafe 边 e_k (timestamp=k) 时：
- 只看到 `timestamp ≤ k` 的边 = 原始图 + e1..ek
- 这等价于按顺序处理到 e_k 时的图状态
- 多个搜索并行执行互不干扰（只读 CSR）

因此 `versioned` 的总匹配数 = Σ(每条 unsafe 边在其时间点的匹配数) = `single` 的结果。

CPU 版本已验证 4/4 正确 ✅

---

## 5. 与 CPU versioned 的差异

| 方面 | CPU versioned | GPU versioned |
|---|---|---|
| timestamp 存储 | `edge_timestamps_` 在 Graph 对象 | `d_csr_timestamps_` 在 GPU global memory |
| 过滤位置 | `FindMatches_versioned()` 循环内 | BFS expand/count kernel 内 |
| 并行粒度 | OMP per-unsafe-edge | CUDA per-partial-match (更细) |
| u_min 选择 | 需要遍历计数 visible 边 | 可用 warp reduction 或近似 |
| 内存开销 | ~2MB (Amazon) | ~17MB (Amazon, GPU HBM) |
