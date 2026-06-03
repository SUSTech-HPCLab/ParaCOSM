# GPU-Accelerated CSM 实验报告

> 基于 ParaCOSM 框架，在 A100 80GB PCIe 上的 CPU-GPU 协同 CSM 探索
> 数据集: Amazon (403K V, 2.2M E), LiveJournal (4.8M V, 43M E)
> 日期: 2026-04-23

---

## 1. 模块概览

| 模块 | 文件 | 命令行 | 策略 |
|---|---|---|---|
| GPU Classify | `gpu_classifier.cu/h` | `-m gpu` | GPU 批量 safe/unsafe 边分类 |
| GPU Candidate Filter | `gpu_candidate_filter.cu/h` | PFM2 内自动触发 | GPU 并行候选过滤（label + joinability + visited） |
| GPU DFS Search | `gpu_search.cu/h` | `-m gpu_all` | 每个 GPU thread 独立 DFS 递归 |
| GPU BFS Search | `gpu_bfs_search.cu/h` | `-m gpu_bfs` | BFS 逐层展开，三阶段 kernel pipeline |

---

## 2. Phase 1: GPU Classify (`-m gpu`)

### 设计
- 查询边 label 三元组存 GPU constant memory（≤15 条）
- 数据图顶点标签 (`vlabels_`) 存 GPU global memory
- 每个 CUDA thread 处理一条数据边，遍历查询边三元组做双向 label 匹配
- 滑动窗口模式：窗口内边批量发 GPU 分类，找到第一个 unsafe 后 CPU 执行 AddEdge + 匹配

### 优化
- **vlabel 延迟更新**：只在有新顶点添加时才 cudaMemcpy vlabels（Amazon 更新流中 0 次顶点添加 → vlabel update 从 8.8s → 0ms）

### 结果

| 查询 | CPU single | GPU classify | 说明 |
|---|---:|---:|---|
| 8v sparse Q11 | 3.3s | 7.1s | GPU 2x 慢 |
| 9v 5.13x | 48.5s | 51.8s | 接近 |
| 10v pf_001 | 118s | 139s | GPU 慢 |

### 分析
- **Classify 计算量太小**（≤30 次比较/边），H2D/D2H 传输开销主导
- nsys profile: kernel 总时间仅 410ms（11.7 万次调用），cudaMemcpy 吃了 15.8s
- **结论**: Classify 不适合 GPU offloading，但基础设施（CUDA 编译链、CPU↔GPU 数据管道）已搭建

---

## 3. Phase 2: GPU Candidate Filter (PFM2 内)

### 设计
- 在 `Parallel_FindMatches2` 的 Layer 1 串行循环前加 GPU 路径
- 当 `u_min_nbrs.size() >= GPU_FILTER_THRESHOLD` 时：
  - 构建 `CandidateFilterTask`（候选列表 + joinability 约束的 flattened neighbor lists）
  - GPU kernel: 每 thread 一个候选，做 label check + binary search joinability + visited check
  - 返回有效候选 indices

### 结果

| 指标 | Amazon | LiveJournal |
|---|---|---|
| avg degree | 11 | 16 |
| `u_min` 典型候选数 | 10-50 | 10-50 |
| GPU 调用占 PFM2 比例 | 0.1% (threshold=1024) | 0.1% |
| 性能影响 | 无（threshold 太高不触发） | — |

Threshold=16 时（LJ 上 81% 调用走 GPU）：
- LJ 6v tree: CPU 41s → GPU 48s（**慢 17%**）

### 分析
- `u_min` 的设计目标就是选最小 degree 邻居 → 候选数天然很小
- 每次 GPU 调用 ~20-50μs（1 kernel launch + 6 cudaMemcpy），CPU 串行处理 10-50 候选仅 1-5μs
- **结论**: per-PFM2-call GPU offloading 不可行，threshold 设为 1024（基本不触发）

---

## 4. Phase 3: GPU DFS Search (`-m gpu_all`)

### 设计
- **Inter-update 并行**：所有 unsafe 边先加入图，再一次性搜索
- 数据图以 **padded CSR** 格式常驻 GPU（每顶点预留 2x 容量，支持增量更新）
- 查询图 + 匹配顺序在 constant memory
- 每个 GPU thread = 一个 `(unsafe_edge × query_edge × direction)` 的完整 DFS
- Device-side 递归 `gpu_dfs()`，stack size 8KB/thread

### CSR 增量更新
- 每次 `AddEdge` 只更新 2 个顶点的邻居数据（cudaMemcpy ~100 bytes × 2）
- 避免 ~35ms 的全量 CSR 重建
- Padded CSR: `capacity[v] = max(degree * 2, degree + 8)`

### 结果

| 查询 | CPU batch_all 8T | GPU DFS | GPU/CPU |
|---|---:|---:|---:|
| 8v sparse Q11 | 827ms | 1,260ms | 1.5x 慢 |
| 8v sparse Q100 | 705ms | 3,848ms | 5.5x 慢 |
| 9v 5.13x | 36.5s | 149.5s | 4.1x 慢 |

### nsys Profile (8v Q11, 30s run)

| 指标 | 值 |
|---|---|
| gpu_search_kernel 总计 | 1.34s（366 次，avg 3.7ms） |
| 最大单次 kernel | 23ms |
| cudaMemcpy | 0.97s |

### 分析
- **SIMT 分支发散**: 2.3M threads，大部分开头就因 label 不匹配 return（空跑），有效 thread < 5%
- **递归开销**: 8v 查询 6 层递归 × 120B/level = 720B/thread, 2.3M threads × 720B = ~1.7GB local memory
- **负载极不均**: 有些 thread 找 0 个 match（1μs），有些找 100 万（100ms+），warp 内等最慢
- **结论**: DFS 递归 + 分支发散是 GPU SIMT 模型最不擅长的场景

---

## 5. Phase 4: GPU BFS Search (`-m gpu_bfs`) ⭐

### 设计
三阶段 kernel pipeline，消除递归和 warp 分歧：

```
Init kernel:   unsafe_edges × query_edges × 2 → depth=2 partial matches
                (每 thread 一个 task，做 label check，通过则写 partial match)
     ↓
Expand kernel: depth d → depth d+1
                (每 thread 一个 partial match，找 u_min，遍历候选，
                 label + joinability + visited check，通过则写新 partial match)
     ↓ (repeat for each BFS level: depth 2→3, 3→4, ..., Q-2→Q-1)
Count kernel:  depth Q-1
                (每 thread 一个 partial match，计数有效候选而不展开)
```

- Ping-pong 双 buffer：`d_buf_a_` ↔ `d_buf_b_`，交替读写
- 每个 partial match = `[order_idx, m[0], ..., m[Q-1]]` = `(Q+1) × 4 bytes`
- Buffer 容量: 200M entries（9v ≈ 8GB per buffer，A100 80GB 充足）
- 溢出处理：分块输入 + 自适应 chunk size

### 为什么 BFS 比 DFS 好

| 问题 | DFS | BFS |
|---|---|---|
| Warp 分歧 | 严重：同 warp 不同搜索路径 | **无**：同层所有 thread 做相同操作 |
| 递归 stack | 720B/thread → 低 occupancy | **零**：无递归 |
| 负载均衡 | 极不均：subtree 大小差 1000x | **自然均衡**：有效候选产生多个 output |

### 结果

| 查询 | CPU batch_all 8T | GPU BFS | 加速比 | 正确性 |
|---|---:|---:|---:|:---:|
| **8v sparse Q11** | 824ms | **397ms** | **2.07x** | ✓ |
| **8v sparse Q100** | 719ms | **386ms** | **1.86x** | ✓ |
| 8v dense Q1 | 110ms | 349ms | 0.31x | ✓ |

GPU BFS 时间分解（8v sparse Q11）：

| 阶段 | 时间 |
|---|---:|
| classify (CPU OMP) | 11ms |
| add edges (CPU) | 63ms |
| build CSR (GPU) | 35ms |
| BFS search (GPU) | 50ms |
| **总计** | **397ms** |

### 溢出问题
中间层 partial matches 可能爆炸（如 9v 5.13x 在 depth 5 产生 25M+），超出 buffer 导致结果不完整。已实现分块处理但大量溢出时仍会丢失 match。

### 适用条件
- ✅ 中间层 partial matches < 200M（buffer 容量）
- ✅ 查询图稀疏（低扇出，中间层增长可控）
- ❌ 查询图产生大量中间匹配（需要混合 BFS+DFS 策略）

---

## 6. 性能全景对比

### 6.1 GPU 模式对比 (parallel_graphflow 各执行模式)

Amazon 数据集（403K V, 2.2M E, A100 80GB），`--time-limit 120`

| Dataset | Query | single | batch_all 8T | gpu (classify) | gpu_all (DFS) | gpu_bfs | 备注 |
|---------|-------|-------:|-------------:|---------------:|--------------:|--------:|------|
| AZ 6v   | Q_1   | 313ms  | 69ms         | 2,128ms        | 355ms         | 347ms   | pos一致 |
| AZ 6v   | Q_5   | 348ms  | 54ms         | 2,151ms        | 365ms         | 374ms   | pos一致 |
| AZ 6v   | Q_10  | 327ms  | 72ms         | 82,314ms       | 417ms         | 353ms   | gpu classify 异常慢 |
| AZ 8v   | Q_1   | 4,147ms| 733ms        | 72,057ms       | 3,771ms       | **459ms** | gpu_bfs 最优 |
| AZ 8v   | Q_11  | 102,844ms | 20,534ms | 125,715ms    | T/O           | 21,559ms* | *buffer溢出: pos=1.28B vs 3.86B |
| AZ 8v   | Q_20  | T/O    | 92,132ms     | T/O            | T/O           | —        | 仅 batch_all 跑出 |

### 6.2 六算法横向对比 (Amazon 8V, `-m single`)

`--time-limit 120`, Amazon 数据集, 8 顶点查询图 (queryset_release/AZ/8v)

| Query | graphflow | p_graphflow 8T | symbi | p_symbi 8T | turboflux | p_turboflux 8T |
|-------|----------:|---------------:|------:|-----------:|----------:|---------------:|
| Q_1   | 4,119ms   | 15,240ms       | 2,316ms | 1,884ms  | 74,607ms  | 251,221ms      |
| Q_2   | 27,399ms  | 35,544ms       | 33,763ms | 23,191ms | T/O       | T/O            |
| Q_3   | 8,539ms   | 16,456ms       | 10,114ms | 4,966ms  | 55,497ms  | 11,559ms       |
| Q_4   | 63,850ms  | T/O            | 56,542ms | 41,616ms | 26,548ms  | 14,516ms       |
| Q_5   | 63,536ms  | T/O            | 55,975ms | 36,725ms | 23,496ms  | 14,757ms       |
| Q_6   | 64,297ms  | T/O            | 55,417ms | 32,359ms | 25,418ms  | 47,084ms       |
| Q_7   | 49,858ms  | T/O            | 7,654ms  | 4,098ms  | 53,009ms  | SEGF           |
| Q_8   | 37,686ms  | T/O            | 8,157ms  | 3,929ms  | 50,904ms  | SEGF           |
| Q_9   | T/O       | T/O            | T/O      | T/O      | 62,451ms  | 54,646ms       |
| Q_10  | T/O       | T/O            | T/O      | T/O      | T/O       | 12,246ms       |
| Q_11  | T/O       | T/O            | 97,351ms | 92,163ms | T/O       | T/O            |

> **注**: v2/v3 两轮跑出，Q_1~Q_3 取 v3（最新），Q_4~Q_11 取 v2。`p_graphflow` 在 v2 全 T/O（旧脚本无 `-m single`），v3 修复后 Q_1~Q_3 跑出但比 serial 慢。

### 6.3 关键发现

**GPU 模式**:
1. **GPU BFS 在 8v Q_1 上取得最优结果** — 459ms vs CPU batch_all 733ms（**1.6x 加速**）
2. **GPU classify (`-m gpu`) 全面慢于 CPU** — kernel 计算量太小，H2D/D2H 传输开销主导
3. **GPU DFS (`gpu_all`) 接近 CPU batch_all 但无优势** — warp 分歧 + 递归开销
4. **BFS buffer 溢出是主要瓶颈** — Q_11 仅恢复 33% 的匹配（1.28B vs 3.86B），需调大 buffer

**六算法对比**:
5. **p_graphflow 在 `-m single` 下反而比 serial 慢** — Q_1~Q_3 慢 2-4x，OMP 并行开销 > 收益
6. **symbi/p_symbi 是 graphflow 系列中最快的** — Q_1 最佳 1,884ms（p_symbi），Q_7/Q_8 也远超其他
7. **turboflux 在某些查询上独占优势** — Q_9/Q_10 只有 turboflux 系列能跑出来
8. **parallel_turboflux 不稳定** — Q_1/Q_3/Q_7/Q_8 SEGFAULT

### 6.4 TODO (下次测试)

- [ ] **调大 BFS buffer**: 当前 200M entries，Q_11 溢出严重。建议调至 500M+ (约 20GB/buffer，A100 80GB 可承受)
- [ ] 补跑 Q_12~Q_20 的六算法对比
- [ ] LiveJournal 数据集测试
- [ ] 修复 parallel_turboflux SEGFAULT

---

## 7. 修改的源文件

| 文件 | 改动 |
|---|---|
| `CMakeLists.txt` | 添加 `CUDA` 语言、4 个 `.cu` 文件、`CUDA::cudart` |
| `core/gpu/gpu_classifier.cu/h` | GPU 批量 Classify |
| `core/gpu/gpu_candidate_filter.cu/h` | GPU 候选过滤 |
| `core/gpu/gpu_search.cu/h` | GPU DFS + padded CSR + 增量更新 |
| `core/gpu/gpu_bfs_search.cu/h` | GPU BFS 三阶段 pipeline |
| `core/inter_executor/inter_executor.h/cpp` | `BatchUpdates_GPU` / `GPU_AllAtOnce` / `GPU_BFS` |
| `matching_executor/matching.h` | `GetQueryGraph()` / `GetDataGraph()` / `GetMatchingOrders()` |
| `matching_executor/Parallel_GraphFlow/parallel_graphflow.h/cpp` | GPU 成员 + PFM2 内 GPU 路径 + `GetMatchingOrders()` |
| `matching_executor/main.cpp` | `gpu` / `gpu_all` / `gpu_bfs` 模式分支 |

---

## 8. 未来方向

1. **调大 BFS buffer 重测**: 当前 200M entries 导致 Q_11 buffer 溢出（仅恢复 33% 匹配）。调至 500M~1B entries（~20-40GB/buffer），充分利用 A100 80GB 显存
2. **混合 BFS + DFS**: BFS 展开到中间层（如 depth 5），如果 partial matches > threshold → 回退到 GPU DFS 计数（每个 partial match 一个 thread）
3. **LiveJournal 测试**: 更大图可能让 GPU 内存带宽优势显现（CSR ~340MB vs A100 2TB/s HBM2e）
4. **Warp-level 协作 BFS**: 一个 warp 处理一个 partial match 的所有候选（而不是 1 thread 1 partial match），减少 atomicAdd 竞争
5. **Multi-stream pipeline**: 将 unsafe edges 分批，用多个 CUDA stream 重叠 CSR 更新 + 搜索
6. **修复 parallel_turboflux**: 多个 8v 查询 SEGFAULT，需排查内存越界
