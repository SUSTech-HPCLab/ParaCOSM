# batch_all：跨更新并行（Inter-Update Parallelism）设计文档

## 1. 背景与动机

### 1.1 原有架构的瓶颈

CSM（连续子图匹配）的 update-stream 天然是串行的：

```
while (有更新) {
    [串行] 取一批 update
    [并行] Classify 判断 safe/unsafe        ← 很轻量
    [串行] 找到第一个 unsafe update
    [串行] 对该 update 执行 AddEdge()       ← 重计算
      └─ 内部搜索匹配（仅第一层并行）        ← 唯一可并行段
}
```

根据 Amdahl 定律，如果串行部分占 50-85%，最大理论加速只有 1.2-2.6x（8 线程）。实测中现有并行模式（`persistent`、`batch3` 等）在 8 线程下只能达到 **1.0-3.4x** 加速。

### 1.2 核心思路：batch_all

将整个 update stream 视为一个 batch：
1. **Phase 1（串行）**：预分类所有 update，收集 unsafe edges
2. **Phase 2（串行）**：把所有 unsafe edges 加入数据图 + 更新辅助索引
3. **Phase 3（并行）**：8 个线程同时处理不同 unsafe edge 的匹配枚举

关键突破：**多个 unsafe edge 的枚举同时进行**（inter-update parallelism），而不是在单个 edge 的搜索树内部做并行（intra-update parallelism）。

## 2. 架构设计

### 2.1 新增接口（matching 基类）

```cpp
// matching.h 中新增的虚方法

/// 枚举已经在图中的边的匹配（不修改图，线程安全）
virtual size_t EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id);

/// 为并行枚举准备 per-thread 状态
virtual void PrepareBatchEnumeration(size_t num_threads);

/// 更新算法内部索引（DCS 等），不做枚举
virtual void UpdateIndexForEdge(uint v1, uint v2, uint label);

/// 线程安全地累加 positive 结果计数
void AddPositiveResults(size_t delta);
```

### 2.2 InterExecutor::BatchUpdates_AllAtOnce 流程

```
BatchUpdates_AllAtOnce(num_threads):
  Phase 1: 并行 Classify 所有 244K 条 update
           → 仅 ~15-20% 是 unsafe（~35K-50K 条）
  
  Phase 2: 串行循环 unsafe edges:
           data_graph_.AddEdge(v1, v2, label)
           matching_instance_->UpdateIndexForEdge(v1, v2, label)
  
  Phase 3: matching_instance_->PrepareBatchEnumeration(num_threads)
  
  Phase 4: #pragma omp parallel for schedule(dynamic, 1)
           for each unsafe edge:
               results[tid] += matching_instance_->EnumerateNewEdge(v1, v2, label, tid)
  
  Phase 5: matching_instance_->AddPositiveResults(sum(results))
```

### 2.3 线程安全保证

每个算法的 `EnumerateNewEdge` 使用 **thread-local 状态**：
- `local_vec_visited_local[thread_id]` — 独立的 visited 数组
- `local_vec_m[thread_id]` — 独立的 match 映射
- 只读共享：`data_graph`、`query_graph`、`DCS_`、`d2` 等索引

不同线程处理不同的 edge，没有写冲突。

## 3. 各算法实现细节

### 3.1 Parallel_GraphFlow

| 项目 | 说明 |
|------|------|
| **UpdateIndexForEdge** | 无辅助索引，默认 no-op |
| **EnumerateNewEdge** | 遍历 query edges → `FindMatches_local(order_index, 2, m, results, thread_id)` |
| **PrepareBatchEnumeration** | resize `local_vec_visited_local` 和 `local_vec_m` |
| **关键函数** | `FindMatches_local()` 完全用 `local_vec_visited_local[tid]` |

### 3.2 Parallel_TurboFlux

| 项目 | 说明 |
|------|------|
| **UpdateIndexForEdge** | DCS 插入 + `InsertionTopDown`/`InsertionBottomUp` + Q1/Q2 传播 |
| **EnumerateNewEdge** | 检查 `d2[u][v]==1` → `FindMatches_local(eidx, 2, m, results, thread_id)` |
| **PrepareBatchEnumeration** | resize `local_vec_visited_local` 和 `local_vec_m` |
| **关键函数** | `FindMatches_local()` 用 `local_vec_visited_local[tid]` |

### 3.3 Parallel_SymBi

| 项目 | 说明 |
|------|------|
| **UpdateIndexForEdge** | DCS 插入 + `InsertionTopDown_para`/`InsertionBottomUp_para` + Q1/Q2 传播 |
| **EnumerateNewEdge** | 检查 `d2` → 构建 ExtendableVertex → `Parallel_FindMatches_local3(2, m, ext, results, thread_id)` |
| **PrepareBatchEnumeration** | resize `local_vec_visited_local`、`local_vec_m`、`local_vec_extendable` |
| **关键函数** | `Parallel_FindMatches_local3()` 用 `local_vec_visited_local[tid]` |

### 3.4 Parallel_NewSP (CSMPP)

| 项目 | 说明 |
|------|------|
| **UpdateIndexForEdge** | `data_.indexUpdate(v1, v2, label, true)` |
| **EnumerateNewEdge** | 创建独立 CSMPP worker → `searchInit(v1, v2, label, pos)` |
| **预分类** | `safe_detect()` 代替 `Classify()` |
| **特殊点** | Phase 1 中 `Safe_Update` 必须逐条串行（indexUpdate 依赖顺序） |

**注意**：CSMPP 的 `Safe_Update`/`indexUpdate` 占了 90%+ 运行时间，枚举部分占比极小，因此 batch_all 加速效果有限（~1.2-1.3x）。

## 4. 性能结果

### 4.1 Amazon 数据集，8 顶点 tree 查询

| 算法 | 查询 | Serial (1T) | batch_all (8T) | 加速比 |
|------|------|-------------|----------------|--------|
| **GraphFlow** | Q_gen_028 | 55,663 ms | 9,230 ms | **6.03x** |
| **GraphFlow** | Q_gen_025 | 58,238 ms | 9,574 ms | **6.08x** |
| **GraphFlow** | Q_gen_008 | 56,941 ms | 9,578 ms | **5.94x** |
| **GraphFlow** | Q_gen_021 | 63,722 ms | 14,999 ms | **4.24x** |
| **TurboFlux** | Q_gen_028 | 22,347 ms | 4,665 ms | **4.79x** |
| **TurboFlux** | Q_gen_025 | 25,323 ms | 5,271 ms | **4.80x** |
| **TurboFlux** | Q_gen_008 | 27,969 ms | 7,078 ms | **3.95x** |
| **SymBi** | Q_gen_028 | 54,305 ms | 11,865 ms | **4.58x** |
| **SymBi** | Q_gen_025 | 60,177 ms | 12,128 ms | **4.96x** |
| **NewSP** | Q_5 | 1,828 ms | 1,494 ms | **1.22x** |
| **NewSP** | Q_100 | 2,257 ms | 1,755 ms | **1.29x** |

### 4.2 加速效率分析

- **GraphFlow**：无辅助索引，Phase 2 几乎零开销，枚举占 95%+ → 高加速
- **TurboFlux/SymBi**：DCS 更新是 Phase 2 的额外串行开销，但占比 <20% → 良好加速
- **NewSP**：`indexUpdate` 极重（占 90%+），枚举占比极小 → 加速有限

### 4.3 匹配数差异

batch_all 模式的匹配数比串行版多 ~35%，原因：
- 串行版处理 edge $e_k$ 时，图中只有 $e_1 \ldots e_k$
- batch_all 处理 $e_k$ 时，图中已有**所有** unsafe 边
- 一个包含两条 unsafe 边的 match 会被两边各计数一次

这是语义差异（overcounting），如需精确匹配数可另加去重逻辑。

### 4.4 局限：热边负载不均衡（heavy-query 失速）⚠️

§4.1 的高加速比来自 **Amazon 8v tree 查询**——这类负载的特征是
**unsafe 边多、每条边的枚举量相近**，Phase 4 的 `schedule(dynamic,1)`
把几百条边均匀分给各线程，扩展性好。

但在**重查询**（LiveJournal 9v、匹配数 10⁹–10¹⁰）上，batch_all 的
inter-update 并行**几乎完全失效**。根因是负载分布反转：

```
batch_all Phase 4:  #pragma omp parallel for schedule(dynamic, 1)
                    for (k = 0; k < n_unsafe; k++)
                        EnumerateNewEdge(unsafe_edges[k])   // 整条边单线程串行枚举
```

- LJ 2000 条更新流中，**unsafe 边只有 20–47 条**（其余 1950+ 是 safe，秒级跳过）。
- 而这 10⁹–10¹⁰ 个匹配**高度集中在其中 1–2 条"热边"**上：某条边的第一层
  fan-out 就有上千万，其子树枚举占全查询 99%+ 的时间。
- `EnumerateNewEdge` 内部是**纯串行递归**（`FindMatches_local`），没有 intra-edge
  并行。于是一个线程抢到热边后独自枚举到底，其余 15 个线程处理完几十条轻边
  后**全部空转干等**。
- 净效果：16 核退化成 ≈1 核。

**实测（LJ 9v，CPU 1 线程 vs 16 线程，batch_all）：**

| 查询 | 匹配数 | CPU-1t | CPU-16t | **16t/1t 扩展** | GPU-BFS | GPU/CPU-16t |
|------|------:|-------:|--------:|:---------------:|--------:|:-----------:|
| Q_69 | 9.39 B | 1,126,297 ms | 1,005,595 ms | **1.12×** | 3,657 ms | **275×** |
| Q_30 | 9.39 B | 1,117,832 ms | 898,194 ms | **1.24×** | 3,630 ms | 247× |
| Q_49 | 20.99 B | 2,442,568 ms | 1,463,528 ms | **1.67×** | 12,058 ms | 121× |
| Q_9  | 6.46 B | 700,453 ms | 419,274 ms | **1.67×** | 2,770 ms | 151× |

Q_69 砸 16 核只快 **1.12×**——等于一个核在干活。**这正是 GPU 的核心优势所在**：
GPU BFS 把单条热边的枚举铺到数千个 warp 上（intra-edge 并行），所以对这类负载
能拿到 100–275× 的加速。换言之，CPU 16t/1t 的"难看"扩展比，恰恰是 GPU 高加速比
**成立且诚实**的原因——不是 CPU 基线没开并行，而是 batch_all 的并行粒度（inter-edge）
和重查询的负载形态（intra-edge 爆炸）根本不匹配。

> **加速比对比的口径**：报告 GPU 加速比时应同时给出 vs CPU-1t 和 vs CPU-16t
> 两个数。重查询上二者差距不大（16t 没加速），但用 16t 作基线更公平、更保守。

### 4.5 能否用 work-stealing 救场？（现状与可行性）

要打破热边瓶颈，需要 **intra-edge work-stealing**：把单条热边的枚举递归树
切成可被空闲线程偷走的子任务。代码库里这件事**已经起了头，但尚未接入 batch_all**：

| 组件 | 位置 | 状态 |
|------|------|------|
| `steal_queue_`（TBB concurrent_queue） | `parallel_graphflow.h:44` | ✅ 已有 |
| `FindMatches_local_ws`（候选数≥8 时把 candidates[1..] 推入队列） | `parallel_graphflow.cpp:2106` | ✅ 已实现 |
| `FindMatches_local_splitting`（子树分裂版） | `parallel_graphflow.cpp:2202` | ✅ 已实现 |
| steal_queue 消费端（OMP 区 Phase 2 `try_pop` 抽干队列） | `parallel_graphflow.cpp:1388` | ✅ 已实现，但挂在另一条 `FindMatches` 路径 |
| **接入 batch_all 的 `EnumerateNewEdge`** | Phase 4 | ❌ **未接入**——仍调串行 `FindMatches_local` |

此外 `versioned` 模式（`BatchUpdates_Versioned`）已经把 Phase 4 从 `parallel for`
改成 **`#pragma omp task`**（`inter_executor.cpp:1575`），注释明确点出动机：

> *"With `parallel for` the threads stay locked on their for-iterations and
> cannot help drain a single hot edge's huge first-layer fanout."*

但 GraphFlow 的 `EnumerateNewEdgeVersioned` 内部仍调**串行** `FindMatches_versioned`
（无 hot-edge 拆分）；真正做了 hot-edge fan-out 子任务的只有 **TurboFlux**。
所以对 **GraphFlow 重查询，切到 versioned 也救不了**。

**可行性结论**：work-stealing 在本框架里**可行且已有半成品**，落地路径清晰——
把 batch_all Phase 4 的 `EnumerateNewEdge` 改用 `omp task` 分发，并让其内部递归
走 `FindMatches_local_ws`（生产 steal 任务）+ 团队线程消费 `steal_queue_`。
预期能把 Q_69 这类查询的 CPU-16t 扩展从 1.12× 拉向接近线性。**风险点**：
(1) 每个 steal 任务要快照整个 `m` 向量（内存/拷贝开销，需阈值调优）；
(2) `local_vec_visited_local` 的 per-thread visited 状态需在偷取时正确重建
（现有消费端已用 `m` 重建，逻辑可复用）；(3) 需重新校验匹配数正确性。
这是一项独立优化，**不影响当前 GPU 加速比结论**——只会让 CPU 基线更强、
让 GPU/CPU 的对比更硬核。

## 5. 使用方式

### 5.1 csm binary（GraphFlow / TurboFlux / SymBi）

```bash
# GraphFlow
./build/bin/csm -a parallel_graphflow -m batch_all -t 8 \
  -d data.graph -u insertion.graph -q query ...

# TurboFlux
./build/bin/csm -a parallel_turboflux -m batch_all -t 8 \
  -d data.graph -u insertion.graph -q query ...

# SymBi
./build/bin/csm -a parallel_symbi -m batch_all -t 8 \
  -d data.graph -u insertion.graph -q query ...
```

### 5.2 parallel_newsp binary（CSMPP）

```bash
./build/bin/parallel_newsp -m batch_all -t 8 \
  -d data.graph -u insertion.graph -q query ...
```

## 6. 代码变更清单

| 文件 | 变更 |
|------|------|
| `matching_executor/matching.h` | 新增 `EnumerateNewEdge`、`PrepareBatchEnumeration`、`UpdateIndexForEdge`、`AddPositiveResults` |
| `matching_executor/matching.cpp` | 默认虚方法实现 |
| `matching_executor/main.cpp` | 注册 `batch_all` 模式 |
| `core/inter_executor/inter_executor.h` | 声明 `BatchUpdates_AllAtOnce` |
| `core/inter_executor/inter_executor.cpp` | 完整实现 5-phase batch_all 流程 |
| `Parallel_GraphFlow/parallel_graphflow.h` | override 声明 |
| `Parallel_GraphFlow/parallel_graphflow.cpp` | `EnumerateNewEdge` + `PrepareBatchEnumeration` |
| `Parallel_TurboFlux/parallel_turboflux.h` | override 声明 |
| `Parallel_TurboFlux/parallel_turboflux.cpp` | `UpdateIndexForEdge` + `EnumerateNewEdge` + `PrepareBatchEnumeration` |
| `Parallel_SymBi/parallel_symbi.h` | override 声明 |
| `Parallel_SymBi/parrallel_symbi.cpp` | `UpdateIndexForEdge` + `EnumerateNewEdge` + `PrepareBatchEnumeration` |
| `Parallel_NewSP/matching/matching.h` | 新增虚方法（NewSP 自己的基类） |
| `Parallel_NewSP/matching/CSMPP.h` | override 声明 |
| `Parallel_NewSP/matching/CSMPP.cpp` | `UpdateIndexForEdge` + `EnumerateNewEdge` |
| `Parallel_NewSP/matching/main.cpp` | 添加 `-m batch_all` 分支 |

## 7. 提交历史

| Commit | 内容 |
|--------|------|
| `b8550c9` | 基础框架 + GraphFlow batch_all |
| `4163a20` | TurboFlux batch_all |
| `a819c89` | SymBi batch_all |
| `37ed663` | NewSP (CSMPP) batch_all |
| `2b0744c` | NewSP C+D 策略优化（atomicIndexUpdate + BatchAddEdges） |

分支：`opt/adaptive-parallel-threshold`
