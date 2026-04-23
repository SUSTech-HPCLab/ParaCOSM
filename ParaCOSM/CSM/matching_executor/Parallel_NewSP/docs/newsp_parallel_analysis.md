# NewSP Parallel Batch_All Strategy — Analysis & Benchmark

## 1. 实现概述 (C+D Strategy)

将原先串行的 `batch_all` Phase 1 拆分为 4 个并行阶段：

```
Phase 0: 扫描更新流 — vertex/remove 立即处理，edge insert 收集延迟
Phase 1: 并行 atomicIndexUpdate (#pragma omp atomic)
Phase 2: 并行 safe_detect (只读 index，筛选 unsafe 边)
Phase 3: 并行 BatchAddEdges (per-vertex sorted merge 代替逐条 vector::insert)
Phase 4: 并行 enumerate unsafe edges (已有)
```

### 修改文件

| 文件 | 修改内容 |
|------|---------|
| `graph/graph.h` | 新增 `atomicIndexUpdate`, `BatchAddEdges` 声明 |
| `graph/graph.cpp` | 实现 `atomicIndexUpdate` (omp atomic) 和 `BatchAddEdges` (per-vertex 并行排序合并) |
| `matching/main.cpp` | `batch_all` 模式重写为 4 阶段流水线 |
| `matching/CSMPP.cpp` | `AddVertex` 增加 index 数组扩容 |

### 关键优化点

1. **indexUpdate 并行化**：`index[v][label]++` 是可交换操作，用 `#pragma omp atomic` 实现无锁并行
2. **safe_detect 并行化**：只读 index 数组，天然并行
3. **BatchAddEdges**：将 N 条 `vector::insert` (O(N×deg)) 替换为一次性 per-vertex sorted merge (O(N log N + deg))
4. **safe_detect 提前于 AddEdge**：indexUpdate 不依赖 neighbors，可先更新 index → 并行 safe_detect → 最后批量修改 neighbors

## 2. 正确性验证 (8_self, data graph 6)

| Query | Initial | Positive | Negative | Status |
|-------|---------|----------|----------|--------|
| 8_self/sparse/Q_3 | 5,021,271 ✓ | 20,505,316 ✓ | 0 ✓ | PASS |
| 8_self/sparse/Q_4 | 4,263,827 ✓ | 12,239,519 ✓ | 0 ✓ | PASS |
| 8_self/sparse/Q_5 | 7,067,421 ✓ | 30,566,143 ✓ | 0 ✓ | PASS |

Serial batch_all (1T) vs Parallel batch_all (8T) 匹配数完全一致。

## 3. Benchmark 结果

### 3.1 8-vertex queries (8_self)

| Type | Query | Serial (ms) | Parallel 8T (ms) | Speedup |
|------|-------|-------------|-------------------|---------|
| sparse | Q_3 | 3,507 | 1,459 | **2.40x** |
| sparse | Q_4 | 3,358 | 929 | **3.61x** |
| sparse | Q_5 | 2,792 | 928 | **3.01x** |

### 3.2 9-vertex queries (9_self)

| Type | Query | Serial (ms) | Parallel 8T (ms) | Speedup | 备注 |
|------|-------|-------------|-------------------|---------|------|
| sparse | Q_1 | 1,945 | 760 | **2.56x** | |
| sparse | Q_10 | 8,296 | 3,143 | **2.64x** | |
| sparse | Q_12 | 60,000 | 53,942 | 1.11x | 双方 enumerate 超时 |
| dense | Q_100 | 952 | 372 | **2.56x** | |
| dense | Q_13 | 1,191 | 412 | **2.89x** | |
| dense | Q_16 | 1,789 | 431 | **4.15x** | |
| tree | Q_11 | timeout | timeout | - | **全卡在 InitialMatching** |
| tree | Q_17 | timeout | timeout | - | **全卡在 InitialMatching** |
| tree | Q_2 | timeout | timeout | - | **全卡在 InitialMatching** |

### 3.3 10-vertex queries (10_self)

| Type | Query | Serial (ms) | Parallel 8T (ms) | Speedup | 备注 |
|------|-------|-------------|-------------------|---------|------|
| sparse | Q_1 | 7,063 | 2,148 | **3.29x** | |
| sparse | Q_100 | 22,294 | 9,135 | **2.44x** | |
| sparse | Q_10 | 60,000 | 60,000 | 1.00x | 双方 enumerate 超时 |
| dense | Q_12 | 805 | 335 | **2.40x** | |
| dense | Q_13 | 9,477 | 3,309 | **2.87x** | |
| dense | Q_15 | 1,744 | 606 | **2.88x** | |
| tree | Q_23 | timeout | timeout | - | **全卡在 InitialMatching** |
| tree | Q_28 | timeout | timeout | - | **全卡在 InitialMatching** |
| tree | Q_36 | timeout | timeout | - | **全卡在 InitialMatching** |

### 3.4 时间分解示例

**10/sparse/Q_1 (serial):**
- Load: 2,259ms
- Preprocessing: 124ms
- InitialMatching: 1,048ms → **4,108,412 matches**
- Incremental (batch_all): 7,063ms → 92,506 unsafe / 244,341 total → **21,580,637 positive**

**10/dense/Q_13 (serial):**
- InitialMatching: 15,951ms → **10,245,232 matches**
- Incremental (batch_all): 9,477ms → 51,031 unsafe / 244,341 total → **124,233,250 positive**

**9/tree/Q_11 (serial):**
- InitialMatching: >60,000ms（超时）→ 10,767,567,948 matches（未完成）
- Incremental: **根本没有运行到**

## 4. 瓶颈分析

### 4.1 当前 Incremental 阶段加速已达 2.4x–4.2x

C+D 策略成功将原先占 90%+ 的串行更新阶段并行化：
- Phase 1 (atomicIndexUpdate): ~10ms（244K 边，几乎零开销）
- Phase 2 (safe_detect): ~50ms (并行只读)
- Phase 3 (BatchAddEdges): ~100ms (per-vertex merge)
- Phase 4 (enumerate): 占剩余 95%+ 时间

### 4.2 剩余瓶颈

#### 瓶颈 1: InitialMatching 完全串行 ★★★★★

```cpp
void CSMPP::InitialMatching() {
    for(int i = 0; i < this->data_.NumEdges(); ++i){
        this->searchInit(v1, v2, elabel, init);  // 串行逐边
    }
}
```

所有 9v/10v tree 查询全部卡死在此处。初始数据图 ~2.2M 边，每边触发一次搜索。
tree 查询产生数十亿匹配，单线程根本跑不完。

**解法：** 可复用 batch_all 策略 — 对已有边并行 searchInit，每个边独立 worker：
```cpp
#pragma omp parallel for schedule(dynamic, 1)
for(int i = 0; i < this->data_.NumEdges(); ++i){
    CSMPP worker = CreateBranchWorker();
    worker.searchInit(v1, v2, elabel, init);
    #pragma omp critical { MergeWorkerCounters(worker); }
}
```

#### 瓶颈 2: Worker 创建开销 ★★★

每次 OMP parallel for 迭代都会 `new CSMPP` + 深拷贝：
- `queryVec` 整个 vector（含 `isolatedVertexTimes` mutable map）
- `updateEdgeFindQuery` / `initEdgeFindQuery` (std::map)
- `match`, `matchCandidate`, `visited_` 数组

**解法：**
- `updateEdgeFindQuery` / `initEdgeFindQuery` 改为 `const` 引用/指针（只读数据）
- 预创建 thread-local worker pool（每线程一个 worker，复用）
- `isolatedVertexTimes` 改为 worker-local 临时分配

#### 瓶颈 3: Phase 4 负载不均 ★★

Dynamic scheduling 已有，但少数高 fan-out unsafe 边可能贡献大部分搜索时间。

**解法：**
- 按预估搜索代价对 unsafe 边排序（例如按两端度数乘积降序），大任务先处理
- `kParallelCandidateMinSize` 从 256 降到 64-128，让更多中等候选也能在 `SearchFreeVertexCandidates` 内并行

#### 瓶颈 4: Enumerate 超时 (大query) ★★

9/sparse/Q_12 和 10/sparse/Q_10 在 60s 内无法完成 enumerate。
这是算法本身的搜索空间过大问题，不是并行度问题。

**解法方向：**
- 更激进的剪枝策略（DFS early termination）
- 搜索空间估计 + adaptive timeout per edge

## 5. 优先级排序

| 优先级 | 优化方向 | 预期收益 | 实现难度 |
|--------|---------|---------|---------|
| P0 | 并行化 InitialMatching | tree 查询从不可跑 → ~4-6x | 低（复用现有 worker 模式）|
| P1 | Worker pool 复用 | 降低 overhead 10-30% | 中 |
| P2 | 负载均衡优化 | 提升 ~10-20%（长尾查询） | 低 |
| P3 | 降低 kParallelCandidateMinSize | 提升中等查询 ~10% | 极低（改常量） |
