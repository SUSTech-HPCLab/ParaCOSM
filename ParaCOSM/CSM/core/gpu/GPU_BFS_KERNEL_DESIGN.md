# GPU BFS Kernel 设计文档：从 warp 协作到按图自适应算子

> 本文档记录 `core/gpu/gpu_bfs_search.cu` 在 2026/06 这轮的 kernel 级优化。
> 硬件 A100 80GB PCIe，CUDA 13.3。Golden = CPU `batch_all` 匹配数。
> 配套数据见 `analysis/gpu_ncu/`（ncu 报告 + `GPU_BFS_OPTIMIZATION_LOG.md`，本地）。
>
> 前置：`GPU_CSM_REPORT.md` 是更早的实验报告（DFS vs BFS 选型）。本文档接续它，
> 解决其中遗留的"溢出丢匹配""warp 分歧""稀疏图慢"等问题。

---

## 0. 一句话结论

GPU BFS 的最优策略**依赖数据图的形态**：Amazon（中密、unsafe 边多、中间层爆炸）
和 LiveJournal（稀疏幂律、unsafe 边少、最终匹配数爆炸）是**两个相反的 regime**，
需要**不同的 kernel**。本轮做了一套 warp 协作优化（对两者都好）+ 一个**按图特征
运行时自动选择**的稀疏图专用 count 算子（LJ 上 4.6×）。

---

## 1. BFS pipeline 回顾

对 Q 个顶点的查询，BFS 逐层展开 partial match：

```
Init   : 每个 (unsafe_edge × query_edge × dir) → depth-2 partial match
Expand : depth d → d+1，每层把 frontier 物化进 ping-pong buffer（d_buf_a/b）
  ...（depth 2→3, 3→4, ..., Q-2→Q-1）
Count  : depth Q-1，只计数不物化
```

每个 partial match 布局 `[order_idx, m[0..Q-1]]`，stride = Q+1 个 uint32。
host 端 `BFSFromDepth` 驱动逐层 launch，层间回 host 读 count、判溢出、swap buffer。

---

## 2. 优化链（按提交顺序）

| 代号 | 改动 | 收益 | 提交 |
|---|---|---|---|
| OPT-A | versioned `u_min` 选枢轴：O(Σdeg) 时间戳重扫 → O(1) 物理度数近似 | versioned −16% | `0ef6d0c` |
| OPT-B | warp 聚合原子（每 warp 1 次 atomicAdd 而非每 lane） | 正确，6v 中性 | `0ef6d0c` |
| **OPT-C** | **expand/count：1 线程/pm → 1 warp/pm** | **expand −51%** | `0ef6d0c` |
| gtid 修复 | 全局线程索引改 64 位 | 修 >2^27 pm 回绕 | `0ef6d0c` |
| **P1a** | emitter 计数器 `out_count` → uint64 | **修 10v/11v 的 2^32 垃圾值** | `b44c084` |
| **P1b** | 融合末层 expand+count，最大层不物化 | **修 9v 溢出丢匹配；9v Q_3 16s→1.3s** | `b44c084` |
| **LJ 算子** | 稀疏图内层并行 count，运行时选择 | **LJ 4.6×** | `3a9ce5d` |
| R2 | 融合末 K=2 层（device 递归） | **净亏损 2.9–5.4×，已弃** | — |

### 2.1 OPT-C：warp-per-partial-match（核心提升）

**问题**：旧 expand 是 1 线程处理 1 个 partial match。一个 warp 里挤了 32 个
**互不相关**的 pm，它们候选列表长度悬殊（度数倾斜）→ warp 等最慢 lane；且多数
候选在第一个 label check 就 `continue` → warp 塌缩到 ~2 个活跃 lane。

ncu 实测（`bfs_expand_kernel`，8v Q_5，depth 5→6）：

| 指标 | 旧(1线程/pm) | 新(1warp/pm) |
|---|---:|---:|
| **活跃线程/warp**（满32） | **2.40** | **15.56** |
| Duration | 10.46ms | 5.10ms (−51%) |
| SM 计算吞吐 | 26.7% | 46.8% |
| DRAM 吞吐 | 2.2% | 4.6% |

**改法**：1 warp 协作扫**同一个** pm 的候选（`for i = lane; i < nbr_count; i += 32`）。
同 warp 走相同代码路径、负载均摊、无 pm 间发散。

> **重要副结论**：DRAM 全程 2–5%、L2 命中 ~97% → 这个 kernel **不是访存瓶颈**。
> 所以 `bfs_bsearch` 的访存优化（uint2 打包、`__restrict__`、循环展开）**无收益**，
> 不做。瓶颈永远是 warp 利用率。

### 2.2 P1a：uint64 计数器（正确性）

稠密查询单层 expand 产出 >2^32 候选（10v Q_3、11v Q_1）。`out_count` 是 uint32 →
`atomicAdd` 回绕 → host 读到垃圾值（10v Q_3 报 4,294,974,963 ≈ 2^32，真值 1.32B）。
`d_count_`、4 个 writer kernel 的 `out_count`、`warp_agg_inc`、host 读写全改 64 位。

### 2.3 P1b：末层融合（防溢出 + 提速）

末层 expand 通常物化 depth Q-1（**最大的一层**，9v Q_1 达 2.77 亿）进 buffer，再
单独 count。`bfs_expand_count_kernel` 让每个 depth-(Q-2) 候选在寄存器里直接展开到
depth Q-1 并计数，**最大层永不物化**。两个收益：① 省一次最大层的全局读写；
② 消除 buffer 溢出的 flush（9v Q_1：1,839,971,202 丢匹配/非确定 → 精确 1,839,971,284）。

---

## 3. 核心洞察：LJ vs Amazon 两个 regime

实测（`gpu_bfs`，golden = CPU）：

| Case | unsafe 边 | 最大中间层 | 溢出? | 匹配数 | GPU 搜索 |
|---|---:|---:|:---:|---:|---:|
| AZ 8v Q_5 | 88K | 3.5亿 | **是(flush)** | 1.40B | 0.47s |
| AZ 10v Q_3 | — | >4亿 | **是** | 1.32B | 2.5s |
| LJ 6v Q_1 | 5,691 | 35万 | 否 | 64.4M | 0.22s |
| LJ 9v Q_1(100k) | 1,888 | 4348万 | **否** | **171亿** | 22.3s |

**两个 regime 完全相反：**

- **Amazon（中密、度均匀）**：unsafe 边多 → 中间层爆炸到 10⁸–10⁹ → **buffer 溢出**。
  瓶颈 = frontier 物化 + flush。OPT-C / P1a / P1b 都为此服务。
- **LiveJournal（稀疏、幂律社交图）**：unsafe 边极少 → 中间层很小（10⁵–10⁷）→ **从不溢出** —
  但**最终匹配数爆炸**（1888 边 → 171 亿匹配）。瓶颈转移到末层 count 的**纯枚举计算**。

---

## 4. LJ 专用 count 算子（按图自适应）

### 4.1 为什么需要单独算子

ncu（LJ 7v Q_1，融合 count kernel，LJ 的 ~100% 热路径）：

| 指标 | LJ count | (Amazon 同 kernel) |
|---|---:|---:|
| **活跃线程/warp** | **1.64 / 32** | 4.52 |
| DRAM 吞吐 | 0.02% | 0.22% |
| 计算吞吐 | 58.9% | 50% |

LJ 上利用率只有 **1.64/32（5%）**，比 Amazon 还惨。**根因**：融合 count kernel 把
并行分在**外层**候选（order[Q-2]）上，内层（order[Q-1]）每 lane 串行。Amazon 中密 →
外层候选多 → lane 有活；**LJ 稀疏 + label 过滤狠 → 外层候选极少 → 多数 lane 闲置，
只剩 1-2 lane 在内层深挖**。

→ **稀疏图上并行维度选错了。**

### 4.2 解法：内层并行

`bfs_expand_count_lj_kernel`：语义与默认 kernel 完全相同，但**并行维度翻转** ——
warp 统一遍历外层候选 v，把**内层**候选 w（LJ 上撞 hub 时很长）分给 32 lane。

```
默认(Amazon): for v in 外层[lane::32]:  for w in 内层: count   // 外层分lane
LJ:           for v in 外层(整warp一起): for w in 内层[lane::32]: count  // 内层分lane
```

### 4.3 运行时选择（算子选择逻辑）★

`BuildCSR` 里按图特征设 `sparse_mode_`：

```cpp
const char* e = getenv("GPU_BFS_SPARSE");      // A/B 测试用，强制 0/1
if (e) sparse_mode_ = (e[0] == '1');
else {
    double avg_deg = NumEdges() * 2.0 / V;
    sparse_mode_ = (V >= 1'000'000 && avg_deg < 32.0);   // LJ-like → 内层并行
}
```

- **`V ≥ 1M 且 平均度 < 32`** → 判为稀疏（LJ 类）→ 用内层并行 LJ 算子
- 否则（Amazon 类）→ 默认外层并行算子
- `GPU_BFS_SPARSE=0/1` 可强制覆盖（用于 benchmark 对照）

> 启发式依据：LJ 4.85M 顶点 / 平均度 15.9；Amazon 0.40M 顶点 / 平均度 ~11。
> 顶点数是更干净的区分信号（LJ 候选更分散、过滤更狠 → 外层短）。阈值保守，
> 误判只影响性能不影响正确性（两个 kernel 语义相同）。

### 4.4 实测收益

A/B（同 case，env 强制）：

| | active threads/warp | kernel 耗时 | 整体搜索 |
|---|---:|---:|---:|
| outer-par（旧） | 1.64 | 38.8ms | LJ 6v 223ms |
| **inner-par（LJ）** | **18.39（11.2×）** | **15.5ms（−60%）** | **LJ 6v 48.5ms（4.6×）** |

跨规模稳健（inner == outer 匹配数逐位一致）：

| LJ 查询 | 匹配数 | outer | inner | 提速 |
|---|---:|---:|---:|---:|
| 6v Q_1 | 64,404,068 | 223ms | 48.5ms | 4.6× |
| 8v Q_1 | 2,734,934,567 | 3437ms | 944ms | 3.6× |
| 9v Q_1 | 17,125,480,863 | 22316ms | 7571ms | 2.9× |

**Amazon 不受影响**：自动选 outer-par，8v Q_5 仍 1.40B @ 470ms。

---

## 5. 负面结论：R2（融合末 K=2 层）—— 不要按原样重试

R2 把 P1b 推广为 device 递归 `gpu_extend_count`，融合末 K 层。目标是消除 10v/11v
（最大层在倒数第 3-4 层）的溢出。**正确性没问题，但 K=2 全面慢 2.9–5.4×**。

| 查询 | P1b(K=1) | R2(K=2) |
|---|---:|---:|
| 7v Q_1 | 384ms | 2091ms |
| 8v Q_5 | 474ms | 2307ms |

**原理（关键教训）**：融合一层 = 用**32-lane 并行物化**换**warp 串行递归**。物化一层
成本 = 一次全局写+读（A100 ~2TB/s，便宜）；串行化一层 = 把该层并行度塌缩。当被消除
的层高扇出（9v Q_1 depth-7 有 2.77 亿），这个交换永远亏。**warp-local 深搜只在
buffer 真不够、且被消层低扇出时才赢** —— A100 80GB + 400M buffer 几乎不进这个区间。

---

## 6. 方法论：为什么 LJ 算子赢了而 R2 输了

两者都是"warp-level 重写"，区别只有一个：**改的是不是 profile 指出的那个瓶颈。**

- **R2**：凭直觉改（"融合更多层应该更好"）→ 把高扇出层串行化 → 亏损。
- **LJ 算子**：先 ncu（1.64/32，外层维度退化）→ 针对性改对的维度 → 11.2× 利用率。

**先 profile，再改 profile 指出的东西。** 这是本轮（含 R2 弯路）最重要的工程纪律。

---

## 7. 当前状态与后续

**已完成**（`report` 分支，commit `0ef6d0c` / `b44c084` / `3a9ce5d`）：
- OPT-A/B/C + gtid 修复 + P1a/b + LJ 自适应算子
- `--query-dir` 批量模式（`68c6179`，加载一次图跑全部查询，全量验证用）

**正确性**：Amazon 6v–11v、LJ 6v–9v 全部对齐 CPU 黄金值，确定。

**未来方向**：
- versioned 模式也套 P1b 融合（目前仅非 versioned）
- `sparse_mode` 启发式可进一步用"中间层增长率"在线判定（比静态阈值更准）
- 大于 11v 的查询若前段层（非末层）爆 buffer，需流式分块
- `__launch_bounds__` + block size 调优（OPT-C 后 occupancy 68→56%）
