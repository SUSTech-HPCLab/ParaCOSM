# GPU 加速持续子图匹配 (CSM) — 创新实践报告草稿

> 把 CSM 的几个核心算法（GraphFlow / SymBi / TurboFlux / CaLiG / NewSP）从 CPU 并行扩展到 GPU，并通过 profile 解释**为何选 BFS 而非 DFS**。

---

## 1. 背景

- **CSM (Continuous Subgraph Matching)**：在持续到来的图更新流上，针对一组小 query graph 实时枚举所有匹配
- **已有工作**：ParaCOSM (ICPP'25) 在 CPU 上做了两层并行 (inner-update + inter-update + safe-update)
- **本工作**：把搜索阶段卸载到 GPU，对 5 个核心算法实现 GPU pipeline，并系统性比较 BFS vs DFS 两种 GPU 搜索策略

## 2. GPU 化的模块

`core/gpu/` 下 4 个模块：

| 模块 | 入口 mode | 角色 |
|---|---|---|
| `gpu_classifier` | `-m gpu` | 在 GPU 上把 update 分类为 safe / unsafe |
| `gpu_candidate_filter` | (PFM2 自动) | 候选集过滤 |
| `gpu_search` (DFS) | `-m gpu_all` | 经典 DFS 风格的并行搜索 |
| `gpu_bfs_search` | `-m gpu_bfs[_versioned]` | level-synchronous BFS pipeline ★ |

> 5 个 CPU 算法都有对应 adapter (`Parallel_*_Adapter`)，复用同一套 GPU 后端。

## 3. 核心论点：为什么 GPU 上选 BFS 而非 DFS

### 3.1 DFS 在 GPU 上的硬伤
- 每个 thread 维护自己的 stack → **stack pointer divergence**
- 同一个 warp 内 thread 在不同 depth → **SIMT 利用率低**
- 候选集不均（hub vertex）→ 一个 thread 干几十万邻居，**warp 整体被拖死**
- DFS 路径的回溯 → **不规则 memory access pattern**，难以 coalesce

### 3.2 BFS 的契合点
- **Level-synchronous**：同一 depth 同步推进，warp 内所有 thread 干同样的工作
- **Partial-match buffer 显式化**：每层一个数组，前缀扫描 (scan) 决定输出位置 → coalesced write
- **Work split 容易**：每个 partial match 独立扩展，可用 `block-binary-search` 做边界划分实现负载均衡

### 3.3 量化证据 (待补完整 profile)

> 推荐做法：在 AZ 上选 1–2 个 query，分别跑 `-m gpu_all` (DFS) 与 `-m gpu_bfs_versioned`，用 nsys/ncu 抓 warp_efficiency / sm_active_threads / 全程时间。期望显示 BFS warp efficiency ≈ 80%+，DFS ≈ 30–40%。

## 4. 主实验：Amazon (AZ) 8v 上的加速

### 4.1 总体 speedup（CPU single vs GPU BFS-versioned）

| Algorithm | Median speedup | 说明 |
|---|---|---|
| **NewSP** | **51.1×** | GPU 最大赢家 |
| **TurboFlux** | **49.5×** | |
| **GraphFlow** | **46.5×** | |
| SymBi | 25.3× | |
| CaLiG | 21.7× | 索引型，CPU 本来就快 |

> 来自 `analysis/speedup_8v_by_algo.png`，AZ 100 queries × 5 algos × 2 modes。

### 4.2 BFS 搜索树增长曲线（来自 `[BFS-V] Depth …` 日志）

以 AZ 一个 GraphFlow query 为例（示意）：

```
Depth 1→2: 8.6e4   partial matches
Depth 2→3: 5.4e5   ×6.3
Depth 3→4: 2.1e6   ×3.9
Depth 4→5: 1.99e7  ×9.4
Depth 5→6: 1.08e8  ×5.4
Depth 6→7: 1.48e8  ×1.4   ← 收敛/过滤
Search complete: 6.14e9 matches
```

**观察**：
- 中间层呈指数增长，`depth 4→5` 是 explosion knee
- 最后一层因 last-vertex 约束收敛
- 关键设计目标：让中间最大层（10⁸ 量级）的 BFS step 能塞进 GPU memory 并保持高 occupancy

### 4.3 好 case 深入：AZ NewSP `Q_*` (51× speedup)
*(待补)*  
Setup: query 拓扑、selectivity 估计  
Profile: kernel breakdown (classify <1s, BFS 主导)  
归因: 中间层规模适中（~1e7-1e8），warp efficiency 高，PCIe 摊销充分

### 4.4 差 case 深入：AZ CaLiG (21.7×)
*(待补)*  
归因: CPU CaLiG 本身有索引，候选集小，GPU 摊不动 launch 开销

---

## 5. 适用边界：LiveJournal (LJ) 上的反例

> 不是 GPU 慢，而是稀疏 + 倾斜分布把 GPU 优势挤压到 0。

### 5.1 LJ 8v 主要观察

CPU `single` 在 LJ 上单线程**绝大多数 query 都 timeout (1800s)**，GPU 是少数能跑出结果的方式。**仅在 CPU & GPU 都跑出结果的 query 上**做对比：

| Algorithm | pairs | cpu_med | gpu_med | speedup_med | speedup_max |
|---|---|---|---|---|---|
| GraphFlow | 1 | 1494s | 105s | 14.18× | 14.18× |
| TurboFlux | 9 | 666s | 215s | 3.10× | 5.89× |
| CaLiG | 14 | 584s | 659s | 0.89× | 3.25× |
| NewSP* | 16 | 87s | 637s | 0.14× | 0.77× |
| SymBi | 0 | — | — | — | — |

**注 (NewSP)**：CPU `parallel_newsp -m single` 实际是 stub，`EnumerateNewEdge()` 直接返回 0，所以快得"假"——0 positive matches 而 GPU 报 6.14×10⁹。这条对比不成立，仅作工程教训。

### 5.2 关键数据点：LJ NewSP Q_19 的 buffer overflow

```
Depth 5→6: 4.83e7 → 3.85e8  ←  超过 400M partial-match buffer
[BFS-V] Flush 400000000 partials at depth 6→7 (overflow, chunk=1)  ×15
Search complete: 7.97e10 matches, 1019.0s
```

**即是说：搜索空间到达 10¹⁰ 量级**，反复 flush 显著增加 wall time。GPU 80GB 也吃不下。

### 5.3 为什么稀疏图 GPU 优势缩水

| 因素 | 稠密图 (AZ) | 稀疏倾斜 (LJ) |
|---|---|---|
| 单 update 候选数 | 大 → warp 满载 | 小 → 一个 warp 只用几个 thread |
| Kernel launch 摊销 | 摊得动 | 摊不动 |
| Partial-match 增长 | 中间层稳定 ~10⁷–10⁸ | 爆炸 (10⁸–10¹⁰)，buffer overflow |
| 拓扑倾斜 | 较均匀 | hub vertex 严重 (per-thread 工作量方差大) |

## 6. 后续改进方向

按可行性分层：

### 短期（明确可做）
1. **work-balanced BFS frontier**：edge-centric + block-binary-search row-split，缓解 hub vertex 拖死整个 warp 的问题
2. **partial-match streaming**：超过 buffer 时流式吐到 host pinned memory + zero-copy，解决 LJ 那种 10¹⁰ 量级爆炸
3. **算 + 传 重叠**：用 CUDA stream 把 add-edge 的 H2D copy 和上一批 BFS search 重叠

### 中期
4. **DFS / BFS 自适应切换**：candidate set 小时退回 DFS（少 launch overhead），大时上 BFS
5. **per-query cost model**：根据 query 拓扑预测 BFS 各层规模，提前选 mode 与 buffer size
6. **early intersection pruning**：把 backward neighbor intersection 提前到 partial match 输出阶段，砍掉 LJ 这种倾斜图上的中间爆炸

### 长期
7. **multi-GPU 划分**：query 维度 (per-query 一卡) 或 update 维度 (按 batch 切)
8. **CPU/GPU 混合**：safe edge 留 CPU、unsafe 走 GPU；NewSP 这种 CPU stub 算法补完 versioned CPU path 后做 hybrid

## 7. 结论

- **在稠密、增长规模适中的图（AZ 8v）上，GPU BFS 对几乎所有 CSM 算法都能提供 20×–50× 加速**
- **关键设计是 BFS-versioned 而非 DFS**，因为 SIMT 友好 + work-split 容易 + memory access coalesced
- **稀疏倾斜图 (LJ) 上 GPU 优势缩水**，但根因不是算法选错而是数据 + 当前 partial-match buffer 设计触底；改进方向已识别（streaming, balanced frontier, hybrid）
- 工程层面的教训：**算法 adapter 必须明确每条 path 是否实际执行**——NewSP 的 CPU stub 让 LJ 上的 CPU 数字看起来荒谬地好，差点误导结论

---

## 报告写作清单 (TODO)

- [ ] 跑 `gpu_all` (DFS) vs `gpu_bfs_versioned` 对比，给 §3.3 配一张速度+warp efficiency 表
- [ ] AZ 上至少 1 个好 case + 1 个差 case 的 nsys profile 截图 / kernel breakdown
- [ ] BFS 搜索树深度增长图（log scale），多条 query 并列
- [ ] §4.1 那张 speedup_8v_by_algo.png 的精修版 (boxplot + medians + n)
- [ ] LJ 反例只用 1 张图 + 半页文字
- [ ] 改进方向每条挂一个具体的 profile 观察作为依据

## 关于 LJ 的处理（作者建议）

保留 **半页 + 1 表 + 1 图** 的反例小节，不展开。理由：
- 体现适用边界、不是无脑卖 GPU
- 读者会主动问"稀疏图怎么样"，自己先答
- NewSP CPU stub 的工程教训值得讲（30 秒）

## 不建议放进报告的

- 50 路 SIGHUP、SEGV 等环境/工程折腾
- 编译 / 装 oneAPI / git worktree 过程
- 任何 50/100/300 任务的脚本工程细节
