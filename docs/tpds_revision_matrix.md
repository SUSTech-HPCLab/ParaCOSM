# ParaCOSM TPDS 扩展对照与证据矩阵

更新日期：2026-09-09
基线论文：*ParaCOSM: A Parallel Framework for Continuous Subgraph Matching*, ICPP 2025
代码基线：`e3e0814`

## 1. 本文档的用途

本文档用于区分以下三类内容：

1. ICPP 论文已经发表的贡献；
2. ICPP 之后新增、可构成 TPDS 实质扩展的内容；
3. 虽然代码中已有入口，但尚未达到论文可声明状态的内容。

状态约定：

- **已有**：ICPP 已经公开，不能作为期刊新增贡献；
- **候选新增**：实现基本存在，但仍需补齐正确性、实验或论证；
- **可用证据**：已有代码和较系统的数据，可以进入稿件准备；
- **阻塞**：存在语义、正确性、复现性或实验覆盖问题，解决前不能作主张。

## 2. ICPP 与 TPDS 内容对照

| 主题 | ICPP 2025 已有内容 | ICPP 后新增内容 | 代码位置/提交 | 当前证据 | TPDS 处理建议 | 状态 |
|---|---|---|---|---|---|---|
| 通用框架 | 将多个单线程 CSM 算法接入统一并行框架 | matching 接口扩展为 batch/versioned/GPU 所需的查询图、匹配顺序和枚举接口 | `matching_executor/matching.h` | 五种 CPU 算法已有接入代码 | 保留为背景，不作为新增贡献 | 已有 |
| Inner-update 并行 | BFS 切分搜索树、并行枚举子树 | 持久线程团队、减少 fork/join；递归映射复用、高度数 hash join 等工程优化 | commits `32d2a26`、`efe4ed1`--`eabe142` | 有局部性能记录，缺统一 ablation | 作为增强实现和 CPU baseline 优化，不单独包装成核心新算法 | 候选新增 |
| Inter-update 并行 | safe/unsafe classifier、窗口/batch executor、三阶段过滤 | `batch_all` 将整个流的 unsafe-edge enumeration 并行化 | `BatchUpdates_AllAtOnce`；commit `b8550c9` | Amazon 上约 4--6x；重查询有明显热边失衡 | 仅作为动机、性能上界或 ablation；不可作为语义正确的主方法 | 阻塞 |
| Batch 正确语义 | 原论文按更新顺序执行 unsafe update | timestamp/versioned visibility：一次加入批内边，枚举更新 `e_k` 时只允许访问时间戳不大于 `k` 的边 | `BatchUpdates_Versioned`；commit `5c23bd7` | 早期仅少量查询与 `single` 对齐 | 作为 TPDS 的第一个核心新增贡献；补形式化定义、定理与大规模验证 | 候选新增 |
| CPU 热边处理 | 原论文有单次更新内部的 load balancing | versioned executor 使用 OpenMP task；TurboFlux/SymBi/CaLiG/NewSP 有不同程度的深层 task split | commit `622d5e2` 等 | TurboFlux 路径最完整；GraphFlow versioned 热边仍主要为串行递归 | 不应统一声称五种算法都解决了热边；逐算法做能力表和扩展性实验 | 阻塞 |
| GPU classifier | 无 | GPU 批量 safe/unsafe 分类 | `gpu_classifier.cu` | kernel 很轻，传输/launch 开销导致总体慢于 CPU | 作为负面实验或设计选择依据，不作贡献 | 可用证据 |
| GPU candidate filter | 无 | 将候选过滤 offload 到 GPU | `gpu_candidate_filter.cu` | 候选集通常只有 10--50，GPU 慢约 17% | 作为负面实验说明为何需要整棵搜索树 GPU 化 | 可用证据 |
| GPU DFS | 无 | 每个 CUDA thread 执行一个递归 DFS 任务 | `gpu_search.cu` | 分支发散、local stack、负载不均，整体慢于 CPU batch | 作为 BFS 设计动机和消融，不作为最终方案 | 可用证据 |
| GPU BFS | 无 | 无递归的逐层 frontier expansion；ping-pong/三缓冲与末层融合 | `gpu_bfs_search.cu`；commits `5590d77`、`364730c`、`b44c084` | Amazon/LJ 已有大量结果和 NCU 数据 | 作为 TPDS 第二个核心新增贡献 | 可用证据 |
| GPU warp 协作 | 无 | warp-per-partial-match、warp 聚合计数、64-bit emitter | commits `0ef6d0c`、`b44c084` | expand kernel 时间下降约 51%；计数溢出修复 | 给出 kernel 映射、复杂度及逐步 ablation | 可用证据 |
| 图特征自适应 | 无 | 根据图规模和平均度选择 outer-par 或 inner-par fused-count kernel | commit `3a9ce5d` | LJ 6v/8v/9v 分别约 4.6x/3.6x/2.9x kernel-path 改进 | 作为 GPU 方法的重要组成部分；需增加更多图验证启发式泛化性 | 候选新增 |
| GPU versioned semantics | 无 | timestamp-aware CSR、partial match 携带可见版本、version-aware joinability | `BatchUpdates_GPU_BFS_Versioned`、`SearchBatchEdgesBFS_Versioned` | 有实现；现有大规模 speedup 数据主要来自非 versioned `gpu_bfs` | 作为 TPDS 的语义桥梁；必须使用它重跑主实验 | 阻塞 |
| 批量实验工具 | 单查询脚本 | `--query-dir`、一次加载数据图、子进程 crash isolation、GPU watchdog、CSV 输出 | commits `68c6179`、`7ed2b8f`、`35ee1b1` | 可运行的采集基础设施 | 继续扩展为可复现实验 harness，不作为算法贡献 | 候选新增 |

## 3. 建议的 TPDS 核心主张

### C1：语义保持的版本化批处理

对插入更新流

\[
\Delta G = (e_1,e_2,\ldots,e_n),
\]

为每条边赋予单调时间戳 `ts(e_k)=k`。在枚举由 `e_k` 触发的新匹配时，只允许搜索访问满足 `ts(e) <= k` 的边。这样可以先物化整个批次，再并行执行多个更新的枚举，同时保持与顺序执行一致的可见图。

建议论文主张限定为：

> 对仅含边插入、且顶点集合固定的更新批次，versioned executor 产生的逐更新新增匹配总数与顺序 CSM executor 一致。

在删除边、插入/删除顶点以及同一边重复更新未完成设计和验证前，不应扩展上述主张范围。

需要补充：

- 可见图定义 `G_k = G_0 union {e_i | i <= k}`；
- soundness：versioned search 不会使用未来边；
- completeness：`G_k` 中包含 `e_k` 的每个新增匹配都会被枚举；
- uniqueness：每个新增匹配只归属于使其首次成立的最新更新；
- 分类器不会漏掉可能产生匹配的 unsafe update 的条件。

### C2：面向不同并行粒度的 CPU/GPU executor

CPU 适合在多个相对均匀的 unsafe updates 之间并行；GPU BFS 适合展开由少量热边产生的巨大搜索 frontier。论文不宜只写“增加 GPU 版本”，而应解释并验证两个后端分别利用不同的并行维度：

| 工作负载 | 主要并行维度 | 推荐 executor |
|---|---|---|
| unsafe updates 多且代价较均匀 | update-level | CPU versioned batch |
| unsafe updates 少但存在热边 | frontier/partial-match-level | GPU versioned BFS |
| 匹配量很小 | 避免 GPU 固定开销 | CPU sequential/parallel |

如果进一步实现运行时选择器，可将系统主张提升为异构自适应框架；否则只能表述为两个互补后端。

### C3：适合不规则 CSM 搜索的 GPU BFS pipeline

可以包含以下技术点：

- 用 level-synchronous BFS 消除 device recursion 和跨深度 warp divergence；
- 一 warp 协作处理一个 partial match；
- 融合最后一次 expand 和 count，避免最大 frontier 物化；
- 64-bit 计数和大规模 frontier 溢出处理；
- 面向稀疏幂律图的 inner-parallel count kernel；
- 根据数据图特征选择 count kernel。

其中 GPU classifier、细粒度 candidate offloading 和 GPU DFS 应作为负面结果/消融，用来说明最终设计为何必要。

## 4. 算法与模式能力矩阵

以下“实现”仅表示存在相应 override/入口，不等于已通过论文级正确性验证。

| 算法 | ICPP CPU inner/inter | `batch_all` | CPU `versioned` 接口 | hot-edge task split | GPU BFS 可取得 matching order | versioned 正确性证据 | TPDS 主表资格 |
|---|---:|---:|---:|---:|---:|---:|---|
| GraphFlow | 是 | 是 | 是 | 弱：versioned 枚举仍主要串行 | 是 | 少量 + 非 versioned 大样本 | 条件具备，需重跑 |
| TurboFlux | 是 | 是 | 是 | 强：有 taskloop/deep-layer split | 是 | 不充分 | 需正确性和稳定性验证 |
| SymBi | 是 | 是 | 是 | 有 chunk/task 路径 | 是 | 不充分 | 需正确性和稳定性验证 |
| CaLiG | 是 | 部分/新适配 | 是 | 有 versioned chunk 路径 | 是 | 不充分 | 需与原 CaLiG 输出逐查询对齐 |
| NewSP | 是 | 是 | 是，但为 adapter 实现 | 有 versioned chunk 路径 | 是 | 不充分；需确认保留 NewSP 语义与剪枝 | 暂不进入主 GPU claim |

建议第一轮只把 GraphFlow 作为 GPU 算法实例，把五种算法保留在 CPU framework evaluation。其他算法只有在 versioned correctness 和性能均完成后再纳入 GPU 横向结果。

## 5. 正确性阻塞项

### P0-1：`batch_all` overcount

已知全流可见会使同一匹配被多个 unsafe edge 重复计数。`batch_all` 只能作为 relaxed-semantics upper bound 或消融，不能作为 correctness ground truth。

### P0-2：删除更新未纳入 versioning

CPU/GPU versioned 主路径目前只对 `type == 'e' && is_add` 建立版本并搜索。边删除被跳过，顶点删除在批次末尾集中处理。因此论文必须限定 insertion-only，或者实现 interval/version visibility：

- 插入边可见区间 `[insert_ts, delete_ts)`；
- 删除 `e_k` 时在 `G_{k-1}` 上枚举 expired matches；
- 随后令该边在 `ts >= k` 不可见。

### P0-3：顶点更新时序

当前 versioned executor 在分类边之前先应用批次内全部 vertex additions，并在末尾处理 vertex removals。对顶点与边混合更新流，这不等价于顺序执行。第一版实验应使用固定顶点集合的纯边插入流。

### P0-4：索引算法的版本一致性

TurboFlux、SymBi、CaLiG 和 NewSP 依赖不同 ADS/DCS。需要确认：

- batch 后构建/更新的索引是否包含未来信息；
- 未来信息即便只用于 pruning，是否可能导致 false negative；
- GPU 通用 BFS 是否绕过算法特有索引，从而不能再称为对应算法的 GPU 版本。

### P0-5：frontier overflow

任何发生丢弃、非确定计数或只报告性能量级的查询都不能进入 correctness/performance 主表。必须记录每次运行的 overflow flag、最大 frontier、flush 次数和最终计数一致性。

## 6. 实验矩阵

### 6.1 E0：正确性门槛（必须最先完成）

| ID | 对比 | 数据/查询 | 更新类型 | 重复 | 通过条件 |
|---|---|---|---|---:|---|
| E0.1 | `single` vs CPU `versioned` | 4 数据集，6v--10v，每类至少 20 queries | edge insertion only | 1 | 每查询 positive/negative count 完全相同 |
| E0.2 | `single` vs `gpu_bfs_versioned` | 同 E0.1 | edge insertion only | 1 | 每查询计数完全相同且 overflow=0 |
| E0.3 | CPU vs GPU versioned | 重查询集合 | edge insertion only | 3 | 三次 GPU 结果确定且等于 CPU |
| E0.4 | 小型手工图逐更新 oracle | 覆盖同标签、多 unsafe edge、重复边、空匹配 | insertion only | 1 | 每个 update 的 delta count 都一致，而非只比总和 |
| E0.5 | ADS 算法交叉验证 | GraphFlow、TurboFlux、SymBi、CaLiG、NewSP | insertion only | 1 | 各自 sequential 与 versioned 输出一致 |

任何算法未通过 E0，就不能进入后续主性能图。

### 6.2 E1：总体性能

推荐统一比较：

1. 原始单线程算法；
2. ICPP ParaCOSM；
3. CPU versioned batch；
4. GPU versioned BFS；
5. 可比的最新 CPU/GPU CSM baseline。

报告：end-to-end runtime、纯 incremental search、throughput、平均/P95/P99 update latency、peak host/GPU memory、timeout success rate、正确计数查询比例。

不能把不同语义的 `batch_all` 和 `gpu_bfs` 混入主表；它们单独进入 relaxed/ablation 表。

### 6.3 E2：CPU 扩展性与 batch size

- Threads：`1, 2, 4, 8, 16, 32, 64`；
- Batch/window：`16, 32, 64, 128, 256, 512, 1024, all`；
- 每点至少 3 次；
- 同时报告 throughput 和 latency，避免只展示 batch 越大吞吐越高；
- 记录 unsafe update 数、每条 unsafe edge 的枚举时间和最大/中位代价比。

### 6.4 E3：CPU/GPU 适用区间

按以下特征分桶，而不只按 query vertex 数分桶：

- match count；
- unsafe-edge count/ratio；
- 最大 frontier；
- 最大单边枚举时间占比；
- query density；
- data graph 顶点数、平均度及 degree skew。

目标是找到 CPU/GPU crossover，并验证“均匀多更新”和“少数热边”两个 regime。

### 6.5 E4：GPU ablation

建议顺序：

1. GPU DFS；
2. level-synchronous BFS；
3. `+` warp-per-partial-match；
4. `+` fused last-level count；
5. `+` 64-bit/overflow-safe buffering；
6. `+` sparse inner-parallel count；
7. `+` automatic kernel selection；
8. non-versioned vs versioned 的额外开销。

应同时给出 end-to-end、kernel-only、CSR build 和 H2D/D2H 时间，避免只报告 kernel speedup。

### 6.6 E5：可扩展性和资源限制

- Query size 至少覆盖 6--11，超过 11 的结果只在无 overflow 时纳入；
- V100 32GB 与 A100 80GB 分别测量，说明显存容量和架构差异；
- 报告 frontier buffer 实际峰值，而不是只写固定 400M 上限；
- CPU 绑定 NUMA/线程亲和性，64-core 结果应说明使用单 socket 还是跨 socket；
- 对 timeout 使用 survival/success-rate 表达，不用外推时间替代实测主结果。

## 7. 建议论文结构

1. **Introduction**：ICPP 局限——窗口内仍有顺序点、热边导致 CPU 利用率低、GPU 细粒度 offload 无效。
2. **Background and ICPP ParaCOSM**：压缩介绍原有 inner/inter-update 设计，明确哪些内容来自会议版。
3. **Versioned Batch Execution**：语义、算法、正确性证明、ADS 兼容条件。
4. **CPU Versioned Executor**：跨更新 tasking、热边分裂、负载均衡。
5. **GPU Versioned BFS Executor**：CSR version、BFS pipeline、warp mapping、buffer 管理。
6. **Workload-aware Kernel/Backend Selection**：有自动后端选择器时作为正式章节；否则改为 analysis/discussion。
7. **Evaluation**：正确性门槛、总体性能、扩展性、crossover、ablation、memory。
8. **Related Work**。
9. **Limitations and Discussion**：插入流范围、显存/frontier 上限、单 GPU、后端选择限制。
10. **Conclusion**。

## 8. 执行优先级

### P0：决定论文是否成立

- [ ] 固定主语义为 insertion-only、fixed-vertex update stream；
- [ ] 建立逐更新 count oracle，而不只比较最终总数；
- [ ] 用 CPU/GPU versioned 跑 E0；
- [ ] 禁止 overflow/非确定结果进入主表；
- [ ] 审计五种算法的 version-aware ADS 行为。

### P1：形成完整实验故事

- [ ] 在目标 64-core CPU 上跑线程扩展和 batch-size 曲线；
- [ ] 在相同 query/update 输入上比较 CPU/GPU versioned；
- [ ] 重跑目前使用 `batch_all`/`gpu_bfs` 得到的主要加速数据；
- [ ] 为 GPU 优化建立可复现的逐步 ablation；
- [ ] 收集方差、尾延迟、内存和 overflow 指标。

### P2：复现与论文整理

- [x] CMake 去除硬编码 oneAPI 路径；
- [x] CUDA architecture 支持 `sm_70` 和 `sm_80` 或由用户配置；
- [ ] CPU-only build 不强制 CUDA；
- [ ] 统一算法名、模式名和 CSV schema；
- [ ] 保存 raw logs、运行命令、commit、硬件和软件版本；
- [ ] 在 TPDS 稿件中明确列出相对 ICPP 新增内容。

## 9. 当前推荐的最小可发表范围

如果不扩展删除语义，建议采用下面的最小闭环：

> 在固定顶点集合、边插入更新流上，ParaCOSM 的 versioned batch execution 保持逐更新 CSM 语义；CPU executor 在更新维度并行，单 GPU BFS executor 在搜索 frontier 维度并行。实验展示两类后端的适用区间、正确性、扩展性和优化贡献。

这个范围比声称“完整支持任意动态图更新和五种算法的 GPU 加速”更窄，但与当前代码最接近，也更容易形成可信、可复现的 TPDS 论文。
