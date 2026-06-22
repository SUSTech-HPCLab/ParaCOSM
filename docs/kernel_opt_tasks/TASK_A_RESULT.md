# Task A 结果：2-Phase Binary Search — 负优化，已回退

## 结论（2026-06）

**在 LiveJournal 上 2-phase binary search 是负优化，重查询慢 25–50%，已回退到 B 版（shared m[]）。计数全部 bit-identical（正确性无问题）。**

## 实测数据（LJ 9v，2000 边流，A/B 同机）

| query | B (shared m[]) | A (2-phase) | 结果 |
|---|---:|---:|---:|
| 9v Q_49 | 9842 ms | 14737 ms | **慢 1.50×** |
| 9v Q_69 | 3193 ms | 3979 ms | **慢 1.25×** |
| 9v Q_9  | 2252 ms | 2951 ms | **慢 1.31×** |
| 9v Q_85 |  866 ms |  860 ms | 持平（轻查询） |

计数：全部逐位一致（Q_49=20990632840, Q_85=349669 golden 等）。

## 为什么在 LJ 上变慢

TASK_A.md 的设计前提是「degree≈1000 的长邻居列表，bsearch 要 ~10 次 global load」。
这个前提来自 GraphMiner 的稠密图场景，**不匹配 LiveJournal 的幂律稀疏特性**：

1. **LJ 度数太小，2-phase 的前提不成立。**
   实测度数分布：99.1% 顶点 ≤32 邻居，99.8% ≤64，中位数仅 2。
   对 degree≤32 的列表，原始 bsearch 只需 ≤5 次比较，**整条列表已在 L1**（命中率 85.6%）。
   2-phase 想省的「长列表 global 访问」在 LJ 上几乎不存在。

2. **pivot 采样 + `__syncwarp()` 是纯增开销。**
   每次进入新搜索列表都要 `my_cache[lane] = csr_neighbors[s + lane*nc/32]` + `__syncwarp()` ——
   这是额外的 32 次 global load + 一次 warp 同步，对短列表远超它省下的几次比较。

3. **`nc_j >= 32` 门槛让 2-phase 几乎不触发，却付了重构代价。**
   LJ 99% 的列表 < 32，几乎全走 else 分支（原始 bsearch），2-phase 路径基本没用上。
   但为了协作加载 cache，INNER 循环被重构成 `k_base` 步进 + `jn` 标志 + `__any_sync` 早退，
   这套机制本身引入了发散控制开销，STACK spill 也从 24 B 涨到 40 B。

## 改动范围（已回退）

A 改动横跨 150 行 / 8 个 kernel（不只 LJ count，还有 expand/count/expand_count/versioned 的
bsearch 调用），并重构了 LJ INNER 循环。回退到 ed5ee8d（B 版 shared m[]）后代码恢复干净。

## 教训

和 v2（outer-par）、早期的 shared-mem 误判一样：**移植的优化假设必须先对照目标数据集的特性验证。**
GraphMiner 的 2-phase 对稠密图有效，但 LJ 是幂律稀疏图（中位度数 2），把并行/缓存投在
「长列表」上是错配——LJ 上根本没有长列表。

下一步真正的瓶颈见 ncu detailed 报告（`analysis/gpu_ncu/lj_count_DETAILED_9vQ49_shmem.ncu-rep`）：
**非合并访存（uncoalesced global access），预期加速 50.97%** —— INNER 循环的 `vlabels[w]`、
`csr_neighbors[ns2+k]` 是散点访问，53% sector 被浪费。这才是下一刀。
