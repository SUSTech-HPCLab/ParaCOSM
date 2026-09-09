# GPU Kernel 优化任务总览

## 目标文件
`/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/core/gpu/gpu_bfs_search.cu`

## 背景
ncu profiling 显示 `bfs_expand_count_lj_kernel` 是 memory-latency-bound：
- 9.15/17 stall cycles 为 `long_scoreboard`（等待 global memory load）
- 热点：joinability 检查中的 `bfs_bsearch` 密集调用 global memory
- 当前 occupancy 被寄存器用量限制（40 regs → 53%，优化到 32 regs 后 +10-12%）

## 已有优化（勿重复）
- OPT(B): warp-aggregated atomic increment
- OPT(C): warp-per-partial-match (32 lanes 协作)
- OPT(P1b): fused expand+count (最后两层融合)
- LJ kernel: inner-parallel 设计（INNER 长列表分 lanes，OUTER 短列表 warp-uniform）
- LJ kernel: `__launch_bounds__(256, 8)` 强制 32 regs
- LJ kernel: shared memory 缓存 m[]

## 子任务分配

| Agent | 任务 | 文档路径 | 改动区域 |
|-------|------|----------|----------|
| A | 2-Phase Binary Search | [TASK_A.md](./TASK_A.md) | `bfs_bsearch` 函数 + joinability 循环 |
| B | Shared Memory m[] + launch_bounds | [TASK_B.md](./TASK_B.md) | kernel 入口声明区域 |
| C | BlockReduce + Early Termination | [TASK_C.md](./TASK_C.md) | kernel 归约出口 + sweep 循环 |

## 冲突规避
- **Agent A** 只修改 `bfs_bsearch` 函数体及其调用方式（添加 cache 参数）
- **Agent B** 只修改 kernel 函数的开头（shared memory 声明、`__launch_bounds__` 注解）
- **Agent C** 只修改 kernel 函数的末尾归约 + sweep 循环中的 early-exit 逻辑

三者改动位置不重叠，可并行开发。

## 优先级
1. Task A (最高) — 直接缓解 profiling 发现的 #1 瓶颈
2. Task B (中) — 低风险移植已验证 pattern
3. Task C (低) — 实验性优化

## 参考资料
- GraphMiner 优化技术文档: `/home/v-haibinlai/haibin/ParaCOSM/docs/gpu_kernel_optimization_techniques.md`
- GraphMiner 源码: `/home/v-haibinlai/haibin/GraphMiner/include/set_intersect.cuh`
- GraphMiner search: `/home/v-haibinlai/haibin/GraphMiner/include/search.cuh`
