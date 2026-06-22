# GPU BFS expand kernel — Nsight Compute profiles

Profiles of `bfs_expand_kernel` on the **Amazon 8v Q_5** workload (50k insertion
stream subset), capturing the **depth 5→6** expand launch (`--launch-skip 3
--launch-count 1`, `--set full`). A100 80GB PCIe, CUDA 13.3, Nsight Compute 2026.2.

The two reports bracket the warp-cooperative rework (`OPT(C)`):

| File | Code state | Scheme |
|---|---|---|
| `expand_BEFORE_baseline_8vQ5.ncu-rep` | commit `0ef6d0c^` (pre-OPT-C) | 1 thread per partial match |
| `expand_AFTER_optc_8vQ5.ncu-rep` | commit `b44c084` (current) | 1 warp (32 lanes) per partial match |

## How to open

```bash
# GUI (best — full sections, source/SASS, roofline)
ncu-ui analysis/gpu_ncu/expand_AFTER_optc_8vQ5.ncu-rep

# or open both and use the GUI "Baseline" feature to diff:
ncu-ui analysis/gpu_ncu/expand_BEFORE_baseline_8vQ5.ncu-rep \
       analysis/gpu_ncu/expand_AFTER_optc_8vQ5.ncu-rep

# CLI summary of any report
ncu --import analysis/gpu_ncu/expand_AFTER_optc_8vQ5.ncu-rep --page details
```

## Headline comparison (depth 5→6 expand)

| Metric | BEFORE (1-thread/pm) | AFTER (1-warp/pm) | Change |
|---|---:|---:|---|
| **Avg. Active Threads / Warp** (of 32) | **2.40** | **15.56** | **+548%** ⬆ |
| Kernel Duration | 10.46 ms | 5.10 ms | **−51%** ⬇ |
| Compute (SM) Throughput | 26.7% | 46.8% | +75% ⬆ |
| Warp Cycles / Issued Instruction | 35.75 | 19.06 | −47% ⬇ (less stalling) |
| Achieved Occupancy | 68.6% | 55.8% | −19% (acceptable: each warp now heavier) |
| DRAM Throughput | 2.20% | 4.59% | both negligible |

## Reading

- **Active-threads/warp was the bottleneck.** The old scheme left ~30 of 32
  lanes idle while a few chewed through high-degree candidate lists. The warp
  rework lifts effective lane utilisation 6.5× (2.40 → 15.56), which is exactly
  why duration halves and SM throughput jumps.
- **Not a memory-bound kernel.** DRAM throughput is ~2–5% in both — the candidate
  data lives in L1/L2 (L2 hit ~97%). This is why optimising `bfs_bsearch`
  memory access (uint2 packing, `__restrict__`, unrolling) was *not* pursued.
- **Occupancy dipped** (more registers per warp) but the trade is a clear net
  win because the active lanes are doing useful work instead of idling.

## Reproduce

```bash
source ~/intel/oneapi/setvars.sh
export PATH="$HOME/.local/bin:/usr/local/cuda/bin:$PATH"
D=~/haibin/az_data/AZ/data_graph/data.graph
Q=~/haibin/az_data/AZ/8v/Q_5
U=/tmp/az_stream_50k.graph   # head -50000 of AZ/data_graph/insertion.graph

sudo ncu --set full --kernel-name bfs_expand_kernel \
  --launch-skip 3 --launch-count 1 --force-overwrite \
  --export analysis/gpu_ncu/expand_AFTER_optc_8vQ5 \
  ParaCOSM/CSM/build/bin/csm -q $Q -d $D -u $U \
  -a parallel_graphflow -m gpu_bfs -t 8 --report-initial 0
```
