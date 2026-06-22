# GPU BFS Search — optimization log

Workload: ParaCOSM `gpu_bfs` mode, Amazon (403K V, 2.2M E) + LiveJournal,
A100 80GB PCIe, CUDA 13.3. Golden = CPU `batch_all` match counts.

This log records what was tried on `core/gpu/gpu_bfs_search.cu`, the measured
result, and **why** — including one negative result (R2) kept on purpose so it
isn't re-attempted.

---

## Summary table

| Opt | What | Result | Committed |
|---|---|---|---|
| OPT-A | versioned `u_min` pivot: O(Σdeg) timestamp scan → O(1) physical-degree proxy | versioned −16% | `0ef6d0c` |
| OPT-B | warp-aggregated atomic for slot counters (1 atomic/warp vs 1/lane) | correct; neutral on 6v | `0ef6d0c` |
| **OPT-C** | **expand/count: 1-thread-per-pm → 1-warp-per-pm** | **expand −51%, ~1.6–1.85×** | `0ef6d0c` |
| gtid fix | 64-bit global thread index (was uint32) | fixes wrap at >2^27 pm | `0ef6d0c` |
| **P1a** | emitter counter `out_count`/`d_count_` → uint64 | **fixes 10v Q_3 / 11v Q_1 (were 2^32 garbage)** | `b44c084` |
| **P1b** | fuse final expand+count: largest layer never materialised | **fixes 9v Q_1 overflow drop; 9v Q_3 ~16s→1.3s** | `b44c084` |
| R2 | extend fusion to last K=2–3 layers (device recursion) | **NET LOSS 2.9–5.4× — reverted** | — |
| LJ-inner | sparse fused-count: outer-par → inner-par (parallelise the hub list) | **4.6× on LJ 6v Q_1; 11.2× warp util** | (sparse_mode_) |
| **A — launch_bounds** | `__launch_bounds__(256,8)` on LJ count kernel: REG 40→32, occupancy 53%→full | **+10–12% on heavy LJ 9v** | `054153a` |
| v2 outer-par | flip LJ inner-par → outer-par (kill 32× outer redundancy) | **NEG: 1.8–3.2× SLOWER — kept off** | `9a026ec` |
| ① shfl pivot | lane0 computes warp-uniform pivot + `__shfl` broadcast | **NEG: REG unchanged (40), ~0% — dropped** | — |
| specialization | `template<bool PADDED>` to fold the compact/padded branch | **NEG: REG unchanged (40), degrees ptr in const bank — dropped** | — |
| ③ label-sel pivot | pick pivot by label-effective candidates, not raw degree | **BLOCKED: LJ has 1 edge label (zero selectivity) — not pursued** | — |

Current HEAD `054153a`. Cross-scale regression (Amazon, gpu_bfs vs CPU):
6v/7v/9v/10v Q_1/Q_3 all exact and deterministic. LJ 9v Q_49/69/9/85 counts
bit-identical across every variant above (correctness never regressed).

---

## Why the wins worked

### OPT-C — warp-per-partial-match (the big one)
ncu, `bfs_expand_kernel`, depth 5→6, 8v Q_5 (reports in this dir):

| Metric | BEFORE (1-thr/pm) | AFTER (1-warp/pm) |
|---|---:|---:|
| Avg active threads / warp (of 32) | **2.40** | **15.56** |
| Duration | 10.46 ms | 5.10 ms (−51%) |
| SM throughput | 26.7% | 46.8% |
| Warp cycles / issued instr | 35.75 | 19.06 |
| DRAM throughput | 2.2% | 4.6% |

Root cause of the old slowness: **1 thread per partial match put 32 *unrelated*
pms in one warp.** Their candidate lists had wildly different lengths (degree
skew), so the warp waited on the slowest lane; and most candidates die at the
first label check (`continue`), collapsing the warp to ~2 active lanes.
The fix makes **32 lanes cooperate on ONE pm's candidate list** (strided
`i = lane; i < nbr_count; i += 32`): identical code path, balanced work, no
inter-pm divergence. Utilisation 2.4→15.6 is exactly why duration halves.

Note DRAM stays at 2–5% throughout — **this kernel is not memory-bound** (L2 hit
~97%). That is why `bfs_bsearch` memory tweaks (uint2 packing, `__restrict__`,
unrolling) were *not* pursued; they optimise the wrong resource.

### P1a — uint64 emitter counter (correctness)
A single dense expand layer can emit >2^32 candidates (10v Q_3, 11v Q_1). The
counter was uint32, so `atomicAdd` wrapped and the host read garbage
(10v Q_3 reported 4,294,974,963 ≈ 2^32 instead of 1,315,908,335). Promoting
`d_count_` + the four writer kernels' `out_count` + `warp_agg_inc` + host
reset/read paths to 64-bit fixed it. No perf change — pure correctness.

### P1b — fuse the final expand+count
The last expand normally materialises depth Q-1 (the single **largest** layer —
277M pms for 9v Q_1) into the ping-pong buffer, then a separate count kernel
reads it. `bfs_expand_count_kernel` extends each depth-(Q-2) candidate
in-register straight to depth Q-1 and counts, so that layer is never written.
Two wins: (1) removes a full global write+read of the biggest layer; (2)
eliminates the buffer-overflow *flush* that was silently dropping matches
(9v Q_1: 1,839,971,202 lossy/non-deterministic → exact 1,839,971,284).
Net positive because it fuses **only one** layer — the smallest tail layer,
cheap to serialise — while removing the largest materialisation.

---

## Why R2 (extend fusion to K=2–3) was a net loss — DO NOT REATTEMPT (as-is)

R2 generalised P1b: a device-recursive `gpu_extend_count(m_local, depth)` that
extends a partial match through the last **K** depths in-register. Goal: kill
the overflow at its source for 10v/11v (whose biggest layer sits at depth
Q-2..Q-3, which K=1 can't reach).

Correctness was fine (all goldens matched, 9v Q_1 overflow eliminated). But
**measured K=2 was 2.9–5.4× slower across every query**:

| Query | P1b (K=1) | R2 (K=2) |
|---|---:|---:|
| 7v Q_1 | 384 ms | 2091 ms |
| 8v Q_5 | 474 ms | 2307 ms |
| 8v Q_1 | 44 ms | 127 ms |
| 9v Q_1 | 1754 ms | 5446 ms |

**The principle (the takeaway):** fusing a layer trades its **32-lane parallel
materialisation** for **serial per-warp recursion**. Materialising a layer costs
one global write + read at ~2 TB/s — cheap. Serialising a high-fanout layer
collapses its parallelism (9v Q_1 depth-7 goes from 277M parallel items to 37.6M
warps each serially deep-searching). With a 400M buffer the materialise is
almost always cheaper. **Warp-local deep search only wins under genuine memory
pressure** (buffer can't hold the tail AND the fused layers are low-fanout) —
a regime A100 80GB + 400M buffer rarely enters.

Secondary finding: even forced to K=1, the generic device-recursion path
(`m_local[20]` copied to local memory + call overhead) ran ~2× slower than
P1b's hand-inlined register-only two-level kernel. Generality has a tax here.

Reverted via `git checkout` (no commit). A full warp-centric BFS (R1:
shared-mem frontier + global spill) would only pay off under the same memory
pressure that R2 needs — not worth it for the current hardware/data.

---

## Candidate next steps (not yet done)

- **Block-level count reduction.** Count/fused kernels currently do
  `warp_reduce_add_u64` + **1 atomicAdd per warp** to the single global
  `result_count`. The next increment is a block-level reduction (shared mem) →
  **1 atomicAdd per block**, cutting global atomics ~8× (256-thread block / 32).
  **MEASURED — not worth it (with P1b fusion on).** ncu of the fused count
  kernel `bfs_expand_count_kernel` (8v Q_5, `count_kernel_profile.csv`):
  DRAM throughput **0.22%**, memory throughput 22% → the global atomic is NOT a
  bottleneck (atomics travel the memory subsystem; a hot atomic shows as high
  memory throughput, which it isn't). The real limiter is **Avg active threads
  / warp = 4.52 / 32** (14%), caused by the kernel's nested candidate loops
  (outer v at depth Q-2, inner w at depth Q-1 with a `break` in joinability) —
  the warp re-diverges in the inner loop. So the lever is *warp utilisation*,
  not atomic count. (And improving the inner loop = R2 territory, already shown
  to lose.) Block-reduce would optimise something already <1% of runtime.
- `__launch_bounds__` + block-size sweep (occupancy dipped 68→56% after OPT-C).
- Apply P1b fusion to the versioned path (currently non-versioned only).

---

## LiveJournal vs Amazon — two different regimes (per-graph behaviour)

LJ (4.85M V, 38.6M E, avg deg 15.9, power-law social graph) behaves very
differently from Amazon (403K V, 2.2M E, co-purchase, more uniform degree).

Measured (gpu_bfs, golden = CPU batch_all where it finishes):

| Case | unsafe edges | max mid layer | overflow? | matches | gpu search |
|---|---:|---:|:---:|---:|---:|
| LJ 6v Q_1 (500k) | 5,691 | 355K | no | 64.4M | 224 ms |
| LJ 7v Q_1 (500k) | 12,962 | 224K | no | 24.5M | 31 ms |
| LJ 8v Q_1 (500k) | — | — | no | 30.9B | 35.7 s (CPU T/O) |
| LJ 9v Q_1 (100k) | 1,888 | 43.5M | **no** | **17.1B** | **22.3 s** |
| AZ 8v Q_5 (full) | 88K | 355M | **yes (flush)** | 1.40B | 0.47 s |
| AZ 10v Q_3 (full) | — | >400M | **yes** | 1.32B | 2.5 s |

**The regimes are opposite:**

- **Amazon**: many unsafe edges, **mid layers explode to 10^8–10^9 → buffer
  overflow**. Bottleneck = frontier materialisation + flush. OPT-C (warp coop),
  P1a (uint64 count), P1b (fuse largest layer) all target *this*.
- **LJ (sparse, power-law)**: very few unsafe edges, **mid layers stay small
  (10^5–10^7), never overflow** — but the **final match count explodes**
  (1888 unsafe edges → 17.1B matches; ~9M matches per edge). Bottleneck moves
  to the **last-level count kernel's raw enumeration** (22 s of pure counting),
  not memory at all.

### Implication for "one BFS kernel per graph"
The data supports specialising, but the cleaner form is **runtime-adaptive**
rather than two hand-maintained kernels:

| Signal (cheap to measure at runtime) | LJ-like | Amazon-like |
|---|---|---|
| unsafe-edge count after Classify | low (10^3) | high (10^5) |
| mid-layer growth rate (out/in per level) | moderate | explosive |
| Best strategy | skip fusion/buffer mgmt; **optimise the count kernel** (the hot path) — e.g. better warp utilisation in the nested count loop, or a count-only path that avoids the fused inner re-divergence | warp-coop + P1b fusion + buffer/uint64 (already done) |

So the next LJ-specific win is **not** more buffer/overflow work (LJ never
overflows) — it's making the **terminal count enumeration** faster. That is the
kernel whose ncu showed Avg active threads/warp = 4.52 (the nested inner loop).
For LJ this kernel is ~100% of runtime, so improving its warp utilisation would
pay off here even though it was <1% atomic on Amazon.

(Caveat: improving the nested count inner loop is adjacent to R2's failed
territory — must be measured, not assumed. But the lever is utilisation of the
inner candidate sweep, not atomics.)

---

## LJ count-kernel profile → the LJ-specific bottleneck (ncu, 7v Q_1)

`bfs_expand_count_kernel` is ~100% of LJ runtime (LJ never overflows; the cost
is terminal enumeration). ncu (`LJ_count_7vQ1.csv`):

| Metric | LJ count | (AZ count for ref) |
|---|---:|---:|
| **Avg active threads / warp (of 32)** | **1.64** | 4.52 |
| Not-predicated-off threads / warp | 1.49 | — |
| DRAM throughput | **0.02%** | 0.22% |
| Compute (SM) throughput | 58.9% | 50% |
| Achieved occupancy | 53% | 50% |
| Block limit (registers) | 6 of 8 | — |

**Bottleneck = warp under-utilisation (1.64/32 = 5%), NOT memory (DRAM 0.02%).**

Why LJ is *worse* than Amazon here: the kernel parallelises the **outer**
candidate list (order[Q-2]'s neighbours) across the 32 lanes, then each lane
serially sweeps the **inner** list (order[Q-1]). On Amazon (mid-density) the
outer list is long → lanes have work (4.52). On **LJ (sparse, heavy label
filtering)** the outer list is short → most lanes fail `i < nbr_count`
immediately and idle; 1–2 lanes do the deep inner sweep → 1.64.

So the wrong dimension is parallelised for sparse graphs. **LJ-specific kernel
directions (to be profiled/measured, not assumed — R2 rule):**
1. **Inner-parallel**: serialise outer v (few on LJ), split the inner w list
   across 32 lanes (longer on LJ, esp. near hubs).
2. **Flattened work-list**: collect valid (v) first, then 32-lane over (v,w)
   pairs — balances regardless of which list is short.
3. **warp-per-edge** instead of warp-per-pm: LJ has very few unsafe edges
   (~1888) but huge per-edge match counts — regranularise.

Register pressure (block limit 6/8) also says the kernel is heavy; an
inner-parallel rewrite that drops the per-lane `m_local` footprint could lift
occupancy too. This is the basis for a **separate LJ (sparse-graph) BFS kernel**
selected at runtime by unsafe-edge count / mid-layer growth.

---

## LJ-specific kernel: inner-parallel fused count — WIN (4.6×)

Implemented `bfs_expand_count_lj_kernel`: same fused expand+count semantics as
`bfs_expand_count_kernel`, but parallelism flipped — the warp walks the OUTER
candidate (order[Q-2]) uniformly and splits the INNER candidate list
(order[Q-1]) across 32 lanes. Selected at runtime in BuildCSR via `sparse_mode_`
(heuristic: V ≥ 1M and avg degree < 32 → LJ-like; override with env
`GPU_BFS_SPARSE=0/1`).

A/B on the SAME case (LJ 6v Q_1, 500k stream):

| Fused-count kernel | BFS search |
|---|---:|
| outer-parallel (Amazon default) | 223.1 ms |
| **inner-parallel (LJ)** | **48.5 ms** |
| | **4.6× faster** |

ncu (LJ 7v Q_1), same kernel slot:

| | outer-par | inner-par (LJ) |
|---|---:|---:|
| Duration | 38.8 ms | **15.5 ms (−60%)** |
| Compute (SM) throughput | 58.9% | 61.7% |
| **Avg active threads / warp** | **1.64** | **18.39 (11.2×)** |
| DRAM | 0.02% | 0.05% |

Correctness: LJ 6v/7v Q_1 match golden exactly (64,404,068 / 24,537,556).
Amazon unaffected — auto-selects outer-par, 8v Q_5 still 1,399,751,223 @ 470ms.

**Why this worked where R2 didn't:** both were "warp-level" reworks, but R2
changed parallelism *blindly* (fused more layers, serialised a high-fanout one →
loss). This change was *profile-directed*: ncu said the bottleneck was the wrong
parallel dimension on sparse graphs (outer list short → 1.64/32 lanes), so we
parallelised the dimension that is actually long on LJ (the inner/hub list).
Measure first, then change the thing the measurement points at.

---

## Round 3 — squeezing the LJ fused-count kernel (2026-06)

Goal: make `bfs_expand_count_lj_kernel` (the inner-parallel LJ count kernel)
faster. Four ideas tried; one won (+11%), three were killed by measurement
before wasting an implementation. The discipline that paid off: **a cheap
diagnostic (ncu, register dump, or a deliberately-wrong short-circuit build)
before each real change.**

### ncu diagnosis: the kernel is memory-latency-bound, not compute/bandwidth

ncu on the LJ count kernel (A100, RmProfilingAdminOnly=1 so `sudo ncu`):

| Metric | Value | Reading |
|---|---:|---|
| DRAM throughput | 0.01% | not bandwidth-bound |
| Compute (SM) | 50% | not compute-bound |
| L1 hit | 81.6% | data mostly cached |
| Achieved occupancy | 53% (theo. 75%) | **occupancy-limited** |
| Registers/thread | 40 → Block Limit 6 | **the occupancy cap** |
| Warp cycles/issued instr | 16.9 | — |
| └ `long_scoreboard` | **9.15 (54%)** | **waiting on global loads** |
| └ `wait` | 3.80 | exec dependency |

So: warps stall waiting on global loads, and occupancy is too low (53%) to
cover that latency. The lever is occupancy, and occupancy is gated by registers.

### A (WIN, `054153a`): `__launch_bounds__(256, 8)` → +10–12%

Forcing ptxas to a register budget that admits 8 blocks×256 = 2048 threads/SM
(full occupancy on A100) drops REG 40→32. Costs a 32 B stack spill, but spilled
slots hit L1 while the extra resident warps hide the load latency.

| query | matches | v1 | A | speedup |
|---|---:|---:|---:|---:|
| 9v Q_49 | 20.99 B | 12041 ms | 10985 ms | 1.10× |
| 9v Q_69 |  9.39 B |  3641 ms |  3244 ms | 1.12× |
| 9v Q_9  |  6.46 B |  2763 ms |  2484 ms | 1.11× |
| 9v Q_85 | 349 K   |   840 ms |   838 ms | ~1.0× |

Counts bit-identical. This is the only Round-3 win.

### ① (NEG, dropped): lane-0 pivot + `__shfl` broadcast

The `u_min`/`u2_min` pivots are warp-uniform (all 32 lanes compute the same
scan). Hypothesis: compute on lane 0, broadcast, save registers + redundant
work. Reality: **REG stayed 40** — the pivot temporaries were never the live-set
peak (that's the INNER candidate path), so occupancy didn't move. ~0% change.
Confirmed the peak is INNER, not the pivot.

### specialization (NEG, dropped): `template<bool PADDED>`

LJ builds a compact CSR (`bfs_use_padded==0`), so `(bfs_use_padded && degrees) ?
degrees[x] : offsets[x+1]-off` is always the offsets branch. Templating on
PADDED=false to fold the branch at compile time: **REG stayed 40**. The `degrees`
pointer lives in the constant bank (not a general register), and both branches
need `off`/`offsets` anyway, so folding freed nothing. Dropped.

### B (NEG by diagnosis, not built): cache m[] neighbour lists in shared memory

Intuition (correct on the data): 99.1% of LJ vertices have degree ≤32, 99.8%
≤64, so the neighbour lists the INNER joinability bsearch hits are tiny and
could sit in shared memory. **But a cheap diagnostic killed it first:** a build
that short-circuits the entire INNER bsearch (deliberately wrong counts, timing
only) measured how much of the kernel that bsearch actually costs:

| query | A (full) | bsearch short-circuited | bsearch share |
|---|---:|---:|---:|
| 9v Q_49 | 10985 ms | 10520 ms | **4%** |
| 9v Q_69 |  3244 ms |  2964 ms | **9%** |
| 9v Q_9  |  2484 ms |  2336 ms | **6%** |

The bsearch is only 4–9% of the kernel (small tables already hit L1). Shared
memory could optimise a fraction of that, while risking the occupancy that A
just bought. Not worth it. **The other 91–96% is the INNER sweep itself —
reading the w candidates (`csr_neighbors[w]`, `vlabels[w]`) over a long list.**
The bottleneck is the *number of candidates*, not the cost per candidate.

### ③ (BLOCKED): label-selectivity pivot

Follow-on from B: reduce the candidate count by choosing the pivot whose
label-effective fanout is smallest, not its raw degree. Blocked by the data —
**LJ has exactly 1 edge label** (all edges label 0; 30 vertex labels, uniform).
The edge-label filter in the kernel never prunes anything, so the `(edge_label,
vertex_label)` joint selectivity collapses to vertex-label only. And the INNER
sweep must still scan all `nc2` neighbours to *count* the label-effective ones —
labels can't skip the traversal without a per-label bucketed neighbour index
(a large data-structure change). Not pursued on LJ; would need a multi-edge-label
dataset to even measure a benefit.

### Round-3 takeaway

The kernel is at a real floor for this approach: latency-bound, occupancy now
full (A), bsearch negligible (B), candidate-count irreducible without changing
the on-device data structure or the dataset's label profile (③). Further GPU
gains need a different lever (cut the total matches to count, or a bucketed
index), not a kernel micro-rewrite. Three of four ideas were retired by a
sub-5-minute diagnostic instead of a full build-test-revert cycle.
