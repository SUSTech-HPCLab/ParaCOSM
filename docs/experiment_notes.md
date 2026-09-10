# ParaCOSM Experiment Notes

## Build record: 2026-09-09

### Repository

- Commit before local build fixes: `e3e0814`
- Working directory: `/home/haibin/tpds/ParaCOSM/ParaCOSM/CSM`
- Build directory: `build-v100-gcc10`
- Main binary: `build-v100-gcc10/bin/csm`

### Host

- CPU: 8 x Intel Xeon Platinum 8160
- Physical cores: 192 (24 per socket)
- Hardware threads: 384
- NUMA nodes: 8
- Memory: approximately 96 GB per NUMA node, approximately 768 GB total
- OS toolchain: GCC/G++ 10.5 for this build
- OpenMP: GNU OpenMP 4.5 (`libgomp`)
- TBB: oneTBB 2021.5 system package; runtime resolves to `libtbb.so.12`

The logical CPU ranges reported by `numactl --hardware` are:

| NUMA node | Logical CPUs |
|---:|---|
| 0 | 0--47 |
| 1 | 48--95 |
| 2 | 96--143 |
| 3 | 144--191 |
| 4 | 192--239 |
| 5 | 240--287 |
| 6 | 288--335 |
| 7 | 336--383 |

Each range contains both hardware threads of the node's 24 physical cores.
Experimental CPU affinity must therefore distinguish physical cores from SMT
siblings; `-t 192` without explicit placement is not a sufficient description.

### GPUs

- 8 x NVIDIA Tesla V100 PCIe 32 GB
- Compute capability: 7.0 (`sm_70`)
- Driver: 580.173.02
- CUDA Toolkit used for compilation: 11.5.119
- Current implementation is single-GPU; select a device with
  `CUDA_VISIBLE_DEVICES=<index>`.

### Configure and build

CUDA 11.5 fails with the system-default GCC 11 because NVCC cannot compile
GCC 11's `serializeintrin.h`. GCC 10 is installed and builds successfully:

```bash
cd /home/haibin/tpds/ParaCOSM/ParaCOSM/CSM
CC=/usr/bin/gcc-10 CXX=/usr/bin/g++-10 cmake -S . -B build-v100-gcc10 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=70 \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-10
cmake --build build-v100-gcc10 -j 16
```

Built targets:

- `csm`
- `calig`
- `parallel_calig`
- `parallel_newsp`

### Smoke test

A three-vertex path query, a data graph containing its first edge, and one edge
insertion completing the path were used. Initial matching was disabled. All
three modes reported exactly two positive matches and zero negative matches:

| Mode | Algorithm | Positive | Negative | Result |
|---|---|---:|---:|---|
| `single` | `parallel_graphflow` | 2 | 0 | pass |
| `versioned` | `parallel_graphflow` | 2 | 0 | pass |
| `gpu_bfs_versioned` | `parallel_graphflow` | 2 | 0 | pass on V100 GPU 0 |

During this test, the GPU versioned path exposed and fixed a 32/64-bit
allocation mismatch for `d_count_`. The counter is declared and copied as 64
bits, but the versioned allocation still reserved only 32 bits. Before the fix,
CUDA returned `invalid argument` and the mode incorrectly reported zero
matches.

This smoke test verifies build, dynamic linking, CUDA architecture selection,
kernel launch, and one minimal count result. It is not the full correctness
validation required by E0 in `docs/tpds_revision_matrix.md`.

### Remaining build/runtime warnings

- Several slot-pool functions trigger `may be used uninitialized` warnings;
  their loops appear to assign the value before returning, but should be made
  warning-clean before the final artifact build.
- `Parallel_NewSP::Timer::countUillisSecond()` has no return statement. It is
  currently unused but is undefined behavior if called.
- The GPU BFS implementation allocates three frontier buffers using a fixed
  `MAX_BUF_MATCHES=400,000,000`. Depending on query size this can reserve most
  of a 32 GB V100, even for a tiny input. Buffer capacity should become a
  runtime/memory-aware setting before large V100 experiments.
- CMake still requires CUDA for the unified `csm` target. A separate CPU-only
  build option remains to be implemented.
