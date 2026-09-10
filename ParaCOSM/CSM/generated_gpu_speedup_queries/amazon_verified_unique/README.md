# Amazon full-stream GPU-speedup queries

These five pairwise non-isomorphic labeled query graphs were verified using
all 244,341 Amazon insertion updates.

- CPU: `versioned`, 16 OpenMP threads, NUMA node 0
- GPU: `gpu_bfs_versioned`, one V100 (GPU 0), NUMA node 0
- Correctness: positive-match counts are exactly equal and negative counts are zero
- Timing: the `Incremental Matching` time reported by CSM

See `manifest.csv` for the measured results. Raw logs and run configurations
remain in the sibling `amazon_full_*` directories.
