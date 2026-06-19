#include "core/gpu/gpu_bfs_search.h"
#include "graph_storage/graph.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <omp.h>

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,  \
                    cudaGetErrorString(err));                                   \
        }                                                                      \
    } while (0)

// ============================================================
// Query graph in constant memory
// ============================================================
static constexpr int BFS_MAX_Q = 20;
// Directed neighbor entries = sum of query-vertex degrees = 2 * (#query edges).
// Must cover the worst case: a complete graph on BFS_MAX_Q vertices has
// BFS_MAX_Q*(BFS_MAX_Q-1) = 20*19 = 380 directed entries. Was 64, which silently
// overflowed constant memory for dense queries (e.g. a 9v/33e near-complete query
// needs 66 > 64), corrupting bfs_q_elabels and breaking joinability pruning ->
// match explosion. SetupQuery now also asserts against this bound.
static constexpr int BFS_MAX_Q_EDGES = 384;

__constant__ uint32_t bfs_q_num_vertices;
__constant__ uint32_t bfs_q_vlabels[BFS_MAX_Q];
__constant__ uint32_t bfs_q_offsets[BFS_MAX_Q + 1];
__constant__ uint32_t bfs_q_neighbors[BFS_MAX_Q_EDGES];
__constant__ uint32_t bfs_q_elabels[BFS_MAX_Q_EDGES];

// Flag: 0=compact CSR (use offsets for degree), 1=padded CSR (use d_degrees)
__constant__ uint32_t bfs_use_padded;

// ============================================================
// Device: binary search
// ============================================================
__device__ __forceinline__ int bfs_bsearch(
    const uint32_t* arr, int size, uint32_t target
) {
    int lo = 0, hi = size - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        uint32_t val = arr[mid];
        if (val == target) return mid;
        if (val < target) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}

// ============================================================
// OPT(B): warp-aggregated atomic increment.
// Each active lane reserves a unique slot, but only ONE atomicAdd is issued
// per warp-group instead of one per lane. Correct inside divergent loops:
// __activemask() captures exactly the lanes converged at this call site.
// The set of returned slots is identical to per-lane atomicAdd(ctr,1) (a
// permutation), and the final counter value is unchanged.
// ============================================================
__device__ __forceinline__ uint64_t warp_agg_inc(unsigned long long* ctr) {
    uint32_t mask = __activemask();
    uint32_t lane = threadIdx.x & 31u;
    uint32_t leader = __ffs(mask) - 1;
    uint32_t rank = __popc(mask & ((1u << lane) - 1u)); // active lanes before me
    unsigned long long base = 0;
    if (lane == leader) {
        base = atomicAdd(ctr, (unsigned long long)__popc(mask)); // one atomic for the group
    }
    base = __shfl_sync(mask, base, leader);             // broadcast group base
    return (uint64_t)base + rank;
}

// ============================================================
// OPT(C): warp-level sum reduction of a per-lane uint64 counter.
// Standard butterfly reduction. Requires the whole warp converged here, which
// the callers guarantee: the only per-warp early-return (`warp_id >= in_count`
// and `u_min == UINT32_MAX`) is uniform across all 32 lanes. Lane 0 returns the
// warp total.
// ============================================================
__device__ __forceinline__ uint64_t warp_reduce_add_u64(uint64_t v) {
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1)
        v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

// ============================================================
// Kernel 1: Init — create depth=2 partial matches
// One thread per (edge_idx * num_q_edges * 2 + qe_idx * 2 + dir)
// ============================================================
__global__ void bfs_init_kernel(
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ edges_v1,
    const uint32_t* __restrict__ edges_v2,
    const uint32_t* __restrict__ edges_label,
    const uint32_t* __restrict__ all_orders,
    uint32_t* __restrict__ out_buf,      // output partial matches
    unsigned long long* __restrict__ out_count,    // atomic counter (uint64)
    uint32_t num_data_edges,
    uint32_t num_q_edges,
    uint32_t Q,
    uint32_t max_out
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t total = num_data_edges * num_q_edges * 2;
    if (idx >= total) return;

    uint32_t edge_idx = idx / (num_q_edges * 2);
    uint32_t remainder = idx % (num_q_edges * 2);
    uint32_t qe_idx = remainder / 2;
    uint32_t dir = remainder % 2;

    uint32_t v1 = edges_v1[edge_idx];
    uint32_t v2 = edges_v2[edge_idx];
    uint32_t elabel = edges_label[edge_idx];

    const uint32_t* order = all_orders + qe_idx * Q;
    uint32_t u1 = order[0];
    uint32_t u2 = order[1];

    // Find edge label between u1 and u2 in query
    uint32_t qs = bfs_q_offsets[u1];
    uint32_t qe = bfs_q_offsets[u1 + 1];
    uint32_t found_el = UINT32_MAX;
    for (uint32_t j = qs; j < qe; j++) {
        if (bfs_q_neighbors[j] == u2) { found_el = bfs_q_elabels[j]; break; }
    }
    if (found_el == UINT32_MAX) return;

    // Direction
    uint32_t mv1 = dir == 0 ? v1 : v2;
    uint32_t mv2 = dir == 0 ? v2 : v1;

    // Label check
    if (vlabels[mv1] != bfs_q_vlabels[u1]) return;
    if (vlabels[mv2] != bfs_q_vlabels[u2]) return;
    if (elabel != found_el) return;

    // Create partial match
    uint64_t slot = warp_agg_inc(out_count);  // OPT(B): warp-aggregated atomic
    if (slot >= max_out) return;  // buffer full

    uint32_t stride = Q + 1;
    uint32_t* out = out_buf + slot * stride;
    out[0] = qe_idx;  // order index
    for (uint32_t j = 0; j < Q; j++) out[1 + j] = UINT32_MAX;
    out[1 + u1] = mv1;
    out[1 + u2] = mv2;
}

// ============================================================
// Kernel 2: Expand — expand one BFS level (depth d → d+1)
// One thread per input partial match.
// For each valid candidate at depth d, writes a new partial match.
// Supports both compact CSR (degrees from offsets) and padded CSR (d_degrees array).
// ============================================================
__global__ void bfs_expand_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,  // NULL for compact CSR
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint32_t* __restrict__ out_buf,
    unsigned long long* __restrict__ out_count,
    uint32_t depth,
    uint32_t Q,
    uint32_t max_out
) {
    // OPT(C): warp-per-partial-match. 32 lanes cooperatively sweep the candidate
    // list of ONE partial match. Profile showed 2.35/32 active threads with the
    // old 1-thread-per-pm scheme — most lanes idled while a few chewed through
    // high-degree candidate lists. Both early-returns below are warp-uniform
    // (warp_id and u_min are identical across the warp), so the warp stays
    // converged into the strided loop where __ballot_sync is safe.
    // 32 lanes cooperatively sweep one partial match's candidate list.
    // gtid MUST be 64-bit: at >2^27 partial matches, cur_count*32 exceeds
    // UINT32_MAX and a 32-bit gtid wraps, aliasing warp_ids → double counts.
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 1;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    const uint32_t* m = pm + 1;  // m[0..Q-1]

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u = order[depth];

    // Find u_min: matched query neighbor of u with smallest degree (warp-uniform)
    uint32_t u_min = UINT32_MAX;
    uint32_t u_min_label = 0;
    uint32_t u_min_deg = UINT32_MAX;

    uint32_t qs = bfs_q_offsets[u];
    uint32_t qe = bfs_q_offsets[u + 1];

    for (uint32_t j = qs; j < qe; j++) {
        uint32_t u_other = bfs_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        uint32_t off = csr_offsets[m[u_other]];
        uint32_t deg = (bfs_use_padded && degrees) ? degrees[m[u_other]] : (csr_offsets[m[u_other] + 1] - off);
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = bfs_q_elabels[j];
        }
    }

    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = (bfs_use_padded && degrees) ? degrees[m[u_min]] : (csr_offsets[m[u_min] + 1] - nbr_start);

    // Strided sweep: every lane participates in each round (i may be out of range,
    // but the lane still reaches __ballot_sync so the warp never deadlocks).
    for (uint32_t base = 0; base < nbr_count; base += 32) {
        uint32_t i = base + lane;
        bool emit = false;
        uint32_t v = 0;

        if (i < nbr_count) {
            v = csr_neighbors[nbr_start + i];
            // 1. Label check
            if (vlabels[v] == bfs_q_vlabels[u] &&
                csr_elabels[nbr_start + i] == u_min_label) {
                // 2. Joinability
                bool joinable = true;
                for (uint32_t j = qs; j < qe; j++) {
                    uint32_t u_other = bfs_q_neighbors[j];
                    if (m[u_other] == UINT32_MAX || u_other == u_min) continue;
                    uint32_t s = csr_offsets[m[u_other]];
                    uint32_t nc = (bfs_use_padded && degrees) ? degrees[m[u_other]] : (csr_offsets[m[u_other] + 1] - s);
                    int pos = bfs_bsearch(csr_neighbors + s, nc, v);
                    if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) {
                        joinable = false; break;
                    }
                }
                // 3. Visited check
                if (joinable) {
                    bool visited = false;
                    for (uint32_t d = 0; d < Q; d++) {
                        if (m[d] == v) { visited = true; break; }
                    }
                    emit = !visited;
                }
            }
        }

        // 4. Warp-cooperative write via __syncwarp-fenced ballot.
        // __syncwarp() forces all 32 lanes to reconverge after the divergent
        // joinability/visited checks above, so the full-mask ballot/shfl below
        // see a consistent active set (required on Volta+ independent scheduling).
        __syncwarp();
        uint32_t emit_mask = __ballot_sync(0xffffffffu, emit);
        uint32_t cnt = __popc(emit_mask);
        if (cnt) {
            unsigned long long base_slot = 0;
            uint32_t leader = __ffs(emit_mask) - 1;
            if (lane == leader) base_slot = atomicAdd(out_count, (unsigned long long)cnt);
            base_slot = __shfl_sync(0xffffffffu, base_slot, leader);
            if (emit) {
                uint32_t rank = __popc(emit_mask & ((1u << lane) - 1u));
                uint64_t slot = base_slot + rank;
                if (slot < max_out) {  // buffer full — skip write but keep counting
                    uint32_t* out = out_buf + (size_t)slot * stride;
                    out[0] = order_idx;
                    for (uint32_t j = 0; j < Q; j++) out[1 + j] = m[j];
                    out[1 + u] = v;
                }
            }
        }
    }
}

// ============================================================
// Kernel 3: Count — at the last level, count matches instead of expanding
// ============================================================
__global__ void bfs_count_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,  // NULL for compact CSR
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint64_t* __restrict__ result_count,
    uint32_t depth,
    uint32_t Q
) {
    // OPT(C): warp-per-partial-match (see expand kernel). 32 lanes sweep one
    // partial match's candidate list; each lane accumulates locally, then a
    // single warp reduction + one atomicAdd folds the warp's total in.
    // gtid is 64-bit to avoid wrap at >2^27 partial matches (cur_count*32 > 2^32).
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 1;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    const uint32_t* m = pm + 1;

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u = order[depth];

    uint32_t u_min = UINT32_MAX;
    uint32_t u_min_label = 0;
    uint32_t u_min_deg = UINT32_MAX;

    uint32_t qs = bfs_q_offsets[u];
    uint32_t qe = bfs_q_offsets[u + 1];

    for (uint32_t j = qs; j < qe; j++) {
        uint32_t u_other = bfs_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        uint32_t off = csr_offsets[m[u_other]];
        uint32_t deg = (bfs_use_padded && degrees) ? degrees[m[u_other]] : (csr_offsets[m[u_other] + 1] - off);
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = bfs_q_elabels[j];
        }
    }

    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = (bfs_use_padded && degrees) ? degrees[m[u_min]] : (csr_offsets[m[u_min] + 1] - nbr_start);
    uint64_t local_count = 0;

    for (uint32_t i = lane; i < nbr_count; i += 32) {
        uint32_t v = csr_neighbors[nbr_start + i];

        if (vlabels[v] != bfs_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;

        bool joinable = true;
        for (uint32_t j = qs; j < qe; j++) {
            uint32_t u_other = bfs_q_neighbors[j];
            if (m[u_other] == UINT32_MAX || u_other == u_min) continue;
            uint32_t s = csr_offsets[m[u_other]];
            uint32_t nc = (bfs_use_padded && degrees) ? degrees[m[u_other]] : (csr_offsets[m[u_other] + 1] - s);
            int pos = bfs_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;

        bool visited = false;
        for (uint32_t d = 0; d < Q; d++) {
            if (m[d] == v) { visited = true; break; }
        }
        if (visited) continue;

        local_count++;
    }

    // Warp reduction: lane 0 issues a single atomicAdd for the whole warp.
    local_count = warp_reduce_add_u64(local_count);
    if (lane == 0 && local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// OPT(P1b): Fused expand+count for the last TWO query depths.
// Input: partial matches at depth Q-2. For each candidate v at depth Q-2 the
// thread does NOT materialise a partial match; instead it extends in-register
// to depth Q-1 and counts directly. This avoids materialising the largest BFS
// layer (depth Q-1), which is both the main buffer-overflow trigger and a full
// global read+write pass. One warp per input partial match; each lane owns a
// slice of the depth-(Q-2) candidate list and keeps a private running count.
// ============================================================
__global__ void bfs_expand_count_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,  // NULL for compact CSR
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint64_t* __restrict__ result_count,
    uint32_t depth,   // = Q-2 (the depth of the input partial matches' frontier)
    uint32_t Q
) {
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 1;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    const uint32_t* m = pm + 1;

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u  = order[depth];      // query vertex to map at depth Q-2
    uint32_t u2 = order[depth + 1];  // query vertex to map at depth Q-1 (final)

    // --- pivot u_min for u (depth Q-2) ---
    uint32_t u_min = UINT32_MAX, u_min_label = 0, u_min_deg = UINT32_MAX;
    uint32_t qs = bfs_q_offsets[u], qe = bfs_q_offsets[u + 1];
    for (uint32_t j = qs; j < qe; j++) {
        uint32_t uo = bfs_q_neighbors[j];
        if (m[uo] == UINT32_MAX) continue;
        uint32_t off = csr_offsets[m[uo]];
        uint32_t deg = (bfs_use_padded && degrees) ? degrees[m[uo]] : (csr_offsets[m[uo] + 1] - off);
        if (deg < u_min_deg) { u_min_deg = deg; u_min = uo; u_min_label = bfs_q_elabels[j]; }
    }
    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = (bfs_use_padded && degrees) ? degrees[m[u_min]] : (csr_offsets[m[u_min] + 1] - nbr_start);

    // query-neighbour ranges for u2 (the final vertex), precomputed once
    uint32_t q2s = bfs_q_offsets[u2], q2e = bfs_q_offsets[u2 + 1];

    uint64_t local_count = 0;

    // Each lane sweeps a slice of u's candidate list (depth Q-2 expansion).
    for (uint32_t i = lane; i < nbr_count; i += 32) {
        uint32_t v = csr_neighbors[nbr_start + i];

        // validate v as the depth-(Q-2) mapping of u
        if (vlabels[v] != bfs_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;
        bool joinable = true;
        for (uint32_t j = qs; j < qe; j++) {
            uint32_t uo = bfs_q_neighbors[j];
            if (m[uo] == UINT32_MAX || uo == u_min) continue;
            uint32_t s = csr_offsets[m[uo]];
            uint32_t nc = (bfs_use_padded && degrees) ? degrees[m[uo]] : (csr_offsets[m[uo] + 1] - s);
            int pos = bfs_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) { joinable = false; break; }
        }
        if (!joinable) continue;
        // visited check for v against current m[]
        bool vis = false;
        for (uint32_t d = 0; d < Q; d++) { if (m[d] == v) { vis = true; break; } }
        if (vis) continue;

        // --- v is a valid depth-(Q-2) mapping; now extend to depth Q-1 (final) ---
        // Tentatively set m[u] = v, then find u2's pivot among matched neighbours
        // (which now includes u). We don't write m; we evaluate u2 inline.
        uint32_t u2_min = UINT32_MAX, u2_min_label = 0, u2_min_deg = UINT32_MAX;
        for (uint32_t j = q2s; j < q2e; j++) {
            uint32_t uo = bfs_q_neighbors[j];
            uint32_t muo = (uo == u) ? v : m[uo];   // u just got mapped to v
            if (muo == UINT32_MAX) continue;
            uint32_t off = csr_offsets[muo];
            uint32_t deg = (bfs_use_padded && degrees) ? degrees[muo] : (csr_offsets[muo + 1] - off);
            if (deg < u2_min_deg) { u2_min_deg = deg; u2_min = uo; u2_min_label = bfs_q_elabels[j]; }
        }
        if (u2_min == UINT32_MAX) continue;

        uint32_t m_u2min = (u2_min == u) ? v : m[u2_min];
        uint32_t ns2 = csr_offsets[m_u2min];
        uint32_t nc2 = (bfs_use_padded && degrees) ? degrees[m_u2min] : (csr_offsets[m_u2min + 1] - ns2);

        for (uint32_t k = 0; k < nc2; k++) {
            uint32_t w = csr_neighbors[ns2 + k];
            if (vlabels[w] != bfs_q_vlabels[u2]) continue;
            if (csr_elabels[ns2 + k] != u2_min_label) continue;
            // joinability of w against all matched neighbours of u2 (incl. u→v)
            bool jn = true;
            for (uint32_t j = q2s; j < q2e; j++) {
                uint32_t uo = bfs_q_neighbors[j];
                if (uo == u2_min) continue;
                uint32_t muo = (uo == u) ? v : m[uo];
                if (muo == UINT32_MAX) continue;
                uint32_t s = csr_offsets[muo];
                uint32_t nc = (bfs_use_padded && degrees) ? degrees[muo] : (csr_offsets[muo + 1] - s);
                int pos = bfs_bsearch(csr_neighbors + s, nc, w);
                if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) { jn = false; break; }
            }
            if (!jn) continue;
            // visited: w must differ from all of m[] and from v
            bool vis2 = (w == v);
            if (!vis2) for (uint32_t d = 0; d < Q; d++) { if (m[d] == w) { vis2 = true; break; } }
            if (vis2) continue;

            local_count++;
        }
    }

    local_count = warp_reduce_add_u64(local_count);
    if (lane == 0 && local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// LJ variant (sparse graphs): inner-parallel fused count.
// Profile on LiveJournal showed bfs_expand_count_kernel hit only 1.64/32 active
// threads: it splits the OUTER candidate list (order[Q-2]) across lanes, but on
// sparse/heavily-filtered graphs that list is short, so most lanes idle while a
// few do the deep inner sweep. This variant flips it: the warp iterates the
// outer candidate v *together* (uniform), and splits the INNER candidate list
// (order[Q-1]) across the 32 lanes. When the inner list is long (hubs on LJ),
// all lanes stay busy. Identical match semantics → same count.
// Selected at runtime for low-unsafe-edge / sparse workloads.
// ============================================================
__global__ void bfs_expand_count_lj_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,  // NULL for compact CSR
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint64_t* __restrict__ result_count,
    uint32_t depth,   // = Q-2
    uint32_t Q
) {
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 1;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    const uint32_t* m = pm + 1;

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u  = order[depth];
    uint32_t u2 = order[depth + 1];

    // pivot for u (warp-uniform: every lane computes the same)
    uint32_t u_min = UINT32_MAX, u_min_label = 0, u_min_deg = UINT32_MAX;
    uint32_t qs = bfs_q_offsets[u], qe = bfs_q_offsets[u + 1];
    for (uint32_t j = qs; j < qe; j++) {
        uint32_t uo = bfs_q_neighbors[j];
        if (m[uo] == UINT32_MAX) continue;
        uint32_t off = csr_offsets[m[uo]];
        uint32_t deg = (bfs_use_padded && degrees) ? degrees[m[uo]] : (csr_offsets[m[uo] + 1] - off);
        if (deg < u_min_deg) { u_min_deg = deg; u_min = uo; u_min_label = bfs_q_elabels[j]; }
    }
    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = (bfs_use_padded && degrees) ? degrees[m[u_min]] : (csr_offsets[m[u_min] + 1] - nbr_start);
    uint32_t q2s = bfs_q_offsets[u2], q2e = bfs_q_offsets[u2 + 1];

    uint64_t local_count = 0;

    // OUTER: whole warp walks the same v sequentially (warp-uniform).
    for (uint32_t i = 0; i < nbr_count; i++) {
        uint32_t v = csr_neighbors[nbr_start + i];
        if (vlabels[v] != bfs_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;
        bool joinable = true;
        for (uint32_t j = qs; j < qe; j++) {
            uint32_t uo = bfs_q_neighbors[j];
            if (m[uo] == UINT32_MAX || uo == u_min) continue;
            uint32_t s = csr_offsets[m[uo]];
            uint32_t nc = (bfs_use_padded && degrees) ? degrees[m[uo]] : (csr_offsets[m[uo] + 1] - s);
            int pos = bfs_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) { joinable = false; break; }
        }
        if (!joinable) continue;
        bool vis = false;
        for (uint32_t d = 0; d < Q; d++) { if (m[d] == v) { vis = true; break; } }
        if (vis) continue;

        // u2 pivot given m[u]=v (warp-uniform)
        uint32_t u2_min = UINT32_MAX, u2_min_label = 0, u2_min_deg = UINT32_MAX;
        for (uint32_t j = q2s; j < q2e; j++) {
            uint32_t uo = bfs_q_neighbors[j];
            uint32_t muo = (uo == u) ? v : m[uo];
            if (muo == UINT32_MAX) continue;
            uint32_t off = csr_offsets[muo];
            uint32_t deg = (bfs_use_padded && degrees) ? degrees[muo] : (csr_offsets[muo + 1] - off);
            if (deg < u2_min_deg) { u2_min_deg = deg; u2_min = uo; u2_min_label = bfs_q_elabels[j]; }
        }
        if (u2_min == UINT32_MAX) continue;

        uint32_t m_u2min = (u2_min == u) ? v : m[u2_min];
        uint32_t ns2 = csr_offsets[m_u2min];
        uint32_t nc2 = (bfs_use_padded && degrees) ? degrees[m_u2min] : (csr_offsets[m_u2min + 1] - ns2);

        // INNER: 32 lanes split the w candidate list (long on hubs → lanes busy).
        for (uint32_t k = lane; k < nc2; k += 32) {
            uint32_t w = csr_neighbors[ns2 + k];
            if (vlabels[w] != bfs_q_vlabels[u2]) continue;
            if (csr_elabels[ns2 + k] != u2_min_label) continue;
            bool jn = true;
            for (uint32_t j = q2s; j < q2e; j++) {
                uint32_t uo = bfs_q_neighbors[j];
                if (uo == u2_min) continue;
                uint32_t muo = (uo == u) ? v : m[uo];
                if (muo == UINT32_MAX) continue;
                uint32_t s = csr_offsets[muo];
                uint32_t nc = (bfs_use_padded && degrees) ? degrees[muo] : (csr_offsets[muo + 1] - s);
                int pos = bfs_bsearch(csr_neighbors + s, nc, w);
                if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j]) { jn = false; break; }
            }
            if (!jn) continue;
            bool vis2 = (w == v);
            if (!vis2) for (uint32_t d = 0; d < Q; d++) { if (m[d] == w) { vis2 = true; break; } }
            if (vis2) continue;
            local_count++;
        }
    }

    local_count = warp_reduce_add_u64(local_count);
    if (lane == 0 && local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// Versioned Kernels: each partial match carries max_ts
// Partial match layout: [order_idx, max_ts, m[0], ..., m[Q-1]]
//                        stride = Q + 2
// ============================================================

__global__ void bfs_init_versioned_kernel(
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ edges_v1,
    const uint32_t* __restrict__ edges_v2,
    const uint32_t* __restrict__ edges_label,
    const uint32_t* __restrict__ edges_max_ts,
    const uint32_t* __restrict__ all_orders,
    uint32_t* __restrict__ out_buf,
    unsigned long long* __restrict__ out_count,
    uint32_t num_data_edges,
    uint32_t num_q_edges,
    uint32_t Q,
    uint32_t max_out
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t total = num_data_edges * num_q_edges * 2;
    if (idx >= total) return;

    uint32_t edge_idx = idx / (num_q_edges * 2);
    uint32_t remainder = idx % (num_q_edges * 2);
    uint32_t qe_idx = remainder / 2;
    uint32_t dir = remainder % 2;

    uint32_t v1 = edges_v1[edge_idx];
    uint32_t v2 = edges_v2[edge_idx];
    uint32_t elabel = edges_label[edge_idx];
    uint32_t max_ts = edges_max_ts[edge_idx];

    const uint32_t* order = all_orders + qe_idx * Q;
    uint32_t u1 = order[0];
    uint32_t u2 = order[1];

    uint32_t qs = bfs_q_offsets[u1];
    uint32_t qe = bfs_q_offsets[u1 + 1];
    uint32_t found_el = UINT32_MAX;
    for (uint32_t j = qs; j < qe; j++) {
        if (bfs_q_neighbors[j] == u2) { found_el = bfs_q_elabels[j]; break; }
    }
    if (found_el == UINT32_MAX) return;

    uint32_t mv1 = dir == 0 ? v1 : v2;
    uint32_t mv2 = dir == 0 ? v2 : v1;

    if (vlabels[mv1] != bfs_q_vlabels[u1]) return;
    if (vlabels[mv2] != bfs_q_vlabels[u2]) return;
    if (elabel != found_el) return;

    uint64_t slot = warp_agg_inc(out_count);  // OPT(B): warp-aggregated atomic
    if (slot >= max_out) return;  // count but skip write

    uint32_t stride = Q + 2;
    uint32_t* out = out_buf + slot * stride;
    out[0] = qe_idx;
    out[1] = max_ts;       // store per-edge timestamp cap
    for (uint32_t j = 0; j < Q; j++) out[2 + j] = UINT32_MAX;
    out[2 + u1] = mv1;
    out[2 + u2] = mv2;
}

__global__ void bfs_expand_versioned_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ csr_timestamps,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint32_t* __restrict__ out_buf,
    unsigned long long* __restrict__ out_count,
    uint32_t depth,
    uint32_t Q,
    uint32_t max_out
) {
    // OPT(C): warp-per-partial-match (see non-versioned expand kernel).
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;  // 64-bit: avoid wrap >2^27 pm
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 2;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    uint32_t max_ts = pm[1];
    const uint32_t* m = pm + 2;

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u = order[depth];

    // Find u_min (warp-uniform)
    uint32_t u_min = UINT32_MAX, u_min_label = 0, u_min_deg = UINT32_MAX;
    uint32_t qs = bfs_q_offsets[u];
    uint32_t qe = bfs_q_offsets[u + 1];

    for (uint32_t j = qs; j < qe; j++) {
        uint32_t u_other = bfs_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        // OPT(A): use O(1) physical degree as a proxy for visible degree.
        // u_min is only the search pivot (smallest candidate set); the per-edge
        // timestamp filter below still runs on every neighbour, so an approximate
        // pivot cannot change the match count — only the work done to find it.
        // This removes an O(sum of degrees) scan per partial match (the hot spot
        // that made versioned ~80% slower than plain BFS on LiveJournal).
        uint32_t deg = csr_offsets[m[u_other] + 1] - csr_offsets[m[u_other]];
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = bfs_q_elabels[j];
        }
    }
    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = csr_offsets[m[u_min] + 1] - nbr_start;

    for (uint32_t base = 0; base < nbr_count; base += 32) {
        uint32_t i = base + lane;
        bool emit = false;
        uint32_t v = 0;

        if (i < nbr_count && csr_timestamps[nbr_start + i] <= max_ts) {  // version filter
            v = csr_neighbors[nbr_start + i];
            // Label check
            if (vlabels[v] == bfs_q_vlabels[u] &&
                csr_elabels[nbr_start + i] == u_min_label) {
                // Joinability (version-aware)
                bool joinable = true;
                for (uint32_t j = qs; j < qe; j++) {
                    uint32_t u_other = bfs_q_neighbors[j];
                    if (m[u_other] == UINT32_MAX || u_other == u_min) continue;
                    uint32_t s = csr_offsets[m[u_other]];
                    uint32_t nc = csr_offsets[m[u_other] + 1] - s;
                    int pos = bfs_bsearch(csr_neighbors + s, nc, v);
                    if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j] ||
                        csr_timestamps[s + pos] > max_ts) {
                        joinable = false; break;
                    }
                }
                if (joinable) {
                    bool visited = false;
                    for (uint32_t d = 0; d < Q; d++) {
                        if (m[d] == v) { visited = true; break; }
                    }
                    emit = !visited;
                }
            }
        }

        // Warp-cooperative write
        uint32_t emit_mask = __ballot_sync(0xffffffffu, emit);
        uint32_t cnt = __popc(emit_mask);
        if (cnt) {
            unsigned long long base_slot = 0;
            uint32_t leader = __ffs(emit_mask) - 1;
            if (lane == leader) base_slot = atomicAdd(out_count, (unsigned long long)cnt);
            base_slot = __shfl_sync(0xffffffffu, base_slot, leader);
            if (emit) {
                uint32_t rank = __popc(emit_mask & ((1u << lane) - 1u));
                uint64_t slot = base_slot + rank;
                if (slot < max_out) {
                    uint32_t* out = out_buf + (size_t)slot * stride;
                    out[0] = order_idx;
                    out[1] = max_ts;
                    for (uint32_t j = 0; j < Q; j++) out[2 + j] = m[j];
                    out[2 + u] = v;
                }
            }
        }
    }
}

__global__ void bfs_count_versioned_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ csr_timestamps,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ all_orders,
    const uint32_t* __restrict__ in_buf,
    uint32_t in_count,
    uint64_t* __restrict__ result_count,
    uint32_t depth,
    uint32_t Q
) {
    // OPT(C): warp-per-partial-match (see non-versioned count kernel).
    uint64_t gtid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;  // 64-bit: avoid wrap >2^27 pm
    uint32_t warp_id = (uint32_t)(gtid >> 5);
    uint32_t lane = (uint32_t)(gtid & 31u);
    if (warp_id >= in_count) return;

    uint32_t stride = Q + 2;
    const uint32_t* pm = in_buf + (size_t)warp_id * stride;
    uint32_t order_idx = pm[0];
    uint32_t max_ts = pm[1];
    const uint32_t* m = pm + 2;

    const uint32_t* order = all_orders + order_idx * Q;
    uint32_t u = order[depth];

    uint32_t u_min = UINT32_MAX, u_min_label = 0, u_min_deg = UINT32_MAX;
    uint32_t qs = bfs_q_offsets[u];
    uint32_t qe = bfs_q_offsets[u + 1];

    for (uint32_t j = qs; j < qe; j++) {
        uint32_t u_other = bfs_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        // OPT(A): O(1) physical degree proxy for visible degree (see expand kernel).
        // Pivot choice does not affect correctness; the timestamp filter below is
        // still applied per neighbour.
        uint32_t deg = csr_offsets[m[u_other] + 1] - csr_offsets[m[u_other]];
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = bfs_q_elabels[j];
        }
    }
    if (u_min == UINT32_MAX) return;

    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = csr_offsets[m[u_min] + 1] - nbr_start;
    uint64_t local_count = 0;

    for (uint32_t i = lane; i < nbr_count; i += 32) {
        if (csr_timestamps[nbr_start + i] > max_ts) continue;

        uint32_t v = csr_neighbors[nbr_start + i];
        if (vlabels[v] != bfs_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;

        bool joinable = true;
        for (uint32_t j = qs; j < qe; j++) {
            uint32_t u_other = bfs_q_neighbors[j];
            if (m[u_other] == UINT32_MAX || u_other == u_min) continue;
            uint32_t s = csr_offsets[m[u_other]];
            uint32_t nc = csr_offsets[m[u_other] + 1] - s;
            int pos = bfs_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != bfs_q_elabels[j] ||
                csr_timestamps[s + pos] > max_ts) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;

        bool visited = false;
        for (uint32_t d = 0; d < Q; d++) {
            if (m[d] == v) { visited = true; break; }
        }
        if (visited) continue;

        local_count++;
    }

    local_count = warp_reduce_add_u64(local_count);
    if (lane == 0 && local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(result_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// Host implementation
// ============================================================

GPUBFSSearch::GPUBFSSearch() = default;
GPUBFSSearch::~GPUBFSSearch() { Destroy(); }

void GPUBFSSearch::BuildCSR(const Graph& data) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (d_csr_offsets_) CUDA_CHECK(cudaFree(d_csr_offsets_));
    if (d_csr_neighbors_) CUDA_CHECK(cudaFree(d_csr_neighbors_));
    if (d_csr_elabels_) CUDA_CHECK(cudaFree(d_csr_elabels_));
    if (d_vlabels_) CUDA_CHECK(cudaFree(d_vlabels_));

    uint32_t V = data.NumVertices();
    num_vertices_ = V;

    // Select the fused-count kernel variant. Default heuristic: large, sparse
    // graphs (many vertices, low avg degree) behave like LiveJournal — short
    // outer candidate lists starve the outer-parallel kernel, so use the
    // inner-parallel LJ kernel. Env GPU_BFS_SPARSE=0/1 forces it for A/B tests.
    {
        const char* e = getenv("GPU_BFS_SPARSE");
        if (e) {
            sparse_mode_ = (e[0] == '1');
        } else {
            double avg_deg = V ? (double)data.NumEdges() * 2.0 / V : 0.0;
            sparse_mode_ = (V >= 1000000 && avg_deg < 32.0);  // LJ-like
        }
    }

    std::vector<uint32_t> offsets(V + 1);
    offsets[0] = 0;
    for (uint32_t v = 0; v < V; v++)
        offsets[v + 1] = offsets[v] + data.GetDegree(v);
    uint32_t total = offsets[V];

    std::vector<uint32_t> neighbors(total), elabels(total);
    for (uint32_t v = 0; v < V; v++) {
        const auto& nbrs = data.GetNeighbors(v);
        const auto& labs = data.GetNeighborLabels(v);
        uint32_t off = offsets[v];
        for (size_t j = 0; j < nbrs.size(); j++) {
            neighbors[off + j] = nbrs[j];
            elabels[off + j] = labs[j];
        }
    }

    CUDA_CHECK(cudaMalloc(&d_csr_offsets_, (V + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_neighbors_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_elabels_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vlabels_, V * sizeof(uint32_t)));

    CUDA_CHECK(cudaMemcpy(d_csr_offsets_, offsets.data(), (V+1)*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_, neighbors.data(), total*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_, elabels.data(), total*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, data.vlabels_.data(), V*sizeof(uint32_t), cudaMemcpyHostToDevice));

    csr_built_ = true;
    padded_csr_ = false;
    uint32_t use_padded = 0;
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_use_padded, &use_padded, sizeof(uint32_t)));
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    printf("[BFS] CSR built: V=%u, E=%u, %.1fms\n", V, total, ms);
}

void GPUBFSSearch::BuildCSR_Versioned(const Graph& data) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (d_csr_offsets_) CUDA_CHECK(cudaFree(d_csr_offsets_));
    if (d_csr_neighbors_) CUDA_CHECK(cudaFree(d_csr_neighbors_));
    if (d_csr_elabels_) CUDA_CHECK(cudaFree(d_csr_elabels_));
    if (d_csr_timestamps_) CUDA_CHECK(cudaFree(d_csr_timestamps_));
    if (d_vlabels_) CUDA_CHECK(cudaFree(d_vlabels_));

    uint32_t V = data.NumVertices();
    num_vertices_ = V;

    std::vector<uint32_t> offsets(V + 1);
    offsets[0] = 0;
    for (uint32_t v = 0; v < V; v++)
        offsets[v + 1] = offsets[v] + data.GetDegree(v);
    uint32_t total = offsets[V];

    std::vector<uint32_t> neighbors(total), elabels(total), timestamps(total);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (uint32_t v = 0; v < V; v++) {
        const auto& nbrs = data.GetNeighbors(v);
        const auto& labs = data.GetNeighborLabels(v);
        const auto& ts = data.GetEdgeTimestamps(v);
        uint32_t off = offsets[v];
        for (size_t j = 0; j < nbrs.size(); j++) {
            neighbors[off + j] = nbrs[j];
            elabels[off + j] = labs[j];
            timestamps[off + j] = ts[j];
        }
    }

    CUDA_CHECK(cudaMalloc(&d_csr_offsets_, (V + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_neighbors_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_elabels_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_timestamps_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vlabels_, V * sizeof(uint32_t)));

    CUDA_CHECK(cudaMemcpy(d_csr_offsets_, offsets.data(), (V+1)*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_, neighbors.data(), total*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_, elabels.data(), total*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_timestamps_, timestamps.data(), total*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, data.vlabels_.data(), V*sizeof(uint32_t), cudaMemcpyHostToDevice));

    csr_built_ = true;
    padded_csr_ = false;
    uint32_t use_padded = 0;
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_use_padded, &use_padded, sizeof(uint32_t)));

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    printf("[BFS-Versioned] CSR built: V=%u, E=%u, %.1fms\n", V, total, ms);
}

void GPUBFSSearch::BuildPaddedCSR(const Graph& data) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (d_csr_offsets_) CUDA_CHECK(cudaFree(d_csr_offsets_));
    if (d_csr_neighbors_) CUDA_CHECK(cudaFree(d_csr_neighbors_));
    if (d_csr_elabels_) CUDA_CHECK(cudaFree(d_csr_elabels_));
    if (d_vlabels_) CUDA_CHECK(cudaFree(d_vlabels_));
    if (d_degrees_) CUDA_CHECK(cudaFree(d_degrees_));

    uint32_t V = data.NumVertices();
    num_vertices_ = V;

    // Build padded CSR: each vertex gets capacity = max(degree * 2, degree + 8)
    h_offsets_.resize(V + 1);
    h_capacities_.resize(V);
    std::vector<uint32_t> h_degrees(V);

    h_offsets_[0] = 0;
    for (uint32_t v = 0; v < V; v++) {
        uint32_t deg = data.GetDegree(v);
        uint32_t cap = std::max(deg * 2, deg + 8u);
        h_capacities_[v] = cap;
        h_degrees[v] = deg;
        h_offsets_[v + 1] = h_offsets_[v] + cap;
    }
    uint32_t total_padded = h_offsets_[V];

    // Build flattened neighbors and elabels (with padding zeros)
    std::vector<uint32_t> neighbors(total_padded, 0);
    std::vector<uint32_t> elabels(total_padded, 0);
    for (uint32_t v = 0; v < V; v++) {
        const auto& nbrs = data.GetNeighbors(v);
        const auto& labs = data.GetNeighborLabels(v);
        uint32_t off = h_offsets_[v];
        for (size_t j = 0; j < nbrs.size(); j++) {
            neighbors[off + j] = nbrs[j];
            elabels[off + j] = labs[j];
        }
    }

    CUDA_CHECK(cudaMalloc(&d_csr_offsets_, (V + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_neighbors_, total_padded * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_elabels_, total_padded * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vlabels_, V * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_degrees_, V * sizeof(uint32_t)));

    CUDA_CHECK(cudaMemcpy(d_csr_offsets_, h_offsets_.data(), (V+1)*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_, neighbors.data(), total_padded*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_, elabels.data(), total_padded*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, data.vlabels_.data(), V*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_, h_degrees.data(), V*sizeof(uint32_t), cudaMemcpyHostToDevice));

    csr_built_ = true;
    padded_csr_ = true;
    uint32_t use_padded = 1;
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_use_padded, &use_padded, sizeof(uint32_t)));

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    printf("[BFS] Padded CSR built: V=%u, padded_edges=%u, %.1fms\n", V, total_padded, ms);
}

bool GPUBFSSearch::IncrementalUpdate(const Graph& data, uint32_t v1, uint32_t v2) {
    if (!csr_built_ || !padded_csr_) return false;

    uint32_t new_deg1 = data.GetDegree(v1);
    uint32_t new_deg2 = data.GetDegree(v2);

    if (v1 >= num_vertices_ || v2 >= num_vertices_ ||
        new_deg1 > h_capacities_[v1] || new_deg2 > h_capacities_[v2]) {
        return false;  // need full rebuild
    }

    // Update v1
    const auto& nbrs1 = data.GetNeighbors(v1);
    const auto& labs1 = data.GetNeighborLabels(v1);
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_ + h_offsets_[v1], nbrs1.data(),
                           new_deg1 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_ + h_offsets_[v1], labs1.data(),
                           new_deg1 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_ + v1, &new_deg1, sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Update v2
    const auto& nbrs2 = data.GetNeighbors(v2);
    const auto& labs2 = data.GetNeighborLabels(v2);
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_ + h_offsets_[v2], nbrs2.data(),
                           new_deg2 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_ + h_offsets_[v2], labs2.data(),
                           new_deg2 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_ + v2, &new_deg2, sizeof(uint32_t), cudaMemcpyHostToDevice));

    return true;
}

void GPUBFSSearch::IncrementalUpdateVertex(uint32_t v, const uint32_t* neighbors, const uint32_t* elabels, uint32_t degree) {
    if (!csr_built_ || !padded_csr_ || v >= num_vertices_) return;
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_ + h_offsets_[v], neighbors,
                           degree * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_ + h_offsets_[v], elabels,
                           degree * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_ + v, &degree, sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void GPUBFSSearch::SetupQuery(const Graph& query) {
    uint32_t Q = query.NumVertices();
    num_query_vertices_ = Q;

    std::vector<uint32_t> vl(Q);
    for (uint32_t u = 0; u < Q; u++) vl[u] = query.GetVertexLabel(u);

    std::vector<uint32_t> qoff(Q + 1);
    qoff[0] = 0;
    for (uint32_t u = 0; u < Q; u++) qoff[u+1] = qoff[u] + query.GetDegree(u);
    uint32_t qt = qoff[Q];

    // Guard against silent constant-memory overflow. qt = sum of degrees =
    // 2 * (#edges). If it exceeds the constant array size, cudaMemcpyToSymbol
    // would clobber adjacent symbols and corrupt joinability pruning.
    if (Q > (uint32_t)BFS_MAX_Q) {
        throw std::runtime_error("GPUBFSSearch::SetupQuery: query has " +
            std::to_string(Q) + " vertices > BFS_MAX_Q=" +
            std::to_string(BFS_MAX_Q));
    }
    if (qt > (uint32_t)BFS_MAX_Q_EDGES) {
        throw std::runtime_error("GPUBFSSearch::SetupQuery: query neighbor "
            "entries (2*edges)=" + std::to_string(qt) +
            " > BFS_MAX_Q_EDGES=" + std::to_string(BFS_MAX_Q_EDGES) +
            "; raise BFS_MAX_Q_EDGES");
    }

    std::vector<uint32_t> qn(qt), qe(qt);
    for (uint32_t u = 0; u < Q; u++) {
        const auto& nbrs = query.GetNeighbors(u);
        const auto& labs = query.GetNeighborLabels(u);
        uint32_t off = qoff[u];
        for (size_t j = 0; j < nbrs.size(); j++) {
            qn[off+j] = nbrs[j]; qe[off+j] = labs[j];
        }
    }

    CUDA_CHECK(cudaMemcpyToSymbol(bfs_q_num_vertices, &Q, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_q_vlabels, vl.data(), Q*sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_q_offsets, qoff.data(), (Q+1)*sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_q_neighbors, qn.data(), qt*sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(bfs_q_elabels, qe.data(), qt*sizeof(uint32_t)));

    query_set_ = true;
}

void GPUBFSSearch::SetAllMatchingOrders(const uint32_t* all_orders, uint32_t num_edges, uint32_t num_vertices) {
    if (d_all_orders_) CUDA_CHECK(cudaFree(d_all_orders_));
    all_orders_num_edges_ = num_edges;
    size_t total = (size_t)num_edges * num_vertices;
    CUDA_CHECK(cudaMalloc(&d_all_orders_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_all_orders_, all_orders, total*sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void GPUBFSSearch::EnsureEdgesCapacity(size_t n) {
    if (n <= edges_capacity_) return;
    if (d_edges_v1_) CUDA_CHECK(cudaFree(d_edges_v1_));
    if (d_edges_v2_) CUDA_CHECK(cudaFree(d_edges_v2_));
    if (d_edges_label_) CUDA_CHECK(cudaFree(d_edges_label_));
    if (d_edges_max_ts_) CUDA_CHECK(cudaFree(d_edges_max_ts_));
    edges_capacity_ = n + n/4;
    CUDA_CHECK(cudaMalloc(&d_edges_v1_, edges_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_edges_v2_, edges_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_edges_label_, edges_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_edges_max_ts_, edges_capacity_ * sizeof(uint32_t)));
}

void GPUBFSSearch::EnsureBufCapacity(uint32_t q) {
    if (d_buf_a_) return;  // already allocated
    size_t stride = q + 1;
    size_t buf_bytes = MAX_BUF_MATCHES * stride * sizeof(uint32_t);
    CUDA_CHECK(cudaMalloc(&d_buf_a_, buf_bytes));
    CUDA_CHECK(cudaMalloc(&d_buf_b_, buf_bytes));
    CUDA_CHECK(cudaMalloc(&d_buf_c_, buf_bytes));  // third buffer for overflow flush
    CUDA_CHECK(cudaMalloc(&d_count_, sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_result_, sizeof(uint64_t)));
    printf("[BFS] Buffers allocated: 3 x %zuMB (max %zuM partial matches)\n",
           buf_bytes / (1024*1024), MAX_BUF_MATCHES / 1000000);
}

// ---------------------------------------------------------------------------
// BFSFromDepth: process partial matches from start_depth through Q-1,
// handling overflow by flushing accumulated results through remaining
// levels before continuing. Results are accumulated into d_result_.
// ---------------------------------------------------------------------------
void GPUBFSSearch::BFSFromDepth(uint32_t* in_buf, uint32_t in_count,
                                 uint32_t* scratch_buf, uint32_t start_depth,
                                 uint32_t Q, uint32_t stride, bool versioned)
{
    if (in_count == 0) return;

    int block = 256;
    int grid;
    uint64_t zero64 = 0;

    uint32_t* cur_buf = in_buf;
    uint32_t* next_buf = scratch_buf;
    uint32_t cur_count = in_count;

    for (uint32_t depth = start_depth; depth < Q; depth++) {
        bool is_last = (depth == Q - 1);

        // OPT(P1b): fuse the last two depths. At depth Q-2, instead of
        // materialising the (largest) depth-(Q-1) layer and then counting it,
        // extend each candidate in-register to depth Q-1 and count directly.
        // Only for the non-versioned compact/padded path; versioned keeps the
        // explicit two-kernel flow (timestamp logic).
        if (!versioned && depth == Q - 2 && Q >= 3) {
            grid = (uint32_t)(((uint64_t)cur_count * 32 + block - 1) / block);
            // Pick fusion kernel by graph density (set in BuildCSR). Sparse
            // graphs (LJ) have short outer candidate lists → the default
            // outer-parallel kernel starves lanes (1.64/32); the LJ variant
            // parallelises the inner list instead.
            if (sparse_mode_) {
                bfs_expand_count_lj_kernel<<<grid, block>>>(
                    d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                    d_degrees_, d_all_orders_, cur_buf, cur_count, d_result_, depth, Q
                );
            } else {
                bfs_expand_count_kernel<<<grid, block>>>(
                    d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                    d_degrees_, d_all_orders_, cur_buf, cur_count, d_result_, depth, Q
                );
            }
            CUDA_CHECK(cudaGetLastError());
            printf("[BFS] Fused depth %u->%u (%s count, no materialise): %u input pm\n",
                   depth, depth + 1, sparse_mode_ ? "LJ inner-par" : "outer-par", cur_count);
            break;  // final level handled; done
        }

        if (is_last) {
            // OPT(C): one warp per partial match → 32 threads each.
            grid = (uint32_t)(((uint64_t)cur_count * 32 + block - 1) / block);
            if (versioned) {
                bfs_count_versioned_kernel<<<grid, block>>>(
                    d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_csr_timestamps_,
                    d_vlabels_, d_all_orders_, cur_buf, cur_count, d_result_, depth, Q
                );
            } else {
                bfs_count_kernel<<<grid, block>>>(
                    d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                    d_degrees_, d_all_orders_, cur_buf, cur_count, d_result_, depth, Q
                );
            }
            CUDA_CHECK(cudaGetLastError());
        } else {
            uint32_t in_processed = 0;
            uint32_t total_out = 0;

            while (in_processed < cur_count) {
                uint32_t chunk = cur_count - in_processed;
                uint32_t max_chunk = static_cast<uint32_t>(MAX_BUF_MATCHES / 2);
                if (chunk > max_chunk) chunk = max_chunk;

                bool overflow;
                uint64_t out_count;  // TRUE logical emitter count (may exceed 2^32 for dense layers)
                do {
                    CUDA_CHECK(cudaMemcpy(d_count_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));
                    grid = (uint32_t)(((uint64_t)chunk * 32 + block - 1) / block);  // OPT(C): warp per pm
                    if (versioned) {
                        bfs_expand_versioned_kernel<<<grid, block>>>(
                            d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_csr_timestamps_,
                            d_vlabels_, d_all_orders_,
                            cur_buf + (size_t)in_processed * stride, chunk,
                            next_buf + (size_t)total_out * stride, d_count_,
                            depth, Q, static_cast<uint32_t>(MAX_BUF_MATCHES - total_out)
                        );
                    } else {
                        bfs_expand_kernel<<<grid, block>>>(
                            d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                            d_degrees_, d_all_orders_,
                            cur_buf + (size_t)in_processed * stride, chunk,
                            next_buf + (size_t)total_out * stride, d_count_,
                            depth, Q, static_cast<uint32_t>(MAX_BUF_MATCHES - total_out)
                        );
                    }
                    CUDA_CHECK(cudaGetLastError());
                    CUDA_CHECK(cudaDeviceSynchronize());

                    CUDA_CHECK(cudaMemcpy(&out_count, d_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
                    overflow = ((uint64_t)total_out + out_count > MAX_BUF_MATCHES);
                    if (overflow && chunk > 1) {
                        chunk = chunk / 2;
                    } else {
                        break;
                    }
                } while (true);

                if (overflow) {
                    if (chunk == 1 && total_out == 0) {
                        // Single partial match overflows — its emitters were written up to
                        // MAX_BUF_MATCHES (rest dropped by the kernel's slot<max_out guard).
                        uint32_t written = static_cast<uint32_t>(std::min(out_count, (uint64_t)MAX_BUF_MATCHES));
                        if (written > 0) {
                            printf("[BFS%s] Single-item overflow at depth %u->%u (%u/%llu outputs), flush\n",
                                   versioned ? "-V" : "", depth, depth+1, written,
                                   (unsigned long long)out_count);
                            BFSFromDepth(next_buf, written, d_buf_c_, depth + 1, Q, stride, versioned);
                        }
                        total_out = 0;
                    } else {
                        // Multi-item overflow — flush the buffer-full prefix and continue.
                        total_out = MAX_BUF_MATCHES;
                        printf("[BFS%s] Flush %u partials at depth %u->%u (overflow, chunk=%u)\n",
                               versioned ? "-V" : "", total_out, depth, depth+1, chunk);
                        BFSFromDepth(next_buf, total_out, d_buf_c_, depth + 1, Q, stride, versioned);
                        total_out = 0;
                    }
                } else {
                    total_out += static_cast<uint32_t>(out_count);
                }
                in_processed += chunk;
            }

            printf("[BFS%s] Depth %u→%u: %u → %u partial matches\n",
                   versioned ? "-V" : "", depth, depth+1, cur_count, total_out);

            if (total_out == 0) break;
            std::swap(cur_buf, next_buf);
            cur_count = total_out;
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
}

uint64_t GPUBFSSearch::SearchBatchEdgesBFS(
    const uint32_t* edges_v1,
    const uint32_t* edges_v2,
    const uint32_t* edges_label,
    uint32_t num_edges_data,
    uint32_t num_query_edges,
    uint32_t num_query_vertices
) {
    if (!csr_built_ || !query_set_ || num_edges_data == 0) return 0;

    uint32_t Q = num_query_vertices;
    uint32_t stride = Q + 1;

    EnsureEdgesCapacity(num_edges_data);
    EnsureBufCapacity(Q);

    // Copy edges
    CUDA_CHECK(cudaMemcpy(d_edges_v1_, edges_v1, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_v2_, edges_v2, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_label_, edges_label, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Reset result
    uint64_t zero64 = 0;
    CUDA_CHECK(cudaMemcpy(d_result_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    auto t_start = std::chrono::high_resolution_clock::now();

    // ---- Phase 1: Init — create depth=2 partial matches ----
    CUDA_CHECK(cudaMemcpy(d_count_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    uint32_t total_init_tasks = num_edges_data * num_query_edges * 2;
    int block = 256;
    int grid = (total_init_tasks + block - 1) / block;

    bfs_init_kernel<<<grid, block>>>(
        d_vlabels_, d_edges_v1_, d_edges_v2_, d_edges_label_,
        d_all_orders_, d_buf_a_, d_count_,
        num_edges_data, num_query_edges, Q,
        static_cast<uint32_t>(MAX_BUF_MATCHES)
    );
    CUDA_CHECK(cudaGetLastError());

    uint64_t cur_count64 = 0;
    CUDA_CHECK(cudaMemcpy(&cur_count64, d_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
    uint32_t cur_count = static_cast<uint32_t>(std::min(cur_count64, (uint64_t)MAX_BUF_MATCHES));

    auto t_init = std::chrono::high_resolution_clock::now();
    double init_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_init - t_start).count() / 1000.0;
    printf("[BFS] Init: %u tasks → %u partial matches (%.1fms)\n",
           total_init_tasks, cur_count, init_ms);

    if (cur_count == 0) return 0;

    // ---- Phase 2: BFS level-by-level with overflow handling ----
    BFSFromDepth(d_buf_a_, cur_count, d_buf_b_, 2, Q, stride, false);

    CUDA_CHECK(cudaDeviceSynchronize());

    uint64_t total_count;
    CUDA_CHECK(cudaMemcpy(&total_count, d_result_, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() / 1000.0;
    printf("[BFS] Search complete: %lu matches, %.1fms\n", total_count, total_ms);

    return total_count;
}

// ============================================================
// Versioned BFS search: inner-update semantics with batch parallelism
// ============================================================
uint64_t GPUBFSSearch::SearchBatchEdgesBFS_Versioned(
    const uint32_t* edges_v1,
    const uint32_t* edges_v2,
    const uint32_t* edges_label,
    const uint32_t* edges_max_ts,
    uint32_t num_edges_data,
    uint32_t num_query_edges,
    uint32_t num_query_vertices
) {
    if (!csr_built_ || !query_set_ || num_edges_data == 0 || !d_csr_timestamps_) return 0;

    uint32_t Q = num_query_vertices;
    uint32_t stride = Q + 2;  // order_idx + max_ts + m[0..Q-1]

    EnsureEdgesCapacity(num_edges_data);

    // Allocate versioned BFS buffers (stride = Q+2, not Q+1)
    size_t buf_bytes = MAX_BUF_MATCHES * stride * sizeof(uint32_t);
    if (!d_buf_a_) {
        CUDA_CHECK(cudaMalloc(&d_buf_a_, buf_bytes));
        CUDA_CHECK(cudaMalloc(&d_buf_b_, buf_bytes));
        CUDA_CHECK(cudaMalloc(&d_buf_c_, buf_bytes));
        CUDA_CHECK(cudaMalloc(&d_count_, sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_result_, sizeof(uint64_t)));
    }

    // Copy edges + timestamps
    CUDA_CHECK(cudaMemcpy(d_edges_v1_, edges_v1, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_v2_, edges_v2, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_label_, edges_label, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));
    if (!d_edges_max_ts_) {
        CUDA_CHECK(cudaMalloc(&d_edges_max_ts_, edges_capacity_ * sizeof(uint32_t)));
    }
    CUDA_CHECK(cudaMemcpy(d_edges_max_ts_, edges_max_ts, num_edges_data*sizeof(uint32_t), cudaMemcpyHostToDevice));

    uint64_t zero64 = 0;
    CUDA_CHECK(cudaMemcpy(d_result_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    auto t_start = std::chrono::high_resolution_clock::now();

    // ---- Init ----
    CUDA_CHECK(cudaMemcpy(d_count_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    uint32_t total_init_tasks = num_edges_data * num_query_edges * 2;
    int block = 256;
    int grid = (total_init_tasks + block - 1) / block;

    bfs_init_versioned_kernel<<<grid, block>>>(
        d_vlabels_, d_edges_v1_, d_edges_v2_, d_edges_label_, d_edges_max_ts_,
        d_all_orders_, d_buf_a_, d_count_,
        num_edges_data, num_query_edges, Q,
        static_cast<uint32_t>(MAX_BUF_MATCHES)
    );
    CUDA_CHECK(cudaGetLastError());

    uint64_t cur_count64 = 0;
    CUDA_CHECK(cudaMemcpy(&cur_count64, d_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
    uint32_t cur_count = static_cast<uint32_t>(std::min(cur_count64, (uint64_t)MAX_BUF_MATCHES));

    auto t_init = std::chrono::high_resolution_clock::now();
    printf("[BFS-V] Init: %u tasks → %u partial matches (%.1fms)\n",
           total_init_tasks, cur_count,
           std::chrono::duration_cast<std::chrono::microseconds>(t_init - t_start).count() / 1000.0);

    if (cur_count == 0) return 0;

    // ---- BFS levels ----
    BFSFromDepth(d_buf_a_, cur_count, d_buf_b_, 2, Q, stride, true);

    CUDA_CHECK(cudaDeviceSynchronize());
    uint64_t total_count;
    CUDA_CHECK(cudaMemcpy(&total_count, d_result_, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    auto t_end = std::chrono::high_resolution_clock::now();
    printf("[BFS-V] Search complete: %lu matches, %.1fms\n", total_count,
           std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() / 1000.0);

    return total_count;
}

uint64_t GPUBFSSearch::SearchSingleEdgeBFS(
    uint32_t v1, uint32_t v2, uint32_t label,
    uint32_t num_query_edges,
    uint32_t num_query_vertices
) {
    if (!csr_built_ || !query_set_) return 0;

    uint32_t Q = num_query_vertices;

    EnsureEdgesCapacity(1);
    EnsureBufCapacity(Q);

    // Copy single edge
    CUDA_CHECK(cudaMemcpy(d_edges_v1_, &v1, sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_v2_, &v2, sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_edges_label_, &label, sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Reset result
    uint64_t zero64 = 0;
    CUDA_CHECK(cudaMemcpy(d_result_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    // ---- Init ----
    CUDA_CHECK(cudaMemcpy(d_count_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));

    uint32_t total_init_tasks = 1 * num_query_edges * 2;
    int block = 256;
    int grid = (total_init_tasks + block - 1) / block;
    if (grid < 1) grid = 1;

    bfs_init_kernel<<<grid, block>>>(
        d_vlabels_, d_edges_v1_, d_edges_v2_, d_edges_label_,
        d_all_orders_, d_buf_a_, d_count_,
        1, num_query_edges, Q,
        static_cast<uint32_t>(MAX_BUF_MATCHES)
    );
    CUDA_CHECK(cudaGetLastError());

    uint64_t cc64 = 0;
    CUDA_CHECK(cudaMemcpy(&cc64, d_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
    uint32_t cur_count = static_cast<uint32_t>(std::min(cc64, (uint64_t)MAX_BUF_MATCHES));
    if (cur_count == 0) return 0;

    // ---- BFS levels ----
    uint32_t* cur_buf = d_buf_a_;
    uint32_t* next_buf = d_buf_b_;

    for (uint32_t depth = 2; depth < Q; depth++) {
        bool is_last = (depth == Q - 1);

        if (is_last) {
            grid = (uint32_t)(((uint64_t)cur_count * 32 + block - 1) / block);  // OPT(C): warp per pm
            bfs_count_kernel<<<grid, block>>>(
                d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                d_degrees_, d_all_orders_, cur_buf, cur_count, d_result_, depth, Q
            );
            CUDA_CHECK(cudaGetLastError());
        } else {
            CUDA_CHECK(cudaMemcpy(d_count_, &zero64, sizeof(uint64_t), cudaMemcpyHostToDevice));
            grid = (uint32_t)(((uint64_t)cur_count * 32 + block - 1) / block);  // OPT(C): warp per pm
            bfs_expand_kernel<<<grid, block>>>(
                d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_,
                d_degrees_, d_all_orders_,
                cur_buf, cur_count,
                next_buf, d_count_,
                depth, Q,
                static_cast<uint32_t>(MAX_BUF_MATCHES)
            );
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaDeviceSynchronize());

            uint64_t oc64 = 0;
            CUDA_CHECK(cudaMemcpy(&oc64, d_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
            uint32_t out_count = static_cast<uint32_t>(std::min(oc64, (uint64_t)MAX_BUF_MATCHES));
            if (out_count == 0) break;

            std::swap(cur_buf, next_buf);
            cur_count = out_count;
        }
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    uint64_t total_count;
    CUDA_CHECK(cudaMemcpy(&total_count, d_result_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
    return total_count;
}

void GPUBFSSearch::Destroy() {
    if (d_csr_offsets_) { CUDA_CHECK(cudaFree(d_csr_offsets_)); d_csr_offsets_ = nullptr; }
    if (d_csr_neighbors_) { CUDA_CHECK(cudaFree(d_csr_neighbors_)); d_csr_neighbors_ = nullptr; }
    if (d_csr_elabels_) { CUDA_CHECK(cudaFree(d_csr_elabels_)); d_csr_elabels_ = nullptr; }
    if (d_csr_timestamps_) { CUDA_CHECK(cudaFree(d_csr_timestamps_)); d_csr_timestamps_ = nullptr; }
    if (d_vlabels_) { CUDA_CHECK(cudaFree(d_vlabels_)); d_vlabels_ = nullptr; }
    if (d_degrees_) { CUDA_CHECK(cudaFree(d_degrees_)); d_degrees_ = nullptr; }
    if (d_all_orders_) { CUDA_CHECK(cudaFree(d_all_orders_)); d_all_orders_ = nullptr; }
    if (d_buf_a_) { CUDA_CHECK(cudaFree(d_buf_a_)); d_buf_a_ = nullptr; }
    if (d_buf_b_) { CUDA_CHECK(cudaFree(d_buf_b_)); d_buf_b_ = nullptr; }
    if (d_count_) { CUDA_CHECK(cudaFree(d_count_)); d_count_ = nullptr; }
    if (d_result_) { CUDA_CHECK(cudaFree(d_result_)); d_result_ = nullptr; }
    if (d_edges_v1_) { CUDA_CHECK(cudaFree(d_edges_v1_)); d_edges_v1_ = nullptr; }
    if (d_edges_v2_) { CUDA_CHECK(cudaFree(d_edges_v2_)); d_edges_v2_ = nullptr; }
    if (d_edges_label_) { CUDA_CHECK(cudaFree(d_edges_label_)); d_edges_label_ = nullptr; }
    if (d_edges_max_ts_) { CUDA_CHECK(cudaFree(d_edges_max_ts_)); d_edges_max_ts_ = nullptr; }
    edges_capacity_ = 0;
    csr_built_ = false;
    query_set_ = false;
    padded_csr_ = false;
}
