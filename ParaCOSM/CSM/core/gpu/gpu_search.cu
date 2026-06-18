#include "core/gpu/gpu_search.h"
#include "graph_storage/graph.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#include <chrono>

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,  \
                    cudaGetErrorString(err));                                   \
        }                                                                      \
    } while (0)

// ============================================================
// Constants for query graph (small, fits in constant memory)
// ============================================================
static constexpr int MAX_Q = 20;        // max query vertices
static constexpr int MAX_Q_EDGES = 64;  // max query edge entries (2 * num_edges)

__constant__ uint32_t c_q_num_vertices;
__constant__ uint32_t c_q_vlabels[MAX_Q];
__constant__ uint32_t c_q_offsets[MAX_Q + 1];     // query CSR offsets
__constant__ uint32_t c_q_neighbors[MAX_Q_EDGES]; // query CSR neighbors
__constant__ uint32_t c_q_elabels[MAX_Q_EDGES];   // query edge labels
__constant__ uint32_t c_order[MAX_Q];              // matching order for current order_index

// For inter-update: all matching orders stored in global memory
// (too large for constant memory: num_query_edges * MAX_Q)
// Set by SetAllMatchingOrders(), accessed by batch kernel

// ============================================================
// Device binary search on sorted array
// ============================================================
__device__ __forceinline__ int device_bsearch(
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
// Device recursive DFS: FindMatches_local equivalent
// ============================================================
__device__ uint64_t gpu_dfs(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,
    uint32_t* m,         // local matching array (in registers/local mem)
    uint32_t depth
) {
    uint32_t Q = c_q_num_vertices;
    if (depth >= Q) return 0;

    uint32_t u = c_order[depth];

    // Find u_min: matched query neighbor of u with smallest data degree
    uint32_t u_min = UINT32_MAX;
    uint32_t u_min_label = 0;
    uint32_t u_min_deg = UINT32_MAX;

    uint32_t q_start = c_q_offsets[u];
    uint32_t q_end = c_q_offsets[u + 1];

    for (uint32_t j = q_start; j < q_end; j++) {
        uint32_t u_other = c_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        uint32_t deg = degrees[m[u_other]];
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = c_q_elabels[j];
        }
    }

    if (u_min == UINT32_MAX) return 0;

    uint64_t count = 0;
    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = degrees[m[u_min]];

    for (uint32_t i = 0; i < nbr_count; i++) {
        uint32_t v = csr_neighbors[nbr_start + i];

        // 1. Label check
        if (vlabels[v] != c_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;

        // 2. Joinability: v must be adjacent to all other matched neighbors of u
        bool joinable = true;
        for (uint32_t j = q_start; j < q_end; j++) {
            uint32_t u_other = c_q_neighbors[j];
            if (m[u_other] == UINT32_MAX || u_other == u_min) continue;

            uint32_t s = csr_offsets[m[u_other]];
            uint32_t nc = degrees[m[u_other]];
            int pos = device_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != c_q_elabels[j]) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. Visited check: v must not already appear in m[]
        bool visited = false;
        for (uint32_t d = 0; d < Q; d++) {
            if (m[d] == v) { visited = true; break; }
        }
        if (visited) continue;

        // 4. Valid candidate — set and recurse
        m[u] = v;
        if (depth == Q - 1) {
            count++;
        } else {
            count += gpu_dfs(csr_offsets, csr_neighbors, csr_elabels,
                             vlabels, degrees, m, depth + 1);
        }
        m[u] = UINT32_MAX;
    }

    return count;
}

// ============================================================
// Device recursive DFS with explicit order (for inter-update batch)
// ============================================================
__device__ uint64_t gpu_dfs_order(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,
    const uint32_t* __restrict__ order,  // matching order for this search
    uint32_t* m,
    uint32_t depth
) {
    uint32_t Q = c_q_num_vertices;
    if (depth >= Q) return 0;

    uint32_t u = order[depth];

    uint32_t u_min = UINT32_MAX;
    uint32_t u_min_label = 0;
    uint32_t u_min_deg = UINT32_MAX;

    uint32_t q_start = c_q_offsets[u];
    uint32_t q_end = c_q_offsets[u + 1];

    for (uint32_t j = q_start; j < q_end; j++) {
        uint32_t u_other = c_q_neighbors[j];
        if (m[u_other] == UINT32_MAX) continue;
        uint32_t deg = degrees[m[u_other]];
        if (deg < u_min_deg) {
            u_min_deg = deg;
            u_min = u_other;
            u_min_label = c_q_elabels[j];
        }
    }

    if (u_min == UINT32_MAX) return 0;

    uint64_t count = 0;
    uint32_t nbr_start = csr_offsets[m[u_min]];
    uint32_t nbr_count = degrees[m[u_min]];

    for (uint32_t i = 0; i < nbr_count; i++) {
        uint32_t v = csr_neighbors[nbr_start + i];

        if (vlabels[v] != c_q_vlabels[u]) continue;
        if (csr_elabels[nbr_start + i] != u_min_label) continue;

        bool joinable = true;
        for (uint32_t j = q_start; j < q_end; j++) {
            uint32_t u_other = c_q_neighbors[j];
            if (m[u_other] == UINT32_MAX || u_other == u_min) continue;
            uint32_t s = csr_offsets[m[u_other]];
            uint32_t nc = degrees[m[u_other]];
            int pos = device_bsearch(csr_neighbors + s, nc, v);
            if (pos < 0 || csr_elabels[s + pos] != c_q_elabels[j]) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;

        bool visited = false;
        for (uint32_t d = 0; d < Q; d++) {
            if (m[d] == v) { visited = true; break; }
        }
        if (visited) continue;

        m[u] = v;
        if (depth == Q - 1) {
            count++;
        } else {
            count += gpu_dfs_order(csr_offsets, csr_neighbors, csr_elabels,
                                   vlabels, degrees, order, m, depth + 1);
        }
        m[u] = UINT32_MAX;
    }

    return count;
}

// ============================================================
// Batch kernel: one thread per (data_edge, query_edge, direction)
// Total threads = num_data_edges * num_query_edges * 2
// ============================================================
__global__ void gpu_batch_search_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,
    const uint32_t* __restrict__ edges_v1,
    const uint32_t* __restrict__ edges_v2,
    const uint32_t* __restrict__ edges_label,
    const uint32_t* __restrict__ all_orders,  // [num_q_edges * Q]
    uint64_t* __restrict__ total_count,
    uint32_t num_data_edges,
    uint32_t num_q_edges
) {
    // Thread layout: idx = edge_idx * (num_q_edges * 2) + qe_idx * 2 + dir
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t total_tasks = num_data_edges * num_q_edges * 2;
    if (idx >= total_tasks) return;

    uint32_t Q = c_q_num_vertices;

    // Decode task
    uint32_t edge_idx = idx / (num_q_edges * 2);
    uint32_t remainder = idx % (num_q_edges * 2);
    uint32_t qe_idx = remainder / 2;
    uint32_t dir = remainder % 2;

    uint32_t v1 = edges_v1[edge_idx];
    uint32_t v2 = edges_v2[edge_idx];
    uint32_t elabel = edges_label[edge_idx];

    // Get matching order for this query edge
    const uint32_t* order = all_orders + qe_idx * Q;
    uint32_t u1 = order[0];
    uint32_t u2 = order[1];

    // Get query edge label triple for (u1, u2)
    // Search in query CSR for edge u1-u2
    uint32_t q_start = c_q_offsets[u1];
    uint32_t q_end = c_q_offsets[u1 + 1];
    uint32_t found_elabel = UINT32_MAX;
    for (uint32_t j = q_start; j < q_end; j++) {
        if (c_q_neighbors[j] == u2) {
            found_elabel = c_q_elabels[j];
            break;
        }
    }
    if (found_elabel == UINT32_MAX) return;

    // Direction: dir=0 → map u1→v1, u2→v2; dir=1 → u1→v2, u2→v1
    uint32_t map_v1, map_v2;
    if (dir == 0) { map_v1 = v1; map_v2 = v2; }
    else          { map_v1 = v2; map_v2 = v1; }

    // Label check: vlabels[map_v1] == q_vlabels[u1], vlabels[map_v2] == q_vlabels[u2], edge label matches
    if (vlabels[map_v1] != c_q_vlabels[u1]) return;
    if (vlabels[map_v2] != c_q_vlabels[u2]) return;
    if (elabel != found_elabel) return;

    // Initialize matching: m[u1] = map_v1, m[u2] = map_v2, rest = UNMATCHED
    uint32_t m[MAX_Q];
    for (uint32_t j = 0; j < Q; j++) m[j] = UINT32_MAX;
    m[u1] = map_v1;
    m[u2] = map_v2;

    // Run DFS from depth=2 using this order
    uint64_t local_count = gpu_dfs_order(
        csr_offsets, csr_neighbors, csr_elabels, vlabels, degrees,
        order, m, 2
    );

    if (local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(total_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// Main kernel: one thread per vertex_vector entry
// ============================================================
__global__ void gpu_search_kernel(
    const uint32_t* __restrict__ csr_offsets,
    const uint32_t* __restrict__ csr_neighbors,
    const uint32_t* __restrict__ csr_elabels,
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ degrees,
    const uint32_t* __restrict__ vv_v,      // Layer 1 candidate v3
    const uint32_t* __restrict__ vv_umin,   // u_min for ProcessNeighbor
    const uint32_t* __restrict__ vv_ulabel, // u_min_label
    const uint32_t* __restrict__ vv_i,      // candidate index
    const uint32_t* __restrict__ vv_m,      // flattened partial matches [N * Q]
    uint64_t* __restrict__ total_count,     // atomic result
    uint32_t num_entries,
    uint32_t start_depth     // depth parameter from PFM2 (typically 2)
) {
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_entries) return;

    uint32_t Q = c_q_num_vertices;

    // Load partial match into local memory
    uint32_t m[MAX_Q];
    for (uint32_t j = 0; j < Q; j++) {
        m[j] = vv_m[idx * Q + j];
    }

    uint32_t u_min3 = vv_umin[idx];
    uint32_t u_min_label3 = vv_ulabel[idx];
    uint32_t i2 = vv_i[idx];

    uint32_t u2 = c_order[start_depth + 1]; // query vertex at depth+1

    // ---- ProcessNeighbor equivalent ----
    // Get the actual candidate at depth+1
    uint32_t nbr_start = csr_offsets[m[u_min3]];
    uint32_t nbr_count = degrees[m[u_min3]];

    // Bounds check
    if (i2 >= nbr_count) return;

    uint32_t v = csr_neighbors[nbr_start + i2];

    // 1. Label check
    if (vlabels[v] != c_q_vlabels[u2]) return;
    if (csr_elabels[nbr_start + i2] != u_min_label3) return;

    // 2. Joinability check for u2's query neighbors
    uint32_t q2_start = c_q_offsets[u2];
    uint32_t q2_end = c_q_offsets[u2 + 1];

    for (uint32_t j = q2_start; j < q2_end; j++) {
        uint32_t u_other = c_q_neighbors[j];
        if (m[u_other] == UINT32_MAX || u_other == u_min3) continue;

        uint32_t s = csr_offsets[m[u_other]];
        uint32_t nc = degrees[m[u_other]];
        int pos = device_bsearch(csr_neighbors + s, nc, v);
        if (pos < 0 || csr_elabels[s + pos] != c_q_elabels[j]) return;
    }

    // 3. Visited check (v must not be in m[], including v3)
    for (uint32_t d = 0; d < Q; d++) {
        if (m[d] == v) return;
    }

    // ---- Valid candidate: set match and do full DFS ----
    m[u2] = v;

    uint64_t local_count;
    if (start_depth + 1 == Q - 1) {
        local_count = 1;
    } else {
        local_count = gpu_dfs(csr_offsets, csr_neighbors, csr_elabels,
                              vlabels, degrees, m, start_depth + 2);
    }

    if (local_count > 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(total_count),
                  static_cast<unsigned long long>(local_count));
    }
}

// ============================================================
// GPUSearchEngine implementation
// ============================================================

GPUSearchEngine::GPUSearchEngine() = default;

GPUSearchEngine::~GPUSearchEngine() {
    Destroy();
}

void GPUSearchEngine::BuildCSR(const Graph& data) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Free old CSR if exists
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
    total_csr_edges_ = total_padded;

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

    // Allocate and copy to GPU
    CUDA_CHECK(cudaMalloc(&d_csr_offsets_, (V + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_neighbors_, total_padded * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_csr_elabels_, total_padded * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vlabels_, V * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_degrees_, V * sizeof(uint32_t)));

    CUDA_CHECK(cudaMemcpy(d_csr_offsets_, h_offsets_.data(),
                           (V + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_, neighbors.data(),
                           total_padded * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_, elabels.data(),
                           total_padded * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, data.vlabels_.data(),
                           V * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_, h_degrees.data(),
                           V * sizeof(uint32_t), cudaMemcpyHostToDevice));

    csr_built_ = true;

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    printf("[GPUSearch] CSR built: V=%u, padded_edges=%u, %.1fms\n", V, total_padded, ms);
}

bool GPUSearchEngine::IncrementalUpdate(const Graph& data, uint32_t v1, uint32_t v2) {
    if (!csr_built_) return false;

    // Check if padding can accommodate the new degrees
    uint32_t new_deg1 = data.GetDegree(v1);
    uint32_t new_deg2 = data.GetDegree(v2);

    if (v1 >= num_vertices_ || v2 >= num_vertices_ ||
        new_deg1 > h_capacities_[v1] || new_deg2 > h_capacities_[v2]) {
        return false;  // need full rebuild
    }

    // Update v1's neighbor list on GPU
    const auto& nbrs1 = data.GetNeighbors(v1);
    const auto& labs1 = data.GetNeighborLabels(v1);
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_ + h_offsets_[v1], nbrs1.data(),
                           new_deg1 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_ + h_offsets_[v1], labs1.data(),
                           new_deg1 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    // Update degree
    CUDA_CHECK(cudaMemcpy(d_degrees_ + v1, &new_deg1, sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Update v2's neighbor list on GPU
    const auto& nbrs2 = data.GetNeighbors(v2);
    const auto& labs2 = data.GetNeighborLabels(v2);
    CUDA_CHECK(cudaMemcpy(d_csr_neighbors_ + h_offsets_[v2], nbrs2.data(),
                           new_deg2 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_csr_elabels_ + h_offsets_[v2], labs2.data(),
                           new_deg2 * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_degrees_ + v2, &new_deg2, sizeof(uint32_t), cudaMemcpyHostToDevice));

    return true;
}

void GPUSearchEngine::SetupQuery(const Graph& query) {
    uint32_t Q = query.NumVertices();
    if (Q > MAX_Q) {
        fprintf(stderr, "[GPUSearch] Query too large: %u > %d\n", Q, MAX_Q);
        return;
    }

    // Query vertex labels
    std::vector<uint32_t> q_vlabels(Q);
    for (uint32_t u = 0; u < Q; u++) {
        q_vlabels[u] = query.GetVertexLabel(u);
    }

    // Query CSR
    std::vector<uint32_t> q_offsets(Q + 1);
    q_offsets[0] = 0;
    for (uint32_t u = 0; u < Q; u++) {
        q_offsets[u + 1] = q_offsets[u] + query.GetDegree(u);
    }
    uint32_t q_total = q_offsets[Q];
    if (q_total > MAX_Q_EDGES) {
        fprintf(stderr, "[GPUSearch] Query edges too many: %u > %d\n", q_total, MAX_Q_EDGES);
        return;
    }

    std::vector<uint32_t> q_neighbors(q_total);
    std::vector<uint32_t> q_elabels(q_total);
    for (uint32_t u = 0; u < Q; u++) {
        const auto& nbrs = query.GetNeighbors(u);
        const auto& labs = query.GetNeighborLabels(u);
        uint32_t off = q_offsets[u];
        for (size_t j = 0; j < nbrs.size(); j++) {
            q_neighbors[off + j] = nbrs[j];
            q_elabels[off + j] = labs[j];
        }
    }

    // Copy to constant memory
    CUDA_CHECK(cudaMemcpyToSymbol(c_q_num_vertices, &Q, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_q_vlabels, q_vlabels.data(), Q * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_q_offsets, q_offsets.data(), (Q + 1) * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_q_neighbors, q_neighbors.data(), q_total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_q_elabels, q_elabels.data(), q_total * sizeof(uint32_t)));

    query_set_ = true;
    printf("[GPUSearch] Query set: %u vertices, %u edge entries\n", Q, q_total);
}

void GPUSearchEngine::SetMatchingOrder(const uint32_t* order, uint32_t num_vertices) {
    CUDA_CHECK(cudaMemcpyToSymbol(c_order, order, num_vertices * sizeof(uint32_t)));
}

void GPUSearchEngine::EnsureVVCapacity(size_t n, uint32_t q) {
    if (n <= vv_capacity_) return;

    // Free old
    if (d_vv_v_) CUDA_CHECK(cudaFree(d_vv_v_));
    if (d_vv_umin_) CUDA_CHECK(cudaFree(d_vv_umin_));
    if (d_vv_ulabel_) CUDA_CHECK(cudaFree(d_vv_ulabel_));
    if (d_vv_i_) CUDA_CHECK(cudaFree(d_vv_i_));
    if (d_vv_m_) CUDA_CHECK(cudaFree(d_vv_m_));

    vv_capacity_ = n + n / 4;
    CUDA_CHECK(cudaMalloc(&d_vv_v_, vv_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vv_umin_, vv_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vv_ulabel_, vv_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vv_i_, vv_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_vv_m_, vv_capacity_ * q * sizeof(uint32_t)));
}

uint64_t GPUSearchEngine::SearchFromVertexVector(
    const std::vector<std::tuple<uint32_t, uint32_t, uint32_t,
                                  std::vector<uint32_t>, uint32_t>>& vv,
    uint32_t depth,
    uint32_t num_query_vertices
) {
    uint32_t N = static_cast<uint32_t>(vv.size());
    uint32_t Q = num_query_vertices;
    if (N == 0) return 0;

    EnsureVVCapacity(N, Q);

    // Flatten vertex_vector to SoA on host
    std::vector<uint32_t> h_v(N), h_umin(N), h_ulabel(N), h_i(N);
    std::vector<uint32_t> h_m(N * Q);

    for (uint32_t idx = 0; idx < N; idx++) {
        const auto& [v, um, ul, m_vec, i] = vv[idx];
        h_v[idx] = v;
        h_umin[idx] = um;
        h_ulabel[idx] = ul;
        h_i[idx] = i;
        for (uint32_t j = 0; j < Q && j < m_vec.size(); j++) {
            h_m[idx * Q + j] = m_vec[j];
        }
    }

    // H2D copy
    CUDA_CHECK(cudaMemcpy(d_vv_v_, h_v.data(), N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vv_umin_, h_umin.data(), N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vv_ulabel_, h_ulabel.data(), N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vv_i_, h_i.data(), N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vv_m_, h_m.data(), N * Q * sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Reset result counter
    if (!d_total_count_) {
        CUDA_CHECK(cudaMalloc(&d_total_count_, sizeof(uint64_t)));
    }
    uint64_t zero = 0;
    CUDA_CHECK(cudaMemcpy(d_total_count_, &zero, sizeof(uint64_t), cudaMemcpyHostToDevice));

    // Set stack size for device-side recursion
    // ~120 bytes per recursion level, max ~14 levels
    CUDA_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 8192));

    // Launch kernel
    int block = 128;
    int grid = (N + block - 1) / block;

    gpu_search_kernel<<<grid, block>>>(
        d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_, d_degrees_,
        d_vv_v_, d_vv_umin_, d_vv_ulabel_, d_vv_i_, d_vv_m_,
        d_total_count_,
        N, depth
    );

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Read result
    uint64_t total_count;
    CUDA_CHECK(cudaMemcpy(&total_count, d_total_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    return total_count;
}

void GPUSearchEngine::Destroy() {
    if (d_csr_offsets_) { CUDA_CHECK(cudaFree(d_csr_offsets_)); d_csr_offsets_ = nullptr; }
    if (d_csr_neighbors_) { CUDA_CHECK(cudaFree(d_csr_neighbors_)); d_csr_neighbors_ = nullptr; }
    if (d_csr_elabels_) { CUDA_CHECK(cudaFree(d_csr_elabels_)); d_csr_elabels_ = nullptr; }
    if (d_vlabels_) { CUDA_CHECK(cudaFree(d_vlabels_)); d_vlabels_ = nullptr; }
    if (d_degrees_) { CUDA_CHECK(cudaFree(d_degrees_)); d_degrees_ = nullptr; }
    if (d_vv_v_) { CUDA_CHECK(cudaFree(d_vv_v_)); d_vv_v_ = nullptr; }
    if (d_vv_umin_) { CUDA_CHECK(cudaFree(d_vv_umin_)); d_vv_umin_ = nullptr; }
    if (d_vv_ulabel_) { CUDA_CHECK(cudaFree(d_vv_ulabel_)); d_vv_ulabel_ = nullptr; }
    if (d_vv_i_) { CUDA_CHECK(cudaFree(d_vv_i_)); d_vv_i_ = nullptr; }
    if (d_vv_m_) { CUDA_CHECK(cudaFree(d_vv_m_)); d_vv_m_ = nullptr; }
    if (d_total_count_) { CUDA_CHECK(cudaFree(d_total_count_)); d_total_count_ = nullptr; }
    if (d_batch_v1_) { CUDA_CHECK(cudaFree(d_batch_v1_)); d_batch_v1_ = nullptr; }
    if (d_batch_v2_) { CUDA_CHECK(cudaFree(d_batch_v2_)); d_batch_v2_ = nullptr; }
    if (d_batch_label_) { CUDA_CHECK(cudaFree(d_batch_label_)); d_batch_label_ = nullptr; }
    if (d_all_orders_) { CUDA_CHECK(cudaFree(d_all_orders_)); d_all_orders_ = nullptr; }
    vv_capacity_ = 0;
    batch_capacity_ = 0;
    csr_built_ = false;
    query_set_ = false;
}

void GPUSearchEngine::SetAllMatchingOrders(const uint32_t* all_orders, uint32_t num_edges, uint32_t num_vertices) {
    if (d_all_orders_) CUDA_CHECK(cudaFree(d_all_orders_));
    all_orders_num_edges_ = num_edges;
    size_t total = static_cast<size_t>(num_edges) * num_vertices;
    CUDA_CHECK(cudaMalloc(&d_all_orders_, total * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_all_orders_, all_orders, total * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void GPUSearchEngine::EnsureBatchCapacity(size_t n) {
    if (n <= batch_capacity_) return;
    if (d_batch_v1_) CUDA_CHECK(cudaFree(d_batch_v1_));
    if (d_batch_v2_) CUDA_CHECK(cudaFree(d_batch_v2_));
    if (d_batch_label_) CUDA_CHECK(cudaFree(d_batch_label_));
    batch_capacity_ = n + n / 4;
    CUDA_CHECK(cudaMalloc(&d_batch_v1_, batch_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_batch_v2_, batch_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_batch_label_, batch_capacity_ * sizeof(uint32_t)));
}

uint64_t GPUSearchEngine::SearchBatchEdges(
    const uint32_t* edges_v1,
    const uint32_t* edges_v2,
    const uint32_t* edges_label,
    uint32_t num_edges_data,
    uint32_t num_query_edges,
    uint32_t num_query_vertices
) {
    if (!csr_built_ || !query_set_ || num_edges_data == 0) return 0;

    EnsureBatchCapacity(num_edges_data);

    // Copy edge data to GPU
    CUDA_CHECK(cudaMemcpy(d_batch_v1_, edges_v1, num_edges_data * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_batch_v2_, edges_v2, num_edges_data * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_batch_label_, edges_label, num_edges_data * sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Reset counter
    if (!d_total_count_) CUDA_CHECK(cudaMalloc(&d_total_count_, sizeof(uint64_t)));
    uint64_t zero = 0;
    CUDA_CHECK(cudaMemcpy(d_total_count_, &zero, sizeof(uint64_t), cudaMemcpyHostToDevice));

    // Set stack size for recursion
    CUDA_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, 8192));

    // Total tasks: num_edges * num_query_edges * 2 (two directions)
    uint32_t total_tasks = num_edges_data * num_query_edges * 2;

    int block = 128;
    int grid = (total_tasks + block - 1) / block;

    gpu_batch_search_kernel<<<grid, block>>>(
        d_csr_offsets_, d_csr_neighbors_, d_csr_elabels_, d_vlabels_, d_degrees_,
        d_batch_v1_, d_batch_v2_, d_batch_label_,
        d_all_orders_, d_total_count_,
        num_edges_data, num_query_edges
    );

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    uint64_t total_count;
    CUDA_CHECK(cudaMemcpy(&total_count, d_total_count_, sizeof(uint64_t), cudaMemcpyDeviceToHost));
    return total_count;
}
