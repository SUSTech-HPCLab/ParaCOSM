#include "core/gpu/gpu_candidate_filter.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,  \
                    cudaGetErrorString(err));                                   \
        }                                                                      \
    } while (0)

// ---------------------------------------------------------------------------
// Device-side joinability constraint (copied to constant memory)
// ---------------------------------------------------------------------------
struct DeviceJoinConstraint {
    uint32_t expected_elabel;
    uint32_t nbr_offset;
    uint32_t nbr_size;
};

static constexpr int MAX_CONSTRAINTS_CONST = 16;
__constant__ DeviceJoinConstraint c_constraints[MAX_CONSTRAINTS_CONST];
__constant__ int c_num_constraints;
__constant__ uint32_t c_expected_vlabel;
__constant__ uint32_t c_expected_elabel;
__constant__ int c_num_visited;
__constant__ int c_check_visited;

// ---------------------------------------------------------------------------
// Device: binary search in sorted array
// ---------------------------------------------------------------------------
__device__ int binary_search_sorted(
    const uint32_t* __restrict__ arr,
    int size,
    uint32_t target
) {
    int lo = 0, hi = size - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        uint32_t val = arr[mid];
        if (val == target) return mid;
        if (val < target) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;  // not found
}

// ---------------------------------------------------------------------------
// Kernel: filter candidates
// Each thread processes one candidate from m[u_min]'s neighbor list.
// Checks: (1) vertex label, (2) edge label to u_min, (3) joinability
// with all constraint vertices, (4) not in visited set.
// ---------------------------------------------------------------------------
__global__ void filter_candidates_kernel(
    const uint32_t* __restrict__ vlabels,
    const uint32_t* __restrict__ candidates,
    const uint32_t* __restrict__ candidate_elabels,
    const uint32_t* __restrict__ flat_nbrs,
    const uint32_t* __restrict__ flat_elabels,
    const uint32_t* __restrict__ visited,
    uint32_t* __restrict__ valid_indices,
    uint32_t* __restrict__ valid_count,
    int num_candidates
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_candidates) return;

    uint32_t v = candidates[idx];

    // 1. Vertex label check
    if (vlabels[v] != c_expected_vlabel) return;

    // 2. Edge label to u_min check
    if (candidate_elabels[idx] != c_expected_elabel) return;

    // 3. Joinability check: for each constraint, v must be adjacent to
    //    the constraint vertex with the correct edge label
    int nc = c_num_constraints;
    for (int c = 0; c < nc; c++) {
        uint32_t nbr_off = c_constraints[c].nbr_offset;
        uint32_t nbr_sz  = c_constraints[c].nbr_size;
        uint32_t exp_el  = c_constraints[c].expected_elabel;

        // Binary search for v in the sorted neighbor list
        int pos = binary_search_sorted(flat_nbrs + nbr_off, nbr_sz, v);
        if (pos < 0) return;  // v not adjacent to constraint vertex

        // Check edge label
        if (flat_elabels[nbr_off + pos] != exp_el) return;
    }

    // 4. Visited check (linear scan, small set)
    if (c_check_visited) {
        int nv = c_num_visited;
        for (int i = 0; i < nv; i++) {
            if (visited[i] == v) return;
        }
    }

    // All checks passed — add to valid list via atomic
    uint32_t slot = atomicAdd(valid_count, 1);
    valid_indices[slot] = idx;
}

// ---------------------------------------------------------------------------
// GPUCandidateFilter implementation
// ---------------------------------------------------------------------------

GPUCandidateFilter::GPUCandidateFilter() = default;

GPUCandidateFilter::~GPUCandidateFilter() {
    Destroy();
}

void GPUCandidateFilter::Init(const uint32_t* vlabels, size_t num_vertices) {
    if (initialized_) Destroy();

    vlabels_capacity_ = num_vertices + num_vertices / 4;
    CUDA_CHECK(cudaMalloc(&d_vlabels_, vlabels_capacity_ * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, vlabels,
                           num_vertices * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));

    // Pre-allocate buffers for typical sizes
    EnsureCapacity(8192, 65536, 32);

    // Allocate valid_count on device and pinned host
    CUDA_CHECK(cudaMalloc(&d_valid_count_, sizeof(uint32_t)));
    CUDA_CHECK(cudaMallocHost(&h_valid_count_, sizeof(uint32_t)));

    initialized_ = true;
}

void GPUCandidateFilter::UpdateVertexLabels(const uint32_t* vlabels, size_t num_vertices) {
    if (!initialized_) return;
    if (num_vertices > vlabels_capacity_) {
        if (d_vlabels_) CUDA_CHECK(cudaFree(d_vlabels_));
        vlabels_capacity_ = num_vertices + num_vertices / 4;
        CUDA_CHECK(cudaMalloc(&d_vlabels_, vlabels_capacity_ * sizeof(uint32_t)));
    }
    CUDA_CHECK(cudaMemcpy(d_vlabels_, vlabels,
                           num_vertices * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));
}

void GPUCandidateFilter::EnsureCapacity(size_t num_candidates, size_t flat_nbrs_total, size_t num_visited) {
    if (num_candidates > candidates_cap_) {
        if (d_candidates_) CUDA_CHECK(cudaFree(d_candidates_));
        if (d_candidate_elabels_) CUDA_CHECK(cudaFree(d_candidate_elabels_));
        if (d_valid_indices_) CUDA_CHECK(cudaFree(d_valid_indices_));
        if (h_valid_indices_) CUDA_CHECK(cudaFreeHost(h_valid_indices_));

        candidates_cap_ = num_candidates + num_candidates / 4;
        CUDA_CHECK(cudaMalloc(&d_candidates_, candidates_cap_ * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_candidate_elabels_, candidates_cap_ * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_valid_indices_, candidates_cap_ * sizeof(uint32_t)));
        CUDA_CHECK(cudaMallocHost(&h_valid_indices_, candidates_cap_ * sizeof(uint32_t)));
    }

    if (flat_nbrs_total > flat_nbrs_cap_) {
        if (d_flat_nbrs_) CUDA_CHECK(cudaFree(d_flat_nbrs_));
        if (d_flat_elabels_) CUDA_CHECK(cudaFree(d_flat_elabels_));

        flat_nbrs_cap_ = flat_nbrs_total + flat_nbrs_total / 4;
        CUDA_CHECK(cudaMalloc(&d_flat_nbrs_, flat_nbrs_cap_ * sizeof(uint32_t)));
        CUDA_CHECK(cudaMalloc(&d_flat_elabels_, flat_nbrs_cap_ * sizeof(uint32_t)));
    }

    if (num_visited > visited_cap_) {
        if (d_visited_) CUDA_CHECK(cudaFree(d_visited_));
        visited_cap_ = num_visited + 16;
        CUDA_CHECK(cudaMalloc(&d_visited_, visited_cap_ * sizeof(uint32_t)));
    }
}

uint32_t GPUCandidateFilter::FilterCandidates(
    const CandidateFilterTask& task,
    uint32_t* valid_indices
) {
    if (!initialized_ || task.num_candidates == 0) return 0;

    EnsureCapacity(task.num_candidates, task.flat_nbrs_total, task.num_visited);

    // Copy candidates to device
    CUDA_CHECK(cudaMemcpy(d_candidates_, task.candidates,
                           task.num_candidates * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_candidate_elabels_, task.candidate_elabels,
                           task.num_candidates * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));

    // Copy flattened neighbor lists
    if (task.flat_nbrs_total > 0) {
        CUDA_CHECK(cudaMemcpy(d_flat_nbrs_, task.flat_nbrs,
                               task.flat_nbrs_total * sizeof(uint32_t),
                               cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_flat_elabels_, task.flat_elabels,
                               task.flat_nbrs_total * sizeof(uint32_t),
                               cudaMemcpyHostToDevice));
    }

    // Copy visited vertices
    if (task.num_visited > 0 && task.check_visited) {
        CUDA_CHECK(cudaMemcpy(d_visited_, task.visited_vertices,
                               task.num_visited * sizeof(uint32_t),
                               cudaMemcpyHostToDevice));
    }

    // Set constant memory
    DeviceJoinConstraint host_constraints[MAX_CONSTRAINTS_CONST];
    for (uint32_t i = 0; i < task.num_constraints; i++) {
        host_constraints[i].expected_elabel = task.constraints[i].expected_elabel;
        host_constraints[i].nbr_offset = task.constraints[i].nbr_offset;
        host_constraints[i].nbr_size = task.constraints[i].nbr_size;
    }
    int nc = static_cast<int>(task.num_constraints);
    int nv = static_cast<int>(task.num_visited);
    int cv = task.check_visited ? 1 : 0;

    CUDA_CHECK(cudaMemcpyToSymbol(c_constraints, host_constraints,
                                   task.num_constraints * sizeof(DeviceJoinConstraint)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_constraints, &nc, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_expected_vlabel, &task.expected_vlabel, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_expected_elabel, &task.expected_elabel, sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_visited, &nv, sizeof(int)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_check_visited, &cv, sizeof(int)));

    // Reset valid count
    uint32_t zero = 0;
    CUDA_CHECK(cudaMemcpy(d_valid_count_, &zero, sizeof(uint32_t), cudaMemcpyHostToDevice));

    // Launch kernel
    int block = 256;
    int grid = (static_cast<int>(task.num_candidates) + block - 1) / block;
    filter_candidates_kernel<<<grid, block>>>(
        d_vlabels_,
        d_candidates_,
        d_candidate_elabels_,
        d_flat_nbrs_,
        d_flat_elabels_,
        d_visited_,
        d_valid_indices_,
        d_valid_count_,
        static_cast<int>(task.num_candidates)
    );

    // Get results
    CUDA_CHECK(cudaMemcpy(h_valid_count_, d_valid_count_, sizeof(uint32_t),
                           cudaMemcpyDeviceToHost));
    uint32_t count = *h_valid_count_;

    if (count > 0) {
        CUDA_CHECK(cudaMemcpy(h_valid_indices_, d_valid_indices_,
                               count * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost));
        memcpy(valid_indices, h_valid_indices_, count * sizeof(uint32_t));
    }

    return count;
}

void GPUCandidateFilter::Destroy() {
    if (!initialized_) return;

    if (d_vlabels_) { CUDA_CHECK(cudaFree(d_vlabels_)); d_vlabels_ = nullptr; }
    if (d_candidates_) { CUDA_CHECK(cudaFree(d_candidates_)); d_candidates_ = nullptr; }
    if (d_candidate_elabels_) { CUDA_CHECK(cudaFree(d_candidate_elabels_)); d_candidate_elabels_ = nullptr; }
    if (d_flat_nbrs_) { CUDA_CHECK(cudaFree(d_flat_nbrs_)); d_flat_nbrs_ = nullptr; }
    if (d_flat_elabels_) { CUDA_CHECK(cudaFree(d_flat_elabels_)); d_flat_elabels_ = nullptr; }
    if (d_visited_) { CUDA_CHECK(cudaFree(d_visited_)); d_visited_ = nullptr; }
    if (d_valid_indices_) { CUDA_CHECK(cudaFree(d_valid_indices_)); d_valid_indices_ = nullptr; }
    if (d_valid_count_) { CUDA_CHECK(cudaFree(d_valid_count_)); d_valid_count_ = nullptr; }
    if (h_valid_indices_) { CUDA_CHECK(cudaFreeHost(h_valid_indices_)); h_valid_indices_ = nullptr; }
    if (h_valid_count_) { CUDA_CHECK(cudaFreeHost(h_valid_count_)); h_valid_count_ = nullptr; }

    candidates_cap_ = flat_nbrs_cap_ = visited_cap_ = 0;
    initialized_ = false;
}
