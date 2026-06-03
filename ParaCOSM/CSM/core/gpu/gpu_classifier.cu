#include "core/gpu/gpu_classifier.h"
#include <cuda_runtime.h>
#include <cstdio>

// Max query edges we support in constant memory (typically ≤ 30)
static constexpr int MAX_QUERY_EDGES = 64;

// Query edge triples in constant memory for fast broadcast reads
__constant__ QueryEdgeTriple c_query_triples[MAX_QUERY_EDGES];
__constant__ int c_num_query_edges;

// ---------------------------------------------------------------------------
// CUDA Kernel: classify edges in parallel
// Each thread handles one data edge, checks against all query edge triples.
// Output: results[i] = true (safe) if no query edge matches, false (unsafe) otherwise.
// ---------------------------------------------------------------------------
__global__ void classify_edges_kernel(
    const uint32_t* __restrict__ vlabels,
    const EdgeToClassify* __restrict__ edges,
    uint8_t* __restrict__ results,
    int num_edges
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_edges) return;

    uint32_t v1 = edges[idx].v1;
    uint32_t v2 = edges[idx].v2;
    uint32_t elabel = edges[idx].label;

    uint32_t v1_label = vlabels[v1];
    uint32_t v2_label = vlabels[v2];

    int nq = c_num_query_edges;

    for (int i = 0; i < nq; i++) {
        uint32_t qs = c_query_triples[i].src_label;
        uint32_t qd = c_query_triples[i].dst_label;
        uint32_t qe = c_query_triples[i].edge_label;

        // Check direction 1: v1 -> v2
        if (qs == v1_label && qd == v2_label && qe == elabel) {
            results[idx] = 0;  // unsafe
            return;
        }
        // Check direction 2: v2 -> v1
        if (qs == v2_label && qd == v1_label && qe == elabel) {
            results[idx] = 0;  // unsafe
            return;
        }
    }

    results[idx] = 1;  // safe
}

// ---------------------------------------------------------------------------
// Helper: check CUDA errors
// ---------------------------------------------------------------------------
#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,  \
                    cudaGetErrorString(err));                                   \
        }                                                                      \
    } while (0)

// ---------------------------------------------------------------------------
// GPUClassifier implementation
// ---------------------------------------------------------------------------

GPUClassifier::GPUClassifier() = default;

GPUClassifier::~GPUClassifier() {
    Destroy();
}

void GPUClassifier::Init(
    const uint32_t* vlabels,
    size_t num_vertices,
    const QueryEdgeTriple* query_triples,
    size_t num_query_edges
) {
    if (initialized_) {
        Destroy();
    }

    if (num_query_edges > MAX_QUERY_EDGES) {
        fprintf(stderr, "GPUClassifier: too many query edges (%zu > %d)\n",
                num_query_edges, MAX_QUERY_EDGES);
        return;
    }

    num_query_edges_ = num_query_edges;

    // Copy query triples to constant memory
    int nq = static_cast<int>(num_query_edges);
    CUDA_CHECK(cudaMemcpyToSymbol(c_query_triples, query_triples,
                                   num_query_edges * sizeof(QueryEdgeTriple)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_num_query_edges, &nq, sizeof(int)));

    // Allocate and copy vertex labels
    vlabels_capacity_ = num_vertices;
    CUDA_CHECK(cudaMalloc(&d_vlabels_, num_vertices * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(d_vlabels_, vlabels,
                           num_vertices * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));

    // Pre-allocate edge batch buffers (initial capacity 4096)
    edges_capacity_ = 4096;
    CUDA_CHECK(cudaMalloc(&d_edges_, edges_capacity_ * sizeof(EdgeToClassify)));
CUDA_CHECK(cudaMalloc(&d_results_, edges_capacity_ * sizeof(uint8_t)));

    // Pinned host buffers
    CUDA_CHECK(cudaMallocHost(&h_edges_pinned_, edges_capacity_ * sizeof(EdgeToClassify)));
    CUDA_CHECK(cudaMallocHost(&h_results_pinned_, edges_capacity_ * sizeof(uint8_t)));

    initialized_ = true;
    fprintf(stdout, "GPUClassifier initialized: %zu vertices, %zu query edges\n",
            num_vertices, num_query_edges);
}

void GPUClassifier::UpdateVertexLabels(const uint32_t* vlabels, size_t new_num_vertices) {
    if (!initialized_) return;

    if (new_num_vertices > vlabels_capacity_) {
        // Reallocate
        if (d_vlabels_) {
            CUDA_CHECK(cudaFree(d_vlabels_));
        }
        vlabels_capacity_ = new_num_vertices + (new_num_vertices / 4);  // grow by 25%
        CUDA_CHECK(cudaMalloc(&d_vlabels_, vlabels_capacity_ * sizeof(uint32_t)));
    }

    CUDA_CHECK(cudaMemcpy(d_vlabels_, vlabels,
                           new_num_vertices * sizeof(uint32_t),
                           cudaMemcpyHostToDevice));
}

void GPUClassifier::ClassifyBatch(
    const EdgeToClassify* edges,
    size_t num_edges,
    uint8_t* results
) {
    if (!initialized_ || num_edges == 0) return;

    // Grow device/pinned buffers if needed
    if (num_edges > edges_capacity_) {
        if (d_edges_) CUDA_CHECK(cudaFree(d_edges_));
        if (d_results_) CUDA_CHECK(cudaFree(d_results_));
        if (h_edges_pinned_) CUDA_CHECK(cudaFreeHost(h_edges_pinned_));
        if (h_results_pinned_) CUDA_CHECK(cudaFreeHost(h_results_pinned_));

        edges_capacity_ = num_edges + (num_edges / 4);
        CUDA_CHECK(cudaMalloc(&d_edges_, edges_capacity_ * sizeof(EdgeToClassify)));
        CUDA_CHECK(cudaMalloc(&d_results_, edges_capacity_ * sizeof(uint8_t)));
        CUDA_CHECK(cudaMallocHost(&h_edges_pinned_, edges_capacity_ * sizeof(EdgeToClassify)));
        CUDA_CHECK(cudaMallocHost(&h_results_pinned_, edges_capacity_ * sizeof(uint8_t)));
    }

    // Copy edges to pinned buffer, then H2D
    memcpy(h_edges_pinned_, edges, num_edges * sizeof(EdgeToClassify));
    CUDA_CHECK(cudaMemcpy(d_edges_, h_edges_pinned_,
                           num_edges * sizeof(EdgeToClassify),
                           cudaMemcpyHostToDevice));

    // Launch kernel
    int block_size = 256;
    int grid_size = (static_cast<int>(num_edges) + block_size - 1) / block_size;
    classify_edges_kernel<<<grid_size, block_size>>>(
        d_vlabels_, d_edges_, d_results_, static_cast<int>(num_edges)
    );

    // D2H results
    CUDA_CHECK(cudaMemcpy(h_results_pinned_, d_results_,
                           num_edges * sizeof(uint8_t),
                           cudaMemcpyDeviceToHost));

    // Copy to caller's buffer
    memcpy(results, h_results_pinned_, num_edges * sizeof(uint8_t));
}

void GPUClassifier::Destroy() {
    if (!initialized_) return;

    if (d_vlabels_) { CUDA_CHECK(cudaFree(d_vlabels_)); d_vlabels_ = nullptr; }
    if (d_query_triples_) { CUDA_CHECK(cudaFree(d_query_triples_)); d_query_triples_ = nullptr; }
    if (d_edges_) { CUDA_CHECK(cudaFree(d_edges_)); d_edges_ = nullptr; }
    if (d_results_) { CUDA_CHECK(cudaFree(d_results_)); d_results_ = nullptr; }
    if (h_edges_pinned_) { CUDA_CHECK(cudaFreeHost(h_edges_pinned_)); h_edges_pinned_ = nullptr; }
    if (h_results_pinned_) { CUDA_CHECK(cudaFreeHost(h_results_pinned_)); h_results_pinned_ = nullptr; }

    initialized_ = false;
}
