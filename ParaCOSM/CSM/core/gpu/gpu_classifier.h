#ifndef CORE_GPU_CLASSIFIER_H
#define CORE_GPU_CLASSIFIER_H

#include <cstddef>
#include <cstdint>
#include <vector>

// A query edge's label triple: (src_vertex_label, dst_vertex_label, edge_label)
struct QueryEdgeTriple {
    uint32_t src_label;
    uint32_t dst_label;
    uint32_t edge_label;
};

// A data edge to classify: (v1, v2, edge_label)
struct EdgeToClassify {
    uint32_t v1;
    uint32_t v2;
    uint32_t label;
};

class GPUClassifier {
public:
    GPUClassifier();
    ~GPUClassifier();

    // Non-copyable
    GPUClassifier(const GPUClassifier&) = delete;
    GPUClassifier& operator=(const GPUClassifier&) = delete;

    /**
     * Initialize GPU resources.
     * @param vlabels  Host array of vertex labels (data graph)
     * @param num_vertices  Number of vertices in data graph
     * @param query_triples  Array of query edge label triples
     * @param num_query_edges  Number of query edges
     */
    void Init(
        const uint32_t* vlabels,
        size_t num_vertices,
        const QueryEdgeTriple* query_triples,
        size_t num_query_edges
    );

    /**
     * Update vertex labels on GPU after new vertices are added.
     * @param vlabels  Full host array of vertex labels
     * @param new_num_vertices  New total number of vertices
     */
    void UpdateVertexLabels(const uint32_t* vlabels, size_t new_num_vertices);

    /**
     * Classify a batch of edges on GPU.
     * @param edges  Host array of edges to classify
     * @param num_edges  Number of edges
     * @param results  Host output array: 1 = safe, 0 = unsafe
     */
    void ClassifyBatch(
        const EdgeToClassify* edges,
        size_t num_edges,
        uint8_t* results
    );

    /**
     * Release all GPU resources.
     */
    void Destroy();

    bool IsInitialized() const { return initialized_; }

private:
    bool initialized_ = false;

    // Device pointers
    uint32_t* d_vlabels_ = nullptr;
    QueryEdgeTriple* d_query_triples_ = nullptr;
    EdgeToClassify* d_edges_ = nullptr;
    uint8_t* d_results_ = nullptr;

    // Pinned host buffers for async transfer
    EdgeToClassify* h_edges_pinned_ = nullptr;
    uint8_t* h_results_pinned_ = nullptr;

    // Capacities
    size_t vlabels_capacity_ = 0;
    size_t num_query_edges_ = 0;
    size_t edges_capacity_ = 0;  // current allocated capacity for edge batch
};

#endif // CORE_GPU_CLASSIFIER_H
