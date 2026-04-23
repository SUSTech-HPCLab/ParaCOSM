#ifndef CORE_GPU_SEARCH_H
#define CORE_GPU_SEARCH_H

#include <cstddef>
#include <cstdint>
#include <vector>
#include <tuple>

class Graph;  // forward

class GPUSearchEngine {
public:
    GPUSearchEngine();
    ~GPUSearchEngine();

    GPUSearchEngine(const GPUSearchEngine&) = delete;
    GPUSearchEngine& operator=(const GPUSearchEngine&) = delete;

    /**
     * Build CSR on GPU from current data graph state.
     * Also copies vertex labels.
     */
    void BuildCSR(const Graph& data);

    /**
     * Set up query graph in constant memory.
     * @param query  The query graph
     */
    void SetupQuery(const Graph& query);

    /**
     * Set the matching order for a given order_index.
     * Must be called before SearchFromVertexVector.
     */
    void SetMatchingOrder(const uint32_t* order, uint32_t num_vertices);

    /**
     * Set all matching orders at once (for inter-update batch search).
     * @param all_orders  Flattened array [num_edges * Q]
     * @param num_edges   Number of query edges (= number of matching orders)
     * @param num_vertices Number of query vertices
     */
    void SetAllMatchingOrders(const uint32_t* all_orders, uint32_t num_edges, uint32_t num_vertices);

    /**
     * Run GPU DFS search from vertex_vector entries (inner-update).
     */
    uint64_t SearchFromVertexVector(
        const std::vector<std::tuple<uint32_t, uint32_t, uint32_t,
                                      std::vector<uint32_t>, uint32_t>>& vv,
        uint32_t depth,
        uint32_t num_query_vertices
    );

    /**
     * Run GPU DFS for a batch of edges (inter-update parallelism).
     * Each edge generates up to 2*num_query_edges search tasks (two directions).
     * The graph CSR must be built before calling this.
     *
     * @param edges  Array of (v1, v2, edge_label)
     * @param num_edges_data  Number of data edges to enumerate
     * @param num_query_edges Number of query edges (matching orders)
     * @param num_query_vertices Number of query vertices
     * @return Total number of matches found across all edges
     */
    uint64_t SearchBatchEdges(
        const uint32_t* edges_v1,
        const uint32_t* edges_v2,
        const uint32_t* edges_label,
        uint32_t num_edges_data,
        uint32_t num_query_edges,
        uint32_t num_query_vertices
    );

    void Destroy();
    bool IsInitialized() const { return csr_built_; }

    /**
     * Incrementally update GPU CSR for two vertices after an AddEdge.
     * Only rebuilds full CSR when padding is exhausted.
     * @return true if incremental update succeeded, false if full rebuild needed
     */
    bool IncrementalUpdate(const Graph& data, uint32_t v1, uint32_t v2);

private:
    bool csr_built_ = false;
    bool query_set_ = false;

    // Device CSR graph (padded: each vertex has capacity >= actual degree)
    uint32_t* d_csr_offsets_ = nullptr;    // [V+1] — offsets into padded arrays
    uint32_t* d_csr_neighbors_ = nullptr;  // [total_padded]
    uint32_t* d_csr_elabels_ = nullptr;    // [total_padded]
    uint32_t* d_vlabels_ = nullptr;        // [V]
    uint32_t* d_degrees_ = nullptr;        // [V] — actual degree (< padded capacity)
    uint32_t num_vertices_ = 0;
    uint32_t total_csr_edges_ = 0;
    std::vector<uint32_t> h_offsets_;      // host mirror of offsets for incremental update
    std::vector<uint32_t> h_capacities_;   // per-vertex capacity in padded CSR

    // Device vertex_vector buffers (grown as needed)
    uint32_t* d_vv_v_ = nullptr;
    uint32_t* d_vv_umin_ = nullptr;
    uint32_t* d_vv_ulabel_ = nullptr;
    uint32_t* d_vv_i_ = nullptr;
    uint32_t* d_vv_m_ = nullptr;       // flattened [N * Q]
    uint64_t* d_total_count_ = nullptr; // single uint64 for result
    size_t vv_capacity_ = 0;

    // Device buffers for inter-update batch search
    uint32_t* d_batch_v1_ = nullptr;
    uint32_t* d_batch_v2_ = nullptr;
    uint32_t* d_batch_label_ = nullptr;
    uint32_t* d_all_orders_ = nullptr;  // [num_query_edges * Q]
    size_t batch_capacity_ = 0;
    uint32_t all_orders_num_edges_ = 0;

    void EnsureVVCapacity(size_t n, uint32_t q);
    void EnsureBatchCapacity(size_t n);
};

#endif // CORE_GPU_SEARCH_H
