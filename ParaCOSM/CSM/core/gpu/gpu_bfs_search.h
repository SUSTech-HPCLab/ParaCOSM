#ifndef CORE_GPU_BFS_SEARCH_H
#define CORE_GPU_BFS_SEARCH_H

#include <cstddef>
#include <cstdint>
#include <vector>

class Graph;

class GPUBFSSearch {
public:
    GPUBFSSearch();
    ~GPUBFSSearch();

    GPUBFSSearch(const GPUBFSSearch&) = delete;
    GPUBFSSearch& operator=(const GPUBFSSearch&) = delete;

    // Compact CSR (for batch/inter-update mode)
    void BuildCSR(const Graph& data);

    // Compact CSR with per-edge timestamps (for versioned mode)
    void BuildCSR_Versioned(const Graph& data);

    // Padded CSR (for single/inner-update mode, supports incremental updates)
    void BuildPaddedCSR(const Graph& data);
    bool IncrementalUpdate(const Graph& data, uint32_t v1, uint32_t v2);
    void IncrementalUpdateVertex(uint32_t v, const uint32_t* neighbors, const uint32_t* elabels, uint32_t degree);

    void SetupQuery(const Graph& query);
    void SetAllMatchingOrders(const uint32_t* all_orders, uint32_t num_edges, uint32_t num_vertices);

    /**
     * BFS-level search for a batch of unsafe edges (inter-update).
     */
    uint64_t SearchBatchEdgesBFS(
        const uint32_t* edges_v1,
        const uint32_t* edges_v2,
        const uint32_t* edges_label,
        uint32_t num_edges_data,
        uint32_t num_query_edges,
        uint32_t num_query_vertices
    );

    /**
     * BFS-level search for a single unsafe edge (inner-update).
     * Uses padded CSR with degrees array.
     */
    uint64_t SearchSingleEdgeBFS(
        uint32_t v1, uint32_t v2, uint32_t label,
        uint32_t num_query_edges,
        uint32_t num_query_vertices
    );

    /**
     * BFS-level versioned search: inner-update semantics with batch parallelism.
     * Each unsafe edge has a max_timestamp; only CSR edges with ts <= max_ts are visible.
     */
    uint64_t SearchBatchEdgesBFS_Versioned(
        const uint32_t* edges_v1,
        const uint32_t* edges_v2,
        const uint32_t* edges_label,
        const uint32_t* edges_max_ts,   // per-edge timestamp cap
        uint32_t num_edges_data,
        uint32_t num_query_edges,
        uint32_t num_query_vertices
    );

    void Destroy();
    bool IsInitialized() const { return csr_built_; }
    bool IsPadded() const { return padded_csr_; }

    // Public for lazy flush access
    uint32_t num_vertices_ = 0;
    std::vector<uint32_t> h_offsets_;
    std::vector<uint32_t> h_capacities_;

private:
    bool csr_built_ = false;
    bool query_set_ = false;
    bool padded_csr_ = false;  // true if using padded CSR

    // Device CSR
    uint32_t* d_csr_offsets_ = nullptr;
    uint32_t* d_csr_neighbors_ = nullptr;
    uint32_t* d_csr_elabels_ = nullptr;
    uint32_t* d_csr_timestamps_ = nullptr;  // per-edge timestamp (versioned mode)
    uint32_t* d_vlabels_ = nullptr;
    uint32_t* d_degrees_ = nullptr;  // actual degrees (for padded CSR)

    // Query graph in constant memory (set via SetupQuery)
    // Matching orders in global memory
    uint32_t* d_all_orders_ = nullptr;
    uint32_t all_orders_num_edges_ = 0;
    uint32_t num_query_vertices_ = 0;

    // BFS ping-pong buffers for partial matches
    // Each partial match: [order_idx, m[0], ..., m[Q-1]] = (Q+1) uint32
    uint32_t* d_buf_a_ = nullptr;
    uint32_t* d_buf_b_ = nullptr;
    uint32_t* d_buf_c_ = nullptr;  // third buffer for recursive overflow flush
    static constexpr size_t MAX_BUF_MATCHES = 400'000'000;

    // Counters
    uint32_t* d_count_ = nullptr;
    uint64_t* d_result_ = nullptr;

    // Edge input buffers
    uint32_t* d_edges_v1_ = nullptr;
    uint32_t* d_edges_v2_ = nullptr;
    uint32_t* d_edges_label_ = nullptr;
    uint32_t* d_edges_max_ts_ = nullptr;  // per-unsafe-edge timestamp cap (versioned)
    size_t edges_capacity_ = 0;

    void EnsureEdgesCapacity(size_t n);
    void EnsureBufCapacity(uint32_t q);

    // Recursive BFS from a given depth with overflow handling
    void BFSFromDepth(uint32_t* in_buf, uint32_t in_count,
                      uint32_t* scratch_buf, uint32_t start_depth,
                      uint32_t Q, uint32_t stride, bool versioned);
};

#endif
