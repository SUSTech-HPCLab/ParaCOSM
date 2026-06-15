#ifndef NEWSP_ADAPTER_H
#define NEWSP_ADAPTER_H

// Adapter: wraps NewSP's CSMPP (which uses its own Graph/matching base classes)
// into the main framework's matching interface.
//
// NewSP is compiled as a separate translation unit with its own namespace-isolated
// types. This adapter bridges the two type systems by:
//   1. Holding a newsp::CSMPP instance internally
//   2. Forwarding calls from framework's matching interface to newsp::CSMPP
//   3. Converting between framework's Graph and newsp's Graph data

#include "graph_storage/graph.h"
#include "matching_executor/matching.h"

#include <tbb/concurrent_queue.h>
#include <atomic>
#include <vector>

class Parallel_NewSP_Adapter : public matching
{
public:
    Parallel_NewSP_Adapter(Graph& query_graph, Graph& data_graph, uint max_num_results,
                           bool print_prep, bool print_enum, bool homo,
                           size_t num_threads, size_t auto_tuning);
    ~Parallel_NewSP_Adapter() override;

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;

    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    bool Classify(uint v1, uint v2, uint label) override;

    // Batch / versioned support
    size_t EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id) override;
    size_t EnumerateNewEdgeVersioned(uint v1, uint v2, uint label,
                                      size_t thread_id, uint max_timestamp) override;
    void PrepareBatchEnumeration(size_t num_threads) override;
    void UpdateIndexForEdge(uint v1, uint v2, uint label) override;

    const std::vector<std::vector<uint>>& GetMatchingOrders() const override {
        return gpu_orders_;
    }

private:
    struct Impl;
    Impl* impl_;  // pimpl: holds newsp::Graph, newsp::CSMPP, etc.

    size_t NUMTHREAD_;
    size_t auto_tuning_;

    // Versioned DAG + matching orders (same pattern as SymBi/CaLiG)
    struct TreeNode {
        std::vector<uint> forwards_;
        std::vector<uint> backwards_;
        std::vector<uint> neighbors_;
    };
    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_ = 0;
    std::vector<uint> serialized_tree_;
    std::vector<uint> tree_parent_;

    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> backward_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_labels_;
    std::vector<std::vector<uint>> gpu_orders_;

    // Slot pool
    std::vector<std::vector<uint>> slot_m_;
    std::vector<std::vector<bool>> slot_visited_;
    tbb::concurrent_queue<size_t> free_slots_;
    size_t slot_pool_size_ = 0;
    size_t hot_split_threshold_ = 64;
    size_t hot_chunk_size_ = 8;

    void BuildDAGForVersioned();
    void GenerateVersionedMatchingOrders();
    size_t AcquireSlot();
    void ReleaseSlot(size_t s);
    void FindMatches_versioned_v2(uint order_index, uint depth,
        std::vector<uint>& m, std::vector<bool>& visited,
        size_t& num_results, uint max_ts);
    void FindMatches_versioned_chunk(uint order_index, uint depth,
        std::vector<uint>& m, std::vector<bool>& visited,
        size_t& num_results, uint max_ts, size_t lo, size_t hi);
};

#endif // NEWSP_ADAPTER_H
