#ifndef MATCHING_PARALLEL_CALIG_H
#define MATCHING_PARALLEL_CALIG_H

#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <atomic>

#include "graph_storage/graph.h"
#include "matching_executor/matching.h"
#include <tbb/concurrent_queue.h>

class Parallel_CaLiG : public matching
{
public:
    // CaLiG-specific types
    using u_set = std::unordered_set<uint>;
    using vec = std::vector<uint>;

    // Core-shell matching order decomposition
    struct Split {
        vec core;
        std::vector<u_set> core_nei;
        std::unordered_map<uint, u_set> c_s_nei;  // core vertex → shell indices
        vec shell;
        std::vector<u_set> shell_nei;
    };

    // CaLiG candidate index per data vertex
    struct CaLiGIndex {
        // cand[ui][uj] = set of data vertices matching (ui, uj) relationship
        std::unordered_map<uint, std::unordered_map<uint, u_set>> cand;
        // LI[ui] = is this vertex a live candidate for query vertex ui
        std::unordered_map<uint, bool> LI;
    };

private:
    // Query graph analysis
    std::unordered_map<uint, vec> labels_;  // label → query vertices with that label
    // rep_nei_[ui][label] = query neighbors of ui with that label (only if count > 1)
    std::vector<std::unordered_map<uint, vec>> rep_nei_;

    // Matching orders: matching_order_[ui][uj] = Split for edge (ui, uj)
    std::unordered_map<uint, std::unordered_map<uint, Split>> matching_order_;

    // CaLiG index: one CaLiGIndex per data vertex
    std::vector<CaLiGIndex> cidx_;

    // Parallelization
    size_t NUMTHREAD_;
    size_t auto_tuning_;

    // Versioned search: slot pool
    std::vector<std::unordered_map<uint, uint>> slot_m_;
    std::vector<vec> slot_used_;
    tbb::concurrent_queue<size_t> free_slots_;
    size_t slot_pool_size_ = 0;

    // Versioned matching orders for GPU
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> backward_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_labels_;
    std::vector<std::vector<uint>> gpu_orders_;

    // DAG tree structure (built for versioned/GPU path)
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

public:
    Parallel_CaLiG(Graph& query_graph, Graph& data_graph, uint max_num_results,
                    bool print_prep, bool print_enum, bool homo,
                    size_t num_threads, size_t auto_tuning);
    ~Parallel_CaLiG() override {};

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
    // CaLiG index construction
    void AnalyzeQuery();
    void GenerateMatchingOrders();
    void ConstructCandidates();
    void StaticFilter();

    // CaLiG index maintenance
    void TurnOff(uint vi);
    void TurnOnProcess(uint v1, uint v2);
    void TurnOffProcess(uint v1, uint v2);
    void AddAndCheck(uint ui, uint vi, vec& temp_v, vec& temp_u);
    void DeleteAndCheck(uint ui, uint vi);
    bool CheckNei(uint vi, uint ui);
    bool TryNei(uint th, uint vi, uint ui, u_set& used, vec& to_check);

    // CaLiG search
    int SearchCore(int th, std::unordered_map<uint, uint>& m, vec& used,
                   const vec& c, const std::vector<u_set>& c_n,
                   const vec& s, const std::vector<u_set>& s_n,
                   std::unordered_map<uint, u_set>& c2check);
    bool ShellCand(std::vector<u_set>& result, std::unordered_map<uint, uint>& m,
                   const vec& s, const std::vector<u_set>& s_n, vec& used);
    int NumAdd(int th, const std::vector<u_set>& cand, u_set& used);
    bool NotExit(uint shell_v, const u_set& nei, std::unordered_map<uint, uint>& m);

    // Safe-edge classification helpers
    bool TurnOnProcessSafe(uint v1, uint v2);
    bool TurnOffProcessSafe(uint v1, uint v2);

    // Versioned search internals
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

#endif // MATCHING_PARALLEL_CALIG_H
