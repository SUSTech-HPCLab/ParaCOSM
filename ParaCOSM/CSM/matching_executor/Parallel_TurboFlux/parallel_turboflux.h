#ifndef MATCHING_Parallel_TURBOFLUX
#define MATCHING_Parallel_TURBOFLUX

#include <queue>
#include <unordered_map>
#include <vector>
#include <atomic>
#include <memory>

#include "graph_storage/graph.h"
#include "matching_executor/matching.h"

#include <tbb/concurrent_queue.h>
#include "graph_storage/storage_hash_map.hpp"
#include "core/inner_executor/FindMatchesKernel.h"

class Parallel_TurboFlux : public matching
{
private:
    struct TreeNode {
        std::vector<uint> forwards_;
        std::vector<uint> forward_labels_;
        std::vector<uint> backwards_;
        std::vector<uint> backward_labels_;
        std::vector<uint> neighbors_;
    };

    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_;
    std::vector<uint> serialized_tree_;

    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> backward_vs_;

    // GPU-friendly BFS matching orders, one per undirected query edge.
    // Generated independently from TF's tree-based plan (mirrors GraphFlow).
    std::vector<std::vector<uint>> gpu_orders_;

    std::atomic<int> active_workers_{0};

    std::vector<std::vector<std::vector<uint>>> join_check_vs_;
    std::vector<std::vector<std::vector<uint>>> join_check_labels_;

    // KV
    #ifdef USE_UNORDERED_MAP
    std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_;

    std::vector<std::unordered_map<uint, bool>> d1;
    std::vector<std::unordered_map<uint, bool>> d2;

    std::vector<std::unordered_map<uint, uint>> n1;
    std::vector<std::unordered_map<uint, uint>> np1;

    std::vector<std::unordered_map<uint, uint>> n2;
    std::vector<std::unordered_map<uint, uint>> nc2;

    #else
    std::vector<ska::flat_hash_map<uint, std::vector<uint>>> DCS_;

    std::vector<ska::flat_hash_map<uint, bool>> d1;
    std::vector<ska::flat_hash_map<uint, bool>> d2;

    std::vector<ska::flat_hash_map<uint, uint>> n1;
    std::vector<ska::flat_hash_map<uint, uint>> np1;

    std::vector<ska::flat_hash_map<uint, uint>> n2;
    std::vector<ska::flat_hash_map<uint, uint>> nc2;
    #endif

    std::vector<std::vector<uint>> local_vec_m;
    std::vector<std::vector<bool>> local_vec_visited_local;

    // ---- Hot-edge sub-task parallelism (slot pool, decoupled from tid) ----
    // Each slot holds a (m, visited) pair. visited starts all-false; tasks set
    // their path entries and clear them before releasing. Pool is sized larger
    // than num_threads so nested tasks (parent + children on same thread) don't
    // collide. Configurable via env var PARACOSM_HOT_SPLIT_THRESHOLD.
    std::vector<std::vector<uint>> slot_m_;
    std::vector<std::vector<bool>> slot_visited_;
    tbb::concurrent_queue<size_t> free_slots_;
    size_t slot_pool_size_ = 0;
    size_t hot_split_threshold_ = 64;
    size_t hot_chunk_size_ = 8;
    // Idea #2: deep-layer fork. When recursive DFS reaches an inner depth
    // with many neighbors while inside an OMP parallel region, fan out again.
    // Defaults are conservative to avoid task-overhead regressions on small
    // queries; env vars PARACOSM_DEEP_FORK_THRESHOLD / PARACOSM_MAX_FORK_DEPTH
    // / PARACOSM_DEEP_FORK_CHUNK can lower them for hot workloads.
    size_t deep_fork_threshold_ = 2048;
    size_t deep_fork_chunk_ = 64;
    uint   max_fork_depth_ = 3;

    size_t NUMTHREAD;
    size_t auto_tuning;
    
    std::queue<std::pair<uint, uint>> Q1;
    std::queue<std::pair<uint, uint>> Q2;

    // Lightweight work item: shares m snapshot via shared_ptr to avoid copies
    struct WorkItem {
        uint u;
        uint u_min;
        size_t v_idx;
        std::shared_ptr<std::vector<uint>> m_snapshot;
        uint parent_v;  // the v that was matched at the parallelization layer
    };
    std::vector<WorkItem> vertex_vector;

    // Job queue also uses shared_ptr to avoid deep-copying m
    struct StealJob {
        uint u;
        uint u_min;
        size_t v_idx;
        std::shared_ptr<std::vector<uint>> m_snapshot;
        uint depth;
        uint parent_v;
    };
    tbb::concurrent_queue<StealJob> job_queue;

    // Pre-allocated per-call result buffer
    std::vector<size_t> local_num_result_;

public:

    Parallel_TurboFlux(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, size_t NUMTHREAD, size_t auto_tuning);
    ~Parallel_TurboFlux() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    
    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;

    // batch_all support
    size_t EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id) override;
    size_t EnumerateNewEdgeVersioned(uint v1, uint v2, uint label,
                                      size_t thread_id, uint max_timestamp) override;
    void PrepareBatchEnumeration(size_t num_threads) override;
    void UpdateIndexForEdge(uint v1, uint v2, uint label) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

    // GPU pipeline: expose BFS matching orders (NumEdges entries, one per query edge).
    const std::vector<std::vector<uint>>& GetMatchingOrders() const override {
        return gpu_orders_;
    }

private:
    void BuildDAG();
    void BuildDCS();
    void GenerateMatchingOrder();
    void BuildGPUMatchingOrders();
    void CountDownwards(uint u, uint v, std::vector<uint>& num_explicit_pathes,
        std::vector<std::unordered_map<uint, uint>>& num_dp);
    
    void InsertionTopDown(uint u, uint u_c, uint v, uint v_c);
    void InsertionBottomUp(uint u, uint u_p, uint v, uint v_p);
    void DeletionTopDown(uint u, uint u_c, uint v, uint v_c);
    void DeletionBottomUp(uint u, uint u_p, uint v, uint v_p);

    void FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    bool Classify(uint v1, uint v2, uint label) override;

    void Parallel_FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    void ProcessVertex(size_t v_idx2, uint u2, uint depth, uint order_index, 
        std::vector<uint>& m, uint64_t& num_results);

    inline bool process_vertex_layer1(uint order_index,
        uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        // std::vector<ExtendableVertex>& extendable, 
        size_t &num_results, 
         uint depth, size_t thread_id);

         inline bool ProcessVertex(uint u, uint u_min, size_t v_idx, 
            std::vector<uint>& m, 
            size_t &num_results, 
             uint depth, size_t order_index, size_t thread_id);

    void FindMatches_local(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results, size_t thread_id);

    void FindMatches_versioned(uint order_index, uint depth, std::vector<uint>& m,
                                size_t &num_results, size_t thread_id, uint max_ts);

    // Slot-pool variants used by EnumerateNewEdgeVersioned for hot-edge fan-out.
    // _v2 takes m and visited explicitly (decoupled from thread_id); _chunk only
    // iterates [lo, hi) at the entry depth and recurses via _v2 for deeper levels.
    size_t AcquireSlot();
    void ReleaseSlot(size_t s);
    void FindMatches_versioned_v2(uint order_index, uint depth,
                                   std::vector<uint>& m,
                                   std::vector<bool>& visited,
                                   size_t& num_results, uint max_ts);
    void FindMatches_versioned_chunk(uint order_index, uint depth,
                                      std::vector<uint>& m,
                                      std::vector<bool>& visited,
                                      size_t& num_results, uint max_ts,
                                      size_t lo, size_t hi);

    inline bool ProcessVertex_queue(uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        size_t &num_results, 
         uint depth, size_t order_index, size_t thread_id);

    inline void Parallel_FindMatches_delete(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    // void Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);
    void Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);
    inline bool ProcessVertex_TurboFlux_ParaCOSM_Kernel(uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        size_t &num_results, 
         uint depth, size_t order_index, size_t thread_id);

};
#endif //MATCHING_TURBOFLUX
