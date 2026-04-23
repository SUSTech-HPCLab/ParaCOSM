#ifndef MATCHING_PARALLEL_GRAPHFLOW
#define MATCHING_PARALLEL_GRAPHFLOW

#include <atomic>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include <taskflow/taskflow.hpp>

#include "utils/types.h"
#include "graph_storage/graph.h"
#include "matching_executor/matching.h"
#include "core/gpu/gpu_candidate_filter.h"
#include "core/gpu/gpu_search.h"



class Parallel_Graphflow : public matching
{
private:
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> order_csrs_;
    std::vector<std::vector<uint>> order_offs_;

    std::vector<std::vector<uint>> local_vec_m;
    std::vector<std::vector<bool>> local_vec_visited_local;

    // Legacy vertex_vector (used by non-optimized paths)
    std::vector< std::tuple<
    uint,                // v
    uint,                // u_min
    uint,               // u_min_label
    std::vector<uint>,   // m
    uint                 // i
    > > vertex_vector;

    tbb::concurrent_queue< std::tuple<uint, uint, uint,  std::vector<uint>,  uint > > job_queue;

    // ---- Optimized Layer1Group structure (avoids per-entry m copy) ----
    struct Layer1Group {
        uint v;              // layer-1 data vertex mapped to u
        uint u;              // layer-1 query vertex
        uint u_min2;         // layer-2 u_min
        uint u_min_label2;   // layer-2 u_min_label
        size_t start;        // start index in layer2_indices_
        size_t count;        // number of layer-2 entries
    };
    std::vector<Layer1Group> layer1_groups_;
    std::vector<uint> layer2_indices_;  // flat array of layer-2 neighbor indices

    // Thread-local pre-allocated m vectors (avoids heap allocation per task)
    std::vector<std::vector<uint>> tl_m_;  // [max_threads][query_vertices]

    // Precomputed u_min info per (order_index, depth) to avoid re-deriving in FindMatches_local.
    // precomp_u_min_[order_index][depth] = (u_min_query_neighbor_index, u_min_query_vertex)
    // Computed during GenerateMatchingOrder; NOT valid when mapping changes candidate sizes.
    // For Graphflow u_min depends on runtime m[] mapping, so we cache query topology only.
    // order_nbr_of_[order_index][depth] = list of query neighbors of order_vs_[order_index][depth]
    //   that appear BEFORE depth in the order (i.e., already matched)
    struct DepthNbrInfo {
        uint u;           // query vertex at this depth
        std::vector<uint> matched_nbrs;       // q_nbrs indices where m[q_nbrs[j]] != UNMATCHED
        std::vector<uint> matched_nbr_labels; // corresponding edge labels
    };
    std::vector<std::vector<DepthNbrInfo>> order_depth_info_; // [order_index][depth]

    struct StealWork {
        std::vector<uint> m;
        std::vector<uint> ancestors;
        uint start_depth;
        uint order_index;
    };
    tbb::concurrent_queue<StealWork> steal_queue_;

    mutable std::mutex enum_result_mutex_;

    std::unique_ptr<tf::Executor> persistent_executor_;

    // ---- Persistent thread pool (avoids OMP fork/join per update) ----
    std::vector<std::thread> pool_workers_;
    std::atomic<bool> pool_shutdown_{false};
    std::atomic<uint64_t> pool_epoch_{0};
    std::atomic<size_t> pool_next_item_{0};
    std::atomic<size_t> pool_items_done_{0};
    size_t pool_total_items_{0};
    std::function<void(size_t, size_t)> pool_work_fn_;

    void pool_worker_loop_(size_t tid);
    void pool_dispatch_(size_t count, std::function<void(size_t, size_t)> fn);

    size_t NUMTHREAD;
    size_t auto_tuning;

    // GPU candidate filter for Layer 1 acceleration
    GPUCandidateFilter gpu_candidate_filter_;
    bool gpu_filter_initialized_ = false;
    static constexpr size_t GPU_FILTER_THRESHOLD = 1024;  // min candidates to use GPU

    // GPU DFS search engine for deep recursion acceleration
    GPUSearchEngine gpu_search_engine_;
    bool gpu_search_initialized_ = false;
    bool gpu_search_csr_dirty_ = true;
    std::vector<std::pair<uint32_t, uint32_t>> gpu_search_dirty_edges_;
    static constexpr size_t GPU_DFS_THRESHOLD = 2000;  // min vertex_vector size to use GPU DFS


public:
    Parallel_Graphflow(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo,  size_t NUMTHREAD, size_t auto_tuning);
    ~Parallel_Graphflow() override {
        pool_shutdown_.store(true, std::memory_order_release);
        pool_epoch_.fetch_add(1, std::memory_order_release);
        for (auto& w : pool_workers_) w.join();
    };

    void Preprocessing() override;
    void InitialMatching() override;

    bool Classify(uint v1, uint v2, uint label) override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;

    const std::vector<std::vector<uint>>& GetMatchingOrders() const override {
        return order_vs_;
    }

    /**
     * @brief Enumerate matches for an edge already in the graph (thread-safe).
     * Uses FindMatches_local with local_vec_visited_local[thread_id].
     */
    size_t EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id) override;

    /**
     * @brief Resize per-thread state for batch enumeration.
     */
    void PrepareBatchEnumeration(size_t num_threads) override;

    // Single OMP parallel region for entire update stream
    void PersistentParallelUpdate(
        Graph& data_graph,
        size_t& num_v_updates, size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last);

    void AddEdgeWithSubflow(uint v1, uint v2, uint label, tf::Subflow& sf);
    void RemoveEdgeWithSubflow(uint v1, uint v2, tf::Subflow& sf);
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void GenerateMatchingOrder();
    void FindMatches(uint order_index, uint depth,
            std::vector<uint> m, size_t &num_results);

        void taskflow_findmatches(uint order_index, uint depth,
            std::vector<uint> m, size_t &num_results);

        // Taskflow-based matching where candidates at a layer are
        // partitioned into groups, and each task is responsible for
        // processing a group of candidates (possibly exploring multiple
        // deeper layers sequentially within the same task).
        void taskflow_findmatches_layer(
            uint order_index,
            uint depth,
            std::vector<uint> m,
            size_t &num_results);

        void taskflow_findmatches_subflow(uint order_index, uint depth,
            std::vector<uint> m, std::vector<bool> visited,
            size_t &num_results, tf::Subflow& sf);

        void FindMatches_taskflow_local(uint order_index, uint depth,
            std::vector<uint> &m, std::vector<bool> &visited_local,
            size_t &num_results, size_t local_limit);

    // ParaCOSM kernel entry (defined in core/FindMatchesKernel.cpp)
    void Parallel_Graphflow_FindMatches_ParaCOSM_Kernel(
        uint order_index,
        uint depth,
        std::vector<uint> m,
        size_t &num_results);

    void Parallel_FindMatches(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);
    inline bool ProcessCandidate(
        uint u, uint v, uint u_min, uint u_min_label, uint i,
        const std::vector<uint>& q_nbrs, 
        const std::vector<uint>& q_nbr_labels,
        const std::vector<uint>& u_min_nbr_labels,
        std::vector<uint>& m, 
        uint depth, uint order_index,
        size_t& num_results);

    void FindMatches_pure(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);

    inline bool ProcessNeighbor(
        // uint v,                       
        uint u,                       
        uint u_min,                   
        uint u_min_label,             
        uint order_index,             
        uint depth,                    
        std::vector<uint>& m,               
        size_t& num_results,       
        // const std::vector<uint>& u_min_nbr_labels, 
        uint i    
        , size_t thread_id                    
    );

    inline bool ProcessNeighbor_local(
        // uint v,                      
        uint u,                        
        uint u_min,                    
        uint u_min_label,             
        uint order_index,             
        uint depth,                   
        std::vector<uint>& m,               
        size_t& num_results,        
        // const std::vector<uint>& u_min_nbr_labels,
        uint i                        
        // , size_t thread_id
    );

    bool ProcessNeighbor_queue(
        // uint v,                       
        uint u,                     
        uint u_min,                    
        uint u_min_label,              
        uint order_index,            
        uint depth,                   
        std::vector<uint>& m,             
        size_t& num_results,        
        // const std::vector<uint>& u_min_nbr_labels, 
        uint i                        
        , size_t thread_id
    );

    void FindMatches_local_m(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    void FindMatches_local(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results, size_t thread_id);

    // Work-splitting version: splits large subtrees into steal_queue_
    void FindMatches_local_splitting(uint order_index, uint depth,
        std::vector<uint>& m, std::vector<uint>& ancestors,
        size_t &num_results, size_t thread_id);

    void Parallel_FindMatches2(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);

    // Inner version: uses #pragma omp for (works inside existing parallel region)
    void Parallel_FindMatches2_inner(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    void Process_vertex_queue(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);


};

#endif //MATCHING_GRAPHFLOW
