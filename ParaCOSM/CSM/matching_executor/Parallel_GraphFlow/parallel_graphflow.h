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



class Parallel_Graphflow : public matching
{
private:
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> order_csrs_;
    std::vector<std::vector<uint>> order_offs_;

    std::vector<std::vector<uint>> local_vec_m;
    std::vector<std::vector<bool>> local_vec_visited_local;

    std::vector< std::tuple<
    uint,                // v 
    uint,                // u_min
    uint,               // u_min_label 
    std::vector<uint>,   // m 
    uint                 // i 
    > > vertex_vector;

    tbb::concurrent_queue< std::tuple<uint, uint, uint,  std::vector<uint>,  uint > > job_queue;

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

    inline bool ProcessNeighbor_queue(
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

    void FindMatches_local_m(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);

    void FindMatches_local(uint order_index, uint depth, std::vector<uint> m, size_t &num_results, size_t thread_id);

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
