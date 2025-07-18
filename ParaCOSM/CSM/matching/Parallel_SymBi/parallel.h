#ifndef MATCHING_PARRALLEL_SYMBI
#define MATCHING_PARRALLEL_SYMBI

#include <queue>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "matching/matching.h"
#include "graph/Threadpool.cpp"

#include <tbb/concurrent_queue.h>
#include "graph/storage_hash_map.hpp"


class Parrllel_SymBi : public matching
{
private:
    struct TreeNode {
        std::vector<uint> forwards_;
        std::vector<uint> forward_labels_;
        std::vector<uint> backwards_;
        std::vector<uint> backward_labels_;
        std::vector<uint> neighbors_;
    };
    struct ExtendableVertex {
        uint E;
        uint matched_nbrs;
        uint u_min;

        ExtendableVertex()
        : E(NOT_EXIST), matched_nbrs(0u), u_min(NOT_EXIST) {}
        ExtendableVertex(uint E_arg, uint matched_nbrs_arg, uint u_min_arg)
        : E(E_arg), matched_nbrs(matched_nbrs_arg), u_min(u_min_arg) {}
    };

    struct SearchState {
        uint depth;
        std::vector<uint> m;
        std::vector<ExtendableVertex> extendable;
        std::vector<bool> visited;
    };

    std::vector<std::vector<uint>> eidx_;
    std::vector<TreeNode> treeNode_;
    uint q_root_;
    std::vector<uint> serialized_tree_;
    std::vector<std::vector<uint>> pre_defined_order_;
    std::vector<std::vector<uint>> pre_defined_backward_nbr_;

    size_t num_counter;

    // KV
    // #ifdef USE_UNORDERED_MAP
    // std::vector<std::unordered_map<uint, std::vector<uint>>> DCS_;

    // std::vector<std::unordered_map<uint, bool>> d1;
    // std::vector<std::unordered_map<uint, bool>> d2;

    // std::vector<std::unordered_map<uint, uint>> n1;
    // std::vector<std::unordered_map<uint, uint>> np1;

    // std::vector<std::unordered_map<uint, uint>> n2;
    // std::vector<std::unordered_map<uint, uint>> nc2;

    // #else
    std::vector<ska::flat_hash_map<uint, std::vector<uint>>> DCS_;

    std::vector<ska::flat_hash_map<uint, bool>> d1;
    std::vector<ska::flat_hash_map<uint, bool>> d2;

    std::vector<ska::flat_hash_map<uint, uint>> n1;
    std::vector<ska::flat_hash_map<uint, uint>> np1;

    std::vector<ska::flat_hash_map<uint, uint>> n2;
    std::vector<ska::flat_hash_map<uint, uint>> nc2;
    // #endif
    
    // 用于存储需要进行 InsertionTopDown 操作的节点对。
    std::queue<std::pair<uint, uint>> Q1; // queue 是比较难并行的东西.  only TopDown Method
    // 用于存储需要进行 InsertionBottomUp 操作的节点对。当某个节点对 (u, v) 满足一定条件时，会将其加入 Q2。
    // 在 Q2 非空时，会从队列中取出节点对 (u_queue, v_queue)，并对其进行 InsertionBottomUp 操作。
    std::queue<std::pair<uint, uint>> Q2;

    ThreadPool thread_pool;

    std::vector<std::vector<uint>> local_vec_m;

    std::vector<std::vector<Parrllel_SymBi::ExtendableVertex>> local_vec_extendable;
    std::vector<std::vector<bool>> local_vec_visited_local;

    size_t NUM_THREAD;

    size_t BIG_THREAD;

    std::vector< std::tuple<uint, uint, size_t, std::vector<uint>,
        std::vector<ExtendableVertex>,  uint , uint> > vertex_vector;
    
    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>,
    std::vector<ExtendableVertex>,  uint , uint> > job_queue;

    size_t NUMTHREAD;
    size_t auto_tuning;

    // 并行idea：把 Q1 做成std::vector<std::queue<std::pair<uint, uint>>>，每个线程处理一个队列。

public:

    Parrllel_SymBi(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, 
            std::vector<std::vector<uint>> orders, size_t NUMTHREAD, size_t auto_tuning);
    ~Parrllel_SymBi() override {};

    void Preprocessing() override;
    void InitialMatching() override;
    
    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

    void Update_Graph(Graph& data_, int v1, int v2, int label);
    void UpdateAndFind(uint v1, uint v2, uint label);
    void AddEdgeAsync(uint v1, uint v2, uint label);

private:
    void BuildDAG();
    void BuildDCS();
    
    void InsertionTopDown(uint u, uint u_c, uint v, uint v_c);
    void InsertionBottomUp(uint u, uint u_p, uint v, uint v_p);
    void DeletionTopDown(uint u, uint u_c, uint v, uint v_c);
    void DeletionBottomUp(uint u, uint u_p, uint v, uint v_p);

    bool Classify(uint v1, uint v2, uint label) override;

    void InsertionTopDown_para(uint u, uint u_c,  uint v_c, 
    std::queue<std::pair<uint, uint>> &Q1_para, std::queue<std::pair<uint, uint>> &Q2_para);
    void InsertionBottomUp_para(uint u, uint u_p, uint v_p, std::queue<std::pair<uint, uint>> &Q2_para);

    void FindMatches(uint depth, std::vector<uint>& m, 
            std::vector<ExtendableVertex>& extendable, size_t &num_results);
    void Parallel_FindMatches(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local);

    void Parallel_FindMatches_stack(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local);

    void FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);
    void Parallel_FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

    // void Parallel_FindMatches_local(uint depth, std::vector<uint>& m, 
    //     std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local);

    void Parallel_FindMatches_local2(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local);

    void Parallel_FindMatches_local(uint initial_depth, std::vector<uint> initial_m,
        std::vector<ExtendableVertex> initial_extendable, size_t &num_results, std::vector<bool> initial_visited);

    void Parallel_FindMatches_local_stack(uint initial_depth, std::vector<uint> initial_m,
        std::vector<ExtendableVertex> initial_extendable, size_t &num_results, std::vector<bool> initial_visited);

    void FindMatches2(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);
    
    bool process_vertex(uint u, uint u_min, size_t v_idx, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, 
        std::vector<bool> &visited_local, uint depth);
    bool process_vertex(uint u, uint u_min, size_t v_idx, std::vector<uint>& m, 
            std::vector<ExtendableVertex>& extendable, size_t &num_results, 
             uint depth);

    void Parallel_FindMatches_local2(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);

    void Parallel_FindMatches3(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);

    void Parallel_FindMatches4(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);

    bool process_vertex2(
        uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, 
        size_t &num_results, 
        uint depth) ;

    inline bool process_vertex_more(
        uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, 
        size_t &num_results, 
        uint depth,
        size_t thread_id
    ) ;

    inline bool process_vertex_layer1(
        uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, 
        size_t &num_results, 
         uint depth, size_t thread_id);

    inline void Parallel_FindMatches_local3(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, size_t thread_id);
             
    void AddEdge_Single(uint v1, uint v2, uint label);

    void Parallel_FindMatches(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);

    void FindMatches_test(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results);


    void processVertexK(
        uint u, uint u_min, uint v, 
        std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, 
        size_t& num_results, 
        uint depth
        // const std::vector<std::vector<uint>>& n2, 
        // const std::vector<std::vector<bool>>& visited_, 
        // const std::vector<TreeNode>& treeNode_, 
        // const std::vector<std::vector<std::unordered_map<uint, std::vector<uint>>>>& DCS_, 
        // const std::vector<std::vector<uint>>& eidx_, 
        // bool homomorphism_
    );

    void FindMatches_NEW(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results,size_t thread_id);

    inline bool process_vertex_visit(
            uint u, uint u_min, size_t v_idx, 
            std::vector<uint>& m, 
            std::vector<ExtendableVertex>& extendable, 
            size_t &num_results, 
             uint depth, std::vector<bool> &visited_local) ;

    void Parallel_FindMatches_local_MMM(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local);

        inline bool FindMatches_task(
            uint u, uint u_min, size_t v_idx, 
            std::vector<uint>& m, 
            std::vector<ExtendableVertex>& extendable, 
            size_t &num_results, 
             uint depth, size_t thread_id) ;

    void Parallel_FindMatches_delete(uint depth, std::vector<uint>& m, 
                std::vector<ExtendableVertex>& extendable, size_t &num_results);

                inline bool process_vertex_layer_local(
                    uint u, uint u_min, size_t v_idx, 
                    std::vector<uint>& m, 
                    std::vector<ExtendableVertex>& extendable, 
                    size_t &num_results, 
                     uint depth, size_t thread_id) ;

};
#endif //MATCHING_PARRALLEL_SYMBI
