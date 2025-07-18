#ifndef MATCHING_Parallel_TURBOFLUX
#define MATCHING_Parallel_TURBOFLUX

#include <queue>
#include <unordered_map>
#include <vector>

#include "graph/graph.h"
#include "matching/matching.h"

#include <tbb/concurrent_queue.h>
#include "graph/storage_hash_map.hpp"

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

    size_t NUMTHREAD;
    size_t auto_tuning;
    
    std::queue<std::pair<uint, uint>> Q1;
    std::queue<std::pair<uint, uint>> Q2;

    std::vector< std::tuple<uint, uint, size_t, std::vector<uint>,
      uint> > vertex_vector;

    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>,
    uint , uint> > job_queue;

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
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void BuildDAG();
    void BuildDCS();
    void GenerateMatchingOrder();
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

    inline bool ProcessVertex_queue(uint u, uint u_min, size_t v_idx, 
        std::vector<uint>& m, 
        size_t &num_results, 
         uint depth, size_t order_index, size_t thread_id);

    inline void Parallel_FindMatches_delete(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results);

};
#endif //MATCHING_TURBOFLUX
