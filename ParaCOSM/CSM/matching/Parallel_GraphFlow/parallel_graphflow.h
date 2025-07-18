#ifndef MATCHING_PARALLEL_GRAPHFLOW
#define MATCHING_PARALLEL_GRAPHFLOW

#include <vector>

#include "utils/types.h"
#include "graph/graph.h"
#include "matching/matching.h"



class Parallel_Graphflow : public matching
{
private:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint>> order_vs_;
    std::vector<std::vector<uint>> order_csrs_;
    std::vector<std::vector<uint>> order_offs_;

    std::vector<std::vector<uint>> local_vec_m;
    std::vector<std::vector<bool>> local_vec_visited_local;

    // std::vector< std::tuple<uint, uint, size_t, std::vector<uint>,
    // uint , uint> > vertex_vector;
    std::vector< std::tuple<
    uint,                // v 
    // uint,                // u 
    uint,                // u_min
    uint,               // u_min_label 
    // uint,                // order_index 
    // uint,                // depth 
    std::vector<uint>,   // m 
    // size_t&,             // num_results
    uint                 // i 
    > > vertex_vector;

    tbb::concurrent_queue< std::tuple<uint, uint, uint,  std::vector<uint>,  uint > > job_queue;

    size_t NUMTHREAD;
    size_t auto_tuning;


public:
    Parallel_Graphflow(Graph& query_graph, Graph& data_graph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo,  size_t NUMTHREAD, size_t auto_tuning);
    ~Parallel_Graphflow() override {};

    void Preprocessing() override;
    void InitialMatching() override;

    bool Classify(uint v1, uint v2, uint label) override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void GenerateMatchingOrder();
    void FindMatches(uint order_index, uint depth,
            std::vector<uint> m, size_t &num_results);

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

    void Parallel_FindMatches2(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);

    void Process_vertex_queue(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);


};

#endif //MATCHING_GRAPHFLOW
