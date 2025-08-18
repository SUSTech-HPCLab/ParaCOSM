#ifndef CSM_MATCHING_H
#define CSM_MATCHING_H

#include <vector>

#include "utils/types.h"
#include "graph/graph.h"

enum searchType{
    init, pos, neg
};

class matching
{
protected:
    Graph& query_;
    Graph& data_;
    std::vector<Graph>queryVec; // for different graph and different query

    // config
    const size_t max_num_results_;
    const bool print_preprocessing_results_;
    const bool print_enumeration_results_;
    const bool homomorphism_;
public:
    // execution info
    std::vector<bool> visited_;
    std::vector<std::vector<bool>> visited_parallel;

    size_t num_initial_results_;
    size_t num_positive_results_;
    size_t num_negative_results_;
    size_t num_intermediate_results_before_index_check_;
    size_t num_intermediate_results_after_index_check_;
    size_t num_intermediate_results_after_joinability_check_;
    size_t num_intermediate_results_after_visit_check_;
    size_t num_intermediate_results_with_empty_candidate_set_;
    size_t num_intermediate_results_without_results_;

    Timer indexCheckTime, searchVertexTime, searchInitTime, indexupdateTime, indexBuildTime, matchOrderBuildTime, findQueryGraphTime, DescListTime, WCTime, LRTime;

public:
    matching(std::vector<Graph> & queryVec,
        Graph& query_graph, Graph& data_graph,
        size_t max_num_results = ULONG_MAX, 
        bool print_preprocessing_results = true,
        bool print_enumeration_results = true, 
        bool homomorphism = false
        );
    virtual ~matching() = default;

    virtual void Preprocessing();
    virtual void InitialMatching();

    virtual void AddEdge(uint v1, uint v2, uint label);
    virtual void RemoveEdge(uint v1, uint v2);
    virtual void AddVertex(uint id, uint label);
    virtual void RemoveVertex(uint id);

    virtual bool safe_detect(uint v1, uint v2, uint label, searchType type);
    virtual void Safe_Update(uint v1, uint v2, uint label);

    virtual void Safe_Update_remove(uint v1, uint v2);
    
    virtual void GetMemoryCost(size_t &num_edges, size_t &num_vertices);
    void TimePrint(bool motif);
    // get execution info
    void GetNumInitialResults(size_t &num_initial_results);
    void GetNumPositiveResults(size_t &num_positive_results);
    void GetNumNegativeResults(size_t &num_negative_results);

    void PrintCounter();
};

#endif //CSM_MATCHING_H
