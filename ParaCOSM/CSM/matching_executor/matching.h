#ifndef CSM_MATCHING_H
#define CSM_MATCHING_H

#include <vector>

#include "utils/types.h"
#include "graph_storage/graph.h"

class matching
{
protected:
    Graph& query_;
    Graph& data_;

    // std::vector<Graph>queryVec;

    // config
    const size_t max_num_results_;
    const bool print_preprocessing_results_;
    const bool print_enumeration_results_;
    const bool homomorphism_;

    // execution info
    std::vector<bool> visited_;
    size_t num_initial_results_;
    size_t num_positive_results_;
    size_t num_negative_results_;
    size_t num_intermediate_results_before_index_check_;
    size_t num_intermediate_results_after_index_check_;
    size_t num_intermediate_results_after_joinability_check_;
    size_t num_intermediate_results_after_visit_check_;
    size_t num_intermediate_results_with_empty_candidate_set_;
    size_t num_intermediate_results_without_results_;

    // Timer indexCheckTime, searchVertexTime, searchInitTime, indexupdateTime, indexBuildTime, matchOrderBuildTime, findQueryGraphTime, DescListTime, WCTime, LRTime;

public:
    matching(Graph& query_graph, Graph& data_graph,
        size_t max_num_results = ULONG_MAX, 
        bool print_preprocessing_results = true,
        bool print_enumeration_results = false, 
        bool homomorphism = false);
    virtual ~matching() = default;

    virtual void Preprocessing();
    virtual void InitialMatching();

    virtual void AddEdge(uint v1, uint v2, uint label);
    virtual void RemoveEdge(uint v1, uint v2);
    virtual void AddVertex(uint id, uint label);
    virtual void RemoveVertex(uint id);
    
    virtual void GetMemoryCost(size_t &num_edges, size_t &num_vertices);

    // void TimePrint(bool motif);

    virtual bool Classify(uint v1, uint v2, uint label);

    /**
     * @brief Enumerate new matches for an edge that is ALREADY in the data graph.
     *
     * Unlike AddEdge(), this does NOT call data_.AddEdge(). It only enumerates
     * subgraph matches that contain the given edge. Must be safe for concurrent
     * calls from different threads when each thread uses a distinct thread_id.
     *
     * @param v1    Source vertex of the edge (already in graph)
     * @param v2    Target vertex of the edge (already in graph)
     * @param label Edge label
     * @param thread_id  Caller's thread ID for accessing thread-local state
     * @return Number of new matches found
     */
    virtual size_t EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id);

    /**
     * @brief Prepare per-thread state for EnumerateNewEdge batch calls.
     * @param num_threads  Number of threads that will call EnumerateNewEdge concurrently
     */
    virtual void PrepareBatchEnumeration(size_t num_threads);

    // get execution info
    void GetNumInitialResults(size_t &num_initial_results);
    void GetNumPositiveResults(size_t &num_positive_results);
    void GetNumNegativeResults(size_t &num_negative_results);

    /// Thread-safe: atomically add to the positive results counter.
    void AddPositiveResults(size_t delta) { num_positive_results_ += delta; }

    void PrintCounter();
    
};

#endif //CSM_MATCHING_H
