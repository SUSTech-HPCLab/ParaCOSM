#include <unordered_set>
#include <vector>

#include "utils/types.h"
#include "graph_storage/graph.h"
#include "matching_executor/matching.h"

matching::matching(Graph& query_graph, Graph& data_graph,
        size_t max_num_results, 
        bool print_prep,
        bool print_enum, 
        bool homo)
: query_(query_graph)
, data_(data_graph)

, max_num_results_(max_num_results)
, print_preprocessing_results_(print_prep)
, print_enumeration_results_(print_enum)
, homomorphism_(homo)

, visited_(data_.NumVertices(), false)
, num_initial_results_(0ul)
, num_positive_results_(0ul)
, num_negative_results_(0ul)
, num_intermediate_results_before_index_check_(0ul)
, num_intermediate_results_after_index_check_(0ul)
, num_intermediate_results_after_joinability_check_(0ul)
, num_intermediate_results_after_visit_check_(0ul)
, num_intermediate_results_with_empty_candidate_set_(0ul)
, num_intermediate_results_without_results_(0ul)
// , NUMTHREAD
{}

void matching::Preprocessing()
{}

void matching::InitialMatching()
{}

void matching::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);
}

void matching::RemoveEdge(uint v1, uint v2)
{
    data_.RemoveEdge(v1, v2);
}

void matching::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
}

void matching::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void matching::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}

bool matching::Classify(uint v1, uint v2, uint label){
    return false;
}

size_t matching::EnumerateNewEdge(uint v1, uint v2, uint label, size_t /*thread_id*/)
{
    size_t before = num_positive_results_;
    AddEdge(v1, v2, label);
    return num_positive_results_ - before;
}

size_t matching::EnumerateNewEdgeVersioned(uint v1, uint v2, uint label,
                                            size_t thread_id, uint /*max_timestamp*/)
{
    return EnumerateNewEdge(v1, v2, label, thread_id);
}

void matching::PrepareBatchEnumeration(size_t /*num_threads*/)
{
    // Default: nothing to prepare
}

void matching::UpdateIndexForEdge(uint /*v1*/, uint /*v2*/, uint /*label*/)
{
    // Default: no-op (algorithms without auxiliary indices)
}

void matching::GetNumInitialResults(size_t &num_initial_results)
{
    num_initial_results = num_initial_results_;
}

void matching::GetNumPositiveResults(size_t &num_positive_results)
{
    num_positive_results = num_positive_results_;
}

void matching::GetNumNegativeResults(size_t &num_negative_results)
{
    num_negative_results = num_negative_results_;
}

void matching::PrintCounter()
{
    // std::cout << num_intermediate_results_before_index_check_ << " intermediate results before index check.\n";
    // std::cout << num_intermediate_results_after_index_check_ << " intermediate results after index check.\n";
    // std::cout << num_intermediate_results_after_joinability_check_ << " intermediate results after joinability check.\n";
    // std::cout << num_intermediate_results_after_visit_check_ << " intermediate results after visit check.\n";
    // std::cout << num_intermediate_results_with_empty_candidate_set_ << " intermediate results with empty candidate set.\n";
    // std::cout << num_intermediate_results_without_results_ << " intermediate results without results.\n";
}

// bool 