#include <algorithm>
#include <iostream>
#include <vector>

#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/Parallel_GraphFlow/parallel_graphflow.h"

Parallel_Graphflow::Parallel_Graphflow(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,  size_t NUMTHREAD, size_t auto_tuning)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, order_vs_(query_.NumEdges())
, order_csrs_(query_.NumEdges())
, order_offs_(query_.NumEdges())
, NUMTHREAD(NUMTHREAD)
, auto_tuning(auto_tuning)
{

    auto BIG_THREAD =  NUMTHREAD + 2;

    vertex_vector.resize(0);

    if(query_.NumVertices() > 32){
        BIG_THREAD = query_.NumVertices() + 2;
    }else{
        BIG_THREAD = NUMTHREAD + 2;
    }


    for (uint i = 0; i < query_.NumEdges(); ++i)
    {
        order_vs_[i].resize(query_.NumVertices());
        order_csrs_[i].resize(query_.NumEdges() + 1);
        order_offs_[i].resize(query_.NumVertices(), 0);
    }

    // vertex_vector.resize(0);

    // std::cout << auto_tuning << std::endl;
    if(this->auto_tuning == 1){
        if(query_.NumVertices() < 7){
            if(NUMTHREAD > query_.NumVertices()){
                // std::cout << NUMTHREAD << std::endl;
                if(data_graph.NumEdges() > BIG_DATASET){
                    NUMTHREAD = query_.NumVertices();
                }else{
                    NUMTHREAD = query_.NumVertices() -2;
                }
                
                std::cout << "NUMTHREAD is set to " << NUMTHREAD << " for auto-tuning" << std::endl;
                std::cout << "use " << NUMTHREAD << " threads for incresemental matching " << std::endl;
            }
        }
    }

    local_vec_m = std::vector<std::vector<uint>>(BIG_THREAD, std::vector<uint>(query_.NumVertices())); 
    local_vec_visited_local = std::vector<std::vector<bool>>(BIG_THREAD, std::vector<bool>(data_.NumVertices(), false));

}

void Parallel_Graphflow::Preprocessing()
{
    GenerateMatchingOrder();
}


/**
 * @brief Generates optimal vertex matching orders for subgraph matching.
 * 
 * This function creates multiple vertex matching orders to optimize the subgraph enumeration process:
 * 1. An initial global matching order based on connectivity and degree
 * 2. Edge-specific matching orders for incremental matching during graph updates
 * 
 * @details
 * The function works in three main phases:
 * 
 * 1. Initial Order Generation:
 *    - Selects the highest-degree query vertex as the starting point
 *    - For subsequent positions, selects vertices with maximum connectivity to already ordered vertices
 *    - Uses degree as a tie-breaker to prioritize more constrained vertices first
 *    - Records both the vertex order and connectivity relationships in CSR format
 * 
 * 2. Incremental Order Generation:
 *    - Creates specialized matching orders for each edge in the query graph
 *    - Each order starts with the vertices forming the specific edge
 *    - Completes the order using the same connectivity-first strategy as the initial order
 *    - These orders enable efficient edge-triggered subgraph matching
 * 
 * 3. Recording Connection Information:
 *    - For each position in each order, tracks which previously matched vertices connect to it
 *    - Stores this information in compressed sparse row (CSR) format using:
 *      * order_vs_ : The sequence of vertices in the matching order
 *      * order_csrs_ : The connecting vertices for each position
 *      * order_offs_ : Offsets into the order_csrs_ array for each position
 * 
 * When preprocessing results output is enabled, the function displays the complete matching 
 * orders along with the connectivity information between vertices in the matching order.
 */
void Parallel_Graphflow::GenerateMatchingOrder()
{
    // generate the initial matching order, order_*s_[0]
    std::vector<bool> visited(query_.NumVertices(), false);
    uint max_degree = 0u;
    for (size_t i = 0; i < query_.NumVertices(); i++)
    {
        if (query_.GetDegree(i) > max_degree)
        {
            max_degree = query_.GetDegree(i);
            order_vs_[0][0] = i;
        }
    }
    visited[order_vs_[0][0]] = true;

    // loop over all remaining positions of the order
    for (uint i = 1; i < query_.NumVertices(); ++i)
    {
        uint max_adjacent = 0;
        uint max_adjacent_u = NOT_EXIST;
        for (size_t j = 0; j < query_.NumVertices(); j++)
        {
            uint cur_adjacent = 0u;
            if (visited[j]) continue;

            auto& q_nbrs = query_.GetNeighbors(j);
            for (auto& other : q_nbrs)
                if (visited[other])
                    cur_adjacent++;

            if (!cur_adjacent) continue;
            if (
                max_adjacent_u == NOT_EXIST ||
                (cur_adjacent == max_adjacent &&
                    query_.GetDegree(j) > query_.GetDegree(max_adjacent_u)) ||
                cur_adjacent > max_adjacent
            ) {
                max_adjacent = cur_adjacent;
                max_adjacent_u = j;
            }
        }
        order_vs_[0][i] = max_adjacent_u;
        visited[max_adjacent_u] = true;

        order_offs_[0][i] = order_offs_[0][i - 1];
        auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
        for (auto &other: q_nbrs)
            if (visited[other])
                order_csrs_[0][order_offs_[0][i]++] = other;
    }

    // generate other incremental matching orders
    for (uint i = 1; i < query_.NumEdges(); ++i)
    {
        std::vector<bool> visited(query_.NumVertices(), false);
        
        // get the first edge
        std::vector<uint>::iterator it = std::lower_bound(
            order_offs_[0].begin(), order_offs_[0].end(), i + 1
        );
        order_vs_[i][0] = *(order_vs_[0].begin() + std::distance(order_offs_[0].begin(), it));
        order_vs_[i][1] = order_csrs_[0][i];
        order_csrs_[i][0] = order_vs_[i][0];

        visited[order_vs_[i][0]] = true;
        visited[order_vs_[i][1]] = true;

        order_offs_[i][2] = order_offs_[i][1] = 1;
        for (uint j = 2; j < query_.NumVertices(); ++j)
        {
            uint max_adjacent = 0;
            uint max_adjacent_u = NOT_EXIST;
            for (size_t k = 0; k < query_.NumVertices(); k++)
            {
                uint cur_adjacent = 0u;
                if (visited[k]) continue;
                
                auto& q_nbrs = query_.GetNeighbors(k);
                for (auto& other : q_nbrs)
                    if (visited[other])
                        cur_adjacent++;

                if (!cur_adjacent) continue;
                if (
                    max_adjacent_u == NOT_EXIST ||
                    (cur_adjacent == max_adjacent &&
                        query_.GetDegree(k) > query_.GetDegree(max_adjacent_u)) ||
                    cur_adjacent > max_adjacent
                ) {
                    max_adjacent = cur_adjacent;
                    max_adjacent_u = k;
                }
            }
            order_vs_[i][j] = max_adjacent_u;
            visited[max_adjacent_u] = true;

            order_offs_[i][j] = order_offs_[i][j - 1];
            auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs)
                if (visited[other])
                    order_csrs_[i][order_offs_[i][j]++] = other;
        }
    }
    if (print_preprocessing_results_)
    {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i)
        {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j)
            {
                std::cout << order_vs_[i][j];
                if (j == 0)
                {
                    std::cout << "-";
                    continue;
                }

                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
                if (j != query_.NumVertices() - 1)
                    std::cout << "-";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::FindMatches(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        num_intermediate_results_before_index_check_++;
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;
            }
        }
        else
        {
            size_t num_results_before_recursion = num_results;
            FindMatches(order_index, depth + 1, m, num_results);
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::FindMatches_pure(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    // bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            FindMatches(order_index, depth + 1, m, num_results);
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
    }
}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::Parallel_FindMatches(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    std::vector<     std::tuple<
    uint,                       
    uint,                         
    uint,                        
    uint,                        
    uint,                        
    const std::vector<uint>&,     
    const std::vector<uint>&,     
    uint,                        
    uint,                        
    std::vector<uint>&,         
    // size_t&,                    
    const std::vector<uint>&,   
    uint                        
> > vertex_vector;

    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            // FindMatches(order_index, depth + 1, m, num_results);
            auto depth2 = depth + 1;
            uint u2 = order_vs_[order_index][depth2];
            uint u_min2 = NOT_EXIST;
            uint u_min_label2 = NOT_EXIST;
            uint u_min_size2 = UINT_MAX;
        
            // find u_min
            const auto& q_nbrs2 = query_.GetNeighbors(u2);
            // const auto& q_nbr_labels2 = query_.GetNeighborLabels(u2);

            for (uint i2 = 0u; i2 < q_nbrs2.size(); i2++)
            {
                const uint u_other2 = q_nbrs2[i2];
                const uint u_other_label2 = q_nbr_labels[i2];
        
                if (m[u_other2] == UNMATCHED) continue;
        
                const uint cur_can_size2 = data_.GetNeighbors(m[u_other2]).size();
                if (cur_can_size2 < u_min_size2)
                {
                    u_min_size2 = cur_can_size2;
                    u_min2 = u_other2;
                    u_min_label2 = u_other_label2;
                }
            }

            const auto& u_min_nbrs2 = data_.GetNeighbors(m[u_min2]);
            // const auto& u_min_nbr_labels2 = data_.GetNeighborLabels(m[u_min2]);

            for (uint i2 = 0u; i2 < u_min_nbrs2.size(); i2++)
            {
                // const uint v2 = u_min_nbrs2[i2];
                auto thread_id = omp_get_thread_num();

                // vertex_vector.emplace_back(v,
                //     v2, u2, u_min2, u_min_label2, 
                //     q_nbrs2, q_nbr_labels2, 
                //     order_index, depth2, m,
                //     u_min_nbr_labels2, i2
                // );

                ProcessNeighbor( u2, u_min2, u_min_label2,
                    order_index, depth2, m, num_results, 
                     i2, thread_id);
            }

        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
    }

// ok
    for(size_t t_1 = 0; t_1 < vertex_vector.size(); t_1 ++){
        // auto& [v_ori, v3, u3, u_min3, u_min_label3, q_nbrs3, 
        //     q_nbr_labels3, order_index, depth3, m2, u_min_nbr_labels3, i] = vertex_vector[t_1];
        
        // visited_[v_ori] = true;
        // visited_[v3] = true;
        // m2[u3] = v3;

        // ProcessNeighbor(v3, u3, u_min3, u_min_label3,
        //     q_nbrs3, q_nbr_labels3, 
        //    order_index, depth3, m, num_results, 
        //    u_min_nbr_labels3, i);
        // ProcessCandidate(u3, v3, u_min3, u_min_label3, i,
        //     q_nbrs3, q_nbr_labels3,
        //     u_min_nbr_labels3,
        //     m2,
        //     depth3,
        //     order_index,
        //     num_results);
        // visited_[v3] = false;
        // visited_[v_ori] = false;
        // m2[u3] = UNMATCHED;
    }

}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::Parallel_FindMatches2(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    size_t NUMT = NUMTHREAD;
    std::vector<size_t> local_num_result(NUMT+2,0);

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    vertex_vector.clear();
    job_queue.clear();

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;
        for(size_t i5 = 0; i5< local_vec_visited_local.size(); i5++){
            local_vec_visited_local[i5][v] = true;
        }

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            // FindMatches(order_index, depth + 1, m, num_results);
            auto depth2 = depth + 1;
            uint u2 = order_vs_[order_index][depth2];
            uint u_min2 = NOT_EXIST;
            uint u_min_label2 = NOT_EXIST;
            uint u_min_size2 = UINT_MAX;
        
            // find u_min
            const auto& q_nbrs2 = query_.GetNeighbors(u2);
            // const auto& q_nbr_labels2 = query_.GetNeighborLabels(u2);

            for (uint i2 = 0u; i2 < q_nbrs2.size(); i2++)
            {
                const uint u_other2 = q_nbrs2[i2];
                const uint u_other_label2 = q_nbr_labels[i2];
        
                if (m[u_other2] == UNMATCHED) continue;
        
                const uint cur_can_size2 = data_.GetNeighbors(m[u_other2]).size();
                if (cur_can_size2 < u_min_size2)
                {
                    u_min_size2 = cur_can_size2;
                    u_min2 = u_other2;
                    u_min_label2 = u_other_label2;
                }
            }

            const auto& u_min_nbrs2 = data_.GetNeighbors(m[u_min2]);
            // const auto& u_min_nbr_labels2 = data_.GetNeighborLabels(m[u_min2]);

            for (uint i2 = 0u; i2 < u_min_nbrs2.size(); i2++)
            {

                vertex_vector.emplace_back(v,  u_min2, u_min_label2,
                  m, i2
                );

                // ProcessNeighbor(u2, u_min2, u_min_label2,
                //     order_index, depth2, m, num_results, 
                //      i2);
            }

        }

        
        for(size_t i5 = 0; i5< local_vec_visited_local.size(); i5++){
            local_vec_visited_local[i5][v] = false;
        }
        visited_[v] = false;
        m[u] = UNMATCHED;
    }

    // if(auto_tuning == 1){
        if(vertex_vector.size() < NUMT){
            if(vertex_vector.size() == 0){
                NUMT = 1;
            }
            else{
                NUMT = vertex_vector.size();
            }
        }
    // }

    // if(auto_tuning == 1){
    //     if(vertex_vector.size() > NUMT){
    //         NUMT = vertex_vector.size();
    //     }
    // }

// ok
    #pragma omp parallel for num_threads(NUMT) schedule(auto)
    for(size_t t_1 = 0; t_1 < vertex_vector.size(); t_1 ++){

        size_t thread_id = omp_get_thread_num();
        // auto m3 = m ;
        auto& [v3,  u_min3, u_min_label3,  m2,  i2] = vertex_vector[t_1];
        //             no
        //  v3 = u_min_nbrs[i2];
        // uint v = u_min_nbrs[i2];
        uint u2 = order_vs_[order_index][depth+1];

        // auto m3 = m;
        local_vec_visited_local[thread_id][v3] = true;
        // visited_[v3] = true;
        m2[u2] = v3;
        // m3[u2] = v3;

        ProcessNeighbor(u2, u_min3, u_min_label3,
            order_index, depth+1, m2,
            local_num_result[thread_id], 
            i2, thread_id);
        // ProcessCandidate(u3, v3, u_min3, u_min_label3, i,
        //     q_nbrs3, q_nbr_labels3,
        //     u_min_nbr_labels3,
        //     m2,
        //     depth3,
        //     order_index,
        //     num_results);
        // visited_[v3] = false;
        local_vec_visited_local[thread_id][v3] = false;
        // visited_[u_min3] = true;
        m2[u2] = UNMATCHED;
        // m3[u2] = UNMATCHED;

        // for task spilt
        if(!job_queue.empty() && (t_1 > vertex_vector.size() - NUMT)){
            std::tuple<uint, uint, uint,  std::vector<uint>,  uint >job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto& [v3,  u_min3, u_min_label3,  m2,  i2] = job;
                local_vec_visited_local[thread_id][v3] = true;
                m2[u] = v3;
                // uint u2 = order_vs_[order_index][depth+1];
                ProcessNeighbor_queue(u2, u_min3, u_min_label3,
                    order_index, depth+1, m2,
                    local_num_result[thread_id], 
                    i2, thread_id);
                // Process_vertex_queue
                local_vec_visited_local[thread_id][v3] = false; 
                m2[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

}



/**
 * @brief Processes a vertex candidate in the subgraph matching process.
 * 
 * This function evaluates whether a vertex from the data graph can match the current query vertex,
 * applies filtering conditions, and if valid, adds it to the partial mapping and continues the
 * matching process recursively.
 * 
 * @param u The query vertex currently being matched
 * @param u_min The previously matched neighbor of u that provides candidates
 * @param u_min_label The edge label between u and u_min
 * @param order_index The index of the matching order being used
 * @param depth The current depth in the matching process
 * @param m The current mapping from query to data vertices (passed by reference)
 * @param num_results Reference to the counter for matches found
 * @param i The index of the candidate vertex in u_min's adjacency list
 * @param thread_id The ID of the thread executing this function
 * 
 * @return True if the candidate was successfully processed, false if rejected by filters
 * 
 * @details
 * The function applies a series of filtering steps to validate the candidate:
 * 1. Vertex and edge label matching: Ensures the vertex and connecting edge have compatible labels
 * 2. Connectivity validation: Verifies the candidate connects properly with all previously mapped vertices
 * 3. Visited check: For isomorphism matching, ensures the vertex hasn't been used in this match path
 * 
 * If all filters are passed, the function:
 * - Adds the vertex to the current mapping
 * - Marks it as visited in the thread-local visited array
 * - For complete matches, increments the result counter
 * - For partial matches, continues matching recursively
 * - Backtracks by restoring the mapping and visited status
 * 
 * Thread safety is maintained by using thread-local visited arrays and mappings.
 */
inline bool Parallel_Graphflow::ProcessNeighbor(
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
) {

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const uint v = u_min_nbrs[i];

    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    // 1. Check vertex and edge labels
    if (data_.GetVertexLabel(v) != query_.GetVertexLabel(u) || 
        u_min_nbr_labels[i] != u_min_label) {
        return false;
    }

    // 2. Check joinability with existing mappings
    bool joinable = true;
    for (uint j = 0; j < q_nbrs.size(); ++j) {
        const uint u_other = q_nbrs[j];
        const uint u_other_label = q_nbr_labels[j];

        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        const auto& nbrs = data_.GetNeighbors(m[u_other]);
        auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
        
        if (it == nbrs.end() || *it != v) {
            joinable = false;
            break;
        }

        uint dis = std::distance(nbrs.begin(), it);
        if (data_.GetNeighborLabels(m[u_other])[dis] != u_other_label) {
            joinable = false;
            break;
        }
    }
    if (!joinable) return false;

    // 3. Check visited (for isomorphism)
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Update mappings and recurse
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true;

    if (depth == query_.NumVertices() - 1) {
        ++num_results;
    } else {
        FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
inline bool Parallel_Graphflow::ProcessNeighbor_local(
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
) {

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const uint v = u_min_nbrs[i];

    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    // 1. Check vertex and edge labels
    if (data_.GetVertexLabel(v) != query_.GetVertexLabel(u) || 
        u_min_nbr_labels[i] != u_min_label) {
        return false;
    }

    // 2. Check joinability with existing mappings
    bool joinable = true;
    for (uint j = 0; j < q_nbrs.size(); ++j) {
        const uint u_other = q_nbrs[j];
        const uint u_other_label = q_nbr_labels[j];

        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        const auto& nbrs = data_.GetNeighbors(m[u_other]);
        auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
        
        if (it == nbrs.end() || *it != v) {
            joinable = false;
            break;
        }

        uint dis = std::distance(nbrs.begin(), it);
        if (data_.GetNeighborLabels(m[u_other])[dis] != u_other_label) {
            joinable = false;
            break;
        }
    }
    if (!joinable) return false;

    // 3. Check visited (for isomorphism)
    if (!homomorphism_ && visited_[v]) return false;

    // 4. Update mappings and recurse
    m[u] = v;
    visited_[v] = true;

    if (depth == query_.NumVertices() - 1) {
        ++num_results;
    } else {
        FindMatches_local_m(order_index, depth + 1, m, num_results);
    }

    // Backtrack
    visited_[v] = false;
    m[u] = UNMATCHED;

    return true;
}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
inline bool Parallel_Graphflow::ProcessCandidate(
    uint u, uint v, uint u_min, uint u_min_label, uint i,
    const std::vector<uint>& q_nbrs, 
    const std::vector<uint>& q_nbr_labels,
    const std::vector<uint>& u_min_nbr_labels,
    std::vector<uint>& m, 
    uint depth, uint order_index,
    size_t& num_results)
{
    // 1. check labels
    if (data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
        u_min_nbr_labels[i] != u_min_label) {
        return false;
    }

    // 2. check if joinable
    bool joinable = true;
    for (uint j = 0u; j < q_nbrs.size(); j++) {
        const uint u_other = q_nbrs[j];
        const uint u_other_labels = q_nbr_labels[j];

        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(
            data_.GetNeighbors(m[u_other]).begin(), 
            data_.GetNeighbors(m[u_other]).end(), 
            v);
        uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
        if (it == data_.GetNeighbors(m[u_other]).end() ||
            *it != v ||
            data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels) {
            joinable = false;
            break;
        }
    }
    if (!joinable) return false;


    // 3. check if visited
    if (!homomorphism_ && visited_[v]) return false;

    // 4. add a vertex mapping
    m[u] = v;
    visited_[v] = true;

    if (depth == query_.NumVertices() - 1) {
        num_results++;
    } else {
        FindMatches(order_index, depth + 1, m, num_results);
    }
    
    visited_[v] = false;
    m[u] = UNMATCHED;

    std::cout << num_results << std::endl;
    
    return true;
}



/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::FindMatches_local(uint order_index, uint depth, std::vector<uint> m, size_t &num_results, size_t thread_id)
{
    // if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    // bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        local_vec_visited_local[thread_id][v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
        }
        
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
    }
}




/**
 * @brief Processes a single vertex in the parallel matching algorithm using job queue approach.
 * 
 * This function evaluates whether a candidate vertex from the data graph can be matched to
 * a query vertex, and if valid, adds subsequent matching tasks to a concurrent job queue
 * for parallel processing by worker threads.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The parent query vertex of u in the matching order.
 * @param v_idx The index of the candidate vertex in the DCS array.
 * @param m A reference to the current partial matching (query vertex -> data vertex).
 * @param num_results A reference to the thread-local counter for valid matches found.
 * @param depth The current depth in the matching process.
 * @param order_index The index of the current matching order being used.
 * @param thread_id The ID of the thread executing this function.
 * 
 * @return Returns true if the candidate vertex was successfully processed, false otherwise.
 * 
 * @details
 * The function follows these steps:
 * 1. Retrieves the candidate vertex from the DCS structure.
 * 2. Applies filtering based on the index constraint (d2[u][v]).
 * 3. Verifies candidate joinability with previously matched vertices.
 * 4. Ensures the vertex hasn't been visited in the current search path (for isomorphism).
 * 5. If all checks pass, adds the vertex to the current mapping and marks it as visited.
 * 6. For complete matches (at maximum depth), increments the match counter.
 * 7. For partial matches, instead of recursive calls, adds new matching tasks to the job queue.
 * 8. Backtracks by unmarking the vertex and removing it from the mapping.
 * 
 * This function enables work-stealing parallelism by scheduling tasks that can be executed
 * by any available thread, improving load balancing in the parallel matching process.
 */
void Parallel_Graphflow::Process_vertex_queue(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    size_t NUMT = NUMTHREAD;
    std::vector<size_t> local_num_result(NUMT+2,0);

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    vertex_vector.clear();
    job_queue.clear();

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;
        for(size_t i5 = 0; i5< local_vec_visited_local.size(); i5++){
            local_vec_visited_local[i5][v] = true;
        }

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            // FindMatches(order_index, depth + 1, m, num_results);
            auto depth2 = depth + 1;
            uint u2 = order_vs_[order_index][depth2];
            uint u_min2 = NOT_EXIST;
            uint u_min_label2 = NOT_EXIST;
            uint u_min_size2 = UINT_MAX;
        
            // find u_min
            const auto& q_nbrs2 = query_.GetNeighbors(u2);
            // const auto& q_nbr_labels2 = query_.GetNeighborLabels(u2);

            for (uint i2 = 0u; i2 < q_nbrs2.size(); i2++)
            {
                const uint u_other2 = q_nbrs2[i2];
                const uint u_other_label2 = q_nbr_labels[i2];
        
                if (m[u_other2] == UNMATCHED) continue;
        
                const uint cur_can_size2 = data_.GetNeighbors(m[u_other2]).size();
                if (cur_can_size2 < u_min_size2)
                {
                    u_min_size2 = cur_can_size2;
                    u_min2 = u_other2;
                    u_min_label2 = u_other_label2;
                }
            }

            const auto& u_min_nbrs2 = data_.GetNeighbors(m[u_min2]);
            // const auto& u_min_nbr_labels2 = data_.GetNeighborLabels(m[u_min2]);

            for (uint i2 = 0u; i2 < u_min_nbrs2.size(); i2++)
            {

                vertex_vector.emplace_back(v,  u_min2, u_min_label2,
                  m, i2
                );

                // ProcessNeighbor(u2, u_min2, u_min_label2,
                //     order_index, depth2, m, num_results, 
                //      i2);
            }

        }

        
        for(size_t i5 = 0; i5< local_vec_visited_local.size(); i5++){
            local_vec_visited_local[i5][v] = false;
        }
        visited_[v] = false;
        m[u] = UNMATCHED;
    }
}



/**
 * @brief Processes a vertex candidate and enqueues subsequent matching tasks for parallel execution.
 * 
 * This function evaluates whether a vertex from the data graph can match a query vertex,
 * and if valid, adds new matching tasks to a concurrent job queue rather than recursively
 * processing them, enabling parallel work distribution.
 * 
 * @param u The current query vertex being matched
 * @param u_min The parent query vertex in the matching order
 * @param u_min_label The edge label between u and u_min
 * @param order_index The index of the current matching order
 * @param depth The current depth in the search tree
 * @param m Reference to the current mapping from query to data vertices
 * @param num_results Reference to the counter for matches found
 * @param i Index of the candidate vertex in u_min's neighborhood
 * @param thread_id ID of the thread executing this function
 * 
 * @return true if the candidate was successfully processed, false otherwise
 * 
 * @details
 * The function implements a modified backtracking approach with task queue integration:
 * 
 * 1. Candidate Retrieval and Validation:
 *    - Gets the candidate vertex from u_min's adjacency list
 *    - Checks vertex and edge label compatibility
 *    - Verifies connectivity with all previously matched vertices
 *    - Ensures the vertex hasn't been visited before (for isomorphism)
 * 
 * 2. Match Processing:
 *    - If valid, adds the vertex to the current mapping
 *    - For complete matches (at maximum depth), increments the match counter
 * 
 * 3. Task Generation:
 *    - For partial matches, instead of recursive calls:
 *      a) Identifies the next query vertex to match
 *      b) Finds the optimal matched neighbor to extend from
 *      c) Creates tasks for each potential match
 *      d) Pushes these tasks to a shared job queue for parallel execution
 * 
 * 4. Backtracking:
 *    - Cleans up by restoring the visited status and mapping
 * 
 * This approach enables dynamic load balancing by allowing threads to steal work
 * from the shared queue when they become idle.
 */
inline bool Parallel_Graphflow::ProcessNeighbor_queue(
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
) {

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const uint v = u_min_nbrs[i];

    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    // 1. Check vertex and edge labels
    if (data_.GetVertexLabel(v) != query_.GetVertexLabel(u) || 
        u_min_nbr_labels[i] != u_min_label) {
        return false;
    }

    // 2. Check joinability with existing mappings
    bool joinable = true;
    for (uint j = 0; j < q_nbrs.size(); ++j) {
        const uint u_other = q_nbrs[j];
        const uint u_other_label = q_nbr_labels[j];

        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        const auto& nbrs = data_.GetNeighbors(m[u_other]);
        auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
        
        if (it == nbrs.end() || *it != v) {
            joinable = false;
            break;
        }

        uint dis = std::distance(nbrs.begin(), it);
        if (data_.GetNeighborLabels(m[u_other])[dis] != u_other_label) {
            joinable = false;
            break;
        }
    }
    if (!joinable) return false;

    // 3. Check visited (for isomorphism)
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Update mappings and recurse
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true;

    if (depth == query_.NumVertices() - 1) {
        ++num_results;
    } else {

  // FindMatches(order_index, depth + 1, m, num_results);
  auto depth2 = depth + 1;
  uint u2 = order_vs_[order_index][depth2];
  uint u_min2 = NOT_EXIST;
  uint u_min_label2 = NOT_EXIST;
  uint u_min_size2 = UINT_MAX;

  // find u_min
  const auto& q_nbrs2 = query_.GetNeighbors(u2);
  // const auto& q_nbr_labels2 = query_.GetNeighborLabels(u2);

  for (uint i2 = 0u; i2 < q_nbrs2.size(); i2++)
  {
      const uint u_other2 = q_nbrs2[i2];
      const uint u_other_label2 = q_nbr_labels[i2];

      if (m[u_other2] == UNMATCHED) continue;

      const uint cur_can_size2 = data_.GetNeighbors(m[u_other2]).size();
      if (cur_can_size2 < u_min_size2)
      {
          u_min_size2 = cur_can_size2;
          u_min2 = u_other2;
          u_min_label2 = u_other_label2;
      }
  }

  const auto& u_min_nbrs2 = data_.GetNeighbors(m[u_min2]);
  // const auto& u_min_nbr_labels2 = data_.GetNeighborLabels(m[u_min2]);

  for (uint i2 = 0u; i2 < u_min_nbrs2.size(); i2++)
  {
        job_queue.push(
            std::make_tuple(v, u_min2, u_min_label2,
                m, i2)
        );

      // ProcessNeighbor(u2, u_min2, u_min_label2,
      //     order_index, depth2, m, num_results, 
      //      i2);
  }
        // FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}




/**
 * @brief Recursively finds all matches for a query graph in the data graph.
 * 
 * This function implements a backtracking algorithm for subgraph matching, exploring
 * all possible mappings between query and data vertices according to a predetermined
 * matching order.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the search tree (number of matched vertices so far).
 * @param m The current vertex mapping from query to data graph vertices.
 * @param num_results A reference to the counter for matches found.
 * 
 * @details
 * The function works in several key steps:
 * 
 * 1. Vertex Selection Strategy:
 *    - Identifies the next query vertex to match according to the matching order
 *    - Finds the optimal already-matched neighbor (u_min) to extend from
 *    - Selects the neighbor with the smallest adjacency list to minimize the search space
 * 
 * 2. Candidate Filtering:
 *    - Examines each neighbor of u_min's mapping in the data graph
 *    - Applies a filtering pipeline with increasingly strict conditions:
 *      a) Label matching: Ensures vertex and edge labels match
 *      b) Structural validation: Verifies connectivity with all previously matched neighbors
 *      c) Isomorphism checking: Prevents reusing vertices when homomorphism is disabled
 * 
 * 3. Recursive Exploration:
 *    - For each valid candidate, adds the mapping and marks the vertex as visited
 *    - For complete matches (maximum depth), increments the result counter
 *    - For partial matches, recursively explores further mappings
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * 4. Statistical Tracking:
 *    - Maintains various metrics for performance analysis
 *    - Tracks empty candidate sets and match attempts without results
 *    - Implements early termination when reaching result limits or time constraints
 * 
 * The function uses efficient neighborhood intersection through binary search to
 * ensure consistency between vertex connections in the match.
 */
void Parallel_Graphflow::FindMatches_local_m(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
{
    // if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;

    // find u_min
    const auto& q_nbrs = query_.GetNeighbors(u);
    const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    for (uint i = 0u; i < q_nbrs.size(); i++)
    {
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size)
        {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    // bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++)
    {
        const uint v = u_min_nbrs[i];

        // 1. check labels
        if (
            data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
            u_min_nbr_labels[i] != u_min_label
        ) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(), v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                it == data_.GetNeighbors(m[u_other]).end() ||
                *it != v ||
                data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            FindMatches_local_m(order_index, depth + 1, m, num_results);
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
    }
}




/**
 * @brief Performs initial subgraph matching to find all matching instances in the data graph.
 * 
 * This function initiates the subgraph matching process by finding all data vertices that match
 * the first vertex in the query matching order. For each candidate vertex, it starts a recursive
 * matching process to discover complete matches.
 * 
 * @details
 * The function works as follows:
 * 
 * 1. Creates an empty match vector to store the mapping between query and data vertices.
 * 
 * 2. Iterates through all vertices in the data graph to find potential matches for the
 *    first vertex in the matching order.
 * 
 * 3. For each data vertex with a matching label:
 *    - Maps the first query vertex to this data vertex
 *    - Marks the data vertex as visited to prevent reuse in the same match
 *    - Calls FindMatches() to recursively extend the partial match
 *    - After returning from the recursive call, backtracks by unmarking the data vertex
 *      and resetting the mapping
 * 
 * 4. The total count of complete matches is accumulated in num_initial_results_.
 * 
 * This function uses the optimal vertex matching order (order_vs_[0]) that was previously
 * generated during preprocessing to guide the matching process efficiently.
 */
void Parallel_Graphflow::InitialMatching()
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);
    
    for (size_t i = 0; i < data_.NumVertices(); i++)
    if (data_.GetVertexLabel(i) != NOT_EXIST)
    if (query_.GetVertexLabel(order_vs_[0][0]) == data_.GetVertexLabel(i))
    {
        m[order_vs_[0][0]] = i;
        visited_[i] = true;

        FindMatches(0, 1, m, num_initial_results_);

        visited_[i] = false;
        m[order_vs_[0][0]] = UNMATCHED;
    }
}


/**
 * @brief Determines if a potential data edge matches any edge in the query graph.
 * 
 * This function checks whether a candidate edge in the data graph could potentially
 * be part of a match for the query graph by comparing its vertex labels and edge label
 * against all edges in the query graph.
 * 
 * @param v1 The ID of the first vertex of the edge in the data graph.
 * @param v2 The ID of the second vertex of the edge in the data graph.
 * @param label The label of the edge between v1 and v2.
 * 
 * @return Returns true if the edge does NOT match any query edge (can be filtered out),
 *         returns false if it matches at least one query edge (should be considered).
 * 
 * @details
 * The function performs a comprehensive check by:
 * 1. Iterating through all edges in the query graph
 * 2. For each query edge, checking both possible orientations of the data edge:
 *    - v1→v2 orientation: Vertex and edge labels of v1,v2 match the query edge
 *    - v2→v1 orientation: Vertex and edge labels of v2,v1 match the query edge
 * 3. If any match is found, returns false (indicating the edge is relevant)
 * 4. If no matches are found after checking all query edges, returns true (edge can be filtered)
 * 
 * This function is used to quickly filter out data graph edges that cannot be part of any
 * valid subgraph match, reducing the search space for subsequent matching operations.
 */
bool Parallel_Graphflow::Classify(uint v1, uint v2, uint label){

    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        uint u1 = order_vs_[i][0], u2 = order_vs_[i][1];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        
        // check if any query edge match (v1 --> v2)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<2>(temp_q_labels) == label
        ) {
            return false;
        }
        // check if any query edge match (v2 --> v1)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<2>(temp_q_labels) == label
        ) {
            return false;
        }
    }
    return true;

}


/**
 * @brief Adds an edge to the data graph and finds all new subgraph matches containing it.
 * 
 * This function adds a new edge to the data graph and then identifies all subgraph matches
 * that are created by this addition. The matches are found using a parallel execution strategy.
 * 
 * @param v1 The ID of the first vertex of the edge to add
 * @param v2 The ID of the second vertex of the edge to add
 * @param label The label to assign to the new edge
 * 
 * @details
 * The function performs the following operations in sequence:
 * 
 * 1. Edge Addition:
 *    - First adds the specified edge to the data graph
 * 
 * 2. Match Enumeration:
 *    - Creates an empty mapping vector to track vertex correspondences
 *    - If maximum result count is 0, returns immediately without finding matches
 *    - Otherwise, loops through all query graph edges to find potential matches:
 * 
 * 3. For Each Query Edge:
 *    - Retrieves the edge labels from the query graph
 *    - Checks if the new edge matches the query edge in either direction:
 *      a) v1→v2 orientation (direct match)
 *      b) v2→v1 orientation (reverse match)
 *    - When a match is found:
 *      a) Maps the corresponding query vertices to data vertices
 *      b) Marks the data vertices as visited in both global and thread-local arrays
 *      c) Uses parallel subgraph matching to find all complete matches containing this edge
 *      d) Restores the visited status and mappings after matching
 *    - Terminates early if maximum result count is reached or time limit is exceeded
 * 
 * 4. Result Tracking:
 *    - Updates the positive results counter with the number of new matches found
 * 
 * The function uses parallel matching to efficiently find all new matches created by the
 * edge addition, leveraging the precomputed matching orders for performance.
 */
void Parallel_Graphflow::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    if (max_num_results_ == 0) return;

    size_t num_results = 0ul;
    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        uint u1 = order_vs_[i][0], u2 = order_vs_[i][1];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);
        
        // check if any query edge match (v1 --> v2)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<2>(temp_q_labels) == label
        ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            
            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; 
            }

            Parallel_FindMatches2(i, 2, m, num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; 
            }

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        // check if any query edge match (v2 --> v1)
        if (
            std::get<0>(temp_q_labels) == data_.GetVertexLabel(v2) &&
            std::get<1>(temp_q_labels) == data_.GetVertexLabel(v1) &&
            std::get<2>(temp_q_labels) == label
        ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; 
            }

            Parallel_FindMatches2(i, 2, m, num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; 
            }

            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:
    num_positive_results_ += num_results;
}



/**
 * @brief Removes an edge from the data graph after finding all affected subgraph matches.
 * 
 * This function first identifies and counts all existing matches that would be affected by
 * the edge removal, then removes the edge from the data graph.
 * 
 * @param v1 The ID of the first vertex of the edge to remove
 * @param v2 The ID of the second vertex of the edge to remove
 * 
 * @details
 * The function performs the following operations in sequence:
 * 
 * 1. Match Identification:
 *    - Retrieves the edge label from the data graph
 *    - For each possible matching order in the query graph:
 *      a) Checks if the edge matches a query edge in either direction (v1→v2 or v2→v1)
 *      b) Creates partial mappings with the matched edge vertices
 *      c) Marks the mapped vertices as visited in both global and thread-local visited arrays
 *      d) Uses parallel matching to find all complete matches containing this edge
 *      e) Restores the visited status and mappings after matching
 *    - Terminates early if the maximum result count is reached or time limit is exceeded
 * 
 * 2. Edge Removal:
 *    - After counting all affected matches, removes the edge from the data graph
 *    - Updates the negative results counter with the number of matches found
 * 
 * This approach allows the algorithm to track which subgraph matches will be invalidated
 * by the edge removal, facilitating efficient maintenance of the match set during graph updates.
 */
void Parallel_Graphflow::RemoveEdge(uint v1, uint v2)
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    std::tuple labels = data_.GetEdgeLabel(v1, v2);
    
    size_t num_results = 0ul;
    if (max_num_results_ == 0) goto END_ENUMERATION;

    for (uint i = 0; i < query_.NumEdges(); i++)
    {
        uint u1 = order_vs_[i][1], u2 = order_csrs_[i][0];
        auto temp_q_labels = query_.GetEdgeLabel(u1, u2);

        if (
            std::get<0>(temp_q_labels) == std::get<0>(labels) &&
            std::get<1>(temp_q_labels) == std::get<1>(labels) &&
            std::get<2>(temp_q_labels) == std::get<2>(labels)
        ) {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; 
            }

            Parallel_FindMatches2(i, 2, m, num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; 
            }

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (
            std::get<1>(temp_q_labels) == std::get<0>(labels) &&
            std::get<0>(temp_q_labels) == std::get<1>(labels) &&
            std::get<2>(temp_q_labels) == std::get<2>(labels)
        ) {
            m[u1] = v2;
            m[u2] = v1;
            visited_[v2] = true;
            visited_[v1] = true;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; 
            }

            Parallel_FindMatches2(i, 2, m, num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; 
            }

            visited_[v2] = false;
            visited_[v1] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
    }
    END_ENUMERATION:

    num_negative_results_ += num_results;
    data_.RemoveEdge(v1, v2);
}


/**
 * @brief Adds a vertex to the data graph in a thread-safe manner.
 * 
 * This function safely adds a new vertex with the specified ID and label to the data graph
 * using OpenMP's critical section to prevent race conditions in a multi-threaded environment.
 * 
 * @param id The ID of the new vertex to add
 * @param label The label to assign to the new vertex
 * 
 * @details
 * The function creates a critical section using OpenMP's #pragma omp critical directive
 * to ensure exclusive access when modifying the graph structure. This prevents potential
 * data corruption that could occur if multiple threads attempted to modify the graph
 * simultaneously.
 * 
 * Within the critical section, the function performs two main operations:
 * 1. Delegates the actual vertex addition to the underlying data graph's AddVertex method
 * 2. Resizes the visited_ array to accommodate the new vertex's ID, initializing the new
 *    entry to false
 * 
 * This approach ensures that both the graph structure and the auxiliary data structures
 * used for matching are updated atomically, maintaining consistency in the data structures.
 */
void Parallel_Graphflow::AddVertex(uint id, uint label)
{
    #pragma omp critical
    {
        data_.AddVertex(id, label);
    
        visited_.resize(id + 1, false);
    }

}


/**
 * @brief Removes a vertex from the data graph in a thread-safe manner.
 * 
 * This function safely removes a vertex and all its associated edges from the data graph
 * using OpenMP's critical section to prevent race conditions in a multi-threaded environment.
 * 
 * @param id The ID of the vertex to remove
 * 
 * @details
 * The function creates a critical section using OpenMP's #pragma omp critical directive
 * to ensure exclusive access when modifying the graph structure. This prevents potential
 * data corruption that could occur if multiple threads attempted to modify the graph
 * simultaneously.
 * 
 * Within the critical section, the function delegates the actual vertex removal to
 * the underlying data graph's RemoveVertex method, which handles:
 * 1. Removing the specified vertex from the graph
 * 2. Removing all edges connected to this vertex
 * 3. Updating any internal data structures or indexes
 * 
 * Note: After calling this function, any references to the removed vertex or its edges
 * will be invalid and should not be accessed.
 */
void Parallel_Graphflow::RemoveVertex(uint id)
{
    #pragma omp critical
    data_.RemoveVertex(id);
}

/**
 * @brief Calculates the memory cost of the data structures used by the matching algorithm.
 * 
 * This function computes the memory usage statistics of the Parallel_Graphflow algorithm,
 * focusing on the number of edges and vertices being maintained in internal data structures.
 * 
 * @param[out] num_edges Reference parameter that will be updated with the total number of 
 *                      edge entries stored in various data structures.
 * @param[out] num_vertices Reference parameter that will be updated with the total number of
 *                         vertex entries stored in various data structures.
 * 
 * @details
 * The function analyzes memory usage of various components:
 * - Matching order data structures (order_vs_, order_csrs_, order_offs_)
 * - Thread-local data structures (local_vec_m, local_vec_visited_local)
 * - Any additional auxiliary structures used for parallel matching
 * 
 * This information can be useful for:
 * 1. Performance profiling and optimization
 * 2. Monitoring memory usage in resource-constrained environments
 * 3. Comparing memory efficiency with other matching algorithms
 * 
 * Note: Currently this function is implemented as a placeholder that sets both
 * parameters to 0. A complete implementation would analyze the actual data structures
 * and return their memory usage.
 */
void Parallel_Graphflow::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}
