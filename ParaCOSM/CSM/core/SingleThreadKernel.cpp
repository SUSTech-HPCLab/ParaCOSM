#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/Parallel_GraphFlow/parallel_graphflow.h"


/**
 * @brief Thread-local version of FindMatches for parallel subgraph matching.
 * 
 * This function performs depth-first search for subgraph matching within a single thread's
 * context. It maintains thread-local data structures to avoid synchronization overhead
 * while ensuring correctness in a multi-threaded environment.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the recursive search process.
 * @param m A vector representing the current partial mapping of query vertices to data vertices.
 * @param num_results A reference to the thread-local counter for matches found.
 * @param thread_id The ID of the thread executing this function instance.
 * 
 * @details
 * The function follows a standard backtracking approach:
 * 
 * 1. Retrieves the next query vertex (u) and its parent (u_min) from the matching order.
 * 
 * 2. Filters candidates through a pipeline of increasingly strict conditions:
 *    a) Index-based filtering: Verifies if the candidate satisfies structural constraints via d2[u][v]
 *    b) Joinability checking: Ensures the candidate connects properly to previously matched vertices
 *    c) Visitation checking: Prevents reusing vertices in the current search path (for isomorphism)
 * 
 * 3. For each valid candidate:
 *    - Adds it to the current mapping and marks it as visited in thread-local structures
 *    - If at maximum depth, increments the match counter
 *    - Otherwise, recursively explores deeper matches
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * This thread-local implementation avoids race conditions by using separate visitation arrays
 * for each thread while maintaining the exploration efficiency of the sequential algorithm.
 */
void Turboflux_FindMatches_SingleThread(uint order_index, uint depth, std::vector<uint>& m, 
     size_t &num_results, size_t thread_id)
{
    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];

    // enumerate each neighbor of m[u_min]
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        if (d2[u][v] == 0) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
        {
            const auto& u_backward = join_check_vs_[order_index][u][i];
            const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
            const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

            if (std::get<2>(d_elabel) != u_backward_elabel)
            {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

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
            Turboflux_FindMatches_SingleThread(order_index, depth + 1, m, num_results, thread_id);
        }
        
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
    }
}




/**
 * @brief Processes a single vertex in the parallel matching algorithm.
 * 
 * This function evaluates whether a candidate vertex from the data graph can be matched to 
 * a query vertex in the current partial matching. It performs filtering based on indexing 
 * constraints, joinability with previously matched vertices, and visitation status.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The parent query vertex of u in the matching order.
 * @param v_idx The index of the candidate vertex in the DCS array.
 * @param m A reference to the current partial matching (query vertex -> data vertex).
 * @param num_results A reference to the counter for valid matches found.
 * @param depth The current depth in the matching process.
 * @param order_index The index of the current matching order being used.
 * @param thread_id The ID of the thread executing this function.
 * 
 * @return Returns true if the candidate vertex was successfully processed, false otherwise.
 * 
 * @details
 * The function follows these steps:
 * 1. Retrieves the candidate vertex from the DCS structure.
 * 2. Checks if it satisfies the index constraint (d2[u][v]).
 * 3. Verifies if the candidate is joinable with all previously matched vertices.
 * 4. Ensures the vertex hasn't been visited in the current search path (for isomorphism).
 * 5. If all checks pass, adds the vertex to the current matching and marks it as visited.
 * 6. Either records a complete match or recursively explores further matches.
 * 7. Backtracks by unmarking the vertex and removing it from the matching.
 * 
 * This thread-safe version uses thread-local data structures to maintain correctness
 * during parallel execution.
 */
inline bool Turboflux_ProcessNode(uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    size_t &num_results, 
     uint depth, size_t order_index, size_t thread_id) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];     

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable2 = true;
    for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++) {
        const auto& u_backward = join_check_vs_[order_index][u][i];
        const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
        const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

        if (std::get<2>(d_elabel) != u_backward_elabel) {
            joinable2 = false;
            break;
        }
    }
    if (!joinable2) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add vertex mapping and mark visited
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true;

    // Process result or recurse
    if (depth == query_.NumVertices() - 1) {
        num_results++;
    } else {
        Turboflux_FindMatches_SingleThread(order_index, depth + 1, m, num_results, thread_id);
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
void Graphflow_FindMatches_SingleThread(uint order_index, uint depth, std::vector<uint> m, size_t &num_results, size_t thread_id)
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
            Graphflow_FindMatches_SingleThread(order_index, depth + 1, m, num_results, thread_id);
        }
        
        local_vec_visited_local[thread_id][v] = false;
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
void _FindMatches_SingleThread(uint order_index, uint depth, std::vector<uint> m, size_t &num_results, size_t thread_id)
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
            Graphflow_FindMatches_SingleThread(order_index, depth + 1, m, num_results, thread_id);
        }
        
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
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
inline bool Graphflow_ProcessNode(
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
        Graphflow_FindMatches_SingleThread(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}


/**
 * @brief Processes a single vertex in the subgraph matching algorithm for layer 1.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * @param thread_id The ID of the thread processing this vertex.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited by the current thread.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Symbi_ProcessNode(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth, size_t thread_id) {

    
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {

        // #ifdef TASK_SPILT
        if(query_.NumVertices() < 7){
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }else{

        if(depth < 5){
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }else{

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                {
                    if (temp_extendable[i].E < isolate_minE2)
                    {
                        isolate_minE2 = temp_extendable[i].E; // local
                        isolate_u2 = i; // local
                    }
                }
                else
                {
                    if (temp_extendable[i].E < non_isolate_minE2)
                    {
                        non_isolate_minE2 = temp_extendable[i].E; // local
                        non_isolate_u2 = i; // local
                    }
                }
            }

            if (non_isolate_minE2 == NOT_EXIST)
                u2 = isolate_u2; // local
            else
                u2 = non_isolate_u2; // local
                
            uint u_min2 = temp_extendable[u2].u_min;
            temp_extendable[u2] = {}; // local
        
            size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            // std::cout << "total_size2: " << total_size2 << std::endl;
            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                auto m3 = m;
                // job
                job_queue.push(std::make_tuple(u2, u_min2, v_idx2, m3, temp_extendable, depth2, v));

            
            }
        }
    }

    }
    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}


uint CaLiG_searchCore(uint th, ska::flat_hash_map<uint, uint>& m, vec& used, const vec& c, const vector<u_set>& c_n, const vec& s, const vector<u_set>& s_n, ska::flat_hash_map<uint, u_set>& c2check){
    uint result = 0;
    if(th == c.size()){
        vector<u_set> candidates;
        if(shellCand(candidates, m, s, s_n, used)) return 0;
        u_set used_v;
        return numAdd(0, candidates, used_v);
    }
    u_set candidates;
    auto p_ui = c_n[th].begin();
    auto ui = *p_ui;
    if(c_n[th].size()>1){
        u_set* temp;
        temp = &G[m[ui]].cand[ui][c[th]];
        candidates = intersection(*temp,G[m[ui]].cand[ui][c[th]]);
        for(++p_ui;p_ui!=c_n[th].end();++p_ui){
            ui = *p_ui;
            candidates = intersection(candidates,G[m[ui]].cand[ui][c[th]]);
            if(candidates.empty()) return 0;
        }
        for(auto& vi : used){
            candidates.erase(vi);
        }
        for(auto& vth : candidates){
            m[c[th]] = vth;
            if(c2check.count(c[th])){
                bool to_continue = false;
                for (auto& shell : c2check[c[th]]){
                    if(notExit(s[shell], s_n[shell], m)) {
                        to_continue = true;
                    }
                }
                if(to_continue) {
                    continue;
                }
            }
            used.emplace_back(vth);
            result += searchCore(th+1, m, used, c, c_n, s, s_n, c2check);
            used.pop_back();
        }
    }
    else{
        candidates = G[m[ui]].cand[ui][c[th]];
        for(auto& vth : candidates){
            if(!isInVec(vth, used)){
                m[c[th]] = vth;
                if(c2check.count(c[th])){
                    bool to_continue = false;
                    for (auto& shell : c2check[c[th]]){
                        if(notExit(s[shell], s_n[shell], m)) {
                            to_continue = true;
                        }
                    }
                    if(to_continue) {
                        continue;
                    }

                }
                used.emplace_back(vth);
                result += searchCore(th+1, m, used, c, c_n, s, s_n, c2check);
                used.pop_back();
            }
        }
    }
    return result;
}