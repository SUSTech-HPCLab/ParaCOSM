#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph_storage/graph.h"
#include "matching_executor/Parallel_GraphFlow/parallel_graphflow.h"
#include "matching_executor/Parallel_TurboFlux/parallel_turboflux.h"
#include "matching_executor/Parallel_SymBi/parallel_symbi.h"
#include "core/inner_executor/FindMatchesKernel.h"

#include "utils/types.h"
#include "utils/globals.h"
#include "graph_storage/graph.h"





inline void Parallel_FindMatches_local3_simplified(uint depth, std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results, size_t thread_id,
    const Graph& query_graph,
    const Graph& data_graph,
    const std::vector<std::vector<uint>>& eidx,
    const std::vector<ska::flat_hash_map<uint, std::vector<uint>>>& DCS,
    const std::vector<ska::flat_hash_map<uint, bool>>& d2,
    const std::vector<ska::flat_hash_map<uint, uint>>& n2,
    std::vector<std::vector<bool>>& local_vec_visited_local,
    bool homomorphism,
    const std::vector<std::vector<uint>>& treeNode_neighbors)
{

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    for (uint i = 0; i < extendable.size(); i++)
    {
        
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_graph.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    // If no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]

    // Main parallel for loop
    auto dcs_it = DCS[eidx[u_min][u]].find(m[u_min]);
    if (dcs_it == DCS[eidx[u_min][u]].end()) {
        return; // No candidates found
    }
    const auto& candidates = dcs_it->second;
    for (size_t v_idx = 0; v_idx < candidates.size(); ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        auto v = candidates[v_idx];

        // 1. Check index
        auto d2_it = d2[u].find(v);
        if (d2_it == d2[u].end() || !d2_it->second) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_neighbors[u])
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto dcs_it = DCS[eidx[u_other][u]].find(m[u_other]);
            if (dcs_it == DCS[eidx[u_other][u]].end()) {
                joinable = false;
                break;
            }
            
            const auto& neighbor_candidates = dcs_it->second;
            auto it = std::lower_bound(neighbor_candidates.begin(), neighbor_candidates.end(), v);
            if (
                it == neighbor_candidates.end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        if (!homomorphism && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        {
            m[u] = v;
            local_vec_visited_local[thread_id][v] = true; // imp:
        }

        std::vector<Parrllel_SymBi::ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_neighbors[u])
        {
            if (m[u_other] != UNMATCHED) continue;

            auto n2_it = n2[eidx[u][u_other]].find(m[u]);
            if (n2_it != n2[eidx[u][u_other]].end() && 
                n2_it->second < temp_extendable[u_other].E)
            {
                temp_extendable[u_other].E = n2_it->second;
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs ++;
        }

        if (depth == query_graph.NumVertices() - 1) // match complete
        {
            num_results++;
        }
        else
        {

                Parallel_FindMatches_local3_simplified(depth + 1, m, temp_extendable, num_results, thread_id,
                    query_graph, data_graph, eidx, DCS, d2, n2, local_vec_visited_local, homomorphism, treeNode_neighbors);
        }

        // for backtrack
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
    }
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
inline bool process_vertex_layer1_simplified(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, 
    size_t &num_results, 
    uint depth, size_t thread_id,
    const Graph& query_graph,
    const Graph& data_graph,
    const std::vector<std::vector<uint>>& eidx,
    const std::vector<ska::flat_hash_map<uint, std::vector<uint>>>& DCS,
    const std::vector<ska::flat_hash_map<uint, bool>>& d2,
    const std::vector<ska::flat_hash_map<uint, uint>>& n2,
    std::vector<std::vector<bool>>& local_vec_visited_local,
    bool homomorphism,
    const std::vector<std::vector<uint>>& treeNode_neighbors,
    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>,
    std::vector<Parrllel_SymBi::ExtendableVertex>,  uint , uint> >& job_queue
    
) {

    
         // __builtin_prefetch(&DCS[eidx[u_min][u]][m[u_min]][v_idx], 1, 0);
     auto dcs_it = DCS[eidx[u_min][u]].find(m[u_min]);
     if (dcs_it == DCS[eidx[u_min][u]].end()) {
         return false; // No candidates found
     }
     auto v = dcs_it->second[v_idx];

         // 1. Check index
     auto d2_it = d2[u].find(v);
     if (d2_it == d2[u].end() || !d2_it->second) return false;

    // 2. Check if joinable
    bool joinable = true;
         for (auto& u_other: treeNode_neighbors[u]) {
         if (m[u_other] == UNMATCHED || u_other == u_min) continue;

         auto dcs_it = DCS[eidx[u_other][u]].find(m[u_other]);
         if (dcs_it == DCS[eidx[u_other][u]].end()) {
             joinable = false;
             break;
         }
         
         const auto& neighbor_candidates = dcs_it->second;
         auto it = std::lower_bound(neighbor_candidates.begin(), neighbor_candidates.end(), v);
         if (it == neighbor_candidates.end() || *it != v) {
             joinable = false;
             break;
         }
     }

    if (!joinable) return false;

         // 3. Check if visited
     if (!homomorphism && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

         std::vector<Parrllel_SymBi::ExtendableVertex> temp_extendable(extendable); // local
     for (auto& u_other: treeNode_neighbors[u]) {
         if (m[u_other] != UNMATCHED) continue;

         auto n2_it = n2[eidx[u][u_other]].find(m[u]);
         if (n2_it != n2[eidx[u][u_other]].end() && 
             n2_it->second < temp_extendable[u_other].E) {
             temp_extendable[u_other].E = n2_it->second;
             temp_extendable[u_other].u_min = u;
         }
         temp_extendable[u_other].matched_nbrs++;
     }

         if (depth == query_graph.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);

        // #ifdef TASK_SPILT
        if(query_graph.NumVertices() < 7){
             Parallel_FindMatches_local3_simplified(depth + 1, m, temp_extendable, num_results, thread_id,
                 query_graph, data_graph, eidx, DCS, d2, n2, local_vec_visited_local, homomorphism, treeNode_neighbors);
         }else{

        if(depth < 5){
            Parallel_FindMatches_local3_simplified(depth + 1, m, temp_extendable, num_results,thread_id,
                query_graph, data_graph, eidx, DCS, d2, n2, local_vec_visited_local, homomorphism, treeNode_neighbors);
        }else{

            // the place for task spilt 

            // if(kongxian){spilt, 逻辑如下:}

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_graph.GetNeighbors(i).size())
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
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            auto dcs_it2 = DCS[eidx[u_min2][u2]].find(m[u_min2]);
            if (dcs_it2 == DCS[eidx[u_min2][u2]].end()) {
                // continue; // No candidates found
            }
             size_t total_size2 = dcs_it2->second.size();

            // std::cout << "total_size2: " << total_size2 << std::endl;
            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                // process_vertex(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                auto m3 = m;
                // queue version
                // vertex_queue.emplace(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                
                // vertex_vector.emplace_back(u2, u_min2, v_idx2, m3, temp_extendable,  depth2, v);
                // job
                job_queue.push(std::make_tuple(u2, u_min2, v_idx2, m3, temp_extendable, depth2, v));

            
            }

                    // this is for backtrack
            // local_vec_visited_local[thread_id][v] = false;
            // for(int i=0; i<local_vec_visited_local.size(); i++){
            //     local_vec_visited_local[i][v] = false;
            // }
            // m[u] = UNMATCHED;
        }
    }



    }
// #else
// Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);}

    // #endif
   

    // }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
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
inline bool Prallel_Graphflow_FindMatches(
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

    // const auto& u_min_nbrs = data_.GetNeighbors(m[u_min]);
    // const uint v = u_min_nbrs[i];

    // const auto& u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    // const auto& q_nbrs = query_.GetNeighbors(u);
    // const auto& q_nbr_labels = query_.GetNeighborLabels(u);

    // // 1. Check vertex and edge labels
    // if (data_.GetVertexLabel(v) != query_.GetVertexLabel(u) || 
    //     u_min_nbr_labels[i] != u_min_label) {
    //     return false;
    // }

    // // 2. Check joinability with existing mappings
    // bool joinable = true;
    // for (uint j = 0; j < q_nbrs.size(); ++j) {
    //     const uint u_other = q_nbrs[j];
    //     const uint u_other_label = q_nbr_labels[j];

    //     if (m[u_other] == UNMATCHED || u_other == u_min) continue;

    //     const auto& nbrs = data_.GetNeighbors(m[u_other]);
    //     auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
        
    //     if (it == nbrs.end() || *it != v) {
    //         joinable = false;
    //         break;
    //     }

    //     uint dis = std::distance(nbrs.begin(), it);
    //     if (data_.GetNeighborLabels(m[u_other])[dis] != u_other_label) {
    //         joinable = false;
    //         break;
    //     }
    // }
    // if (!joinable) return false;

    // // 3. Check visited (for isomorphism)
    // if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // // 4. Update mappings and recurse
    // m[u] = v;
    // local_vec_visited_local[thread_id][v] = true;

    // if (depth == query_.NumVertices() - 1) {
    //     ++num_results;
    // } else {
    //     FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    // }

    // // Backtrack
    // local_vec_visited_local[thread_id][v] = false;
    // m[u] = UNMATCHED;

    return true;
}




/**
 * @brief Parallel version of the FindMatches function for handling deletions.
 * 
 * This function is designed to handle subgraph matching in parallel when edges are deleted
 * from the data graph. It clears the vertex and job queues, identifies the next query vertex
 * to process, and uses parallel threads to explore the search space.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. Clears the `vertex_vector` and `job_queue` to reset the state for processing deletions.
 * 2. Identifies the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 3. Iterates over candidate data vertices for the selected query vertex.
 * 4. Uses OpenMP to parallelize the exploration of the search space.
 * 5. Updates the match and recursively explores further matches.
 * 6. Aggregates results from all threads and updates the global match count.
 */
// Function implementation for Parallel_Symbi_FindMatches
void Parallel_Symbi_FindMatches(
    uint depth, 
    std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, 
    size_t& num_results,
    Graph& query_graph,
    Graph& data_graph,
    std::vector<std::vector<uint>>& eidx,
    std::vector<ska::flat_hash_map<uint, std::vector<uint>>>& DCS,
    std::vector<ska::flat_hash_map<uint, bool>>& d2,
    std::vector<ska::flat_hash_map<uint, uint>>& n2,
    std::vector<std::vector<bool>>& local_vec_visited_local,
    bool homomorphism,
    size_t NUMTHREAD,
    std::vector< std::tuple<uint, uint, size_t, std::vector<uint>,
        std::vector<Parrllel_SymBi::ExtendableVertex>,  uint , uint> >& vertex_vector,
    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>,
    std::vector<Parrllel_SymBi::ExtendableVertex>,  uint , uint> >& job_queue,
    std::vector<std::vector<uint>>& treeNode_neighbors


)
{

    vertex_vector.clear();
    job_queue.clear();

    size_t NUMT = NUMTHREAD;

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    // extendable.size() == query_graph.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        // 获取最小邻居的节点
        if (extendable[i].matched_nbrs == query_graph.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]
    // bool candidate_empty = true;

    // Use the flat_hash_map access pattern
    auto it = DCS[eidx[u_min][u]].find(m[u_min]);
    if (it == DCS[eidx[u_min][u]].end()) {
        return; // No candidates found
    }
    
    const auto& candidates = it->second;
    size_t total_size = candidates.size();
     
    std::vector<size_t> local_num_result(NUMT,0);


    for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        auto v = candidates[v_idx];

        // 1. Check index
        auto d2_it = d2[u].find(v);
        if (d2_it == d2[u].end() || !d2_it->second) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_neighbors[u]) {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto dcs_it = DCS[eidx[u_other][u]].find(m[u_other]);
            if (dcs_it == DCS[eidx[u_other][u]].end()) {
                joinable = false;
                break;
            }
            
            const auto& neighbor_candidates = dcs_it->second;
            auto it = std::lower_bound(neighbor_candidates.begin(), 
                            neighbor_candidates.end(), v);
            if (it == neighbor_candidates.end() || *it != v) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        size_t thread_id = omp_get_thread_num();
        if (!homomorphism && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        m[u] = v;
        // for(int i=0; i<local_vec_visited_local.size(); i++){
        //     local_vec_visited_local[i][v] = true;
        // }
        // local_vec_visited_local[thread_id][v] = true; // imp:
        #ifdef DEBUG
        std::cout << "thread_id: " << thread_id << std::endl;
        std::cout << "u: " << u << " v: " << v << std::endl;
        #endif

        std::vector<Parrllel_SymBi::ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_neighbors[u]) {
            if (m[u_other] != UNMATCHED) continue;

            auto n2_it = n2[eidx[u][u_other]].find(m[u]);
            if (n2_it != n2[eidx[u][u_other]].end() && 
                n2_it->second < temp_extendable[u_other].E) {
                temp_extendable[u_other].E = n2_it->second;
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs++;
        }

        if (depth == query_graph.NumVertices() - 1) { // match complete
            num_results++;
        } else {

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_graph.GetNeighbors(i).size())
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
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            auto dcs_it2 = DCS[eidx[u_min2][u2]].find(m[u_min2]);
            if (dcs_it2 == DCS[eidx[u_min2][u2]].end()) {
                continue; // No candidates found
            }
            size_t total_size2 = dcs_it2->second.size();

            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                // process_vertex(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);

                // queue version
                // vertex_queue.emplace(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                
                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, temp_extendable,  depth2, v);

            
            }
        }

        // this is for backtrack
        // local_vec_visited_local[thread_id][v] = false;
        // for(int i=0; i<local_vec_visited_local.size(); i++){
        //     local_vec_visited_local[i][v] = false;
        // }
        m[u] = UNMATCHED;

    }


    if(vertex_vector.size() < NUMT){
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }
    // std::cout << "wtf" << NUMT << std::endl;

    #pragma omp parallel for num_threads(NUMT) schedule(dynamic, 1)
    for (size_t i = 0; i < vertex_vector.size(); ++i) {
        auto [u, u_min, v_idx, m, extendable,  depth, v] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();

        local_vec_visited_local[thread_id][v] = true;
        m[u] = v;
        process_vertex_layer1_simplified(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id,
            query_graph, data_graph, eidx, DCS, d2, n2, local_vec_visited_local, homomorphism, treeNode_neighbors, job_queue);

        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;

        if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
            std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                 std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto [u, u_min, v_idx, m, extendable,  depth, v] = job;
                local_vec_visited_local[thread_id][v] = true;
                m[u] = v;
                // FindMatches_task(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                process_vertex_layer1_simplified(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id,
                    query_graph, data_graph, eidx, DCS, d2, n2, local_vec_visited_local, homomorphism, treeNode_neighbors, job_queue);
                local_vec_visited_local[thread_id][v] = false; 
                m[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

}



/**
 * @brief Parallel version of the FindMatches function for handling deletions.
 * 
 * This function is designed to handle subgraph matching in parallel when edges are deleted
 * from the data graph. It clears the vertex and job queues, identifies the next query vertex
 * to process, and uses parallel threads to explore the search space.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. Clears the `vertex_vector` and `job_queue` to reset the state for processing deletions.
 * 2. Identifies the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 3. Iterates over candidate data vertices for the selected query vertex.
 * 4. Uses OpenMP to parallelize the exploration of the search space.
 * 5. Updates the match and recursively explores further matches.
 * 6. Aggregates results from all threads and updates the global match count.
 */
inline void Parallel_Turboflux_FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{

    // size_t NUMT = NUMTHREAD;

    // vertex_vector.clear();
    // job_queue.clear();

    // uint u = order_vs_[order_index][depth];
    // uint u_min = backward_vs_[order_index][depth];


    // std::vector<size_t> local_num_result(NUMT+2,0);

    // // enumerate each neighbor of m[u_min]
    // for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    // {
    //     // 1. check index
    //     if (d2[u][v] == 0) continue;

    //     // 2. check if joinable
    //     bool joinable = true;
    //     for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
    //     {
    //         const auto& u_backward = join_check_vs_[order_index][u][i];
    //         const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
    //         const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

    //         if (std::get<2>(d_elabel) != u_backward_elabel)
    //         {
    //             joinable = false;
    //             break;
    //         }
    //     }
    //     if (!joinable) continue;

    //     // 3. check if visited
    //     if (!homomorphism_ && visited_[v]) continue;
        
    //     // 4. add a vertex mapping
    //     m[u] = v;
    //     for(size_t i = 0; i< local_vec_visited_local.size(); i++){
    //         local_vec_visited_local[i][v] = true;
    //     }
    //     visited_[v] = true;

    //     if (depth == query_.NumVertices() - 1)
    //     {
    //         num_results++;
    //     }
    //     else
    //     {
    //         // FindMatches(order_index, depth + 1, m, num_results);
    //         auto depth2 = depth + 1;
    //         uint u2 = order_vs_[order_index][depth2];
    //         uint u_min2 = backward_vs_[order_index][depth2];

    //         auto siez_DCS = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();


    //         for (auto v_idx2 = 0; v_idx2 < siez_DCS; ++v_idx2)
    //         {

    //             vertex_vector.emplace_back(u2, u_min2, v_idx2, m, 
    //                  v);

    //         }
    //     }

    //     for(size_t i = 0; i< local_vec_visited_local.size(); i++){
    //         local_vec_visited_local[i][v] = false;
    //     }
        
    //     visited_[v] = false;
    //     m[u] = UNMATCHED;

    // }

    // if(vertex_vector.size() < NUMT){
    //     if(vertex_vector.size() == 0){
    //         NUMT = 1;
    //     }
    // }

    // #pragma omp parallel for num_threads(NUMT)
    // for(size_t i =0; i < vertex_vector.size(); ++i){
    //     auto& [u3, u_min2, v_idx2, m2, v3] = vertex_vector[i];
    //     size_t thread_id = omp_get_thread_num();
    //     // visited_[v3] = true;
    //     // m2[u3] = v3;
    //     local_vec_visited_local[thread_id][v3] = true;
    //     m2[u3] = v3;

    //     ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
    //         depth + 1, order_index, thread_id);

    //     // visited_[v3] = false;
    //     // m2[u3] = UNMATCHED;
    //     local_vec_visited_local[thread_id][v3] = false;
    //     m2[u3] = UNMATCHED;

    //     if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
    //         std::tuple<uint, uint, size_t, std::vector<uint>,
    //              uint , uint> job;
    //         if(job_queue.try_pop(job)){
    //             size_t thread_id = omp_get_thread_num();
    //             auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
    //             local_vec_visited_local[thread_id][v3] = true;
    //             m2[u] = v3;
    //             ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
    //                 depth3, order_index, thread_id);
    //             local_vec_visited_local[thread_id][v3] = false; 
    //             m2[u] = UNMATCHED;
    //         }
    //     }
    // }

    // for (size_t i = 0; i < local_num_result.size(); ++i) {
    //     num_results += local_num_result[i];
    // }

}



// uint Parallel_CaLiG_FindMatches(uint v1, uint v2){

    // uint update_result = 0;
    
    // // Map
    // std::vector<std::tuple<uint, uint, uint>> tasks; // (u1, u2, v2)
    
    // for(auto& li : G[v1].LI) {
    //     if(li.second) {
    //         uint u1 = li.first;
    //         for(auto& candi : G[v1].cand[u1]) {
    //             uint u2 = candi.first;
    //             if(u2 >= 0 && candi.second.find(v2) != candi.second.end()) {
    //                 tasks.push_back(std::make_tuple(u1, u2, v2));
    //             }
    //         }
    //     }
    // }
    
    // // reduce
    // #pragma omp parallel for reduction(+:update_result)
    // for(int i = 0; i < tasks.size(); i++) {
    //     uint u1 = std::get<0>(tasks[i]);
    //     uint u2 = std::get<1>(tasks[i]);
    //     uint v2_val = std::get<2>(tasks[i]);
        
    //     ska::flat_hash_map<uint, uint> matching;
    //     matching[u1] = v1;
    //     matching[u2] = v2;
    //     vec core_v;
    //     core_v.emplace_back(v1);
    //     core_v.emplace_back(v2);
        
    //     uint local_result = searchCore(0, matching, core_v, matching_order[u1][u2].core,
    //                                   matching_order[u1][u2].core_nei, matching_order[u1][u2].shell,
    //                                   matching_order[u1][u2].shell_nei, matching_order[u1][u2].c_s_nei);
        
    //     update_result += local_result;
    // }
    
    // return update_result;

// }


// here
void Parrllel_SymBi::Parallel_Symbi_FindMatches_ParaCOSM_Kernel(uint depth, std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results){

        vertex_vector.clear();
        job_queue.clear();
    
        size_t NUMT = NUM_THREAD;
    
        uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
        uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;
    
        uint u;
    
        // extendable.size() == query_.NumVertices()
        for (uint i = 0; i < extendable.size(); i++)
        {
            if (m[i] != UNMATCHED) continue;
    
            // Check if an extendable query vertex is isolated or not
            if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
            {
                if (extendable[i].E < isolate_minE)
                {
                    isolate_minE = extendable[i].E; // local
                    isolate_u = i; // local
                }
            }
            else
            {
                if (extendable[i].E < non_isolate_minE)
                {
                    non_isolate_minE = extendable[i].E; // local
                    non_isolate_u = i; // local
                }
            }
        }
    
        if (non_isolate_minE == NOT_EXIST)
            u = isolate_u; // local
        else
            u = non_isolate_u; // local
            
        uint u_min = extendable[u].u_min;
        extendable[u] = {}; // local
    
        size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();
         
        std::vector<size_t> local_num_result(NUMT,0);
    
        for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
        {
            __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
            auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];
    
            // 1. Check index
            if (d2[u][v] == 0) continue;
    
            // 2. Check if joinable
            bool joinable = true;
            for (auto& u_other: treeNode_[u].neighbors_) {
                if (m[u_other] == UNMATCHED || u_other == u_min) continue;
    
                auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), 
                                DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
                if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
                    joinable = false;
                    break;
                }
            }
    
            if (!joinable) continue;
    
            // 3. Check if visited
            size_t thread_id = omp_get_thread_num();
            if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;
    
            // 4. Add a vertex mapping
            m[u] = v;
            for(int i=0; i<local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v] = true;
            }
    
            #ifdef DEBUG
            std::cout << "thread_id: " << thread_id << std::endl;
            std::cout << "u: " << u << " v: " << v << std::endl;
            #endif
    
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
            
                // Enumerate each neighbor of m[u_min]
                // bool candidate_empty = true;
            
                size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();
    
                // std::cout << "total_size2: " << total_size2 << std::endl;
                for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
                {
    
                    vertex_vector.emplace_back(u2, u_min2, v_idx2, m, temp_extendable,  depth2, v);

                }
            }
    
            // this is for backtrack
            // local_vec_visited_local[thread_id][v] = false;
            for(int i=0; i<local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v] = false;
            }
            m[u] = UNMATCHED;
    
        }
    
    
        if(vertex_vector.size() < NUMT){
            if(vertex_vector.size() == 0){
                NUMT = 1;
            }
        }
    
        // std::cout << "NUMT"<< NUMT << std::endl;
    // SOLUTION 1
        #pragma omp parallel for num_threads(NUMT) schedule(auto)
        for (size_t i = 0; i < vertex_vector.size() ; ++i) {
            auto [u3, u_min, v_idx, m, extendable,  depth3, v3] = vertex_vector[i];
            size_t thread_id = omp_get_thread_num();
    
            local_vec_visited_local[thread_id][v3] = true;
            m[u3] = v3;
            process_vertex_layer_local(u3, u_min, v_idx, m, extendable, local_num_result[thread_id], depth3, thread_id);
    
            local_vec_visited_local[thread_id][v3] = false;
            m[u3] = UNMATCHED;
    
            if(!job_queue.empty() 
            // && (i > vertex_vector.size() - NUMT)
            ){
                std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                     std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
                if(job_queue.try_pop(job)){
                    size_t thread_id = omp_get_thread_num();
                    auto [u, u_min, v_idx, m, extendable,  depth, v] = job;
                    local_vec_visited_local[thread_id][v] = true;
                    m[u] = v;
                    process_vertex_layer_local(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                    local_vec_visited_local[thread_id][v] = false; 
                    m[u] = UNMATCHED;
                }
            }

        }
    
        #ifdef TASK_SPILT
        #pragma omp parallel num_threads(NUMT)
        {
            size_t thread_id = omp_get_thread_num();
            std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                       std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
        
            while (true) {
                // take a job from the queue
                if (!job_queue.try_pop(job)) {
                    // if empty, go out
                    if (job_queue.empty()) {
                        break;
                    }
                    continue; // find a job, go
                }
        
                // deal with the job
                auto [u, u_min, v_idx, m, extendable, depth, v] = job;
                local_vec_visited_local[thread_id][v] = true;
                m[u] = v;
                process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                local_vec_visited_local[thread_id][v] = false;
                m[u] = UNMATCHED;
            }
        }
        #endif
    
        for (size_t i = 0; i < local_num_result.size(); ++i) {
            num_results += local_num_result[i];
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
 void Parallel_Graphflow::Parallel_Graphflow_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)
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
 inline bool Parallel_TurboFlux::ProcessVertex_TurboFlux_ParaCOSM_Kernel(uint u, uint u_min, size_t v_idx, 
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
        FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}




 /**
 * @brief Parallel version of the FindMatches function for handling deletions.
 * 
 * This function is designed to handle subgraph matching in parallel when edges are deleted
 * from the data graph. It clears the vertex and job queues, identifies the next query vertex
 * to process, and uses parallel threads to explore the search space.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. Clears the `vertex_vector` and `job_queue` to reset the state for processing deletions.
 * 2. Identifies the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 3. Iterates over candidate data vertices for the selected query vertex.
 * 4. Uses OpenMP to parallelize the exploration of the search space.
 * 5. Updates the match and recursively explores further matches.
 * 6. Aggregates results from all threads and updates the global match count.
 */
void Parallel_TurboFlux::Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{

    size_t NUMT = NUMTHREAD;

    vertex_vector.clear();
    job_queue.clear();

    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];


    std::vector<size_t> local_num_result(NUMT+2,0);

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
        if (!homomorphism_ && visited_[v]) continue;
        
        // 4. add a vertex mapping
        m[u] = v;
        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = true;
        }
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
            uint u_min2 = backward_vs_[order_index][depth2];

            auto siez_DCS = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();


            for (auto v_idx2 = 0; v_idx2 < siez_DCS; ++v_idx2)
            {

                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, 
                     v);

            }
        }

        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = false;
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;

    }

    if(vertex_vector.size() < NUMT){
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }

    #pragma omp parallel
    {

        bool all_finished = false;
        size_t thread_id = omp_get_thread_num();

        size_t base = vertex_vector.size() / NUMT;
        size_t rem  = vertex_vector.size() % NUMT;

        size_t start, end;
        if (thread_id < rem) {
            start = thread_id * (base + 1);
            end   = start + (base + 1);
        } else {
            start = rem * (base + 1) + (thread_id - rem) * base;
            end   = start + base;
        }

        if (start > vertex_vector.size()) start = vertex_vector.size();
        if (end   > vertex_vector.size()) end   = vertex_vector.size();

        for(size_t i = start; i < end; ++i){ 
            auto& [u3, u_min2, v_idx2, m2,  v3] = vertex_vector[i];
            local_vec_visited_local[thread_id][v3] = true;
            m2[u3] = v3;
            // search + share jobs if possible
            ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                depth + 1, order_index, thread_id);
        }

        // finish its job
        workers_free[thread_id] = true;

        while (!all_finished)
        {
            std::tuple<uint, uint, size_t, std::vector<uint>,
                    uint , uint> job;
            while (!job_queue.empty())
            {
                if(job_queue.try_pop(job)){
                    size_t thread_id = omp_get_thread_num();
                    auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
                    local_vec_visited_local[thread_id][v3] = true;
                    m2[u3] = v3;
                    ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                            depth3, order_index, thread_id);
                    local_vec_visited_local[thread_id][v3] = false; 
                    m2[u3] = UNMATCHED;
                }
            }
            
            int finished_count = 0;
            for(int i = 0; i < NUMT; ++i){
                if(workers_free[i] == true){
                    finished_count++;
                }
            }
            if(finished_count == NUMT){
                all_finished = true;
            }
        }
        
        std::tuple<uint, uint, size_t, std::vector<uint>,
                    uint , uint> job;
        while (!job_queue.empty())
        {
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
                local_vec_visited_local[thread_id][v3] = true;
                m2[u3] = v3;
                ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                            depth3, order_index, thread_id);
                local_vec_visited_local[thread_id][v3] = false; 
                m2[u3] = UNMATCHED;
            }
        }
        
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

    // #pragma omp parallel for num_threads(NUMT)
    // for(size_t i =0; i < vertex_vector.size(); ++i){
    //     auto& [u3, u_min2, v_idx2, m2, v3] = vertex_vector[i];
    //     size_t thread_id = omp_get_thread_num();
    //     // visited_[v3] = true;
    //     // m2[u3] = v3;
    //     local_vec_visited_local[thread_id][v3] = true;
    //     m2[u3] = v3;

    //     ProcessVertex_TurboFlux_ParaCOSM_Kernel(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
    //         depth + 1, order_index, thread_id);

    //     // visited_[v3] = false;
    //     // m2[u3] = UNMATCHED;
    //     local_vec_visited_local[thread_id][v3] = false;
    //     m2[u3] = UNMATCHED;

    //     if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
    //         std::tuple<uint, uint, size_t, std::vector<uint>,
    //              uint , uint> job;
    //         if(job_queue.try_pop(job)){
    //             size_t thread_id = omp_get_thread_num();
    //             auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
    //             local_vec_visited_local[thread_id][v3] = true;
    //             m2[u] = v3;
    //             ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
    //                 depth3, order_index, thread_id);
    //             local_vec_visited_local[thread_id][v3] = false; 
    //             m2[u] = UNMATCHED;
    //         }
    //     }
    // }

    // for (size_t i = 0; i < local_num_result.size(); ++i) {
    //     num_results += local_num_result[i];
    // }

}
















// void CSMPP::Parallel_searchVertex(uint queryIndex, uint edgeIndex, searchType type, 
//     uint depth){

//     auto & QUERY_GRAPH = this->queryVec[queryIndex]; // take one query graph
//     const uint cacheStatu = QUERY_GRAPH.getCacheStatu(edgeIndex, depth);
//     const auto & MATCHORDER = QUERY_GRAPH.GetMatchOrder(edgeIndex);
//     const uint queryVexter = MATCHORDER[depth];
//     vertexType currentSearchVertexType = QUERY_GRAPH.getVertexType(edgeIndex, depth);
//     const auto & desItem = QUERY_GRAPH.GetDescList(edgeIndex, depth);//<v1Index,v1Label,eLabel>
//     const auto & freezeIndex = QUERY_GRAPH.getUnfreezeList(edgeIndex, depth);

//     //  ///////////////////////////////////////////////////////////////
    
//     uint vertexLabel = QUERY_GRAPH.GetVertexLabel(queryVexter);

//     //1.intersection worst case
//     //1.1 find min
//     std::vector<uint> candidate;
//     LRAndIndexCheckType decision = QUERY_GRAPH.getDecision(edgeIndex, depth); // adaptive label distribute check

//     uint min_u_index = INT_MAX;
//     uint min_u_neighborSize = INT_MAX;

//     for(int i = 0; i < desItem.size(); i++){
//         const auto & item = desItem[i];

//         uint v = this->match[std::get<0>(item)]; // gaigaigai

//         uint Size = this->data_.getIndexValue(v, vertexLabel);
//         if(Size < min_u_neighborSize){
//             min_u_index = i;
//             min_u_neighborSize = Size;
//         }
//         else if(Size == min_u_neighborSize){
//             uint preID = this->match[std::get<0>(desItem[min_u_index])]; // gaigaigai

//             if(this->data_.GetNeighbors(preID).size() > this->data_.GetNeighbors(v).size()){
//                 min_u_index = i;
//                 min_u_neighborSize = Size;
//             }
//         }
//     }





//     uint min_u = this->match[std::get<0>(desItem[min_u_index])];


//     uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel

//     const auto& MinNeighbor = data_.GetNeighbors(min_u);
//     const auto& q_nbr_labels = data_.GetNeighborLabels(min_u);
//     // std::cout << NUM_L << std::endl;
//     // NUM_L++;8

    
//     for(int i = 0; i < MinNeighbor.size(); i++){
//         const uint v = MinNeighbor[i];
//         //1. check label
//         if( this->data_.GetVertexLabel(v) != vertexLabel || q_nbr_labels[i] != min_u_elabel)continue;

//         //2. check visit
//         if(this->visited_[v] == true && !homomorphism_ ) continue;
        
//         // std::cout << NUM_L << std::endl;
//         // NUM_L++; // 369

//         //3.check if joinable
//         bool joinable = true;
//         for(int k = 0; k < desItem.size(); k++){
//             if(k == min_u_index) continue;

//             const uint data_V = this->match[std::get<0>(desItem[k])];

//             const uint elabel = std::get<2>(desItem[k]);
//             const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
//             auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
//             if(it == dataVNeighbor.end() || *it != v){
//                 joinable = false;
//                 break;
//             }
//             else
//             {
//                 uint dis = std::distance(dataVNeighbor.begin(), it);
//                 if(this->data_.GetNeighborLabels(data_V)[dis] != elabel){
//                     joinable = false;
//                     break;
//                 }
//             }
//         }

//         // std::cout << NUM_L << std::endl;
//         // NUM_L++; // 369

//         if(!joinable){
//             continue;
//         }
//         // std::cout << NUM_L << std::endl;
//         // NUM_L++; 369
//         if(decision == Part1Check && !this->indexCheck(v, queryVexter, queryIndex)){
//             continue;
//         }

//         // std::cout << NUM_L << std::endl;
//         // NUM_L++; 367
        
//         candidate.emplace_back(v);
//     }
        
//     if(candidate.size() == 0){
//         return;
//     }
    
//     // std::cout << NUM_L << std::endl;
//     // NUM_L++; 7



//     if(currentSearchVertexType == freeVertex){
//             // if(candidate.size() > 32){
//             //     std::cout << candidate.size()  << std::endl;    
//             // }
//         size_t NUMT = 64;
//         if(candidate.size() < Thread_MAX){
//             NUMT = candidate.size();
//         }else{
//             NUMT = Thread_MAX;
//         }

//         #pragma omp parallel for num_threads(NUMT) shcedule(dynamic, 1)
//         for(int i =0; i < candidate.size(); i++){
//             size_t thread_id = omp_get_thread_num();
//             matchVertexlocal(candidate[i], depth, match_parallel[thread_id], visited_parallel[thread_id]);
//             auto dataV = candidate[i];
//             // this->matchVertex(dataV, depth);

//             // this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
//             Parallel_searchVertex_local(queryIndex, edgeIndex, type, depth + 1, thread_id);

//             // this->popVertex(dataV, depth);
//             popVertex_local(candidate[i], depth, match_parallel[thread_id], visited_parallel[thread_id]);
//             // ma
//         }

//         // for(auto & dataV : candidate){
//         //     this->matchVertex(dataV, depth);
//         //     this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
//         //     this->popVertex(dataV, depth);
//         // }
//     }
//     else{
        
//         this->matchVertex(candidate, depth);
//         this->matchVertexAll(candidate, depth);
//         // this->match_parallel[0] = this->match;
        

//         if(currentSearchVertexType == isolatedVertex){
//             this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[depth]);
//         }
//         if(depth == this->queryVec[queryIndex].NumVertices() - 1){
//             this->addMatchResult(queryIndex, edgeIndex, type);
//         }
//         else{

//             // Parallel_searchVertex(queryIndex, edgeIndex, type, depth + 1);
            
//             size_t thread_id = omp_get_thread_num();
//             Parallel_searchVertex_local(queryIndex, edgeIndex, type, depth + 1, thread_id);
            
//             // this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
//         }
//         if(currentSearchVertexType == isolatedVertex){
//             this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[depth]);
//         }
//         this->popVertex(depth);
//         this->popVertexAll(depth);
//     }
// }








// void CSMPP::Parallel_searchVertex_local(uint queryIndex, uint edgeIndex, searchType type, 
//     uint depth, size_t thread_id){

//     auto & QUERY_GRAPH = this->queryVec[queryIndex]; // take one query graph
//     const uint cacheStatu = QUERY_GRAPH.getCacheStatu(edgeIndex, depth);
//     const auto & MATCHORDER = QUERY_GRAPH.GetMatchOrder(edgeIndex);
//     const uint queryVexter = MATCHORDER[depth];
//     vertexType currentSearchVertexType = QUERY_GRAPH.getVertexType(edgeIndex, depth);
//     const auto & desItem = QUERY_GRAPH.GetDescList(edgeIndex, depth);//<v1Index,v1Label,eLabel>
//     const auto & freezeIndex = QUERY_GRAPH.getUnfreezeList(edgeIndex, depth);

//     //  ///////////////////////////////////////////////////////////////
    
//     uint vertexLabel = QUERY_GRAPH.GetVertexLabel(queryVexter);


//     //1.intersection worst case
//     //1.1 find min
//     std::vector<uint> candidate;
//     LRAndIndexCheckType decision = QUERY_GRAPH.getDecision(edgeIndex, depth); // adaptive label distribute check

//     uint min_u_index = INT_MAX;
//     uint min_u_neighborSize = INT_MAX;



//     for(int i = 0; i < desItem.size(); i++){
//         const auto & item = desItem[i];

//         // uint v = this->match[std::get<0>(item)]; // gaigaigai
//         uint v = this->match_parallel[thread_id][std::get<0>(item)]; // gaigaigai

//         uint Size = this->data_.getIndexValue(v, vertexLabel);
//         if(Size < min_u_neighborSize){
//             min_u_index = i;
//             min_u_neighborSize = Size;
//         }
//         else if(Size == min_u_neighborSize){
//             // uint preID = this->match[std::get<0>(desItem[min_u_index])]; // gaigaigai
//             uint preID = this->match_parallel[thread_id][std::get<0>(desItem[min_u_index])];

//             if(this->data_.GetNeighbors(preID).size() > this->data_.GetNeighbors(v).size()){
//                 min_u_index = i;
//                 min_u_neighborSize = Size;
//             }
//         }
//     }

//     // std::cout << NUM_L << std::endl;
//     // NUM_L++; //201
//     // std::cout << "min_u_index: " << min_u_index << std::endl;

//     // uint min_u = this->match[std::get<0>(desItem[min_u_index])];
//     if(desItem.size()==0){
//         return;
//     }

//     uint min_u = this->match_parallel[thread_id][std::get<0>(desItem[min_u_index])];

//     uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel

//     const auto& MinNeighbor = data_.GetNeighbors(min_u);
//     const auto& q_nbr_labels = data_.GetNeighborLabels(min_u);

//     // std::cout << NUM_L << std::endl;
//     // NUM_L++; //200
    
//     for(int i = 0; i < MinNeighbor.size(); i++){
//         const uint v = MinNeighbor[i];
//         //1. check label
//         if( this->data_.GetVertexLabel(v) != vertexLabel || q_nbr_labels[i] != min_u_elabel)continue;

//         //2. check visit
//         // if(this->visited_[v] == true && !homomorphism_ ) continue;
//         if(this->visited_parallel[thread_id][v] == true && !homomorphism_){continue;}

//         //3.check if joinable
//         bool joinable = true;
//         for(int k = 0; k < desItem.size(); k++){
//             if(k == min_u_index) continue;

//             // const uint data_V = this->match[std::get<0>(desItem[k])];
//             const uint data_V = this->match_parallel[thread_id][std::get<0>(desItem[k])];

//             const uint elabel = std::get<2>(desItem[k]);
//             const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
//             auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
//             if(it == dataVNeighbor.end() || *it != v){
//                 joinable = false;
//                 break;
//             }
//             else
//             {
//                 uint dis = std::distance(dataVNeighbor.begin(), it);
//                 if(this->data_.GetNeighborLabels(data_V)[dis] != elabel){
//                     joinable = false;
//                     break;
//                 }
//             }
//         }
//         if(!joinable){
//             continue;
//         }
//         if(decision == Part1Check && !this->indexCheck(v, queryVexter, queryIndex)){
//             continue;
//         }
        
//         candidate.emplace_back(v);
//     }
        
//     if(candidate.size() == 0){
//         return;
//     }

            
// // std::cout << "min_u_index: " << min_u_index << std::endl;

 
//     if(currentSearchVertexType == freeVertex){

//         for(int i =0; i < candidate.size(); i++){
            
//             // size_t thread_id = omp_get_thread_num();
//             matchVertexlocal(candidate[i], depth, match_parallel[thread_id], visited_parallel[thread_id]);
//             auto dataV = candidate[i];
            
//             this->Parallel_searchVertex_local(queryIndex, edgeIndex, type, depth + 1, thread_id);

//             popVertex_local(candidate[i], depth, match_parallel[thread_id], visited_parallel[thread_id]);
//         }

//         // for(auto & dataV : candidate){
//         //     this->matchVertex(dataV, depth);
//         //     this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
//         //     this->popVertex(dataV, depth);
//         // }
//     }
//     else{
//         // this->matchVertex(candidate, depth);
//         this->matchVertexlocal(candidate, depth, match_parallel[thread_id], matchCandidate_parallel[thread_id], thread_id);
        
//         if(currentSearchVertexType == isolatedVertex){
//             // this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[depth]);
//             // std::cout << NUM_L << std::endl;
//             // NUM_L++;
//             this->queryVec[queryIndex].isolatedVertexTimesAdd(matchCandidate_parallel[thread_id][depth]);
//             // std::cout << NUM_L << std::endl;
//             // NUM_L++;
//         }
//         if(depth == this->queryVec[queryIndex].NumVertices() - 1){
//             this->addMatchResult(queryIndex, edgeIndex, type);
//         }
//         else{

//             this->Parallel_searchVertex_local(queryIndex, edgeIndex, type, depth + 1, thread_id);
//             // std::cout << NUM_L << std::endl;
//             // NUM_L++;
//             // this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
//         }
//         if(currentSearchVertexType == isolatedVertex){
//             // this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[depth]);
//             this->queryVec[queryIndex].isolatedVertexTimesMinus(matchCandidate_parallel[thread_id][depth]);
//         }
//         // this->popVertex(depth);
//         this->popVertex_local(depth, match_parallel[thread_id], matchCandidate_parallel[thread_id], thread_id);
//     }
// }





// void turnOffProcess_parallel(uint v1, uint v2){

//     G[v1].nei.erase(v2);
//     G[v2].nei.erase(v1);
    
//     // v1 
//     auto cand_v1 = G[v1].cand;
    
//     #pragma omp parallel
//     {

//         std::vector<std::tuple<uint, uint, uint>> local_updates; // (ui, uj, v2)
//         std::vector<uint> local_erase_ui; // ui
        
//         #pragma omp for nowait
//         for (int i = 0; i < cand_v1.size(); i++) {
//             // 使用迭代器访问 map 的第 i 个元素
//             auto it = cand_v1.begin();
//             std::advance(it, i);
//             uint ui = it->first;
//             auto& candi = it->second;
            
//             if(G[v1].label != Q[ui].label){
//                 local_erase_ui.push_back(ui);
//                 continue;
//             }
            
//             for(auto& ui_nei : candi){
//                 uint uj = ui_nei.first;
//                 if(Q[uj].label == G[v2].label){
//                     local_updates.push_back(std::make_tuple(ui, uj, v2));
//                 }
//             }
//         }

//         #pragma omp critical
//         {
//             for(auto ui : local_erase_ui){
//                 G[v1].cand.erase(ui);
//                 G[v1].LI.erase(ui);
//             }
            
//             for(auto& update : local_updates){
//                 uint ui = std::get<0>(update);
//                 uint uj = std::get<1>(update);
//                 uint v_val = std::get<2>(update);
                
//                 G[v1].cand[ui][uj].erase(v_val);
//                 G[v2].cand[uj][ui].erase(v1);
//             }
//         }
//     }
    
//     turnOff(v1);
//     turnOff(v2);
// }


