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
void Parallel_Symbi_FindMatches(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{

    vertex_vector.clear();
    job_queue.clear();

    size_t NUMT = NUMTHREAD;

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    // extendable.size() == query_.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        // 获取最小邻居的节点
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

    // Enumerate each neighbor of m[u_min]
    // bool candidate_empty = true;

    size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();
     
    std::vector<size_t> local_num_result(NUMT,0);


    for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
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
        // for(int i=0; i<local_vec_visited_local.size(); i++){
        //     local_vec_visited_local[i][v] = true;
        // }
        // local_vec_visited_local[thread_id][v] = true; // imp:
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
        process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);

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
                process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
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

    #pragma omp parallel for num_threads(NUMT)
    for(size_t i =0; i < vertex_vector.size(); ++i){
        auto& [u3, u_min2, v_idx2, m2, v3] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();
        // visited_[v3] = true;
        // m2[u3] = v3;
        local_vec_visited_local[thread_id][v3] = true;
        m2[u3] = v3;

        ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
            depth + 1, order_index, thread_id);

        // visited_[v3] = false;
        // m2[u3] = UNMATCHED;
        local_vec_visited_local[thread_id][v3] = false;
        m2[u3] = UNMATCHED;

        if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
            std::tuple<uint, uint, size_t, std::vector<uint>,
                 uint , uint> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
                local_vec_visited_local[thread_id][v3] = true;
                m2[u] = v3;
                ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                    depth3, order_index, thread_id);
                local_vec_visited_local[thread_id][v3] = false; 
                m2[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

}



uint Parallel_CaLiG_FindMatches(uint v1, uint v2){

    uint update_result = 0;
    
    // Map
    std::vector<std::tuple<uint, uint, uint>> tasks; // (u1, u2, v2)
    
    for(auto& li : G[v1].LI) {
        if(li.second) {
            uint u1 = li.first;
            for(auto& candi : G[v1].cand[u1]) {
                uint u2 = candi.first;
                if(u2 >= 0 && candi.second.find(v2) != candi.second.end()) {
                    tasks.push_back(std::make_tuple(u1, u2, v2));
                }
            }
        }
    }
    
    // reduce
    #pragma omp parallel for reduction(+:update_result)
    for(int i = 0; i < tasks.size(); i++) {
        uint u1 = std::get<0>(tasks[i]);
        uint u2 = std::get<1>(tasks[i]);
        uint v2_val = std::get<2>(tasks[i]);
        
        ska::flat_hash_map<uint, uint> matching;
        matching[u1] = v1;
        matching[u2] = v2;
        vec core_v;
        core_v.emplace_back(v1);
        core_v.emplace_back(v2);
        
        uint local_result = searchCore(0, matching, core_v, matching_order[u1][u2].core,
                                      matching_order[u1][u2].core_nei, matching_order[u1][u2].shell,
                                      matching_order[u1][u2].shell_nei, matching_order[u1][u2].c_s_nei);
        
        update_result += local_result;
    }
    
    return update_result;

}