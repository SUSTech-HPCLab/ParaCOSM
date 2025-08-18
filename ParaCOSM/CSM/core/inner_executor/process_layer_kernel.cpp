
















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
 inline bool Parallel_TurboFlux::ProcessVertex(uint u, uint u_min, size_t v_idx, 
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
 inline bool Parallel_TurboFlux::ProcessVertex_queue(uint u, uint u_min, size_t v_idx, 
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
        size_t total_size2 = DCS_[eidx_[u_min][u]][m[u_min]].size();
        for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
        {            
            job_queue.push(std::make_tuple(u, u_min, v_idx2, m, depth+1, v));
        }
        // FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}