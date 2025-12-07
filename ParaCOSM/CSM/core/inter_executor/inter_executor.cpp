#include <vector>
#include <atomic>
#include <iostream>
#include "inter_executor.h"
#include "matching_executor/matching.h"
#include "graph_storage/graph.h"

// Constructor implementation
/**
 * Construct an InterExecutor bound to a data graph and a matching engine.
 * - Does not take ownership of the matching instance.
 * - Thread-safety: construction is single-threaded; the instance is not thread-safe by itself.
 */
InterExecutor::InterExecutor(Graph& data_graph, matching* matching_instance)
    : data_graph_(data_graph), matching_instance_(matching_instance) {}

/**
 * BatchUpdates
 *
 * Windowed classification over a large batch popped from the queue `data_graph_.updates_`.
 * The method classifies edges in parallel using OpenMP, finds the earliest unsafe index
 * within each sliding window, then applies the corresponding window of edge insertions.
 *
 * Parameters:
 * - num_v_updates: cumulative counter of vertex updates (not incremented here; only edges handled)
 * - num_e_updates: cumulative counter of edge updates (incremented by per-item counters)
 * - unsafe_updates: increments when applying a window changes the result counters
 * - count: external total update counter (not modified here)
 * - positive_num_results_last / negative_num_results_last: previous result counters to detect changes
 * - reach_time_limit: stop early if set to true
 *
 * Concurrency:
 * - Uses OpenMP parallel regions for classification and window application.
 * - Shares `safe` and `idx_update` across threads; uses critical sections when updating them.
 */
void InterExecutor::BatchUpdates(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    const size_t batch_size = 80000;
    const size_t sliding_widow_size = 200;
    const size_t classify_threads = 8;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);

    std::vector<bool> is_safes(batch_size, false);
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);

    bool safe = true;
    size_t idx_update = batch_size + 5;
    size_t main_size = batch_size / sliding_widow_size;

    while (!data_graph_.updates_.empty() && !reach_time_limit) {
        batch_updates.clear();
        is_safes.assign(batch_size, false);
        num_e_update_vec.assign(batch_size, 0ul);

        size_t actual_threads = std::min(classify_threads, sliding_widow_size);

        #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
        {
            #pragma omp single
            {
                for (size_t i = 0; i < batch_size && !data_graph_.updates_.empty(); ++i) {
                    batch_updates.push_back(data_graph_.updates_.front());
                    data_graph_.updates_.pop();
                }
            }

            for (size_t local_batch = 0; local_batch < main_size; ++local_batch) {
                #pragma omp for
                for (size_t local_window = 0; local_window < sliding_widow_size; ++local_window) {
                    size_t cnt = local_batch * sliding_widow_size + local_window;
                    if (cnt >= batch_updates.size()) {
                        continue;
                    }

                    const auto& insert = batch_updates[cnt];
                    if (insert.type == 'e' && insert.is_add) {
                        is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                        num_e_update_vec[cnt]++;

                        if (!is_safes[cnt]) {
                            #pragma omp critical
                            {
                                safe = false;
                                if (idx_update >= cnt) {
                                    idx_update = cnt;
                                }
                            }
                            #pragma omp cancel for
                        }
                    }
                }

                #pragma omp barrier

                #pragma omp single
                {
                    if (!safe) {
                        for (size_t i = 0; i < sliding_widow_size; i++) {
                            if (i < batch_updates.size()) {
                                const auto& insert = batch_updates[local_batch * sliding_widow_size + i];
                                matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
                            }
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i < num_e_update_vec.size(); i++) {
            num_e_updates += num_e_update_vec[i];
        }
    }
}

/**
 * BatchUpdates2
 *
 * Queue-based smaller batch processing. Pulls up to `batch_size` updates from
 * `data_graph_.updates_`, classifies edges in parallel, and if an unsafe update is found,
 * applies all subsequent edge insertions in order. Tracks changes in result counters
 * to attribute "unsafe" updates.
 *
 * Parameters:
 * - num_v_updates: cumulative counter of vertex updates (not incremented here)
 * - num_e_updates: cumulative counter of edge updates (incremented when popping from queue)
 * - unsafe_updates: increments when results change after applying updates
 * - count: external total update counter (not modified here)
 * - positive_num_results_last / negative_num_results_last: previous result counters to detect changes
 * - reach_time_limit: stop early if set to true
 *
 * Concurrency:
 * - Uses OpenMP with a single section to pull from the queue and a for-section to classify.
 * - `safe` and `idx_update` are shared and guarded as needed.
 */
void InterExecutor::BatchUpdates2(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    const size_t batch_size = 400;
    const size_t classify_threads = 16;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);
    size_t idx_update = batch_size + 5;
    std::vector<bool> is_safes(batch_size, false);
    bool safe = true;
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);

    while (!data_graph_.updates_.empty() && !reach_time_limit) {
        batch_updates.clear();
        safe = true;
        is_safes.assign(batch_size, false);
        num_e_update_vec.assign(batch_size, 0ul);

        size_t actual_threads = std::min(classify_threads, batch_updates.size());

        #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
        {
            #pragma omp single
            {
                // get from queue
                for (size_t i = 0; i < batch_size && !data_graph_.updates_.empty(); ++i) {
                    batch_updates.push_back(data_graph_.updates_.front());
                    data_graph_.updates_.pop();
                    num_e_updates++;
                }
            }

            #pragma omp for
            for (size_t cnt = 0; cnt < batch_updates.size(); cnt++) {
                if (!safe) {
                    continue;
                }

                const auto& insert = batch_updates[cnt];
                if (insert.type == 'e' && insert.is_add) {
                    is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                    if (!is_safes[cnt]) {
                        {
                            safe = false;
                            if (idx_update >= cnt) {
                                #pragma omp critical
                                idx_update = cnt;
                            }
                        }
                    }
                }
            }
        }

        if (!safe) {
            for (size_t i = idx_update; i < batch_updates.size(); i++) {
                const auto& insert = batch_updates[i];
                matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);

                size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
                matching_instance_->GetNumPositiveResults(positive_num_results_cur);
                matching_instance_->GetNumNegativeResults(negative_num_results_cur);
                if (positive_num_results_cur != positive_num_results_last ||
                    negative_num_results_cur != negative_num_results_last) {
                    positive_num_results_last = positive_num_results_cur;
                    negative_num_results_last = negative_num_results_cur;

                    unsafe_updates++;
                }
            }
        }
    }
}

/**
 * BatchUpdates_OpenMP
 *
 * OpenMP-based sliding window over updates_vec_ (vector form). Classifies updates within a
 * fixed-size window; if an unsafe update is found, apply it immediately; otherwise skip the window.
 */
void InterExecutor::BatchUpdates_OpenMP(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    size_t sliding_window_base = 0;
    const size_t window_size = 160;
    const size_t update_size = data_graph_.updates_vec_.size();

    std::vector<bool> update_safes(window_size + 1, false);
    size_t min_safe = window_size + 10;
    bool is_safes = true;

    while (sliding_window_base < update_size) {
        is_safes = true;
        min_safe = window_size + 10;

        #pragma omp parallel for num_threads(4)
        for (size_t i = 0; i < window_size; i++) {
            if (i + sliding_window_base >= update_size) {
                continue;
            }
            const auto& insert = data_graph_.updates_vec_[i + sliding_window_base];

            if (insert.type == 'e' && insert.is_add) {
                update_safes[i] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);

                if (!update_safes[i]) {
                    {
                        is_safes = false;
                        if (min_safe > i) {
                            min_safe = i;
                        }
                    }
                }
            } else if (insert.type == 'v' && insert.is_add) {
                matching_instance_->AddVertex(insert.id1, insert.label);
            } else if (insert.type == 'v' && !insert.is_add) {
                #pragma omp critical
                {
                    is_safes = false;
                    if (min_safe > i) {
                        min_safe = i;
                    }
                }
            } else if (insert.type == 'e' && !insert.is_add) {
                update_safes[i] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                if (!update_safes[i]) {
                    #pragma omp critical
                    {
                        is_safes = false;
                        if (min_safe > i) {
                            min_safe = i;
                        }
                    }
                }
            }
        }

        if (!is_safes) {
            const auto& insert_unsafe = data_graph_.updates_vec_[min_safe + sliding_window_base];

            if (insert_unsafe.type == 'e' && insert_unsafe.is_add) {
                matching_instance_->AddEdge(insert_unsafe.id1, insert_unsafe.id2, insert_unsafe.label);
            } else if (insert_unsafe.type == 'v' && !insert_unsafe.is_add) {
                matching_instance_->RemoveVertex(insert_unsafe.id1);
            }

            sliding_window_base += min_safe + 1;
            num_e_updates += min_safe + 1;

            size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
            matching_instance_->GetNumPositiveResults(positive_num_results_cur);
            matching_instance_->GetNumNegativeResults(negative_num_results_cur);
            if (positive_num_results_cur != positive_num_results_last ||
                negative_num_results_cur != negative_num_results_last) {
                positive_num_results_last = positive_num_results_cur;
                negative_num_results_last = negative_num_results_cur;

                unsafe_updates++;
            }
        } else {
            sliding_window_base += window_size;
            num_e_updates += window_size;
        }

        // if ((update_size - sliding_window_base) % (update_size / 10) < (window_size - 1)) {
        //     std::cout << "update progress: " << (sliding_window_base) * 100 / update_size << "%" << std::endl;
        // }
    }
}

/**
 * ProcessBatchUpdatesQueue
 *
 * Port of the free function ProcessBatchUpdates: pull a medium-size batch from the queue,
 * classify updates in parallel, and if any unsafe edge is found, apply the whole batch of
 * edge insertions in order. Also applies vertex add/remove immediately.
 */
void InterExecutor::ProcessBatchUpdatesQueue(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    const size_t batch_size = 10000;
    const size_t classify_threads = 8;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);

    std::vector<bool> is_safes(batch_size, false);
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);

    bool safe = true;
    size_t idx_update = batch_size + 5;

    size_t actual_threads = std::min(classify_threads, batch_updates.size());

    #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
    {
        #pragma omp single
        {
            for (size_t i = 0; i < batch_size && !data_graph_.updates_.empty(); ++i) {
                batch_updates.push_back(data_graph_.updates_.front());
                data_graph_.updates_.pop();
            }
        }

        #pragma omp for
        for (size_t cnt = 0; cnt < batch_updates.size(); cnt++) {
            const auto& insert = batch_updates[cnt];
            if (insert.type == 'e' && insert.is_add) {
                is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                num_e_update_vec[cnt]++;
                if (!is_safes[cnt]) {
                    #pragma omp critical
                    {
                        safe = false;
                        if (idx_update >= cnt) {
                            idx_update = cnt;
                        }
                    }
                }
            }
            if (insert.type == 'v' && insert.is_add) {
                matching_instance_->AddVertex(insert.id1, insert.label);
            } else if (insert.type == 'v' && !insert.is_add) {
                matching_instance_->RemoveVertex(insert.id1);
            } else if (insert.type == 'e' && !insert.is_add) {
                is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                num_e_update_vec[cnt]++;
                if (!is_safes[cnt]) {
                    #pragma omp critical
                    {
                        safe = false;
                        if (idx_update >= cnt) {
                            idx_update = cnt;
                        }
                    }
                }
            }
        }
    }

    if (!safe) {
        for (size_t i = 0; i < batch_updates.size(); i++) {
            const auto& insert = batch_updates[i];
            matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
        }
    }

    for (size_t i = 0; i < num_e_update_vec.size(); i++) {
        num_e_updates += num_e_update_vec[i];
    }
}

/**
 * SingleThreadUpdate
 *
 * Port of free function SingleThreadUpdate: sequentially pop all updates from the queue and
 * apply them one by one, updating counters and tracking changes to result counts to measure
 * unsafe updates.
 */
void InterExecutor::SingleThreadUpdate(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    size_t num_updates = data_graph_.updates_.size();
    std::cout << "Start with " << num_updates << " Update" << std::endl;

    for (size_t i = 0; i < num_updates; i++) {
        InsertUnit insert = data_graph_.updates_.front();
        data_graph_.updates_.pop();

        if (insert.type == 'v' && insert.is_add) {
            matching_instance_->AddVertex(insert.id1, insert.label);
            num_v_updates++;
        } else if (insert.type == 'v' && !insert.is_add) {
            matching_instance_->RemoveVertex(insert.id1);
            num_v_updates++;
        } else if (insert.type == 'e' && insert.is_add) {
            matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
            num_e_updates++;
        } else if (insert.type == 'e' && !insert.is_add) {
            matching_instance_->RemoveEdge(insert.id1, insert.id2);
            num_e_updates++;
        }
        if (reach_time_limit) break;

        size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
        matching_instance_->GetNumPositiveResults(positive_num_results_cur);
        matching_instance_->GetNumNegativeResults(negative_num_results_cur);
        if (positive_num_results_cur != positive_num_results_last ||
            negative_num_results_cur != negative_num_results_last) {
            positive_num_results_last = positive_num_results_cur;
            negative_num_results_last = negative_num_results_cur;
            unsafe_updates++;
        }
    }
}

/**
 * ProcessBatchUpdates (a.k.a. BatchUpdates3)
 *
 * Sliding-window scan over `data_graph_.updates_vec_` (vector form). The method classifies
 * updates within a fixed-size window and, upon encountering the first unsafe update in the
 * window, applies that unsafe update immediately. If a whole window is safe, the window is
 * skipped in one step. Tracks changes in result counters to count unsafe updates.
 *
 * Parameters:
 * - num_v_updates: cumulative counter of vertex updates (not incremented here; only safe vertex adds applied)
 * - num_e_updates: cumulative counter of edge updates (advanced by window/found index)
 * - unsafe_updates: increments when results change after applying an unsafe update
 * - count: external total update counter (not modified here)
 * - positive_num_results_last / negative_num_results_last: previous result counters to detect changes
 * - reach_time_limit: stop early if set to true (read-only here)
 * - window_size: number of updates to inspect per iteration
 */
void InterExecutor::ProcessBatchUpdates(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t window_size
) {
    size_t sliding_window_base = 0;
    const size_t update_size = data_graph_.updates_vec_.size();
    
    std::vector<bool> update_safes(window_size + 1, false);
    
    while (sliding_window_base < update_size) {
        bool is_safes = true;
        size_t min_safe = window_size + 10;
        
        // Process updates in the current window
        for (size_t i = 0; i < window_size; i++) {
            if (i + sliding_window_base >= update_size) {
                continue;
            }
            
            const auto& update = data_graph_.updates_vec_[i + sliding_window_base];
            
            if (IsUpdateSafe(update)) {
                update_safes[i] = true;
                // Apply safe vertex additions immediately
                if (update.type == 'v' && update.is_add) {
                    matching_instance_->AddVertex(update.id1, update.label);
                }
            } else {
                update_safes[i] = false;
                is_safes = false;
                if (min_safe > i) {
                    min_safe = i;
                }
                break; // Stop processing this window
            }
        }
        
        if (!is_safes) {
            // Process the first unsafe update
            const auto& unsafe_update = data_graph_.updates_vec_[min_safe + sliding_window_base];
            ApplyUnsafeUpdate(unsafe_update);
            
            sliding_window_base += min_safe + 1;
            num_e_updates += min_safe + 1;
            
            // Check if results changed
            size_t positive_num_results_cur = 0, negative_num_results_cur = 0;
            matching_instance_->GetNumPositiveResults(positive_num_results_cur);
            matching_instance_->GetNumNegativeResults(negative_num_results_cur);
            
            if (positive_num_results_cur != positive_num_results_last || 
                negative_num_results_cur != negative_num_results_last) {
                positive_num_results_last = positive_num_results_cur;
                negative_num_results_last = negative_num_results_cur;
                unsafe_updates++;
            }
        } else {
            // All updates in this window are safe
            sliding_window_base += window_size;
            num_e_updates += window_size;
        }
    }
}

/**
 * ProcessSingleUpdate
 *
 * Classify a single update and either apply it immediately if unsafe or apply safe
 * vertex additions directly. Returns whether the update was safe according to the
 * classification rules used by IsUpdateSafe.
 */
bool InterExecutor::ProcessSingleUpdate(const InsertUnit& update) {
    if (IsUpdateSafe(update)) {
        if (update.type == 'v' && update.is_add) {
            matching_instance_->AddVertex(update.id1, update.label);
        }
        return true;
    } else {
        ApplyUnsafeUpdate(update);
        return false;
    }
}

/**
 * ApplyUnsafeUpdate
 *
 * Apply a single update that has been deemed unsafe for deferred processing.
 * Currently handles edge insertions and vertex deletions. For edge deletions,
 * the base matching engine may not implement RemoveEdge; a warning is printed.
 */
void InterExecutor::ApplyUnsafeUpdate(const InsertUnit& update) {
    if (update.type == 'e' && update.is_add) {
        matching_instance_->AddEdge(update.id1, update.id2, update.label);
    } else if (update.type == 'v' && !update.is_add) {
        matching_instance_->RemoveVertex(update.id1);
    } else if (update.type == 'e' && !update.is_add) {
        // Note: RemoveEdge method might not exist in the base matching class
        // This would need to be implemented or handled differently
        std::cerr << "Warning: Edge removal not implemented in base matching class" << std::endl;
    }
}

/**
 * IsUpdateSafe
 *
 * Determine whether an update can be considered "safe" to defer/apply in bulk:
 * - Edge updates are safe if matching_instance_->Classify(...) returns true.
 * - Vertex additions are treated as safe.
 * - Vertex deletions are always unsafe.
 */
bool InterExecutor::IsUpdateSafe(const InsertUnit& update) {
    if (update.type == 'v' && !update.is_add) {
        // Vertex removal is always unsafe
        return false;
    } else if (update.type == 'e') {
        // Edge updates need classification
        return matching_instance_->Classify(update.id1, update.id2, update.label); // self defined filtering 
    } else if (update.type == 'v' && update.is_add) {
        // Vertex addition is always safe
        return true;
    }
    
    return false;
}

/**
 * DegreePruning
 *
 * Evaluate a simple degree-based rule between candidate data vertices (v1, v2)
 * and query vertices (u1, u2). Returns true if the degree condition holds for
 * the current rule; returns false otherwise.
 */
bool InterExecutor::DegreePruning(uint v1, uint v2, uint u1, uint u2, Graph& data_, Graph& query_){
    if (data_.GetDegree(v1) < query_.GetDegree(u1) && data_.GetDegree(v2) < query_.GetDegree(u2)){
        return true;
    }
    return false;
}

/**
 * QueryGraphPruning
 *
 * Check a structural/edge-label rule derived from the query graph for mapping
 * (u1, u2) to (v1, v2). Returns true if the condition holds; false otherwise.
 */
bool InterExecutor::QueryGraphPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_){
    if (data_.GetEdgeLabel(v1, v2) == query_.GetEdgeLabel(u1, u2)){
        return true;
    }
    return false;
}

/**
 * LabelPruning
 *
 * Verify vertex/edge label compatibility for mapping (u1, u2) -> (v1, v2).
 * Returns true if the combined label condition holds; false otherwise.
 */
bool InterExecutor::LabelPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_){
    if (data_.GetEdgeLabel(v1, v2) == query_.GetEdgeLabel(u1, u2)){
        if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1) && data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){
            return true;
        }
    }
    return false;
}

/**
 * DCSPruning
 *
 * Placeholder for index-based pruning (e.g., DCS). Returns false by default
 * and contains commented guidance for a richer implementation.
 */
bool InterExecutor::DCSPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_){
    // auto v1_label = data_.GetVertexLabel(v1);
    // for (uint u1 = 0; u1 < query_.NumVertices(); u1++){ // take one u1 from query
    //     //  Get v1 that label matched to query label
    //     if (v1_label == query_.GetVertexLabel(u1)){
    //         // Get v2 that label matched to query vertex
    //         for (uint u2 = 0; u2 < query_.NumVertices(); u2++){ // take one u2 from query
    //             if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){           
    //                 // Prune 1: detect if the edge is in the query graph
    //                 // if (label == std::get<2>(query_.GetEdgeLabel(u1, u2))){return false;}                                      
    //                 bool reversed = false;
    //                 if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
    //            s     {
    //                     std::swap(u1, u2); // 
    //                     std::swap(v1, v2); //
    //                     reversed = true;
    //                 }
    //                 if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1)
    //                      != treeNode_[u2].backwards_.end())
    //                 {
    //                     // Prune 2: detect if the edge is in the DCS graph
    //                     // for more complex query graph
    //                     // detect if DCS has update. If so, that means new FM will occur
    //                     return false;
    //                     // Prune 3: detect if the edge is in the FM path
    //                     // for more complex query graph
    //                     // If so, that means new FM will occur
    //                     bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];
    //                     if(old_p_d1 || old_p_d2 || old_c_d2){
    //                         return false;
    //                     }
    //                 }
    //                 if (reversed)
    //                 {
    //                     std::swap(u1, u2);
    //                     std::swap(v1, v2);
    //                 }
    //             }
    //         }
    //     }
    // }

    return false;
}

/**
 * FMPathPruning
 *
 * Placeholder for feature/path-based pruning between (u1, u2) and (v1, v2).
 * Returns false by default.
 */
bool InterExecutor::FMPathPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_){
    (void)v1; (void)v2; (void)u1; (void)u2; (void)query_; (void)data_;
    return false;
}
