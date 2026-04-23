#include <vector>
#include <atomic>
#include <iostream>
#include <thread>
#include <functional>
#include <unordered_set>
#include <immintrin.h>
#include <omp.h>
#include "inter_executor.h"
#include "matching_executor/matching.h"
#include "matching_executor/Parallel_GraphFlow/parallel_graphflow.h"
#include "graph_storage/graph.h"
#include "taskflow/taskflow.hpp"

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
                        #pragma omp critical
                        {
                            safe = false;
                            if (idx_update > cnt) {
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
                    #pragma omp critical
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
            } else if (insert_unsafe.type == 'e' && !insert_unsafe.is_add) {
                matching_instance_->RemoveEdge(insert_unsafe.id1, insert_unsafe.id2);
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
        size_t unsafe_start = idx_update;
        for (size_t i = unsafe_start; i < batch_updates.size(); i++) {
            const auto& insert = batch_updates[i];
            if (insert.type == 'e' && insert.is_add) {
                matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
            } else if (insert.type == 'e' && !insert.is_add) {
                matching_instance_->RemoveEdge(insert.id1, insert.id2);
            } else if (insert.type == 'v' && insert.is_add) {
                matching_instance_->AddVertex(insert.id1, insert.label);
            } else if (insert.type == 'v' && !insert.is_add) {
                matching_instance_->RemoveVertex(insert.id1);
            }
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
 * BatchUpdates_Taskflow
 *
 * Taskflow-based parallel batch processing. Similar to BatchUpdates2 but uses Taskflow
 * for better parallel task scheduling. Pulls a batch from the queue, classifies updates
 * in parallel using Taskflow tasks, and applies unsafe updates sequentially.
 *
 * Parameters:
 * - num_v_updates: cumulative counter of vertex updates
 * - num_e_updates: cumulative counter of edge updates
 * - unsafe_updates: increments when results change after applying updates
 * - count: external total update counter (not modified here)
 * - positive_num_results_last / negative_num_results_last: previous result counters to detect changes
 * - reach_time_limit: stop early if set to true
 */
void InterExecutor::BatchUpdates_Taskflow(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    const size_t batch_size = 400;
    const size_t num_threads = 16;

    // Create Taskflow executor
    tf::Executor executor(num_threads);
    tf::Taskflow taskflow;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);
    std::vector<bool> is_safes(batch_size, false);
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);
    
    std::atomic<bool> safe(true);
    std::atomic<size_t> idx_update(batch_size + 5);
    std::atomic<size_t> local_v_updates(0);

    while (!data_graph_.updates_.empty() && !reach_time_limit) {
        batch_updates.clear();
        safe.store(true);
        is_safes.assign(batch_size, false);
        num_e_update_vec.assign(batch_size, 0ul);
        idx_update.store(batch_size + 5);
        local_v_updates.store(0);

        // Pop batch from queue (single-threaded)
        for (size_t i = 0; i < batch_size && !data_graph_.updates_.empty(); ++i) {
            batch_updates.push_back(data_graph_.updates_.front());
            data_graph_.updates_.pop();
        }

        if (batch_updates.empty()) {
            break;
        }

        // Create Taskflow tasks for parallel classification
        taskflow.clear();
        
        // Parallel classification tasks
        for (size_t cnt = 0; cnt < batch_updates.size(); cnt++) {
            taskflow.emplace([this, &batch_updates, &is_safes, &num_e_update_vec, 
                             &safe, &idx_update, &local_v_updates, cnt]() {
                if (!safe.load()) {
                    return;
                }

                const auto& insert = batch_updates[cnt];
                if (insert.type == 'e' && insert.is_add) {
                    is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                    num_e_update_vec[cnt]++;
                    
                    if (!is_safes[cnt]) {
                        safe.store(false);
                        size_t current_idx = idx_update.load();
                        while (current_idx > cnt && 
                               !idx_update.compare_exchange_weak(current_idx, cnt)) {
                            current_idx = idx_update.load();
                        }
                    }
                } else if (insert.type == 'v' && insert.is_add) {
                    // Vertex additions can be handled immediately
                    matching_instance_->AddVertex(insert.id1, insert.label);
                    local_v_updates++;
                } else if (insert.type == 'v' && !insert.is_add) {
                    // Vertex removals are unsafe
                    safe.store(false);
                    size_t current_idx = idx_update.load();
                    while (current_idx > cnt && 
                           !idx_update.compare_exchange_weak(current_idx, cnt)) {
                        current_idx = idx_update.load();
                    }
                } else if (insert.type == 'e' && !insert.is_add) {
                    is_safes[cnt] = matching_instance_->Classify(insert.id1, insert.id2, insert.label);
                    num_e_update_vec[cnt]++;
                    
                    if (!is_safes[cnt]) {
                        safe.store(false);
                        size_t current_idx = idx_update.load();
                        while (current_idx > cnt && 
                               !idx_update.compare_exchange_weak(current_idx, cnt)) {
                            current_idx = idx_update.load();
                        }
                    }
                }
            });
        }

        // Execute all classification tasks in parallel
        executor.run(taskflow).wait();

        // Accumulate vertex updates
        num_v_updates += local_v_updates.load();

        // Apply unsafe updates sequentially
        if (!safe.load()) {
            size_t unsafe_idx = idx_update.load();
            for (size_t i = unsafe_idx; i < batch_updates.size(); i++) {
                const auto& insert = batch_updates[i];
                if (insert.type == 'e' && insert.is_add) {
                    matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
                } else if (insert.type == 'e' && !insert.is_add) {
                    matching_instance_->RemoveEdge(insert.id1, insert.id2);
                } else if (insert.type == 'v' && !insert.is_add) {
                    matching_instance_->RemoveVertex(insert.id1);
                }

                // Check if results changed
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

        // Accumulate edge update counts
        for (size_t i = 0; i < num_e_update_vec.size(); i++) {
            num_e_updates += num_e_update_vec[i];
        }
    }
}

/**
 * BatchUpdates4
 *
 * 串行地遍历 `data_graph_.updates_vec_` 中的所有更新：
 *  - 顶点增加/删除：直接调用 AddVertex / RemoveVertex
 *  - 边增加/删除：直接调用 AddEdge / RemoveEdge
 * 每次更新后通过全局结果计数器判断是否发生结果变化，从而统计
 * 「unsafe updates」。该实现完全单线程，方便与现有策略做结果对比，
 * 后续也可以在「每个 update 为一个任务单元」的基础上引入 Taskflow 并行。
 */
void InterExecutor::BatchUpdates4(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    (void)count;           // 当前策略不使用该计数器
    (void)reach_time_limit; // 暂时忽略时间限制

    const size_t update_size = data_graph_.updates_vec_.size();
    std::cout << "Start BatchUpdates4 with " << update_size << " updates" << std::endl;

    // 使用单线程 Taskflow，将「每个 update 一个 task」串行执行。
    tf::Executor executor(1);
    tf::Taskflow taskflow;

    std::vector<tf::Task> tasks;
    tasks.reserve(update_size);

    for (size_t i = 0; i < update_size; ++i) {
        const InsertUnit& insert = data_graph_.updates_vec_[i];

        // 为每个 update 创建一个 task，接收 Subflow 以支持内层动态展开
        auto t = taskflow.emplace([this,
                                   &num_v_updates,
                                   &num_e_updates,
                                   &unsafe_updates,
                                   &positive_num_results_last,
                                   &negative_num_results_last,
                                   &insert](tf::Subflow& sf) {
            if (insert.type == 'v' && insert.is_add) {
                matching_instance_->AddVertex(insert.id1, insert.label);
                num_v_updates++;
            } else if (insert.type == 'v' && !insert.is_add) {
                matching_instance_->RemoveVertex(insert.id1);
                num_v_updates++;
            } else if (insert.type == 'e' && insert.is_add) {
                if (auto* pg = dynamic_cast<Parallel_Graphflow*>(matching_instance_)) {
                    pg->AddEdgeWithSubflow(insert.id1, insert.id2, insert.label, sf);
                } else {
                    matching_instance_->AddEdge(insert.id1, insert.id2, insert.label);
                }
                num_e_updates++;
            } else if (insert.type == 'e' && !insert.is_add) {
                if (auto* pg = dynamic_cast<Parallel_Graphflow*>(matching_instance_)) {
                    pg->RemoveEdgeWithSubflow(insert.id1, insert.id2, sf);
                } else {
                    matching_instance_->RemoveEdge(insert.id1, insert.id2);
                }
                num_e_updates++;
            }

            // 参考 SingleThreadUpdate 的统计方式，判断本次更新是否为「unsafe」
            size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
            matching_instance_->GetNumPositiveResults(positive_num_results_cur);
            matching_instance_->GetNumNegativeResults(negative_num_results_cur);
            if (positive_num_results_cur != positive_num_results_last ||
                negative_num_results_cur != negative_num_results_last) {
                positive_num_results_last = positive_num_results_cur;
                negative_num_results_last = negative_num_results_cur;
                unsafe_updates++;
            }
        });

        tasks.push_back(t);
        if (i > 0) {
            // 保证更新按原来的顺序执行，即第 i-1 个任务先于第 i 个任务
            tasks[i - 1].precede(tasks[i]);
        }
    }

    executor.run(taskflow).wait();
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
        matching_instance_->RemoveEdge(update.id1, update.id2);
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

/**
 * BatchUpdates_Persistent
 *
 * Parallel classification + serial matching.
 *
 * Key design:
 * 1. Large window (256): first apply all safe vertex-adds in the window (serial, fast)
 * 2. Parallel classify all edge updates in the window via `#pragma omp for`
 * 3. Find earliest unsafe, apply it serially (AddEdge uses internal OMP parallelism)
 * 4. Advance past unsafe, repeat
 *
 * This parallelizes the classification phase (~30% of total time) while
 * preserving correctness: vertex-adds are applied before edge classification
 * so that Classify sees the correct graph state.
 */
void InterExecutor::BatchUpdates_Persistent(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t num_threads
) {
    const size_t WINDOW = 128;
    const size_t update_size = data_graph_.updates_vec_.size();
    if (num_threads < 2) num_threads = 2;
    const size_t classify_threads = num_threads - 1; // reserve main thread

    // ---- Build a persistent std::thread pool for Classify (separate from OMP) ----
    std::vector<uint8_t> slot_safe(WINDOW, 0);
    std::atomic<size_t> classify_next{0};
    std::atomic<size_t> classify_done{0};
    size_t classify_count = 0;
    std::atomic<uint64_t> classify_epoch{0};
    std::atomic<bool> pool_stop{false};
    size_t classify_base = 0; // window base for classify workers

    auto classify_worker = [&](size_t /*tid*/) {
        uint64_t last_epoch = 0;
        while (!pool_stop.load(std::memory_order_relaxed)) {
            uint64_t cur = classify_epoch.load(std::memory_order_acquire);
            if (cur == last_epoch) { _mm_pause(); continue; }
            last_epoch = cur;
            while (true) {
                size_t idx = classify_next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= classify_count) break;
                if (slot_safe[idx] != 0) {
                    classify_done.fetch_add(1, std::memory_order_release);
                    continue;
                }
                const auto& update = data_graph_.updates_vec_[idx + classify_base];
                bool safe = matching_instance_->Classify(update.id1, update.id2, update.label);
                slot_safe[idx] = safe ? 1 : 2;
                classify_done.fetch_add(1, std::memory_order_release);
            }
        }
    };

    // Launch classify workers
    std::vector<std::thread> classify_pool;
    for (size_t i = 0; i < classify_threads; i++)
        classify_pool.emplace_back(classify_worker, i);

    // ---- Main update loop (main thread) ----
    size_t sliding_window_base = 0;

    while (sliding_window_base < update_size) {
        const size_t actual = std::min(WINDOW, update_size - sliding_window_base);

        // Phase 1: vertex-adds + mark edges (main thread, fast)
        for (size_t i = 0; i < actual; i++) {
            const auto& update = data_graph_.updates_vec_[i + sliding_window_base];
            if (update.type == 'v' && update.is_add) {
                matching_instance_->AddVertex(update.id1, update.label);
                slot_safe[i] = 1;
            } else if (update.type == 'v' && !update.is_add) {
                slot_safe[i] = 2;
            } else {
                slot_safe[i] = 0; // needs classify
            }
        }

        // Phase 2: dispatch classify to std::thread pool
        classify_base = sliding_window_base;
        classify_count = actual;
        classify_next.store(0, std::memory_order_relaxed);
        classify_done.store(0, std::memory_order_relaxed);
        classify_epoch.fetch_add(1, std::memory_order_release); // wake workers

        // Main thread also helps classify
        while (true) {
            size_t idx = classify_next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= actual) break;
            if (slot_safe[idx] != 0) {
                classify_done.fetch_add(1, std::memory_order_release);
                continue;
            }
            const auto& update = data_graph_.updates_vec_[idx + sliding_window_base];
            bool safe = matching_instance_->Classify(update.id1, update.id2, update.label);
            slot_safe[idx] = safe ? 1 : 2;
            classify_done.fetch_add(1, std::memory_order_release);
        }

        // Wait for all classify done
        while (classify_done.load(std::memory_order_acquire) < actual) { _mm_pause(); }

        // Ensure workers have drained (set next past end)
        classify_next.store(actual + classify_threads + 1, std::memory_order_release);

        // Phase 3: find first unsafe
        size_t first_unsafe = actual;
        for (size_t i = 0; i < actual; i++) {
            if (slot_safe[i] == 2) { first_unsafe = i; break; }
        }

        if (first_unsafe < actual) {
            // Phase 4: apply unsafe update (AddEdge uses OMP internally — no conflict!)
            const auto& unsafe_update = data_graph_.updates_vec_[first_unsafe + sliding_window_base];
            ApplyUnsafeUpdate(unsafe_update);

            sliding_window_base += first_unsafe + 1;
            num_e_updates += first_unsafe + 1;

            size_t pos_cur = 0, neg_cur = 0;
            matching_instance_->GetNumPositiveResults(pos_cur);
            matching_instance_->GetNumNegativeResults(neg_cur);
            if (pos_cur != positive_num_results_last ||
                neg_cur != negative_num_results_last) {
                positive_num_results_last = pos_cur;
                negative_num_results_last = neg_cur;
                unsafe_updates++;
            }
        } else {
            sliding_window_base += actual;
            num_e_updates += actual;
        }
    }

    // Shutdown classify pool
    pool_stop.store(true, std::memory_order_release);
    classify_epoch.fetch_add(1, std::memory_order_release);
    for (auto& t : classify_pool) t.join();
}

// ===========================================================================
// BatchUpdates_Pipelined
//
// Pipelined classify caching + prefetch + adaptive window.
// ===========================================================================

void InterExecutor::BatchUpdates_Pipelined(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t num_threads
) {
    const size_t MIN_WINDOW = 64;
    const size_t MAX_WINDOW = 1024;
    size_t window_size = 128; // initial, will adapt
    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;
    if (num_threads < 2) num_threads = 2;
    const size_t classify_threads = num_threads - 1;

    // ---- Global classify cache across windows ----
    // 0 = not classified, 1 = safe, 2 = unsafe
    std::vector<uint8_t> classify_cache(update_size, 0);

    // ---- Persistent classify thread pool (same design as BatchUpdates_Persistent) ----
    std::atomic<size_t> pool_next{0};
    std::atomic<size_t> pool_done{0};
    size_t pool_count = 0;
    size_t pool_base = 0;
    std::atomic<uint64_t> pool_epoch{0};
    std::atomic<bool> pool_stop{false};

    // Pointer to cache so workers can write to it
    uint8_t* cache_ptr = classify_cache.data();

    auto worker_fn = [&](size_t /*tid*/) {
        uint64_t last_epoch = 0;
        while (!pool_stop.load(std::memory_order_relaxed)) {
            uint64_t cur = pool_epoch.load(std::memory_order_acquire);
            if (cur == last_epoch) { _mm_pause(); continue; }
            last_epoch = cur;
            while (true) {
                size_t idx = pool_next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= pool_count) break;
                size_t global_idx = idx + pool_base;
                // Skip if already cached (from prefetch or previous window)
                if (cache_ptr[global_idx] != 0) {
                    pool_done.fetch_add(1, std::memory_order_release);
                    continue;
                }
                const auto& update = data_graph_.updates_vec_[global_idx];
                bool safe = matching_instance_->Classify(update.id1, update.id2, update.label);
                cache_ptr[global_idx] = safe ? 1 : 2;
                pool_done.fetch_add(1, std::memory_order_release);
            }
        }
    };

    std::vector<std::thread> workers;
    for (size_t i = 0; i < classify_threads; i++)
        workers.emplace_back(worker_fn, i);

    // Helper: dispatch classify for [base, base+count) and wait
    auto dispatch_and_wait = [&](size_t base, size_t cnt) {
        pool_base = base;
        pool_count = cnt;
        pool_next.store(0, std::memory_order_relaxed);
        pool_done.store(0, std::memory_order_relaxed);
        pool_epoch.fetch_add(1, std::memory_order_release);

        // Main thread also helps
        while (true) {
            size_t idx = pool_next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= cnt) break;
            size_t global_idx = idx + base;
            if (cache_ptr[global_idx] != 0) {
                pool_done.fetch_add(1, std::memory_order_release);
                continue;
            }
            const auto& update = data_graph_.updates_vec_[global_idx];
            bool safe = matching_instance_->Classify(update.id1, update.id2, update.label);
            cache_ptr[global_idx] = safe ? 1 : 2;
            pool_done.fetch_add(1, std::memory_order_release);
        }
        while (pool_done.load(std::memory_order_acquire) < cnt) { _mm_pause(); }
        // Drain workers
        pool_next.store(cnt + classify_threads + 1, std::memory_order_release);
    };

    // Helper: dispatch prefetch (non-blocking, returns immediately)
    // Call wait_prefetch() later to ensure completion.
    size_t prefetch_base = 0, prefetch_count = 0;
    bool prefetch_active = false;

    auto start_prefetch = [&](size_t base, size_t cnt) {
        prefetch_base = base;
        prefetch_count = cnt;
        prefetch_active = true;
        pool_base = base;
        pool_count = cnt;
        pool_next.store(0, std::memory_order_relaxed);
        pool_done.store(0, std::memory_order_relaxed);
        pool_epoch.fetch_add(1, std::memory_order_release);
        // Don't wait — return immediately so main thread can do matching
    };

    auto wait_prefetch = [&]() {
        if (!prefetch_active) return;
        // Main thread helps finish remaining work
        while (true) {
            size_t idx = pool_next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= prefetch_count) break;
            size_t global_idx = idx + prefetch_base;
            if (cache_ptr[global_idx] != 0) {
                pool_done.fetch_add(1, std::memory_order_release);
                continue;
            }
            const auto& update = data_graph_.updates_vec_[global_idx];
            bool safe = matching_instance_->Classify(update.id1, update.id2, update.label);
            cache_ptr[global_idx] = safe ? 1 : 2;
            pool_done.fetch_add(1, std::memory_order_release);
        }
        while (pool_done.load(std::memory_order_acquire) < prefetch_count) { _mm_pause(); }
        pool_next.store(prefetch_count + classify_threads + 1, std::memory_order_release);
        prefetch_active = false;
    };

    // ---- Main update loop ----
    size_t sliding_window_base = 0;
    size_t consecutive_safe_windows = 0;

    while (sliding_window_base < update_size && !reach_time_limit) {
        // Wait for any outstanding prefetch before we start this window
        wait_prefetch();

        const size_t actual = std::min(window_size, update_size - sliding_window_base);

        // Phase 1: Handle vertex updates + mark in cache (main thread, fast)
        bool has_vertex_change = false;
        for (size_t i = 0; i < actual; i++) {
            size_t gi = i + sliding_window_base;
            const auto& update = data_graph_.updates_vec_[gi];
            if (update.type == 'v' && update.is_add) {
                matching_instance_->AddVertex(update.id1, update.label);
                classify_cache[gi] = 1; // safe
                num_v_updates++;
                has_vertex_change = true;
            } else if (update.type == 'v' && !update.is_add) {
                classify_cache[gi] = 2; // unsafe
                has_vertex_change = true;
            }
            // edge updates: leave cache as-is (may already be prefetched)
        }

        // If vertex label changed, invalidate uncached edges after this point
        // (conservative: only invalidate entries that were prefetched BEFORE the vertex change)
        // In practice, Classify only checks label match which is unchanged by AddEdge,
        // but vertex add/remove can change the label set. Be safe: invalidate if vertex changed.
        if (has_vertex_change) {
            for (size_t i = 0; i < actual; i++) {
                size_t gi = i + sliding_window_base;
                const auto& update = data_graph_.updates_vec_[gi];
                if (update.type == 'e' && classify_cache[gi] != 0) {
                    // This was prefetched before vertex change — might be stale
                    // Only invalidate if there was a vertex change before this edge
                    // For simplicity, invalidate all prefetched edges in this window
                    classify_cache[gi] = 0;
                }
            }
        }

        // Phase 2: Parallel classify edges not yet in cache
        dispatch_and_wait(sliding_window_base, actual);

        // Phase 3: Find first unsafe
        size_t first_unsafe = actual;
        for (size_t i = 0; i < actual; i++) {
            if (classify_cache[i + sliding_window_base] == 2) {
                first_unsafe = i;
                break;
            }
        }

        if (first_unsafe >= actual) {
            // All safe — skip entire window
            sliding_window_base += actual;
            num_e_updates += actual;
            consecutive_safe_windows++;
            // Adaptive: grow window
            if (consecutive_safe_windows >= 2) {
                window_size = std::min(window_size * 2, MAX_WINDOW);
            }
            continue;
        }

        // Phase 4: Apply unsafe update
        // Start prefetching next window in background (unless vertex removal)
        consecutive_safe_windows = 0;
        // Adaptive: shrink window
        window_size = std::max(window_size / 2, MIN_WINDOW);

        size_t next_base = sliding_window_base + first_unsafe + 1;
        size_t next_actual = (next_base < update_size)
            ? std::min(window_size, update_size - next_base) : 0;

        const auto& unsafe_update = data_graph_.updates_vec_[first_unsafe + sliding_window_base];
        bool is_vertex_removal = (unsafe_update.type == 'v' && !unsafe_update.is_add);

        // Pipeline: prefetch next window while we apply unsafe update
        if (next_actual > 0 && !is_vertex_removal) {
            start_prefetch(next_base, next_actual);
        }

        // Apply the unsafe update (this is the expensive part — AddEdge + FindMatches)
        ApplyUnsafeUpdate(unsafe_update);

        sliding_window_base = next_base;
        num_e_updates += first_unsafe + 1;

        // Check results change
        size_t pos_cur = 0, neg_cur = 0;
        matching_instance_->GetNumPositiveResults(pos_cur);
        matching_instance_->GetNumNegativeResults(neg_cur);
        if (pos_cur != positive_num_results_last ||
            neg_cur != negative_num_results_last) {
            positive_num_results_last = pos_cur;
            negative_num_results_last = neg_cur;
            unsafe_updates++;
        }
    }

    // Wait for any outstanding prefetch
    wait_prefetch();

    // Shutdown pool
    pool_stop.store(true, std::memory_order_release);
    pool_epoch.fetch_add(1, std::memory_order_release);
    for (auto& t : workers) t.join();
}

// ===========================================================================
// BatchUpdates_AllAtOnce
//
// Strategy: pre-classify ALL updates, add ALL unsafe edges to graph, then
// enumerate matches for ALL unsafe edges in parallel (inter-update parallelism).
// ===========================================================================

void InterExecutor::BatchUpdates_AllAtOnce(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t num_threads
) {
    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    // ---- Phase 1: Apply vertex updates & classify edge updates ----
    // Vertex updates must be applied in order (they change vertex labels that
    // Classify depends on). Edge classification is read-only w.r.t. vertex
    // labels and can be parallelised.

    struct UnsafeEdge {
        uint v1, v2, label;
    };

    std::vector<uint8_t> update_type(update_size, 0); // 0=safe-edge, 1=unsafe-edge, 2=vertex
    // First pass (serial): apply vertex updates
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v') {
            update_type[i] = 2;
            if (u.is_add) {
                matching_instance_->AddVertex(u.id1, u.label);
            }
            num_v_updates++;
        }
    }

    // Second pass (parallel): classify edge updates
    #pragma omp parallel for num_threads(num_threads) schedule(static, 512)
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] != 0) continue; // already handled (vertex)
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type != 'e' || !u.is_add) {
            update_type[i] = 2; // edge deletion or other → treat as vertex-like
            continue;
        }
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        update_type[i] = safe ? 0 : 1;
    }

    // Collect unsafe edges (preserve order)
    std::vector<UnsafeEdge> unsafe_edges;
    unsafe_edges.reserve(update_size / 10); // heuristic
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (update_type[i] == 1) {
            unsafe_edges.push_back({u.id1, u.id2, u.label});
        }
        if (u.type == 'e') num_e_updates++;
    }

    std::cout << "[batch_all] total updates: " << update_size
              << "  unsafe edges: " << unsafe_edges.size()
              << "  (safe skipped: " << (update_size - unsafe_edges.size() - num_v_updates)
              << ")" << std::endl;

    // ---- Phase 2: Add ALL unsafe edges to data graph + update indices (serial) ----
    for (const auto& e : unsafe_edges) {
        data_graph_.AddEdge(e.v1, e.v2, e.label);
        // Update algorithm-specific index structures (DCS for TurboFlux, etc.)
        // No-op for GraphFlow which has no auxiliary indices.
        matching_instance_->UpdateIndexForEdge(e.v1, e.v2, e.label);
    }

    // ---- Phase 3: Prepare per-thread state for parallel enumeration ----
    matching_instance_->PrepareBatchEnumeration(num_threads);

    // ---- Phase 4: Enumerate matches in parallel across unsafe edges ----
    const size_t n_unsafe = unsafe_edges.size();
    std::vector<size_t> thread_results(num_threads, 0);

    #pragma omp parallel num_threads(num_threads)
    {
        size_t tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 1)
        for (size_t k = 0; k < n_unsafe; k++) {
            if (reach_time_limit) continue;
            const auto& e = unsafe_edges[k];
            size_t r = matching_instance_->EnumerateNewEdge(e.v1, e.v2, e.label, tid);
            thread_results[tid] += r;
        }
    }

    // ---- Phase 5: Aggregate results ----
    size_t total_positive = 0;
    for (size_t t = 0; t < num_threads; t++) {
        total_positive += thread_results[t];
    }

    // Bump the matching instance's positive result counter
    matching_instance_->AddPositiveResults(total_positive);

    size_t pos_cur = 0, neg_cur = 0;
    matching_instance_->GetNumPositiveResults(pos_cur);
    matching_instance_->GetNumNegativeResults(neg_cur);
    positive_num_results_last = pos_cur;
    negative_num_results_last = neg_cur;
    unsafe_updates = unsafe_edges.size();

    // Handle vertex removals (deferred — must be in order after everything)
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v' && !u.is_add) {
            matching_instance_->RemoveVertex(u.id1);
        }
    }
}

// =====================================================================
// Versioned batch update: inner-update semantics + full parallelism
// =====================================================================

void InterExecutor::BatchUpdates_Versioned(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t num_threads
) {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    // ---- Phase 1: Initialize timestamps on existing graph ----
    auto t1 = std::chrono::high_resolution_clock::now();
    data_graph_.InitTimestamps();

    // ---- Phase 2: Classify all edges (parallel) + handle vertex updates ----
    struct VersionedEdge { uint v1, v2, label, timestamp; };
    std::vector<uint8_t> update_type(update_size, 0); // 0=safe_edge, 1=unsafe_edge, 2=vertex

    // Sequential pass: vertex additions (must be in order)
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v') {
            update_type[i] = 2;
            if (u.is_add) matching_instance_->AddVertex(u.id1, u.label);
            num_v_updates++;
        }
    }

    // Parallel pass: classify edge updates
    #pragma omp parallel for num_threads(num_threads) schedule(static, 512)
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] != 0) continue;
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type != 'e' || !u.is_add) {
            update_type[i] = 2;
            continue;
        }
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        update_type[i] = safe ? 0 : 1;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;

    // ---- Phase 3: Add ALL edges with timestamps, collect unsafe edges ----
    auto t3 = std::chrono::high_resolution_clock::now();
    std::vector<VersionedEdge> unsafe_edges;
    unsafe_edges.reserve(update_size / 10);

    uint ts = 1;  // timestamp 0 = original graph edges
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'e' && u.is_add) {
            num_e_updates++;
            data_graph_.AddEdgeVersioned(u.id1, u.id2, u.label, ts);
            if (update_type[i] == 1) {
                unsafe_edges.push_back({u.id1, u.id2, u.label, ts});
            }
            ts++;
        }
    }

    auto t4 = std::chrono::high_resolution_clock::now();
    double add_edges_ms = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0;

    std::cout << "[versioned] classify: " << classify_ms << " ms, "
              << "add edges: " << add_edges_ms << " ms ("
              << unsafe_edges.size() << " unsafe / " << (ts - 1) << " total)" << std::endl;

    // ---- Phase 4: Prepare per-thread state ----
    matching_instance_->PrepareBatchEnumeration(num_threads);

    // ---- Phase 5: Search all unsafe edges in parallel (version-aware) ----
    auto t5 = std::chrono::high_resolution_clock::now();
    const size_t n_unsafe = unsafe_edges.size();
    std::vector<size_t> thread_results(num_threads, 0);

    #pragma omp parallel num_threads(num_threads)
    {
        size_t tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic, 1)
        for (size_t k = 0; k < n_unsafe; k++) {
            if (reach_time_limit) continue;
            const auto& e = unsafe_edges[k];
            size_t r = matching_instance_->EnumerateNewEdgeVersioned(
                e.v1, e.v2, e.label, tid, e.timestamp);
            thread_results[tid] += r;
        }
    }

    auto t6 = std::chrono::high_resolution_clock::now();
    double search_ms = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() / 1000.0;

    // ---- Phase 6: Aggregate results ----
    size_t total_positive = 0;
    for (size_t t = 0; t < num_threads; t++) {
        total_positive += thread_results[t];
    }
    matching_instance_->AddPositiveResults(total_positive);

    size_t pos_cur = 0, neg_cur = 0;
    matching_instance_->GetNumPositiveResults(pos_cur);
    matching_instance_->GetNumNegativeResults(neg_cur);
    positive_num_results_last = pos_cur;
    negative_num_results_last = neg_cur;
    unsafe_updates = unsafe_edges.size();

    // Handle vertex removals
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v' && !u.is_add) {
            matching_instance_->RemoveVertex(u.id1);
        }
    }

    // Clean up timestamps
    data_graph_.ClearTimestamps();

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_total_start).count() / 1000.0;

    std::cout << "[versioned] search: " << search_ms << " ms, "
              << "total: " << total_ms << " ms, "
              << "positive: " << total_positive << std::endl;
}

// =====================================================================
// GPU BFS Versioned: inner-update semantics with GPU batch parallelism
// =====================================================================

void InterExecutor::BatchUpdates_GPU_BFS_Versioned(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    Graph& query = matching_instance_->GetQueryGraph();
    Graph& data = matching_instance_->GetDataGraph();

    // ---- Phase 1: Vertex updates + classify edges (CPU OMP) ----
    auto t1 = std::chrono::high_resolution_clock::now();

    data_graph_.InitTimestamps();

    struct VersionedEdge { uint v1, v2, label; uint32_t timestamp; };
    std::vector<uint8_t> update_type(update_size, 0);

    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v') {
            update_type[i] = 2;
            if (u.is_add) matching_instance_->AddVertex(u.id1, u.label);
            num_v_updates++;
        }
    }

    #pragma omp parallel for schedule(static, 512)
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] != 0) continue;
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type != 'e' || !u.is_add) { update_type[i] = 2; continue; }
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        update_type[i] = safe ? 0 : 1;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;

    // ---- Phase 2: Add ALL edges with timestamps ----
    auto t3 = std::chrono::high_resolution_clock::now();

    std::vector<VersionedEdge> unsafe_edges;
    unsafe_edges.reserve(update_size / 10);

    uint32_t ts = 1;
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'e' && u.is_add) {
            num_e_updates++;
            data_graph_.AddEdgeVersioned(u.id1, u.id2, u.label, ts);
            if (update_type[i] == 1) {
                unsafe_edges.push_back({u.id1, u.id2, u.label, ts});
            }
            ts++;
        }
    }

    auto t4 = std::chrono::high_resolution_clock::now();
    double add_edges_ms = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0;

    // ---- Phase 3: Build versioned CSR on GPU ----
    auto t5 = std::chrono::high_resolution_clock::now();

    if (!gpu_bfs_ready_) {
        gpu_bfs_search_.SetupQuery(query);
    }
    gpu_bfs_search_.BuildCSR_Versioned(data);

    const auto& orders = matching_instance_->GetMatchingOrders();
    uint32_t Q = query.NumVertices();
    uint32_t num_query_edges = static_cast<uint32_t>(orders.size());

    std::vector<uint32_t> flat_orders(num_query_edges * Q);
    for (uint32_t i = 0; i < num_query_edges; i++) {
        for (uint32_t j = 0; j < Q && j < orders[i].size(); j++) {
            flat_orders[i * Q + j] = orders[i][j];
        }
    }
    gpu_bfs_search_.SetAllMatchingOrders(flat_orders.data(), num_query_edges, Q);
    gpu_bfs_ready_ = true;

    auto t6 = std::chrono::high_resolution_clock::now();
    double csr_ms = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() / 1000.0;

    // ---- Phase 4: GPU BFS versioned search ----
    auto t7 = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> ev1(unsafe_edges.size()), ev2(unsafe_edges.size());
    std::vector<uint32_t> elab(unsafe_edges.size()), ets(unsafe_edges.size());
    for (size_t i = 0; i < unsafe_edges.size(); i++) {
        ev1[i] = unsafe_edges[i].v1;
        ev2[i] = unsafe_edges[i].v2;
        elab[i] = unsafe_edges[i].label;
        ets[i] = unsafe_edges[i].timestamp;
    }

    uint64_t total_matches = gpu_bfs_search_.SearchBatchEdgesBFS_Versioned(
        ev1.data(), ev2.data(), elab.data(), ets.data(),
        static_cast<uint32_t>(unsafe_edges.size()),
        num_query_edges, Q
    );

    auto t8 = std::chrono::high_resolution_clock::now();
    double search_ms = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count() / 1000.0;

    // ---- Phase 5: Aggregate results ----
    matching_instance_->AddPositiveResults(static_cast<size_t>(total_matches));

    size_t pos_cur = 0, neg_cur = 0;
    matching_instance_->GetNumPositiveResults(pos_cur);
    matching_instance_->GetNumNegativeResults(neg_cur);
    positive_num_results_last = pos_cur;
    negative_num_results_last = neg_cur;
    unsafe_updates = unsafe_edges.size();

    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v' && !u.is_add) {
            matching_instance_->RemoveVertex(u.id1);
        }
    }

    data_graph_.ClearTimestamps();

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_total_start).count() / 1000.0;

    std::cout << "[GPU-BFS-Versioned] Complete:" << std::endl;
    std::cout << "  classify:    " << classify_ms << " ms" << std::endl;
    std::cout << "  add edges:   " << add_edges_ms << " ms (" << unsafe_edges.size() << " unsafe)" << std::endl;
    std::cout << "  build CSR:   " << csr_ms << " ms" << std::endl;
    std::cout << "  BFS search:  " << search_ms << " ms (" << total_matches << " matches)" << std::endl;
    std::cout << "  total time:  " << total_ms << " ms" << std::endl;
}

// =====================================================================
// GPU-accelerated batch classification
// =====================================================================

/**
 * InitGPUClassifier
 *
 * Extract query edge label triples from the query graph and initialize
 * the GPUClassifier with current data graph vertex labels.
 */
void InterExecutor::InitGPUClassifier() {
    Graph& query = matching_instance_->GetQueryGraph();
    Graph& data = matching_instance_->GetDataGraph();

    // Collect all unique query edge label triples
    // Iterate all query vertices and their neighbors (each edge appears twice)
    std::vector<QueryEdgeTriple> triples;
    for (uint u = 0; u < query.NumVertices(); u++) {
        const auto& nbrs = query.GetNeighbors(u);
        const auto& elabs = query.GetNeighborLabels(u);
        for (size_t j = 0; j < nbrs.size(); j++) {
            uint v = nbrs[j];
            if (u < v) {  // avoid duplicates
                QueryEdgeTriple t;
                t.src_label = query.GetVertexLabel(u);
                t.dst_label = query.GetVertexLabel(v);
                t.edge_label = elabs[j];
                triples.push_back(t);
            }
        }
    }

    std::cout << "InitGPUClassifier: " << triples.size() << " query edge triples, "
              << data.NumVertices() << " data vertices" << std::endl;

    gpu_classifier_.Init(
        data.vlabels_.data(),
        data.NumVertices(),
        triples.data(),
        triples.size()
    );
    gpu_initialized_ = true;
}

/**
 * BatchUpdates_GPU
 *
 * GPU-accelerated sliding window batch processing over `data_graph_.updates_vec_`.
 * Uses GPUClassifier to classify edges in bulk, then applies the first unsafe edge on CPU.
 */
void InterExecutor::BatchUpdates_GPU(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit,
    size_t window_size
) {
    // Initialize GPU classifier on first call
    if (!gpu_initialized_) {
        InitGPUClassifier();
    }

    const size_t update_size = data_graph_.updates_vec_.size();
    size_t sliding_window_base = 0;

    // Buffers for edge batch
    std::vector<EdgeToClassify> edge_batch;
    std::vector<uint8_t> classify_results;
    // Map from position in edge_batch back to position in window
    std::vector<size_t> edge_to_window_idx;

    edge_batch.reserve(window_size);
    classify_results.reserve(window_size);
    edge_to_window_idx.reserve(window_size);

    size_t total_gpu_edges = 0;
    size_t total_unsafe = 0;

    // Timing accumulators (microseconds)
    double time_collect_us = 0, time_gpu_classify_us = 0, time_vlabel_update_us = 0;
    double time_apply_unsafe_us = 0, time_scan_result_us = 0;
    size_t num_windows = 0;
    size_t vlabel_update_count = 0;  // how many times we actually updated vlabels

    // Track whether GPU vlabels need refresh
    bool vlabels_dirty = false;
    size_t last_gpu_num_vertices = matching_instance_->GetDataGraph().NumVertices();

    auto t_start = std::chrono::high_resolution_clock::now();

    while (sliding_window_base < update_size && !reach_time_limit) {
        size_t window_end = std::min(sliding_window_base + window_size, update_size);
        size_t actual_window = window_end - sliding_window_base;
        num_windows++;

        // Phase 1: Collect all edge updates in this window; apply safe vertex ops immediately
        auto t_collect_start = std::chrono::high_resolution_clock::now();
        edge_batch.clear();
        edge_to_window_idx.clear();

        bool has_vertex_removal = false;
        size_t first_vertex_removal = actual_window;

        for (size_t i = 0; i < actual_window; i++) {
            const auto& update = data_graph_.updates_vec_[sliding_window_base + i];

            if (update.type == 'v' && update.is_add) {
                // Vertex addition is always safe — apply immediately
                matching_instance_->AddVertex(update.id1, update.label);
                num_v_updates++;
                vlabels_dirty = true;  // mark that GPU vlabels are stale
            } else if (update.type == 'v' && !update.is_add) {
                // Vertex removal is always unsafe
                has_vertex_removal = true;
                if (i < first_vertex_removal) {
                    first_vertex_removal = i;
                }
                break;  // Stop at first unsafe
            } else if (update.type == 'e') {
                // Edge update — collect for GPU classification
                EdgeToClassify e;
                e.v1 = update.id1;
                e.v2 = update.id2;
                e.label = update.label;
                edge_batch.push_back(e);
                edge_to_window_idx.push_back(i);
            }
        }

        // Phase 2: GPU batch classify all collected edges
        size_t first_unsafe_window_idx = actual_window;  // sentinel
        time_collect_us += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - t_collect_start).count();

        if (!edge_batch.empty() && !has_vertex_removal) {
            // Only update GPU vertex labels if they actually changed
            if (vlabels_dirty) {
                auto t_vl = std::chrono::high_resolution_clock::now();
                Graph& data = matching_instance_->GetDataGraph();
                gpu_classifier_.UpdateVertexLabels(data.vlabels_.data(), data.NumVertices());
                last_gpu_num_vertices = data.NumVertices();
                vlabels_dirty = false;
                vlabel_update_count++;
                time_vlabel_update_us += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - t_vl).count();
            }

            auto t_gpu = std::chrono::high_resolution_clock::now();
            classify_results.resize(edge_batch.size());
            gpu_classifier_.ClassifyBatch(
                edge_batch.data(),
                edge_batch.size(),
                classify_results.data()
            );
            time_gpu_classify_us += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - t_gpu).count();
            total_gpu_edges += edge_batch.size();

            // Find first unsafe edge
            auto t_scan = std::chrono::high_resolution_clock::now();
            for (size_t j = 0; j < edge_batch.size(); j++) {
                if (!classify_results[j]) {
                    // Unsafe
                    first_unsafe_window_idx = edge_to_window_idx[j];
                    break;
                }
            }
            time_scan_result_us += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - t_scan).count();
        } else if (has_vertex_removal) {
            first_unsafe_window_idx = first_vertex_removal;
        }

        // Phase 3: Apply first unsafe update on CPU
        if (first_unsafe_window_idx < actual_window) {
            auto t_apply = std::chrono::high_resolution_clock::now();
            const auto& unsafe_update = data_graph_.updates_vec_[sliding_window_base + first_unsafe_window_idx];
            ApplyUnsafeUpdate(unsafe_update);
            time_apply_unsafe_us += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - t_apply).count();
            total_unsafe++;

            sliding_window_base += first_unsafe_window_idx + 1;
            num_e_updates += first_unsafe_window_idx + 1;

            // Check if results changed
            size_t pos_cur = 0, neg_cur = 0;
            matching_instance_->GetNumPositiveResults(pos_cur);
            matching_instance_->GetNumNegativeResults(neg_cur);
            if (pos_cur != positive_num_results_last ||
                neg_cur != negative_num_results_last) {
                positive_num_results_last = pos_cur;
                negative_num_results_last = neg_cur;
                unsafe_updates++;
            }
        } else {
            // All safe — skip entire window
            sliding_window_base += actual_window;
            num_e_updates += actual_window;
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() / 1000.0;

    std::cout << "[GPU] BatchUpdates_GPU complete: "
              << total_gpu_edges << " edges classified on GPU, "
              << total_unsafe << " unsafe updates, "
              << elapsed_ms << "ms total" << std::endl;
    std::cout << "[GPU] Timing breakdown:" << std::endl;
    std::cout << "  collect edges:    " << time_collect_us / 1000.0 << " ms" << std::endl;
    std::cout << "  vlabel update:    " << time_vlabel_update_us / 1000.0 << " ms" << std::endl;
    std::cout << "  GPU classify:     " << time_gpu_classify_us / 1000.0 << " ms" << std::endl;
    std::cout << "  scan results:     " << time_scan_result_us / 1000.0 << " ms" << std::endl;
    std::cout << "  apply unsafe:     " << time_apply_unsafe_us / 1000.0 << " ms" << std::endl;
    double overhead = elapsed_ms - (time_collect_us + time_vlabel_update_us + time_gpu_classify_us + time_scan_result_us + time_apply_unsafe_us) / 1000.0;
    std::cout << "  other overhead:   " << overhead << " ms" << std::endl;
    std::cout << "  windows:          " << num_windows << " (avg batch " 
              << (num_windows > 0 ? total_gpu_edges / num_windows : 0) << " edges/window)" << std::endl;
    std::cout << "  vlabel syncs:     " << vlabel_update_count << " / " << num_windows << " windows" << std::endl;
}

// =====================================================================
// GPU inter-update: all data on GPU, batch enumeration on GPU
// =====================================================================

void InterExecutor::BatchUpdates_GPU_AllAtOnce(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    Graph& query = matching_instance_->GetQueryGraph();
    Graph& data = matching_instance_->GetDataGraph();

    // ---- Phase 1: Apply vertex updates & classify edge updates (CPU) ----
    auto t1 = std::chrono::high_resolution_clock::now();

    struct UnsafeEdge { uint v1, v2, label; };

    std::vector<uint8_t> update_type(update_size, 0);
    // Apply vertex updates
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v') {
            update_type[i] = 2;
            if (u.is_add) matching_instance_->AddVertex(u.id1, u.label);
            num_v_updates++;
        }
    }

    // Classify edges (parallel on CPU for now)
    #pragma omp parallel for schedule(static, 512)
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] != 0) continue;
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type != 'e' || !u.is_add) { update_type[i] = 2; continue; }
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        update_type[i] = safe ? 0 : 1;
    }

    // Collect unsafe edges
    std::vector<UnsafeEdge> unsafe_edges;
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] == 1) {
            const auto& u = data_graph_.updates_vec_[i];
            unsafe_edges.push_back({u.id1, u.id2, u.label});
        }
        if (data_graph_.updates_vec_[i].type == 'e') num_e_updates++;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;

    // ---- Phase 2: Add ALL unsafe edges to data graph (CPU serial) ----
    auto t3 = std::chrono::high_resolution_clock::now();
    for (const auto& e : unsafe_edges) {
        data_graph_.AddEdge(e.v1, e.v2, e.label);
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    double add_edges_ms = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0;

    // ---- Phase 3: Build CSR on GPU from final graph state (ONCE) ----
    auto t5 = std::chrono::high_resolution_clock::now();

    if (!gpu_search_ready_) {
        gpu_search_engine_.SetupQuery(query);
    }
    gpu_search_engine_.BuildCSR(data);

    // Set all matching orders
    const auto& orders = matching_instance_->GetMatchingOrders();
    uint32_t Q = query.NumVertices();
    uint32_t num_query_edges = static_cast<uint32_t>(orders.size());

    std::vector<uint32_t> flat_orders(num_query_edges * Q);
    for (uint32_t i = 0; i < num_query_edges; i++) {
        for (uint32_t j = 0; j < Q && j < orders[i].size(); j++) {
            flat_orders[i * Q + j] = orders[i][j];
        }
    }
    gpu_search_engine_.SetAllMatchingOrders(flat_orders.data(), num_query_edges, Q);
    gpu_search_ready_ = true;

    auto t6 = std::chrono::high_resolution_clock::now();
    double csr_ms = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() / 1000.0;

    // ---- Phase 4: GPU batch enumeration ----
    auto t7 = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> ev1(unsafe_edges.size()), ev2(unsafe_edges.size()), elab(unsafe_edges.size());
    for (size_t i = 0; i < unsafe_edges.size(); i++) {
        ev1[i] = unsafe_edges[i].v1;
        ev2[i] = unsafe_edges[i].v2;
        elab[i] = unsafe_edges[i].label;
    }

    uint64_t total_matches = gpu_search_engine_.SearchBatchEdges(
        ev1.data(), ev2.data(), elab.data(),
        static_cast<uint32_t>(unsafe_edges.size()),
        num_query_edges, Q
    );

    auto t8 = std::chrono::high_resolution_clock::now();
    double search_ms = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count() / 1000.0;

    // ---- Phase 5: Aggregate results ----
    matching_instance_->AddPositiveResults(static_cast<size_t>(total_matches));

    size_t pos_cur = 0, neg_cur = 0;
    matching_instance_->GetNumPositiveResults(pos_cur);
    matching_instance_->GetNumNegativeResults(neg_cur);
    positive_num_results_last = pos_cur;
    negative_num_results_last = neg_cur;
    unsafe_updates = unsafe_edges.size();

    // Handle vertex removals
    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v' && !u.is_add) {
            matching_instance_->RemoveVertex(u.id1);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_total_start).count() / 1000.0;

    std::cout << "[GPU-AllAtOnce] Complete:" << std::endl;
    std::cout << "  classify:    " << classify_ms << " ms" << std::endl;
    std::cout << "  add edges:   " << add_edges_ms << " ms (" << unsafe_edges.size() << " unsafe)" << std::endl;
    std::cout << "  build CSR:   " << csr_ms << " ms" << std::endl;
    std::cout << "  GPU search:  " << search_ms << " ms (" << total_matches << " matches)" << std::endl;
    std::cout << "  total tasks: " << unsafe_edges.size() * num_query_edges * 2 << std::endl;
    std::cout << "  total time:  " << total_ms << " ms" << std::endl;
}

// =====================================================================
// GPU BFS-level parallel inter-update
// =====================================================================

void InterExecutor::BatchUpdates_GPU_BFS(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    Graph& query = matching_instance_->GetQueryGraph();
    Graph& data = matching_instance_->GetDataGraph();

    // ---- Phase 1: Apply vertex updates & classify edge updates (CPU) ----
    auto t1 = std::chrono::high_resolution_clock::now();

    struct UnsafeEdge { uint v1, v2, label; };
    std::vector<uint8_t> update_type(update_size, 0);

    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v') {
            update_type[i] = 2;
            if (u.is_add) matching_instance_->AddVertex(u.id1, u.label);
            num_v_updates++;
        }
    }

    #pragma omp parallel for schedule(static, 512)
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] != 0) continue;
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type != 'e' || !u.is_add) { update_type[i] = 2; continue; }
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        update_type[i] = safe ? 0 : 1;
    }

    std::vector<UnsafeEdge> unsafe_edges;
    for (size_t i = 0; i < update_size; i++) {
        if (update_type[i] == 1) {
            const auto& u = data_graph_.updates_vec_[i];
            unsafe_edges.push_back({u.id1, u.id2, u.label});
        }
        if (data_graph_.updates_vec_[i].type == 'e') num_e_updates++;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double classify_ms = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;

    // ---- Phase 2: Add ALL unsafe edges to data graph (CPU serial) ----
    auto t3 = std::chrono::high_resolution_clock::now();
    for (const auto& e : unsafe_edges) {
        data_graph_.AddEdge(e.v1, e.v2, e.label);
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    double add_edges_ms = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0;

    // ---- Phase 3: Build CSR + setup query + orders on GPU (ONCE) ----
    auto t5 = std::chrono::high_resolution_clock::now();

    if (!gpu_bfs_ready_) {
        gpu_bfs_search_.SetupQuery(query);
    }
    gpu_bfs_search_.BuildCSR(data);

    const auto& orders = matching_instance_->GetMatchingOrders();
    uint32_t Q = query.NumVertices();
    uint32_t num_query_edges = static_cast<uint32_t>(orders.size());

    std::vector<uint32_t> flat_orders(num_query_edges * Q);
    for (uint32_t i = 0; i < num_query_edges; i++) {
        for (uint32_t j = 0; j < Q && j < orders[i].size(); j++) {
            flat_orders[i * Q + j] = orders[i][j];
        }
    }
    gpu_bfs_search_.SetAllMatchingOrders(flat_orders.data(), num_query_edges, Q);
    gpu_bfs_ready_ = true;

    auto t6 = std::chrono::high_resolution_clock::now();
    double csr_ms = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() / 1000.0;

    // ---- Phase 4: GPU BFS-level search ----
    auto t7 = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> ev1(unsafe_edges.size()), ev2(unsafe_edges.size()), elab(unsafe_edges.size());
    for (size_t i = 0; i < unsafe_edges.size(); i++) {
        ev1[i] = unsafe_edges[i].v1;
        ev2[i] = unsafe_edges[i].v2;
        elab[i] = unsafe_edges[i].label;
    }

    uint64_t total_matches = gpu_bfs_search_.SearchBatchEdgesBFS(
        ev1.data(), ev2.data(), elab.data(),
        static_cast<uint32_t>(unsafe_edges.size()),
        num_query_edges, Q
    );

    auto t8 = std::chrono::high_resolution_clock::now();
    double search_ms = std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count() / 1000.0;

    // ---- Phase 5: Aggregate results ----
    matching_instance_->AddPositiveResults(static_cast<size_t>(total_matches));

    size_t pos_cur = 0, neg_cur = 0;
    matching_instance_->GetNumPositiveResults(pos_cur);
    matching_instance_->GetNumNegativeResults(neg_cur);
    positive_num_results_last = pos_cur;
    negative_num_results_last = neg_cur;
    unsafe_updates = unsafe_edges.size();

    for (size_t i = 0; i < update_size; i++) {
        const auto& u = data_graph_.updates_vec_[i];
        if (u.type == 'v' && !u.is_add) {
            matching_instance_->RemoveVertex(u.id1);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_total_start).count() / 1000.0;

    std::cout << "[GPU-BFS] Complete:" << std::endl;
    std::cout << "  classify:    " << classify_ms << " ms" << std::endl;
    std::cout << "  add edges:   " << add_edges_ms << " ms (" << unsafe_edges.size() << " unsafe)" << std::endl;
    std::cout << "  build CSR:   " << csr_ms << " ms" << std::endl;
    std::cout << "  BFS search:  " << search_ms << " ms (" << total_matches << " matches)" << std::endl;
    std::cout << "  total time:  " << total_ms << " ms" << std::endl;
}

// =====================================================================
// GPU BFS inner-update: padded CSR + per-edge BFS search
// =====================================================================

void InterExecutor::BatchUpdates_GPU_BFS_Single(
    size_t& num_v_updates,
    size_t& num_e_updates,
    size_t& unsafe_updates,
    size_t& count,
    size_t& positive_num_results_last,
    size_t& negative_num_results_last,
    std::atomic_bool& reach_time_limit
) {
    auto t_total_start = std::chrono::high_resolution_clock::now();

    const size_t update_size = data_graph_.updates_vec_.size();
    if (update_size == 0) return;

    Graph& query = matching_instance_->GetQueryGraph();
    Graph& data = matching_instance_->GetDataGraph();

    // ---- One-time setup: build padded CSR + query + orders on GPU ----
    if (!gpu_bfs_ready_) {
        gpu_bfs_search_.SetupQuery(query);
        gpu_bfs_search_.BuildPaddedCSR(data);

        const auto& orders = matching_instance_->GetMatchingOrders();
        uint32_t Q = query.NumVertices();
        uint32_t num_query_edges = static_cast<uint32_t>(orders.size());

        std::vector<uint32_t> flat_orders(num_query_edges * Q);
        for (uint32_t i = 0; i < num_query_edges; i++) {
            for (uint32_t j = 0; j < Q && j < orders[i].size(); j++) {
                flat_orders[i * Q + j] = orders[i][j];
            }
        }
        gpu_bfs_search_.SetAllMatchingOrders(flat_orders.data(), num_query_edges, Q);
        gpu_bfs_ready_ = true;
    }

    uint32_t Q = query.NumVertices();
    const auto& orders = matching_instance_->GetMatchingOrders();
    uint32_t num_query_edges = static_cast<uint32_t>(orders.size());

    size_t total_gpu_searches = 0;
    size_t total_csr_rebuilds = 0;
    size_t total_safe_edges = 0;
    size_t total_flushes = 0;
    double time_classify_us = 0, time_add_us = 0, time_gpu_search_us = 0, time_csr_update_us = 0;

    // Mini-batch window: accumulate unsafe edges, search in one kernel launch
    static constexpr size_t WINDOW_SIZE = 500;
    std::vector<uint32_t> unsafe_v1, unsafe_v2, unsafe_label;
    unsafe_v1.reserve(WINDOW_SIZE);
    unsafe_v2.reserve(WINDOW_SIZE);
    unsafe_label.reserve(WINDOW_SIZE);

    printf("[GPU-BFS-Window] Starting: %zu updates, Q=%u, num_query_edges=%u, window=%zu\n",
           update_size, Q, num_query_edges, WINDOW_SIZE);
    fflush(stdout);

    // Track dirty vertices whose GPU CSR is stale (safe edges modify CPU graph but not GPU)
    std::unordered_set<uint32_t> dirty_vertices;
    size_t total_lazy_skips = 0;

    // Flush dirty vertices to GPU CSR
    auto flush_dirty = [&]() {
        if (dirty_vertices.empty()) return;
        auto tf0 = std::chrono::high_resolution_clock::now();
        for (uint32_t v : dirty_vertices) {
            uint32_t deg = data.GetDegree(v);
            // Check padding capacity
            if (v < gpu_bfs_search_.num_vertices_ && !gpu_bfs_search_.h_capacities_.empty()
                && deg <= gpu_bfs_search_.h_capacities_[v]) {
                // Incremental update single vertex
                const auto& nbrs = data.GetNeighbors(v);
                const auto& labs = data.GetNeighborLabels(v);
                gpu_bfs_search_.IncrementalUpdateVertex(v, nbrs.data(), labs.data(), deg);
            } else {
                // Need full rebuild
                gpu_bfs_search_.BuildPaddedCSR(data);
                total_csr_rebuilds++;
                dirty_vertices.clear();
                auto tf1 = std::chrono::high_resolution_clock::now();
                time_csr_update_us += std::chrono::duration_cast<std::chrono::microseconds>(tf1 - tf0).count();
                return;
            }
        }
        auto tf1 = std::chrono::high_resolution_clock::now();
        time_csr_update_us += std::chrono::duration_cast<std::chrono::microseconds>(tf1 - tf0).count();
        total_lazy_skips += dirty_vertices.size();
        dirty_vertices.clear();
    };

    // ---- Process updates one by one ----
    for (size_t i = 0; i < update_size && !reach_time_limit; i++) {
        const auto& u = data_graph_.updates_vec_[i];

        if (u.type == 'v') {
            if (u.is_add) {
                matching_instance_->AddVertex(u.id1, u.label);
            } else {
                matching_instance_->RemoveVertex(u.id1);
            }
            num_v_updates++;
            continue;
        }

        if (u.type != 'e' || !u.is_add) continue;
        num_e_updates++;

        // Classify
        auto tc0 = std::chrono::high_resolution_clock::now();
        bool safe = matching_instance_->Classify(u.id1, u.id2, u.label);
        auto tc1 = std::chrono::high_resolution_clock::now();
        time_classify_us += std::chrono::duration_cast<std::chrono::microseconds>(tc1 - tc0).count();

        // Add edge to CPU data graph
        auto ta0 = std::chrono::high_resolution_clock::now();
        data_graph_.AddEdge(u.id1, u.id2, u.label);
        auto ta1 = std::chrono::high_resolution_clock::now();
        time_add_us += std::chrono::duration_cast<std::chrono::microseconds>(ta1 - ta0).count();

        if (safe) {
            // Safe edge: only update CPU graph, mark vertices as dirty for lazy GPU sync
            dirty_vertices.insert(u.id1);
            dirty_vertices.insert(u.id2);
            total_safe_edges++;
            continue;
        }

        // Unsafe edge: accumulate for mini-batch GPU BFS
        dirty_vertices.insert(u.id1);
        dirty_vertices.insert(u.id2);
        unsafe_v1.push_back(u.id1);
        unsafe_v2.push_back(u.id2);
        unsafe_label.push_back(u.label);
        unsafe_updates++;

        // When window is full, flush and search
        if (unsafe_v1.size() >= WINDOW_SIZE) {
            total_flushes++;
            // Flush all dirty vertices to GPU CSR
            flush_dirty();

            // GPU BFS search for all accumulated unsafe edges at once
            auto ts0 = std::chrono::high_resolution_clock::now();
            uint64_t matches = gpu_bfs_search_.SearchBatchEdgesBFS(
                unsafe_v1.data(), unsafe_v2.data(), unsafe_label.data(),
                static_cast<uint32_t>(unsafe_v1.size()),
                num_query_edges, Q
            );
            auto ts1 = std::chrono::high_resolution_clock::now();
            time_gpu_search_us += std::chrono::duration_cast<std::chrono::microseconds>(ts1 - ts0).count();

            if (matches > 0) {
                matching_instance_->AddPositiveResults(static_cast<size_t>(matches));
            }
            total_gpu_searches++;

            if (total_gpu_searches <= 5 || total_gpu_searches % 20 == 0) {
                printf("[GPU-BFS-Window] batch #%zu: %zu edges → %lu matches, gpu_time=%.1fms\n",
                       total_gpu_searches, unsafe_v1.size(), matches,
                       std::chrono::duration_cast<std::chrono::microseconds>(ts1 - ts0).count() / 1000.0);
                fflush(stdout);
            }

            unsafe_v1.clear();
            unsafe_v2.clear();
            unsafe_label.clear();
        }

        // Update result counters
        size_t pos_cur = 0, neg_cur = 0;
        matching_instance_->GetNumPositiveResults(pos_cur);
        matching_instance_->GetNumNegativeResults(neg_cur);
        positive_num_results_last = pos_cur;
        negative_num_results_last = neg_cur;
    }

    // Flush remaining unsafe edges
    if (!unsafe_v1.empty()) {
        flush_dirty();
        auto ts0 = std::chrono::high_resolution_clock::now();
        uint64_t matches = gpu_bfs_search_.SearchBatchEdgesBFS(
            unsafe_v1.data(), unsafe_v2.data(), unsafe_label.data(),
            static_cast<uint32_t>(unsafe_v1.size()),
            num_query_edges, Q
        );
        auto ts1 = std::chrono::high_resolution_clock::now();
        time_gpu_search_us += std::chrono::duration_cast<std::chrono::microseconds>(ts1 - ts0).count();
        if (matches > 0) {
            matching_instance_->AddPositiveResults(static_cast<size_t>(matches));
        }
        total_gpu_searches++;
        total_flushes++;
        printf("[GPU-BFS-Window] final batch: %zu edges → %lu matches, gpu_time=%.1fms\n",
               unsafe_v1.size(), matches,
               std::chrono::duration_cast<std::chrono::microseconds>(ts1 - ts0).count() / 1000.0);
        fflush(stdout);

        size_t pos_cur = 0, neg_cur = 0;
        matching_instance_->GetNumPositiveResults(pos_cur);
        matching_instance_->GetNumNegativeResults(neg_cur);
        positive_num_results_last = pos_cur;
        negative_num_results_last = neg_cur;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_total_start).count() / 1000.0;

    std::cout << "[GPU-BFS-Window] Complete:" << std::endl;
    std::cout << "  classify:      " << time_classify_us / 1000.0 << " ms" << std::endl;
    std::cout << "  add edges:     " << time_add_us / 1000.0 << " ms" << std::endl;
    std::cout << "  CSR update:    " << time_csr_update_us / 1000.0 << " ms (" << total_csr_rebuilds << " rebuilds)" << std::endl;
    std::cout << "  GPU BFS:       " << time_gpu_search_us / 1000.0 << " ms (" << total_gpu_searches << " batches, " << unsafe_updates << " unsafe edges)" << std::endl;
    std::cout << "  safe edges:    " << total_safe_edges << " (lazy, " << total_lazy_skips << " dirty vertices flushed in " << total_flushes << " flushes)" << std::endl;
    std::cout << "  window size:   " << WINDOW_SIZE << std::endl;
    std::cout << "  total time:    " << total_ms << " ms" << std::endl;
    std::cout << "  pos matches:   " << positive_num_results_last << std::endl;
}
