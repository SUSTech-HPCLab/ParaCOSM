#ifndef CORE_INTER_EXECUTOR_H
#define CORE_INTER_EXECUTOR_H

#include <vector>
#include <atomic>
#include "matching_executor/matching.h"
#include "graph_storage/graph.h"

// Forward declaration
struct InsertUnit;

/**
 * @brief InterExecutor class for handling batch updates with sliding window approach
 * 
 * This class provides efficient batch processing of graph updates using a sliding window
 * mechanism to identify safe and unsafe updates, optimizing the update process.
 */
class InterExecutor {
private:
    Graph& data_graph_;
    matching* matching_instance_;
    
    // Configuration constants
    static constexpr size_t DEFAULT_WINDOW_SIZE = 16;
    
public:
    /**
     * @brief Constructor for InterExecutor
     * 
     * @param data_graph Reference to the data graph
     * @param matching_instance Pointer to the matching instance
     */
    InterExecutor(Graph& data_graph, matching* matching_instance);
    
    /**
     * @brief Destructor
     */
    ~InterExecutor() = default;
    
    /**
     * @brief Process batch updates using sliding window approach
     * 
     * This method processes graph updates in batches using a sliding window mechanism.
     * It identifies safe updates that can be processed in parallel and unsafe updates
     * that require immediate processing to maintain consistency.
     * 
     * @param num_v_updates Reference to count vertex updates
     * @param num_e_updates Reference to count edge updates
     * @param unsafe_updates Reference to count unsafe updates
     * @param count Reference to total update count
     * @param positive_num_results_last Reference to last positive results count
     * @param negative_num_results_last Reference to last negative results count
     * @param reach_time_limit Reference to time limit flag
     * @param window_size Size of the sliding window (default: 16)
     */
    void ProcessBatchUpdates(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit,
        size_t window_size = DEFAULT_WINDOW_SIZE
    );
    
    /**
     * @brief Process a single update safely
     * 
     * @param update The update to process
     * @return true if the update is safe, false otherwise
     */
    bool ProcessSingleUpdate(const InsertUnit& update);
    
    /**
     * @brief Apply an unsafe update immediately
     * 
     * @param update The unsafe update to apply
     */
    void ApplyUnsafeUpdate(const InsertUnit& update);
    
    /**
     * @brief Check if an update is safe
     * 
     * @param update The update to check
     * @return true if safe, false otherwise
     */
    bool IsUpdateSafe(const InsertUnit& update);
    
    /**
     * @brief Legacy function name for backward compatibility
     * 
     * This is an alias for ProcessBatchUpdates to maintain compatibility
     * with existing code that might call BatchUpdates3.
     */
    void BatchUpdates3(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    ) {
        ProcessBatchUpdates(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last,
            reach_time_limit, DEFAULT_WINDOW_SIZE
        );
    }

    /**
     * @brief Batch processing strategy 2 (queue-based), ported from free function BatchUpdates2
     */
    void BatchUpdates2(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Batch processing strategy 1 (windowed classification), ported from free function BatchUpdates
     */
    void BatchUpdates(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief OpenMP windowed processing over `updates_vec_` (port of BatchUpdates_OpenMP)
     */
    void BatchUpdates_OpenMP(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Queue-based batch processing (port of free function ProcessBatchUpdates)
     */
    void ProcessBatchUpdatesQueue(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Single-threaded update application (port of SingleThreadUpdate)
     */
    void SingleThreadUpdate(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Taskflow-based parallel batch processing
     * 
     * Uses Taskflow library to parallelize the classification and processing
     * of graph updates. Similar to BatchUpdates2 but uses Taskflow's task
     * graph for better parallel scheduling.
     */
    void BatchUpdates_Taskflow(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Sequential batch processing over updates_vec_ (BatchUpdates4).
     *
     * 这个策略按顺序遍历 `data_graph_.updates_vec_` 中的所有更新，
     * 对每条更新直接调用匹配器的 AddEdge / RemoveEdge /
     * AddVertex / RemoveVertex 接口，并在每次更新后统计结果变化，
     * 以统计「unsafe updates」。该实现是单线程的，适合作为
     * 在引入 Taskflow 并行之前的基线版本。
     */
    void BatchUpdates4(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit
    );

    /**
     * @brief Persistent OMP pool mode: single parallel region wrapping entire update loop.
     *
     * Avoids per-update OMP fork/join by keeping threads alive across the entire
     * update stream. The master thread runs the sliding-window classification loop,
     * and inner Parallel_FindMatches2 reuses the existing thread team.
     */
    void BatchUpdates_Persistent(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit,
        size_t num_threads = 8
    );

    /**
     * @brief Batch-all mode: pre-classify ALL updates, add unsafe edges to graph,
     * then enumerate matches for ALL unsafe edges in parallel across threads.
     *
     * This achieves inter-update parallelism: multiple unsafe edges' match
     * enumeration runs concurrently on different threads, each using thread-local
     * visited arrays via EnumerateNewEdge(). Note: match counts may differ
     * slightly from serial processing because all unsafe edges are visible in
     * the graph simultaneously.
     */
    void BatchUpdates_AllAtOnce(
        size_t& num_v_updates,
        size_t& num_e_updates,
        size_t& unsafe_updates,
        size_t& count,
        size_t& positive_num_results_last,
        size_t& negative_num_results_last,
        std::atomic_bool& reach_time_limit,
        size_t num_threads = 8
    );

    /**
     * DegreePruning
     *
     * Check degree-based pruning constraints for mapping query vertices (u1, u2)
     * to data vertices (v1, v2). Uses degrees in the provided graphs to quickly
     * rule out infeasible mappings.
     *
     * Parameters:
     * - v1, v2: candidate data vertices
     * - u1, u2: corresponding query vertices
     * - data_: reference to the data graph
     * - query_: reference to the query graph
     */
    bool DegreePruning(uint v1, uint v2, uint u1, uint u2, Graph& data_, Graph& query_);

    /**
     * QueryGraphPruning
     *
     * Perform structural pruning using query-graph constraints for mapping
     * (u1, u2) -> (v1, v2). Typically considers adjacency/edge existence and
     * other structural relationships present in the query graph.
     *
     * Parameters:
     * - v1, v2: candidate data vertices
     * - u1, u2: corresponding query vertices
     * - query_: reference to the query graph
     * - data_: reference to the data graph
     */
    bool QueryGraphPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_);

    /**
     * LabelPruning
     *
     * Check label compatibility for mapping (u1, u2) -> (v1, v2) by verifying
     * vertex and/or edge labels satisfy query label constraints.
     *
     * Parameters:
     * - v1, v2: candidate data vertices
     * - u1, u2: corresponding query vertices
     * - query_: reference to the query graph
     * - data_: reference to the data graph
     */
    bool LabelPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_);

    /**
     * DCSPruning
     *
     * Apply Degree-Constraint Subgraph (DCS) or similar index-based pruning to
     * eliminate infeasible mappings early for (u1, u2) -> (v1, v2).
     *
     * Parameters:
     * - v1, v2: candidate data vertices
     * - u1, u2: corresponding query vertices
     * - query_: reference to the query graph
     * - data_: reference to the data graph
     */
    bool DCSPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_);

    /**
     * FMPathPruning
     *
     * Prune using feature/path-based constraints (e.g., frequent motif or
     * path-consistency checks) between (u1, u2) and (v1, v2).
     *
     * Parameters:
     * - v1, v2: candidate data vertices
     * - u1, u2: corresponding query vertices
     * - query_: reference to the query graph
     * - data_: reference to the data graph
     */
    bool FMPathPruning(uint v1, uint v2, uint u1, uint u2, Graph& query_, Graph& data_);
};

#endif // CORE_INTER_EXECUTOR_H


