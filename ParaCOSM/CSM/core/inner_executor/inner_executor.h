#ifndef CORE_INNER_EXECUTOR_H
#define CORE_INNER_EXECUTOR_H

#include <vector>
#include <memory>
#include <queue>
#include <unordered_map>
#include "matching_executor/matching.h"
#include "graph_storage/Threadpool.cpp"
#include <tbb/concurrent_queue.h>
#include "graph_storage/storage_hash_map.hpp"

/**
 * @brief InnerExecutor class that inherits from matching base class
 * 
 * This class provides a simplified interface for graph matching operations
 * while maintaining the inheritance hierarchy from the matching base class.
 * It can be used as a drop-in replacement for other matching algorithms.
 */
class InnerExecutor : public matching {
private:
    // Inherited from matching base class:
    // - query_, data_, max_num_results_, print_preprocessing_results_, print_enumeration_results_, homomorphism_
    // - visited_, num_initial_results_, num_positive_results_, num_negative_results_, etc.

    // Additional configuration
    std::vector<std::vector<uint>> orders_;
    size_t num_threads_;
    size_t auto_tuning_;

    // Data structures for the matching algorithm
    std::vector<std::vector<uint>> eidx_;
    std::vector<std::vector<ska::flat_hash_map<uint, std::vector<uint>>>> DCS_;
    std::vector<ska::flat_hash_map<uint, bool>> d1, d2;
    std::vector<ska::flat_hash_map<uint, uint>> n1, np1, n2, nc2;
    
    // Threading support
    ThreadPool thread_pool_;
    std::vector<std::vector<uint>> local_vec_m_;
    std::vector<std::vector<bool>> local_vec_visited_local_;
    
    // Job queue for parallel processing
    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>, uint, uint>> job_queue_;
    
    // Vertex processing structures
    std::vector<std::tuple<uint, uint, size_t, std::vector<uint>, uint, uint>> vertex_vector_;

public:
    /**
     * @brief Constructor for InnerExecutor
     * 
     * @param query_graph Reference to the query graph
     * @param data_graph Reference to the data graph
     * @param max_num_results Maximum number of results to find
     * @param print_prep Whether to print preprocessing information
     * @param print_enum Whether to print enumeration information
     * @param homomorphism Whether to use homomorphic matching
     * @param orders Predefined vertex processing orders
     * @param num_threads Number of threads to use
     * @param auto_tuning Whether to enable auto-tuning
     */
    InnerExecutor(
        Graph& query_graph,
        Graph& data_graph,
        uint max_num_results = 10000000,
        bool print_prep = false,
        bool print_enum = false,
        bool homomorphism = false,
        std::vector<std::vector<uint>> orders = {},
        size_t num_threads = 8,
        size_t auto_tuning = 0
    );

    /**
     * @brief Destructor
     */
    ~InnerExecutor() override = default;

    // Override base class methods
    void Preprocessing() override;
    void InitialMatching() override;
    
    // Override graph modification methods
    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    // Override utility methods
    void GetMemoryCost(size_t& num_edges, size_t& num_vertices) override;
    bool Classify(uint v1, uint v2, uint label) override;

    /**
     * @brief Execute the main matching algorithm
     * 
     * @return Number of matches found
     */
    size_t Execute();

    /**
     * @brief Execute the ParaCOSM Kernel2 algorithm directly
     * 
     * @param depth Current depth in the search tree
     * @param m Current matching vector
     * @param extendable Extendable vertices information
     * @return Number of results found
     */
    size_t ExecuteParaCOSMKernel2(
        uint depth,
        std::vector<uint>& m,
        std::vector<uint>& extendable
    );

    /**
     * @brief Update graph and find matches
     * 
     * @param v1 First vertex of the edge
     * @param v2 Second vertex of the edge
     * @param label Label of the edge
     * @return Number of matches found after update
     */
    size_t UpdateAndFind(uint v1, uint v2, uint label);

    /**
     * @brief Add edge asynchronously
     * 
     * @param v1 First vertex of the edge
     * @param v2 Second vertex of the edge
     * @param label Label of the edge
     */
    void AddEdgeAsync(uint v1, uint v2, uint label);

    /**
     * @brief Get the number of threads used
     * 
     * @return Number of threads
     */
    size_t GetNumThreads() const { return num_threads_; }

    /**
     * @brief Set the number of threads
     * 
     * @param num_threads New number of threads
     */
    void SetNumThreads(size_t num_threads) { num_threads_ = num_threads; }

private:
    /**
     * @brief Initialize data structures for the matching algorithm
     */
    void InitializeDataStructures();

    /**
     * @brief Build edge index mapping
     */
    void BuildEdgeIndex();

    /**
     * @brief Build DCS (Dynamic Candidate Set) structures
     */
    void BuildDCS();

    /**
     * @brief Process vertex layer for parallel execution
     * 
     * @param u Query vertex
     * @param u_min Minimum query vertex
     * @param v_idx Data vertex index
     * @param m Matching vector
     * @param depth Current depth
     * @param thread_id Thread ID
     * @return Number of results found
     */
    size_t ProcessVertexLayer(
        uint u, uint u_min, size_t v_idx,
        std::vector<uint>& m, uint depth, size_t thread_id
    );


};

#endif // CORE_INNER_EXECUTOR_H
