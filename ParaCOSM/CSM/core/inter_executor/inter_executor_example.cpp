#include "inter_executor.h"
#include <iostream>

/**
 * @brief Example usage of InterExecutor class
 * 
 * This file demonstrates how to use the InterExecutor class
 * to process batch updates efficiently.
 */

void ExampleUsage() {
    // Example: How to use InterExecutor
    
    // 1. Create a data graph and matching instance
    // Graph data_graph;
    // matching* matching_instance = new SomeMatchingClass(query_graph, data_graph);
    
    // 2. Create InterExecutor instance
    // InterExecutor executor(data_graph, matching_instance);
    
    // 3. Initialize counters
    // size_t num_v_updates = 0, num_e_updates = 0, unsafe_updates = 0;
    // size_t count = 0, positive_num_results_last = 0, negative_num_results_last = 0;
    // std::atomic_bool reach_time_limit(false);
    
    // 4. Process batch updates
    // executor.ProcessBatchUpdates(
    //     num_v_updates, num_e_updates, unsafe_updates, count,
    //     positive_num_results_last, negative_num_results_last,
    //     reach_time_limit
    // );
    
    // 5. Or use the legacy function name for backward compatibility
    // executor.BatchUpdates3(
    //     num_v_updates, num_e_updates, unsafe_updates, count,
    //     positive_num_results_last, negative_num_results_last,
    //     reach_time_limit
    // );
    
    std::cout << "InterExecutor example usage demonstrated." << std::endl;
    std::cout << "Please uncomment the code above and provide actual Graph and matching instances." << std::endl;
}

// Legacy function wrapper for backward compatibility
// This maintains the exact same interface as the original BatchUpdates3 function
inline void BatchUpdates3(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {
    
    InterExecutor executor(data_graph, mm);
    executor.BatchUpdates3(
        num_v_updates, num_e_updates, unsafe_updates, count,
        positive_num_results_last, negative_num_results_last,
        reach_time_limit
    );
}
