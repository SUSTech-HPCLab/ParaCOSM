# InnerExecutor Class

## Overview

The `InnerExecutor` class is a **direct implementation** that inherits from the `matching` base class, providing a simplified interface for graph matching operations. It maintains the inheritance hierarchy from the existing matching architecture while offering enhanced functionality.

## Features

- **Inheritance-based Design**: Directly inherits from `matching` base class, maintaining architectural consistency
- **Simplified Interface**: Easy-to-use methods for common graph matching operations
- **Error Handling**: Comprehensive exception handling with meaningful error messages
- **Resource Management**: Automatic management of data structures and thread pools
- **Flexible Configuration**: Configurable parameters for threads, auto-tuning, and output options
- **Parallel Processing**: Built-in support for multi-threaded graph operations

## Class Hierarchy

```
matching (base class)
    ↓
InnerExecutor (derived class)
```

### Inherited Members from `matching`

- **Graph References**: `query_`, `data_`
- **Configuration**: `max_num_results_`, `print_preprocessing_results_`, `print_enumeration_results_`, `homomorphism_`
- **Execution Info**: `visited_`, `num_initial_results_`, `num_positive_results_`, `num_negative_results_`
- **Base Methods**: `Preprocessing()`, `InitialMatching()`, `AddEdge()`, `RemoveEdge()`, `AddVertex()`, `RemoveVertex()`

### Additional Members in `InnerExecutor`

- **Configuration**: `orders_`, `num_threads_`, `auto_tuning_`
- **Data Structures**: `eidx_`, `DCS_`, `d1`, `d2`, `n1`, `np1`, `n2`, `nc2`
- **Threading**: `thread_pool_`, `local_vec_m_`, `local_vec_visited_local_`
- **Job Management**: `job_queue_`, `vertex_vector_`

## Constructor Parameters

```cpp
InnerExecutor(
    Graph& query_graph,           // Query graph to match
    Graph& data_graph,            // Data graph to search in
    uint max_num_results = 10000000, // Maximum results to find
    bool print_prep = false,      // Print preprocessing info
    bool print_enum = false,      // Print enumeration info
    bool homomorphism = false,    // Use homomorphic matching
    std::vector<std::vector<uint>> orders = {}, // Vertex processing orders
    size_t num_threads = 8,       // Number of threads
    size_t auto_tuning = 0        // Auto-tuning flag
);
```

## Public Methods

### Overridden Base Class Methods
- `void Preprocessing() override`: Initialize data structures and build indices
- `void InitialMatching() override`: Set up initial matching state
- `void AddEdge(uint v1, uint v2, uint label) override`: Add edge with additional logic
- `void RemoveEdge(uint v1, uint v2) override`: Remove edge with additional logic
- `void AddVertex(uint id, uint label) override`: Add vertex with additional logic
- `void RemoveVertex(uint id) override`: Remove vertex with additional logic
- `void GetMemoryCost(size_t& num_edges, size_t& num_vertices) override`: Get memory usage with InnerExecutor costs
- `bool Classify(uint v1, uint v2, uint label) override`: Implement classification logic

### New Methods
- `size_t Execute()`: Execute the main matching algorithm
- `size_t ExecuteParaCOSMKernel2(...)`: Execute the ParaCOSM Kernel2 algorithm
- `size_t UpdateAndFind(uint v1, uint v2, uint label)`: Update graph and find matches
- `void AddEdgeAsync(uint v1, uint v2, uint label)`: Add edge asynchronously
- `size_t GetNumThreads() const`: Get current thread count
- `void SetNumThreads(size_t num_threads)`: Set thread count

### Private Helper Methods
- `void InitializeDataStructures()`: Initialize all data structures
- `void BuildEdgeIndex()`: Build edge index mapping
- `void BuildDCS()`: Build Dynamic Candidate Set structures
- `size_t ProcessVertexLayer(...)`: Process vertex layer for parallel execution

## Usage Example

```cpp
#include "core/inner_executor.h"
#include "graph_storage/graph.h"

int main() {
    // Create your graphs
    Graph query_graph, data_graph;
    
    // Create executor (inherits from matching)
    InnerExecutor executor(query_graph, data_graph, 10000000, false, false, false, {}, 8, 0);
    
    // Call inherited base class methods
    executor.Preprocessing();
    executor.InitialMatching();
    
    // Execute main algorithm
    size_t results = executor.Execute();
    
    // Use ParaCOSM Kernel2 directly
    uint depth = 0;
    std::vector<uint> m, extendable;
    size_t kernel_results = executor.ExecuteParaCOSMKernel2(depth, m, extendable);
    
    // Use inherited graph modification methods
    executor.AddEdge(1, 2, 1);
    
    // Get execution statistics (inherited from matching)
    size_t initial_results, positive_results, negative_results;
    executor.GetNumInitialResults(initial_results);
    executor.GetNumPositiveResults(positive_results);
    executor.GetNumNegativeResults(negative_results);
    
    return 0;
}
```

## Benefits of Inheritance-based Design

1. **Architectural Consistency**: Follows the existing codebase structure
2. **Code Reuse**: Inherits all base functionality from `matching`
3. **Polymorphism**: Can be used anywhere a `matching` object is expected
4. **Extensibility**: Easy to add new matching algorithms by inheriting from `matching`
5. **Maintainability**: Changes to base class automatically propagate to derived classes

## Dependencies

- **Base Class**: `matching` from `matching_executor/matching.h`
- **Graph**: `Graph` class from `graph_storage/graph.h`
- **Threading**: `ThreadPool` from `graph/Threadpool.cpp`
- **Concurrency**: `tbb::concurrent_queue` from Intel TBB
- **Data Structures**: `ska::flat_hash_map` from `graph_storage/storage_hash_map.hpp`
- **Standard C++**: vector, queue, unordered_map, algorithm, omp.h

## Implementation Notes

- **Data Structure Initialization**: All data structures are initialized in the constructor
- **Thread Safety**: Uses thread-local storage for parallel operations
- **Memory Management**: Automatic cleanup through RAII and smart pointers
- **Error Handling**: Comprehensive exception handling throughout all methods
- **Placeholder Methods**: Some methods (like `Execute()`) contain placeholder implementations that you can customize

## Future Enhancements

- Implement full ParaCOSM Kernel2 algorithm
- Add more sophisticated candidate filtering
- Optimize data structure access patterns
- Add performance profiling capabilities
- Implement adaptive thread management
