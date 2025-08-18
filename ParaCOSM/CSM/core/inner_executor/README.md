
# InnerExecutor and Kernel Functions

### Overview

The `core/inner_executor` module provides the parallelized subgraph matching executor `InnerExecutor`, along with several directly callable matching kernels (SymBi / GraphFlow / TurboFlux / ParaCOSM Kernel). This module can be used either as a unified external interface (inheriting from the `matching` base class) or at the kernel level for customized parallel enumeration.

### Key Features

* **Inheritance-based design**: `InnerExecutor` directly inherits from the `matching` base class, reusing existing execution statistics, graph interfaces, and configurations.
* **Multi-kernel support**: Supports parallelized versions of SymBi, GraphFlow, TurboFlux, and ParaCOSM Kernel.
* **Thread safety & parallelism**: Uses OpenMP for task-level parallelism and `tbb::concurrent_queue` for task splitting and dynamic load balancing.
* **Configurable**: Number of threads, homomorphism mode, auto-tuning, vertex ordering, and more are configurable.
* **Engineering-friendly**: Provides out-of-the-box data structure initialization, resource management, and exception handling.

### Directory Structure

* `inner_executor.h/.cpp`: External interface and implementation of `InnerExecutor`.
* `FindMatchesKernel.h/.cpp`: Parallel matching kernels and auxiliary routines.
* `FindMatchesKernel_simplified.cpp`: Simplified/sample implementation for readability and validation.

### Classes and Functions

#### Class: `InnerExecutor` (in `core/inner_executor/inner_executor.h`)

* Constructor parameters:

  * `Graph& query_graph, Graph& data_graph`
  * `uint max_num_results`
  * `bool print_prep, print_enum, homomorphism`
  * `std::vector<std::vector<uint>> orders`
  * `size_t num_threads, auto_tuning`

* Core methods (override or add-on):

  * `void Preprocessing()`: Build indices and candidate structures, initialize thread-local data.
  * `void InitialMatching()`: Initialize match vectors and statistics.
  * `size_t Execute()`: Main execution entry (extendable with the real solving logic).
  * `size_t ExecuteParaCOSMKernel2(uint depth, std::vector<uint>& m, std::vector<uint>& extendable)`: Example entry to directly call ParaCOSM Kernel2.
  * Graph modification interface (override from base class): `AddEdge/RemoveEdge/AddVertex/RemoveVertex`.
  * Statistics/utilities: `GetMemoryCost(...)`, `Classify(...)`.

* Key data structures (excerpt):

  * `eidx_`: Query graph edge index mapping.
  * `DCS_ / d1 / d2 / n1 / np1 / n2 / nc2`: Candidate and constraint-related indices.
  * `local_vec_m_ / local_vec_visited_local_`: Thread-local match/visit states.
  * `vertex_vector_ / job_queue_`: Used for task splitting and parallel execution.

#### Kernel Functions (in `core/inner_executor/FindMatchesKernel.h/.cpp`)

* **SymBi Kernel**

  * `void Parallel_Symbi_FindMatches(...)`
  * `void Parallel_Symbi_FindMatches_Simplified(...)`
  * `void Parallel_Symbi_FindMatches_ParaCOSM_Kernel(uint depth, std::vector<uint>& m, std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results)`

* **GraphFlow Kernel**

  * `void Parallel_Graphflow_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)`

* **TurboFlux Kernel**

  * `void Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)`

* **Typical auxiliary routines (excerpt)**

  * `Parallel_FindMatches_local3_simplified(...)`
  * `process_vertex_layer1_simplified(...)`

All kernel flows follow a similar structure: **candidate retrieval → constraint filtering (indices/connectivity/visited) → recursive expansion / task splitting → backtracking**, with task queue splitting enabled at key branching points to improve parallel throughput.

### Usage Examples

#### 1) Using `InnerExecutor` as a unified executor

```cpp
#include "core/inner_executor/inner_executor.h"
#include "graph_storage/graph.h"

int main() {
    Graph query_graph, data_graph;

    InnerExecutor executor(
        query_graph,
        data_graph,
        /* max_num_results = */ 10000000,
        /* print_prep = */ false,
        /* print_enum = */ false,
        /* homomorphism = */ false,
        /* orders = */ {},
        /* num_threads = */ 8,
        /* auto_tuning = */ 0
    );

    executor.Preprocessing();
    executor.InitialMatching();
    size_t results = executor.Execute();

    return 0;
}
```

#### 2) Directly calling a simplified kernel (SymBi Simplified)

```cpp
#include "core/inner_executor/FindMatchesKernel.h"
#include "graph_storage/graph.h"

void run_kernel(const Graph& query, const Graph& data) {
    std::vector<uint> m(query.NumVertices(), UNMATCHED);
    std::vector<Parrllel_SymBi::ExtendableVertex> extendable(query.NumVertices());

    size_t num_results = 0;
    size_t num_threads = 8;
    bool homomorphism = false;

    // These structures should be constructed and filled during preprocessing
    std::vector<std::vector<uint>> eidx(query.NumVertices(), std::vector<uint>(query.NumVertices()));
    std::vector<ska::flat_hash_map<uint, std::vector<uint>>> DCS(query.NumEdges() * 2);
    std::vector<ska::flat_hash_map<uint, bool>> d2(query.NumVertices());
    std::vector<ska::flat_hash_map<uint, uint>> n2(query.NumEdges() * 2);
    std::vector<std::vector<bool>> local_visited(num_threads, std::vector<bool>(data.NumVertices(), false));
    std::vector<std::vector<uint>> treeNode_neighbors(query.NumVertices());

    Parallel_Symbi_FindMatches_Simplified(
        /* depth */ 0,
        m,
        extendable,
        num_results,
        query,
        data,
        eidx,
        DCS,
        d2,
        n2,
        local_visited,
        homomorphism,
        num_threads,
        treeNode_neighbors
    );
}
```

### Threads and Task Scheduling

* OpenMP is used for candidate-level and layer-level parallelism: `#pragma omp parallel for ...`.
* `tbb::concurrent_queue` is used for task splitting (branches beyond a depth threshold are enqueued, tail phases are consumed opportunistically).
* Thread-local state (`local_vec_m_`, `local_vec_visited_local_`) reduces lock contention.
* It is recommended to tune `num_threads` according to data size, and enable task splitting if tail utilization is low.

### Configuration and Recommendations

* **Threads**: `InnerExecutor(num_threads=...)` or kernel’s `NUMTHREAD` argument.
* **Homomorphism**: `homomorphism = true/false` (when `false`, duplicate-visit checks are enforced).
* **Vertex order**: Controlled via `orders_` or kernel’s `order_index`.
* **Preprocessing**: Ensure `eidx / DCS / d2 / n2 / treeNode_neighbors` are consistent with the graph and sorted to support `lower_bound`-based filtering.

### Dependencies

* OpenMP (parallel programming).
* Intel TBB: `tbb::concurrent_queue`.
* `ska::flat_hash_map` (see `graph_storage/storage_hash_map.hpp`).
* Graph & matching base class: `graph_storage/graph.h`, `matching_executor/matching.h`.

### Build

* Enable this module in the top-level or `CSM/` `CMakeLists.txt` (ensure OpenMP/TBB are available).
* Typical flags: `-fopenmp`, link with TBB; see project CMake configuration for platform-specific details.

### Implementation Notes

* Candidate retrieval primarily uses pre-built `DCS` and `eidx`; adjacency-based `lower_bound` intersections are used as fallback.
* Visit deduplication uses thread-local boolean arrays, restored in pairs during backtracking.
* Task splitting only occurs beyond a specified depth to avoid oversplitting into excessively fine-grained tasks.


* Complete the actual solving and tuning logic in `Execute()`.
* Enrich/optimize candidate index construction (`BuildDCS()`) and edge indices (`BuildEdgeIndex()`).
* Add performance profiling tools and refine adaptive thread/task-splitting strategies.

