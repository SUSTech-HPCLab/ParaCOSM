# Parallel Utils: Practical Utilities for Parallel Programming

This module provides a collection of practical utilities, test cases, and micro-benchmarks for parallel programming in C++. The focus is on building, testing, and understanding core parallel components such as pipelines, hash maps, vectors, and cache-aware algorithms. The codebase is organized to help learners and practitioners experiment with and evaluate parallel programming techniques and data structures.

---

## 📂 Directory Structure

- **test_component/**: Component-level utilities and micro-benchmarks
  - **pipeline/**: Pipeline processing and buffer management
  - **hashmap/**: Parallel hash map implementations and tests
  - **vector/**: Vector operations and hash functions for vectors
  - `test_DCS.cpp`, `test_omp_getthread.cpp`, `test_openmp.cpp`: Miscellaneous component tests, including OpenMP thread queries and custom data structures

- **test_parallel/**: Parallel algorithm implementations and benchmarks
  - `DFS_omp.cpp`, `DFS_omp2.cpp`, `DFS_omp3.cpp`, `DFS_omp4.cpp`: Variants of parallel depth-first search using OpenMP
  - `C17backtrack.cpp`, `backtrack.cpp`: Backtracking algorithms with and without parallelism

- **test_cache/**: Cache-aware algorithm tests
  - `cache_bfs.cpp`: Breadth-first search optimized for cache performance

---

## 🧩 Component Highlights

### test_component

- **pipeline/**:  
  - `test_pipeline.cpp`, `test_pipeline_buffer.cpp`: Demonstrate pipeline parallelism and buffer management, useful for streaming and staged computation scenarios.
- **hashmap/**:  
  - `test_hashmap.cpp`, `test_hash_vect.cpp`: Test and benchmark parallel hash map implementations.
  - `hashmap_result.md`: Documents results and observations from hash map experiments.
- **vector/**:  
  - `test_vector.cpp`, `test_vec_hash.cpp`: Explore vector operations and custom hash functions for vector data.
- **Miscellaneous**:  
  - `test_DCS.cpp`: Custom data structure or algorithm test.
  - `test_omp_getthread.cpp`, `test_openmp.cpp`: OpenMP thread management and parallel region tests.

### test_parallel

- **Parallel DFS**:  
  - Multiple implementations (`DFS_omp.cpp` to `DFS_omp4.cpp`) show different strategies for parallelizing depth-first search using OpenMP.
- **Backtracking**:  
  - `C17backtrack.cpp`, `backtrack.cpp`: Compare serial and parallel backtracking approaches.

### test_cache

- **Cache-Aware BFS**:  
  - `cache_bfs.cpp`: Implements and tests a cache-optimized breadth-first search, illustrating the impact of memory access patterns on performance.

---

## 🚀 Getting Started

### Prerequisites

- C++17 or later
- OpenMP support (for parallel and thread-related tests)

### Compilation Example

```bash
# Compile a parallel DFS example
g++ DFS_omp.cpp -o DFS_omp -fopenmp

# Compile a pipeline test
g++ test_pipeline.cpp -o test_pipeline -std=c++11

# Compile cache-aware BFS
g++ cache_bfs.cpp -o cache_bfs -fopenmp
```

---

## 📊 Learning Goals

- Understand and experiment with parallel data structures and algorithms
- Benchmark and compare different parallelization strategies
- Explore the impact of cache and memory access patterns on performance
- Learn practical OpenMP usage for real-world problems

---

## 🤝 Contributing

Contributions, bug reports, and suggestions are welcome! Please open an issue or submit a pull request to help improve this utility suite.

---

If you need more detailed explanations or want to see example outputs for any component, feel free to ask!
