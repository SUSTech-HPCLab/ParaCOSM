
# Basic Component Tutorial

This tutorial collection provides a set of C++ code examples covering essential components for high-performance and parallel programming. The examples include sparse matrix operations, parallel and multithreaded vector processing, memory read speed benchmarking, and precise time measurement. These are designed to help beginners quickly understand and experiment with fundamental techniques in modern C++.

---

## Table of Contents

- [eigen_sparse_matrix.cpp](#eigen_sparse_matrixcpp)
- [openmp_vector.cpp](#openmp_vectorcpp)
- [pthread_vector.cpp](#pthread_vectorcpp)
- [pthread.cpp](#pthreadcpp)
- [test_read_speed.cpp](#test_read_speedcpp)
- [time_measure.cpp](#time_measurecpp)

---

## 🧩 Component Overview

### `eigen_sparse_matrix.cpp`

- **Purpose:** Demonstrates how to use the Eigen library to construct a sparse matrix and perform BFS (Breadth-First Search) traversal on a graph represented by a sparse adjacency matrix.
- **Use Case:** Useful for handling large-scale sparse graphs, such as those found in social networks or recommendation systems.

### `openmp_vector.cpp`

- **Purpose:** Shows how to use OpenMP to parallelize the initialization and processing of a `std::vector`.
- **Use Case:** Suitable for simple data-parallel tasks, such as batch data initialization or numerical computations.

### `pthread_vector.cpp`

- **Purpose:** Utilizes C++11 standard threads (`std::thread`) to fill a large vector in parallel, with manual control over thread workload.
- **Use Case:** Ideal for scenarios requiring explicit thread management, such as chunked data processing.

### `pthread.cpp`

- **Purpose:** Demonstrates concurrent execution of graph updates and query operations using multiple threads, with precise timing of operations.
- **Use Case:** Helps understand data structure manipulation and performance analysis in a multithreaded environment.

### `test_read_speed.cpp`

- **Purpose:** Benchmarks memory read speed for different data sizes, simulating the effects of cache and main memory access.
- **Use Case:** Useful for analyzing how hardware cache hierarchies impact program performance.

### `time_measure.cpp`

- **Purpose:** Provides macros for measuring code execution time at nanosecond and microsecond precision.
- **Use Case:** Essential for profiling and optimizing code segments that require accurate timing.

---

## 🚀 Getting Started

1. **Dependencies**
   - C++11 or later
   - [Eigen](https://eigen.tuxfamily.org/) (required for `eigen_sparse_matrix.cpp`)
   - OpenMP support (for `openmp_vector.cpp`)

2. **Compilation Examples**
   ```bash
   # Compile eigen_sparse_matrix.cpp
   g++ eigen_sparse_matrix.cpp -o eigen_sparse_matrix -I /path/to/eigen

   # Compile openmp_vector.cpp
   g++ openmp_vector.cpp -o openmp_vector -fopenmp

   # Compile pthread_vector.cpp
   g++ pthread_vector.cpp -o pthread_vector -std=c++11

   # Compile pthread.cpp
   g++ pthread.cpp -o pthread -std=c++11

   # Compile test_read_speed.cpp
   g++ test_read_speed.cpp -o test_read_speed

   # Compile time_measure.cpp
   g++ time_measure.cpp -o time_measure
   ```

3. **Run the Examples**
   ```bash
   ./eigen_sparse_matrix
   ./openmp_vector
   ./pthread_vector
   ./pthread
   ./test_read_speed
   ./time_measure
   ```

---

## 📌 Contributions & Feedback

Contributions and suggestions for more basic component examples are welcome! Help us make high-performance and parallel programming more accessible to everyone.

---

If you need more detailed code explanations or further documentation, feel free to ask!
