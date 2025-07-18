# OpenMP Tutorial: High-Performance Parallel Programming in C++

This tutorial repository provides a hands-on, in-depth exploration of OpenMP and related parallel programming techniques in C++. It is designed for learners and practitioners who want to master parallelism, thread management, and performance optimization on modern multi-core systems.

---

## 🌟 Overview

OpenMP is a widely-used API for shared-memory parallel programming in C and C++. This tutorial covers:

- The basics of OpenMP: syntax, configuration, and runtime API
- Advanced parallel patterns: fork-join, tasks, and thread pools
- Performance measurement and best practices
- Comparison with Intel TBB (Threading Building Blocks)
- Real-world code examples and practical tips

---

## 📂 Directory Structure

- **Readme.md**: Main documentation and conceptual explanations
- **ThreadPoolDFS/**: Thread pool and parallel DFS (Depth-First Search) implementations
- **openmp_fork_join/**: Classic fork-join parallelism patterns and OpenMP constructs
- **openmp_task/**: Task-based parallelism and OpenMP task API usage
- **tbb/**: Examples using Intel TBB for comparison with OpenMP

---

## 🚀 Getting Started

### Prerequisites

- C++11 or later
- GCC or Clang with OpenMP support
- (Optional) Intel TBB library for the `tbb/` section

### Compilation Example

```bash
# Compile a basic OpenMP example
g++ your_example.cpp -o your_example -fopenmp

# Compile with TBB (if needed)
g++ tbb_test.cpp -o tbb_test -ltbb
```

---

## 🧠 Tutorial Modules

### 1. OpenMP Basics

- **Hello World**: Learn how to write your first OpenMP program and understand parallel regions.
- **Configuration API**: Explore functions like `omp_get_max_threads()`, `omp_set_num_threads()`, and `omp_get_thread_num()` to control and query the OpenMP runtime.
- **Reference Table**: Comprehensive list of OpenMP runtime functions for thread management, scheduling, locks, and timing.

### 2. ThreadPoolDFS

- **ThreadPool.cpp / ThreadPool.h**: Implementation of a simple thread pool for managing worker threads.
- **DFS_omp4.cpp**: Parallel depth-first search using OpenMP and thread pools.
- **Readme.md**: Detailed explanation of thread pool design, DFS parallelization strategies, and performance considerations.

### 3. openmp_fork_join

- **test_parallel_for.cpp**: Demonstrates the use of `#pragma omp parallel for` for data-parallel loops.
- **test_section.cpp**: Shows how to use OpenMP sections for task parallelism.
- **test_single.cpp**: Explains the `single` directive for single-threaded execution within a parallel region.
- **test_recursive.cpp**: Recursive parallelism patterns.
- **get_thread_num.cpp**: How to query thread IDs in OpenMP.
- **test_time.cpp**: Measuring execution time in parallel regions.
- **Readme.md**: In-depth guide to fork-join parallelism and OpenMP constructs.

### 4. openmp_task

- **task.cpp / taskgroup.cpp / test_task.cpp**: Examples of OpenMP task creation, task groups, and dynamic task scheduling.
- **Readme.md**: Overview of task-based parallelism, when to use tasks vs. parallel for, and best practices.

### 5. tbb

- **tbb_test.cpp / test_diff_thread2.cpp**: Simple examples using Intel TBB to illustrate differences and similarities with OpenMP.

---

## 📊 Performance and Best Practices

- Learn how to measure and analyze the performance of parallel code.
- Understand the impact of thread scheduling, load balancing, and synchronization.
- Discover common pitfalls and how to avoid them (e.g., false sharing, oversubscription).

---

## 📚 Further Reading

- [OpenMP Official Documentation](https://www.openmp.org/specifications/)
- [Intel TBB Documentation](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)
- [OpenMP Examples by muatik](https://github.com/muatik/openmp-examples)
- [OpenMP Validation & Verification Suite](https://crpl.cis.udel.edu/ompvv/)

---

## 🤝 Contributing

Contributions, bug reports, and suggestions are welcome! Please open an issue or submit a pull request to help improve this tutorial.

---

## 📝 License

This tutorial is provided for educational purposes. Please check individual files for license information if you plan to reuse code in your own projects.

---

If you need more detailed explanations for any module or want to see example outputs, feel free to ask!
