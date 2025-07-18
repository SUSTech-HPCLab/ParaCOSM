# CSM-Benchmark Tutorial

Welcome to the CSM-Benchmark (Continuous Subgraph Matching Benchmark) tutorial! This tutorial directory contains example code and guides to help you quickly get started and master the key technologies of the project.

## 📁 Directory Structure

```
tutorial/
├── README.md                    # This file - Tutorial overview
├── basic_component_tutorial/    # Basic component tutorials
├── openmp_tutorial/            # OpenMP parallel programming 
├── parallel_utils/             # Parallelization tools and 
└── plot_figure_tutorial/       # Data visualization tutorial
```

## 🚀 Quick Start

### Requirements

- **Compiler**: GCC 7.0+ (C++17 support) / intel icpx
- **Parallel Library**: OpenMP 4.0+
- **Optional Dependencies**: 
  - Eigen3 (sparse matrix operations)
  - Python 3.7+ (data analysis and visualization)
  - Jupyter Notebook (interactive analysis)

### Compilation Examples

Most C++ examples can be compiled using the following commands:

```bash
# Basic compilation
g++ -std=c++17 example.cpp -o example

# Enable OpenMP
g++ -std=c++17 -fopenmp example.cpp -o example

# Optimized compilation
g++ -std=c++17 -fopenmp -O3 example.cpp -o example
```

## 📚 Tutorial Contents

### 1. [Basic Component Tutorial](basic_component_tutorial/)

Learn the core C++ components and technologies used in the project:

- **Time Measurement** (`time_measure.cpp`): High-precision performance timing
- **Sparse Matrix** (`eigen_sparse_matrix.cpp`): Eigen library sparse matrix operations
- **Vector Operations** (`openmp_vector.cpp`, `pthread_vector.cpp`): Parallel vector computation
- **Multi-threading** (`pthread.cpp`): Basic pthread usage
- **File Reading Performance** (`test_read_speed.cpp`): I/O performance testing

**Recommended Learning Order**:
1. Start with `time_measure.cpp` to learn performance measurement
2. Learn `openmp_vector.cpp` to understand OpenMP basics
3. Compare with `pthread_vector.cpp` to understand different parallel approaches

### 2. [OpenMP Parallel Programming Tutorial](openmp_tutorial/)

Deep dive into OpenMP parallel programming techniques:

- **Basic Concepts**: Hello World examples and API configuration
- **Fork-Join Model** (`openmp_fork_join/`): Parallel region creation and synchronization
- **Task Parallelism** (`openmp_task/`): Dynamic task scheduling
- **Thread Pool** (`ThreadPoolDFS/`): Thread pool implementation for depth-first search
- **TBB Integration** (`tbb/`): Intel Threading Building Blocks

**Key Topics**:
- Thread number configuration and management
- Load balancing
- Synchronization and race condition avoidance
- Performance tuning techniques

### 3. [Parallelization Utils](parallel_utils/)

Practical parallelization analysis and testing tools:

- **Hash Table Performance** (`test_hashmap.cpp`, `hashmap_result.md`): Parallel hash table performance analysis
- **Vector Operations** (`test_vector.cpp`): Parallel vector operation performance testing
- **Pipeline** (`test_pipeline.cpp`): Producer-consumer pattern implementation
- **DCS Algorithm** (`test_DCS.cpp`): Dynamic Community Search testing
- **Cache Performance** (`test_cache/`): Cache-friendly data structure

**Performance Analysis**:
- Use `dcs_size.ipynb` to analyze DCS algorithm performance
- Check `exp_result.md` for parallel overhead analysis

### 4. [Data Visualization Tutorial](plot_figure_tutorial/)

Learn how to generate publication-quality experimental figures:

- **Runtime Analysis** (`fig1_runtime_qs8.pdf`): Query size vs runtime
- **Memory Usage** (`fig3_memory_usage.pdf`): Memory consumption visualization
- **Parallel Performance** (`fig4_parallel_runtime.pdf`): Parallel speedup analysis
- **Load Balancing** (`fig5_load_balance.pdf`): Thread load distribution
- **Technical Overhead** (`fig6_tech_loss.pdf`): Technical overhead analysis

**Jupyter Notebook**:
- Use `plot.ipynb` to learn complete data visualization workflow
- Includes data preprocessing, chart generation, and styling

## 🔧 Hands-on Exercises

### Exercise 1: Performance Benchmarking

1. Compile and run `basic_component_tutorial/time_measure.cpp`
2. Modify the code to test different workload performances
3. Compare nanosecond vs microsecond timing precision

### Exercise 2: OpenMP Parallel Optimization

1. Start learning basic concepts from `openmp_tutorial/`
2. Implement a simple parallel search algorithm
3. Test performance scalability with different thread counts

### Exercise 3: Performance Analysis

1. Run performance tests in `parallel_utils/`
2. Use profiling tools to analyze bottlenecks
3. Try optimization and compare results

### Exercise 4: Data Visualization

1. Install Jupyter environment: `pip install jupyter matplotlib pandas`
2. Run `plot_figure_tutorial/plot.ipynb`
3. Modify code to generate your own experimental charts

## 📖 Further Reading

### Recommended Resources

- C++ learning 【慢速学习C++】 https://www.bilibili.com/video/BV1PFf6YGEyW/?share_source=copy_web&vd_source=72eac555730ba7e7a64f9fa1d7f2b2d4
- **OpenMP Official Documentation**: [openmp.org](https://www.openmp.org/)
- **Performance Analysis Tools**: Intel VTune, perf, gprof
- **Related Papers**: Check `docs/` directory for algorithm background

### Frequently Asked Questions

**Q: Getting OpenMP compilation errors?**
A: Make sure to use the `-fopenmp` flag and check compiler version support

**Q: Unstable performance test results?**
A: Run multiple times and take average, close background programs, use CPU binding

**Q: High memory usage?**
A: Check data structure sizes, use sparse representations, consider chunked processing

## 🤝 Contribution Guidelines

If you find issues in the tutorial or have improvement suggestions:

1. Submit an Issue in the project root directory
2. Provide specific error information or suggestions
3. Pull Requests for tutorial improvements are welcome

## 📄 License

This tutorial follows the project's main license and is for academic research use only.

---

**Start your CSM-Benchmark learning journey!** 🎯

Begin with the [Basic Component Tutorial](basic_component_tutorial/) and gradually master the core technologies of continuous subgraph matching.
