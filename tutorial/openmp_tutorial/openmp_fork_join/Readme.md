# OpenMP Quick Start Guide

This guide provides a concise introduction to OpenMP, a widely-used API for parallel programming in C/C++ and Fortran. OpenMP simplifies multi-threaded programming by using compiler directives, library routines, and environment variables to parallelize code. Below, we cover the basics and provide 4-5 practical examples with explanations to help you get started.

## What is OpenMP?
OpenMP (Open Multi-Processing) enables parallel programming on shared-memory systems. It uses pragmas (directives) to specify parallel regions, work-sharing constructs, and synchronization, making it easy to parallelize loops and tasks without manually managing threads.

## Setting Up OpenMP
To use OpenMP in C/C++:
1. Include the OpenMP header: `#include <omp.h>`
2. Compile with an OpenMP-enabled compiler (e.g., `gcc` with `-fopenmp` flag).

## Key Concepts
- **Parallel Regions**: Defined with `#pragma omp parallel`, where multiple threads execute the enclosed code.
- **Work-Sharing Constructs**: Distribute tasks among threads (e.g., `#pragma omp for` for loops).
- **Synchronization**: Ensures thread safety (e.g., `#pragma omp critical` or `omp_set_lock`).
- **Library Functions**: Control threads and query runtime information (e.g., `omp_get_num_threads()`).

## Examples and Explanations

### Example 1: Basic Parallel Region
This example demonstrates a simple parallel region where each thread prints its ID.

```c
#include <stdio.h>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        printf("Hello from thread %d\n", thread_id);
    }
    return 0;
}
```

**Explanation**:
- `#pragma omp parallel` creates a team of threads.
- `omp_get_thread_num()` returns the unique ID of each thread.
- Each thread executes the `printf` statement concurrently.
- Compile with `gcc -fopenmp example1.c -o example1` and run to see output from multiple threads (number depends on system).

### Example 2: Parallelizing a For Loop
This example parallelizes a loop to compute the sum of an array.

```c
#include <stdio.h>
#include <omp.h>

int main() {
    int n = 1000;
    int arr[1000], sum = 0;
    for (int i = 0; i < n; i++) arr[i] = i + 1;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }

    printf("Sum: %d\n", sum);
    return 0;
}
```

**Explanation**:
- `#pragma omp parallel for` divides the loop iterations among threads.
- The `reduction(+:sum)` clause ensures thread-safe accumulation of `sum` by maintaining private copies for each thread and combining them at the end.
- Expected output: `Sum: 500500` (sum of numbers 1 to 1000).

### Example 3: Critical Section for Synchronization
This example shows how to protect a shared resource using a critical section.

```c
#include <stdio.h>
#include <omp.h>

int main() {
    int counter = 0;
    #pragma omp parallel num_threads(4)
    {
        #pragma omp critical
        {
            counter++;
            printf("Thread %d incremented counter to %d\n", omp_get_thread_num(), counter);
        }
    }
    printf("Final counter: %d\n", counter);
    return 0;
}
```

**Explanation**:
- `#pragma omp critical` ensures only one thread at a time can execute the enclosed block, preventing race conditions on `counter`.
- `num_threads(4)` specifies 4 threads for the parallel region.
- The output shows each thread incrementing `counter` sequentially, with the final value being 4.

### Example 4: Dynamic Scheduling for Load Balancing
This example parallelizes a loop with uneven workload using dynamic scheduling.

```c
#include <stdio.h>
#include <omp.h>
#include <math.h>

int main() {
    int n = 100;
    double result[100];

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        result[i] = sin(i * 0.1) * cos(i * 0.1); // Uneven computation
        printf("Thread %d computed result[%d]\n", omp_get_thread_num(), i);
    }

    printf("First result: %f\n", result[0]);
    return 0;
}
```

**Explanation**:
- `#pragma omp parallel for schedule(dynamic)` dynamically assigns loop iterations to threads, ideal for tasks with varying computational costs.
- Each thread computes a portion of `result` based on available iterations.
- The `sin` and `cos` functions simulate uneven workloads.
- Output shows which thread computed each index, with dynamic scheduling balancing the load.

### Example 5: Using Locks for Synchronization
This example uses OpenMP locks to manage access to a shared resource.

```c
#include <stdio.h>
#include <omp.h>

int main() {
    omp_lock_t lock;
    omp_init_lock(&lock);
    int shared_value = 0;

    #pragma omp parallel num_threads(4)
    {
        omp_set_lock(&lock);
        shared_value++;
        printf("Thread %d updated shared_value to %d\n", omp_get_thread_num(), shared_value);
        omp_unset_lock(&lock);
    }

    omp_destroy_lock(&lock);
    printf("Final shared_value: %d\n", shared_value);
    return 0;
}
```

**Explanation**:
- `omp_init_lock` initializes a lock, and `omp_destroy_lock` releases it.
- `omp_set_lock` and `omp_unset_lock` control access to `shared_value`, ensuring only one thread modifies it at a time.
- Each thread increments `shared_value`, and the output shows sequential updates, with the final value being 4.

## Next Steps
- **Experiment**: Modify the examples (e.g., change `num_threads`, scheduling types, or array sizes).
- **Explore**: Learn advanced OpenMP features like task parallelism (`#pragma omp task`) or nested parallelism (`omp_set_nested`).
- **Resources**: Refer to the [OpenMP official documentation](https://www.openmp.org/) for detailed specifications.

This guide provides a foundation for using OpenMP to parallelize your code efficiently. Practice with these examples to build confidence in parallel programming!