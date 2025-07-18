#include <iostream>
#include <chrono>
#include <omp.h>
#include <tbb/tbb.h>

void parallelProgram(int num) {
    tbb::parallel_invoke(
        [num]() {
            #pragma omp parallel for
            for (int i = 0; i < num / 2; i++) {
                printf("i=%d the current thread id: %d\n", i, omp_get_thread_num());
            }
        },
        [num]() {
            #pragma omp parallel for
            for (int i = num / 2; i < num; i++) {
                printf("i=%d the current thread id: %d\n", i, omp_get_thread_num());
            }
        }
    );
}

int main() {
    int num = omp_get_num_procs();
    auto start_time = std::chrono::steady_clock::now();
    parallelProgram(num);
    auto end_time = std::chrono::steady_clock::now();
    std::cout << "parallelProgram elapse time: " << std::chrono::duration<double>(end_time - start_time).count() << " seconds" << std::endl;

    return 0;
}