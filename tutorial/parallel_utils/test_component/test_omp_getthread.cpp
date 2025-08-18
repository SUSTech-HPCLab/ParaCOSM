#include <iostream>
#include <omp.h>
#include <vector>
#include <cstdlib>

int main() {
    const int N = 100000000;
    std::vector<int> data(N, 1);  // 初始化一个大数组
    double start_time, end_time;
    
    // 设置使用5个线程
    omp_set_num_threads(5);  
    
    // 测量omp_get_thread_num() 和 omp_get_num_threads() 的时间
    start_time = omp_get_wtime();  // 记录起始时间

    #pragma omp parallel
    {
        // 测量获取线程信息的时间
        double thread_start_time = omp_get_wtime();
        
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        
        double thread_end_time = omp_get_wtime();  // 记录结束时间

        // 计算调用omp_get_thread_num() 和 omp_get_num_threads() 的耗时
        double thread_time = thread_end_time - thread_start_time;
        
        if (thread_id == 0) {
            std::cout << "Time taken for omp_get_thread_num() and omp_get_num_threads() on thread " << thread_id 
                      << ": " << thread_time << " seconds." << std::endl;
        }
    }

    end_time = omp_get_wtime();  // 记录结束时间
    std::cout << "Total time for parallel region: " << end_time - start_time << " seconds." << std::endl;

    return 0;
    // Time taken for omp_get_thread_num() and omp_get_num_threads() on thread 0: 8.89413e-08 seconds.
}
