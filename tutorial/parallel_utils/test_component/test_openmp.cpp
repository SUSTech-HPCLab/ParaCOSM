#include <iostream>
#include <omp.h>
#include <vector>
#include <cstdlib>

int main() {
    // 设置数组大小
    const int N = 100;
    std::vector<int> data(N, 1);  // 初始化一个大数组
    double start_time, end_time;

    // 设置使用5个线程
    omp_set_num_threads(4);  

    // 测试并行区域的启动时间
    start_time = omp_get_wtime();  // 记录起始时间

    #pragma omp parallel
    {
        // 并行执行的任务 (这里我们只是进行一个简单的加法操作)
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        for (int i = thread_id; i < N; i += num_threads) {
            data[i] += 1;  // 简单的操作，模拟一些计算
        }
    }

    end_time = omp_get_wtime();  // 记录结束时间
    std::cout << "Time taken to run parallel region: " << end_time - start_time << " seconds." << std::endl;

    // 现在可以对比串行和并行区域的耗时，串行部分:
    start_time = omp_get_wtime();  // 记录起始时间

    // 串行执行的任务
    for (int i = 0; i < N; ++i) {
        data[i] += 1;
    }

    end_time = omp_get_wtime();  // 记录结束时间
    std::cout << "Time taken for serial execution: " << end_time - start_time << " seconds." << std::endl;

    return 0;

    // 9.31609e-05 seconds. 这是并行所需时间，其纳秒数是 9.31609e-05 * 1e9 = 98,736 ns
    // 因此，需要提前创建好并行域
}
