#include <omp.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <vector>

void recursive_function(int n, int &sum) {
    if (n <= 0) return;

    //  sleep:
     std::this_thread::sleep_for(std::chrono::milliseconds(400));
    auto K = (omp_get_thread_num());
    std::cout << K*n <<std::endl;

    sum += n;
    
    // 在递归调用时创建任务
    // #pragma omp task
    recursive_function(n - 1, sum);

    
}

int main() {
    int N = 8;

    omp_set_num_threads(4);

    int global_sum = 0;
    // vector sum for opoenmp
    std::vector<int> sum(4, 0);

    #pragma omp parallel
    {
        #pragma omp for
            for (int i = 0; i < 4; ++i) {
                for(int j = 0; j< 2; ++j) {
                auto su1 = sum[omp_get_thread_num()];
                recursive_function(3, sum[omp_get_thread_num()]);  // 每次循环调用递归函数
                }
            }
        
    }

    for (int i = 0; i < 4; i++) {
        std::cout << "sum[" << i << "] = " << sum[i] << std::endl;
    }

    return 0;
}
