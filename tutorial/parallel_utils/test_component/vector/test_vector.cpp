#include <iostream>
#include <vector>
#include <mutex>
#include <omp.h>

// 定义向量的向量和互斥锁
std::vector<std::vector<int>> vec_of_vecs;
std::mutex vec_mutex;

void parallel_write(int num_threads, int num_outer_elements, int num_inner_elements) {
    // #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_outer_elements; ++i) {
        std::vector<int> inner_vec(num_inner_elements, i); // 创建一个长度为num_inner_elements的数组，所有元素都为i
        // std::lock_guard<std::mutex> lock(vec_mutex);
        vec_of_vecs.push_back(inner_vec);
    }
}

void parallel_read(int num_threads, int num_outer_elements) {

    // 并行的读是可以的，但是并行的写就不是很行
    #pragma omp parallel for num_threads(num_threads) 
    for (int i = 0; i < num_outer_elements; ++i) {
        // std::lock_guard<std::mutex> lock(vec_mutex);
        if (i < vec_of_vecs.size()) {
            std::cout << "Outer Element " << i << ": ";
            for (const auto& val : vec_of_vecs[i]) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main() {
    int num_threads = 4;
    int num_outer_elements = 10;
    int num_inner_elements = 5;

    // 并行写入向量的向量
    parallel_write(num_threads, num_outer_elements, num_inner_elements);

    // 并行读取向量的向量
    parallel_read(num_threads, num_outer_elements);

    return 0;
}