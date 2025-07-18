#include <iostream>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <omp.h>

// 定义向量的unordered_map和互斥锁
std::vector<std::unordered_map<int, int>> vec_of_maps;
std::mutex vec_mutex;

void parallel_write(int num_threads, int num_outer_elements, int num_inner_elements) {

    // 并行的写不可以，但是并行读可以

    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_outer_elements; ++i) {
        std::unordered_map<int, int> inner_map;
        for (int j = 0; j < num_inner_elements; ++j) {
            inner_map[j] = i * j; // 填充unordered_map
        }
        std::lock_guard<std::mutex> lock(vec_mutex);
        vec_of_maps.push_back(inner_map);
    }
}

void parallel_read(int num_threads, int num_outer_elements) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_outer_elements; ++i) {
        // std::lock_guard<std::mutex> lock(vec_mutex);
        if (i < vec_of_maps.size()) {
            std::cout << "Outer Element " << i << ": ";
            for (const auto& pair : vec_of_maps[i]) {
                std::cout << "{" << pair.first << ": " << pair.second << "} ";
            }
            std::cout << std::endl;
        }
    }
}

int main() {
    int num_threads = 4;
    int num_outer_elements = 10;
    int num_inner_elements = 5;

    // 并行写入向量的unordered_map
    parallel_write(num_threads, num_outer_elements, num_inner_elements);

    // 并行读取向量的unordered_map
    parallel_read(num_threads, num_outer_elements);

    return 0;
}