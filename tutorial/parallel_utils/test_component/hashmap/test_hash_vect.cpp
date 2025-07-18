#include <iostream>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <omp.h>
#include <chrono>

// 定义哈希表和互斥锁
std::unordered_map<int, std::vector<int>> hashmap;
std::mutex hashmap_mutex;

void parallel_write(int num_threads, int num_elements, bool use_mutex) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_elements; ++i) {
        std::vector<int> vec(i, i); // 创建一个长度为i的数组，所有元素都为i
        if (use_mutex) {
            std::lock_guard<std::mutex> lock(hashmap_mutex);
            hashmap[i%3] = vec;
        } else {
            hashmap[i%3] = vec;
        }
    }
}

void parallel_read(int num_threads, int num_elements, bool use_mutex) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_elements; ++i) {
        if (use_mutex) {
            std::lock_guard<std::mutex> lock(hashmap_mutex);
            if (hashmap.find(i) != hashmap.end()) {
                std::cout << "Key: " << i << ", Values: ";
                for (const auto& val : hashmap[i]) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        } else {
            if (hashmap.find(i) != hashmap.end()) {
                std::cout << "Key: " << i << ", Values: ";
                for (const auto& val : hashmap[i]) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        }
    }
}

int main() {
    int num_threads = 4;
    int num_elements = 20;

    // 实验1：有互斥锁，无数据竞争，读操作
    hashmap.clear();
    auto start = std::chrono::high_resolution_clock::now();
    parallel_write(num_threads, num_elements, false);
    parallel_read(num_threads, num_elements, true);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "实验1：有互斥锁，无数据竞争，读操作，运行时间：" << duration.count() << "秒" << std::endl;

    
    return 0;
}