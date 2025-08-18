#include <iostream>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <omp.h>

// 定义哈希表和互斥锁
std::unordered_map<int, std::vector<int>> hashmap;
std::mutex hashmap_mutex;

void parallel_write(int num_threads, int num_elements) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_elements; ++i) {
        std::lock_guard<std::mutex> lock(hashmap_mutex);
            hashmap[i].push_back(i);
        
        
    }
}

void parallel_read(int num_threads, int num_elements) {
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < num_elements; ++i) {
        // std::lock_guard<std::mutex> lock(hashmap_mutex);
        if (hashmap.find(i) != hashmap.end()) {
            std::cout << "Key: " << i << ", Value: " << hashmap[i%4][0] << std::endl;
        }
        // keep the data for 1s after reading
        auto a = hashmap[i%4][0];
        
        
        // std::this_thread::sleep_for(std::chrono::seconds(1));

    }
}

int main() {
    int num_threads = 8;
    int num_elements = 20;

    // 并行写入哈希表
    parallel_write(num_threads, num_elements);

    // 并行读取哈希表
    parallel_read(num_threads, num_elements);

    return 0;
}