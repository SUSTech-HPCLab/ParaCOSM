#include <iostream>
#include <vector>
#include <thread>

void fill_vector(std::vector<int>& vec, int start_idx, int end_idx) {
    for (int i = start_idx; i < end_idx; ++i) {
        vec[i] = i;
    }
}

int main() {
    const int size = 10000000; // Vector的大小
    std::vector<int> vec(size, -1); // 初始化vector并填充为-1

    const int num_threads = 4; // 使用4个线程
    int chunk_size = size / num_threads; // 每个线程处理的元素数量

    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        int start_idx = i * chunk_size;
        int end_idx = (i == num_threads - 1) ? size : (i + 1) * chunk_size; // 最后一个线程处理余下的元素
        threads.push_back(std::thread(fill_vector, std::ref(vec), start_idx, end_idx));
    }

    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }

    // 输出结果
    // for (int i = 0; i < size; ++i) {
    //     std::cout << vec[i] << " ";
    // }
    std::cout << "Finish" << std::endl;

    return 0;
}
