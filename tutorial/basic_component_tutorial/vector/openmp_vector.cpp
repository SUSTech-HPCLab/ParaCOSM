#include <iostream>
#include <vector>
#include <omp.h>

int main() {
    const int size = 100; // Vector的大小
    std::vector<int> vec(size, -1); // 初始化vector并填充为-1

    // 启用OpenMP进行并行填充
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        vec[i] = i;
    }

    // 输出结果
    for (int i = 0; i < size; ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
