#include <iostream>
#include <vector>
#include <omp.h>

// 阶段1：生成数据
void generate_data(std::vector<int>& buffer, int num_elements) {
    for (int i = 0; i < num_elements; ++i) {
        buffer.push_back(i);
    }
}

// 阶段2：处理数据
void process_data(std::vector<int>& buffer) {
    for (auto& elem : buffer) {
        elem *= 2; // 简单的处理：将每个元素乘以2
    }
}

// 阶段3：输出数据
void output_data(const std::vector<int>& buffer) {
    for (const auto& elem : buffer) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

int main() {
    int num_elements = 10000;
    std::vector<int> buffer1, buffer2;
    std::vector<int>* current_buffer = &buffer1;
    std::vector<int>* next_buffer = &buffer2;

    omp_set_num_threads(4);

    // 使用OpenMP实现Pipeline并行
    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < 2; ++i) { // 运行两次以演示双缓冲
                // 阶段1：生成数据
                #pragma omp task firstprivate(current_buffer)
                {
                    generate_data(*current_buffer, num_elements);
                    std::cout << "Data generated." << std::endl;
                }

                // 阶段2：处理数据
                #pragma omp task depend(in: current_buffer)
                {
                    process_data(*current_buffer);
                    std::cout << "Data processed." << std::endl;
                }

                // 阶段3：输出数据
                #pragma omp task depend(in: current_buffer)
                {
                    output_data(*current_buffer);
                    std::cout << "Data outputted." << std::endl;
                }

                // 交换缓冲区
                std::swap(current_buffer, next_buffer);
                current_buffer->clear();
            }
        }
    }

    return 0;
}