#include <iostream>
#include <vector>
#include <omp.h>

// 阶段1：生成数据
void generate_data(std::vector<int>& data, int num_elements) {
    for (int i = 0; i < num_elements; ++i) {
        data.push_back(i);
    }
}

// 阶段2：处理数据
void process_data(std::vector<int>& data) {
    for (auto& elem : data) {
        elem *= 2; // 简单的处理：将每个元素乘以2
    }
}

// 阶段3：输出数据
void output_data(const std::vector<int>& data) {
    for (const auto& elem : data) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

int main() {
    int num_elements = 4000;
    std::vector<int> data;

    omp_set_num_threads(4);

    // 使用OpenMP实现Pipeline并行
    #pragma omp parallel
    {
        #pragma omp single
        {
            // 阶段1：生成数据
            #pragma omp task
            {
                generate_data(data, num_elements);
                std::cout << "Data generated." << std::endl;
            }

            // 阶段2：处理数据
            #pragma omp task depend(in: data)
            {
                process_data(data);
                std::cout << "Data processed." << std::endl;
            }

            // 阶段3：输出数据
            #pragma omp task depend(in: data)
            {
                output_data(data);
                std::cout << "Data outputted." << std::endl;
            }
        }
    }

    return 0;
}