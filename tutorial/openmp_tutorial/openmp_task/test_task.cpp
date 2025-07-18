#include <iostream>
#include <omp.h>

using namespace std;

// 将并行化部分放在 main 函数中，避免在 print_hi 中创建新的并行区域
inline void print_hi() {
    // 使用 OpenMP 任务并行化
    for (int i = 0; i < 5; ++i) {
        #pragma omp task
        {
            cout << "hi!" << endl;
        }
    }
}

int main() {
    // 使用 OpenMP 并行区域
    #pragma omp parallel
    {
        // 确保只有一个线程执行下面的块
        #pragma omp single
        {
            // 调用包装好的函数
            print_hi();
        }
    }

    return 0;
}
