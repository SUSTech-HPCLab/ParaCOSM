#include <omp.h>
#include <iostream>

void perform_task(int id) {
    std::cout << "Task " << id << " is being executed by thread " << omp_get_thread_num() << std::endl;
}

int main() {
    omp_set_num_threads(4); // 设置并行线程数为4

    #pragma omp parallel
    {
        #pragma omp single
        {
            // 创建多个任务
            #pragma omp task
            for (int i = 0; i < 5; ++i) {
                // #pragma omp task
                {
                    perform_task(i);
                }
            }
            #pragma omp task
            for (int i = 0; i < 5; ++i) {
                #pragma omp task
                {
                    perform_task(i);
                }
            }
        }
    }
    return 0;
}
