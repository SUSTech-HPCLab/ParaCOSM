#include <omp.h>
#include <iostream>

void perform_task(int i) {
    // 模拟任务处理
    std::cout << "Task " << i << " is being processed by thread " << omp_get_thread_num() << std::endl;
}

int main() {
    omp_set_num_threads(4); // 设置并行线程数为4

    #pragma omp parallel
    {
        #pragma omp single
        {
            // 创建多个任务
            for (int i = 0; i < 15; ++i) {
                #pragma omp task
                {
                    // 使用 taskgroup 来确保这些任务由3个线程处理
                    #pragma omp taskgroup
                    {
                        #pragma omp parallel num_threads(3)
                        {
                            #pragma omp single nowait
                            {
                                perform_task(i);
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < 5; ++i) {
                #pragma omp task
                {
                    // 使用 taskgroup 来确保这些任务由1个线程处理
                    #pragma omp taskgroup
                    {
                        #pragma omp parallel num_threads(1)
                        {
                            #pragma omp single nowait
                            {
                                perform_task(i);
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}