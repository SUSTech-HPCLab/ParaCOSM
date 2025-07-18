#include <omp.h>
#include <stdio.h>

void parallel_function() {
    // omp_set_num_threads(4);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        
        printf("并行区域：线程 %d 执行\n", id);
    }
}

int main() {
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        printf("并行区域外：线程 %d 执行\n", id);
        #pragma omp single
        {
            printf("单线程执行部分\n");
            parallel_function();  // 调用并行函数
        }

        parallel_function();
    }

    return 0;
}
