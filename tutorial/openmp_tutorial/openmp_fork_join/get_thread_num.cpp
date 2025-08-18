#include <omp.h>
#include <stdio.h>

int main() {
// 在并行区域外部调用，返回主线程的编号，即0
printf("Outside Parallel Region: Thread Num = %d\n", omp_get_thread_num());

// 设置并行区域内的线程数为4
omp_set_num_threads(4);

// 定义并行区域
#pragma omp parallel
{
// 在并行区域内部调用，返回当前线程的编号
    printf("Inside Parallel Region: Thread Num = %d\n", omp_get_thread_num());
}

return 0;
}