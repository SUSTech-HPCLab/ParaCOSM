#include <stdio.h>
#include <omp.h>

int main() {
#pragma omp parallel sections
{
#pragma omp section
{
printf("Section 1, Thread ID: %d\n", omp_get_thread_num());
}
#pragma omp section
{
printf("Section 2, Thread ID: %d\n", omp_get_thread_num());
}
#pragma omp section
{
printf("Section 3, Thread ID: %d\n", omp_get_thread_num());
}
#pragma omp section
{
printf("Section 4, Thread ID: %d\n", omp_get_thread_num());
}
}
return 0;
}