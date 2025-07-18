#include <iostream>
#include <omp.h>

using namespace std;

void matrix_multiply(const double* A, const double* B, double* C, int N) {
    #pragma omp parallel for collapse(2) // 使用并行化，collapse(2) 合并两个循环
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i * N + j] = 0.0;
            for (int k = 0; k < N; ++k) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

int main() {
    int N = 10;  // 矩阵的大小
    double *A = new double[N * N];
    double *B = new double[N * N];
    double *C = new double[N * N];

    // 初始化矩阵 A 和 B
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i * N + j] = 1.0;  // 假设矩阵 A 的元素都是 1
            B[i * N + j] = 2.0;  // 假设矩阵 B 的元素都是 2
        }
    }

    // 记录开始时间
    double start_time = omp_get_wtime();

    // 调用矩阵乘法函数
    matrix_multiply(A, B, C, N);

    // 记录结束时间
    double end_time = omp_get_wtime();

    // 输出矩阵乘法的耗时
    cout << "Matrix multiplication time: " << end_time - start_time << " seconds" << endl;

    // 清理内存
    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
