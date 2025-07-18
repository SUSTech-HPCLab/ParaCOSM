#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

// 矩阵乘法函数
void matrixMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& result) {
    int m = A.size();    // A 矩阵的行数
    int n = A[0].size(); // A 矩阵的列数 / B 矩阵的行数
    int p = B[0].size(); // B 矩阵的列数

    bool isParallel = true;
    omp_set_num_threads(4);

    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            result[i][j] = 0; // 初始化结果矩阵
            for (int k = 0; k < n; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main() {

    // 设置矩阵大小 N
    int N = 400; // 你可以根据需要修改矩阵大小

    // 初始化随机数生成器
    srand(time(0));

    // 初始化矩阵 A 和 B
    vector<vector<int>> A(N, vector<int>(N));
    vector<vector<int>> B(N, vector<int>(N));

    // 随机生成矩阵 A 和 B 的数值
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = rand() % 100; // 生成 0 到 99 之间的随机数
            B[i][j] = rand() % 100; // 生成 0 到 99 之间的随机数
        }
    }


    int m = A.size();
    int n = B[0].size();
    vector<vector<int>> result_AB(m, vector<int>(n));
    vector<vector<int>> result_BA(m, vector<int>(n));

    std::cout << "Init Complete" << std::endl;
    omp_set_num_threads(4);

    // 同时计算 A * B 和 B * A
    #pragma omp parallel
    {
    #pragma omp  sections // sections 并行
    {
        
        #pragma omp section 
        {
            matrixMultiply(A, B, result_AB);
            std::cout << "Thread Num = " << omp_get_num_threads() << std::endl;
            std::cout << "End" << std::endl;

                // 输出结果
            cout << "Matrix A * B (前5个结果):" << endl;
            for (int i = 0; i < min(5, static_cast<int>(result_AB.size())); ++i) {
                for (int j = 0; j < min(5, static_cast<int>(result_AB[i].size())); ++j) {
                    cout << result_AB[i][j] << " ";
                }
            cout << endl;
            }
        }

        #pragma omp section 
        {
            matrixMultiply(B, A, result_BA);
            std::cout << "Thread Num = " << omp_get_num_threads() << std::endl;
            std::cout << "End" << std::endl;

            cout << "\nMatrix B * A (前5个结果):" << endl;
            for (int i = 0; i < min(5, static_cast<int>(result_BA.size())); ++i) {
                for (int j = 0; j < min(5, static_cast<int>(result_BA[i].size())); ++j) {
                    cout << result_BA[i][j] << " ";
                }
            cout << endl;
            }
        }
    }

    #pragma omp  sections // sections 并行
    {
        // do something
        std::cout << "Thread Num = " << omp_get_num_threads() << std::endl;
    }
    
    }





    return 0;
}
