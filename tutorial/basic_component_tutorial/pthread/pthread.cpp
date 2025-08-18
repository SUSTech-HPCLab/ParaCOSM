#include <iostream>
#include <chrono>
#include <thread>
#include <vector>

// 假设有以下类的定义
class DataGraph {
public:
    void AddEdge(int v1, int v2, int label) {
        // 执行图的更新逻辑，假设是一些计算和内存访问
    }
};

class Query {
public:
    int NumVertices() const { return 100; } // 假设有100个顶点
    // 其他与查询相关的成员函数
};

void Print_Time_Nano(const std::string& message, const std::chrono::high_resolution_clock::time_point& start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout << message << duration << " ns\n";
}

void MethodA(DataGraph& data_, int v1, int v2, int label) {
    // 方法A: 更新图
    data_.AddEdge(v1, v2, label);
}

void MethodB(Query& query_) {
    // 方法B: 遍历查询边
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++) {
        // 假设这是查询逻辑
        // 这里可以进一步优化为异步或并行计算的逻辑
    }
}

int main() {
    DataGraph data_;
    Query query_;

    // 获取当前时间
    auto start_nano = std::chrono::high_resolution_clock::now();

    // 方法A：更新图
    int v1 = 1, v2 = 2, label = 10;
    
    // 创建两个线程
    std::thread threadA(MethodA, std::ref(data_), v1, v2, label);
    std::thread threadB(MethodB, std::ref(query_));

    // 等待两个线程完成
    threadA.join();
    threadB.join();

    // 打印时间
    Print_Time_Nano("My_Duration on DataGraph Updation ", start_nano);
    
    return 0;
}
