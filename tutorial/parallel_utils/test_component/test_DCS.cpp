#include <iostream>
#include <queue>
#include <utility> // for std::pair

void InsertionTopDown_para(uint u1, uint u2, uint v1, uint v2, std::queue<std::pair<uint, uint>>& Q1_para, std::queue<std::pair<uint, uint>>& Q2_para) {
    // 示例实现
    Q1_para.emplace(u1, v1);
    Q2_para.emplace(u2, v2);
}

int main() {
    std::queue<std::pair<uint, uint>> Q1_para, Q2_para;

    uint u1 = 1, u2 = 2, v1 = 3, v2 = 4;
    bool old_p_d1 = true, old_c_d2 = true;

    if (old_p_d1)
        InsertionTopDown_para(u1, u2, v1, v2, Q1_para, Q2_para);
    if (old_c_d2)
        InsertionTopDown_para(u2, u1, v2, v1, Q2_para, Q1_para);

    // 输出队列中的元素
    while (!Q1_para.empty()) {
        auto [u, v] = Q1_para.front();
        Q1_para.pop();
        std::cout << "Q1_para: (" << u << ", " << v << ")" << std::endl;
    }

    while (!Q2_para.empty()) {
        auto [u, v] = Q2_para.front();
        Q2_para.pop();
        std::cout << "Q2_para: (" << u << ", " << v << ")" << std::endl;
    }

    return 0;
}
