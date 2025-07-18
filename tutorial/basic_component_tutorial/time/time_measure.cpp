#include <iostream>
#include <chrono>

#define Print_Time_Nano(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() << " ns" << std::endl
#define Print_Time_Micro(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " us" << std::endl

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // 模拟一些工作
    for (volatile int i = 0; i < 10; ++i);

    Print_Time_Nano("My_Duration: ", start);
    Print_Time_Micro("My_Duration: ", start);

    return 0;
}