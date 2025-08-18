#include <tbb/tbb.h>
#include <iostream>

void print_hello_world(int i) {
    std::cout << "Hello, World! from thread " << i << std::endl;
    // sleep
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
}

int main() {
    // 创建一个TBB并行循环
    tbb::parallel_for(0, 5, 1, [](int i) {
        print_hello_world(i);
    });

    std::cout << "All thread complete" << std::endl;
    
    return 0;
}
