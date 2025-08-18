#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>

const size_t KB = 1024;
const size_t MB = 1024 * KB;
const size_t GB = 1024 * MB;

void measure_read_speed(size_t size) {
    std::vector<char> data(size);
    for (size_t i = 0; i < size; ++i) {
        data[i] = rand() % 256;
    }

    volatile char sink;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < size; ++i) {
        sink = data[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double speed = size / duration.count() / (1024 * 1024); // MB/s
    std::cout << "Read speed for " << size / (1024 * 1024) << " MB: " << speed << " MB/s" << std::endl;
}

int main() {
    std::cout << "Measuring read speed for different memory sizes..." << std::endl;

    measure_read_speed(256 * KB); // L1 Cache size
    measure_read_speed(1 * MB);   // L2 Cache size
    measure_read_speed(8 * MB);   // L3 Cache size
    measure_read_speed(256 * MB); // Main memory size

    return 0;
    // 0.36175872 B / ns
}