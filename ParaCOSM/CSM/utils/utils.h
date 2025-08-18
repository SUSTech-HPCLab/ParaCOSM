#ifndef UTILS_UTILS
#define UTILS_UTILS

#include <algorithm>
#include <future>
#include <iostream>
#include <vector>
#include <string.h>
#include <sys/stat.h> /* For stat() */
#include "utils/types.h"

namespace io {
    inline size_t file_exists(const char *path) {
        struct stat st;
        return stat(path, &st) == 0;
    }
}

inline void execute_with_time_limit(std::function<void()> fun, uint time_limit, std::atomic<bool>& reach_time_limit)
{
    std::future future = std::async(std::launch::async, fun);
    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred)
        {
            std::cout << "Deferred" << std::endl;
            exit(-1);
        }
        else if (status == std::future_status::timeout)
        {
            reach_time_limit = true;
            std::cout << "Timeout " << time_limit << "s\n";
        }
    } while (status != std::future_status::ready);
}

namespace mem {
    /**
     * get peak virtual memory space of the current process
     * https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process#answer-64166
     */
    inline int parseLine(char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    }

    inline int getValue(){ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmPeak:", 7) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
}

namespace randomUtil
{
    inline int getRandomIntBetweenAandB(int a, int b){
        return (rand()%(b-a+1))+ a;
    }

    inline bool randomProbability(double probability){
        int randInt = rand();
        double randRatio = (randInt/(RAND_MAX+0.0));
        return randRatio < probability;
    }
}

#define BIG_DATASET 3000000

        // void sum_local_num_results(const std::vector<size_t>& local_num_results, size_t& num_results) {
        //     size_t size = local_num_results.size();
        //     size_t i = 0;
    
        //     // 使用AVX2的256位寄存器进行并行化
        //     __m256i sum_vec = _mm256_setzero_si256();  // 初始化为0的256位向量
    
        //     // 对齐后的数据处理
        //     for (; i + 4 <= size; i += 4) {
        //         // 加载4个size_t数据
        //         __m256i vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&local_num_results[i]));
                
        //         // 对这些4个size_t进行加法
        //         sum_vec = _mm256_add_epi64(sum_vec, vec);
        //     }
    
        //     // 将sum_vec中的结果存入临时数组
        //     alignas(32) size_t temp[4];
        //     _mm256_store_si256(reinterpret_cast<__m256i*>(temp), sum_vec);
    
        //     // 将sum_vec中的四个元素相加，得到最终的结果
        //     num_results += temp[0] + temp[1] + temp[2] + temp[3];
    
        //     // 处理剩余的元素（如果数据的大小不是4的倍数）
        //     for (; i < size; ++i) {
        //         num_results += local_num_results[i];
        //     }
        // }

#endif //UTILS_UTILS
