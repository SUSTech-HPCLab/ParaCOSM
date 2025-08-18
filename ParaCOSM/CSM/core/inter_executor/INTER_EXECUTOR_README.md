# InterExecutor 重构说明

## 概述

`InterExecutor` 类是对原始 `BatchUpdates` 函数的重构，采用了面向对象的设计模式，提供了更好的代码组织、可维护性和可扩展性。

## 重构内容

### 1. 原始函数分析

原始的 `BatchUpdates` 函数是一个内联函数，具有以下特点：
- 使用滑动窗口机制处理批量图更新
- 区分安全更新和不安全更新
- 维护各种计数器（顶点更新、边更新、不安全更新等）
- 处理顶点和边的添加/删除操作

### 2. 重构后的类结构

#### 核心类：`InterExecutor`

```cpp
class InterExecutor {
private:
    Graph& data_graph_;           // 数据图引用
    matching* matching_instance_;  // 匹配算法实例
    
    static constexpr size_t DEFAULT_WINDOW_SIZE = 16;  // 默认窗口大小
    
public:
    // 构造函数
    InterExecutor(Graph& data_graph, matching* matching_instance);
    
    // 主要方法
    void ProcessBatchUpdates(...);     // 批量处理更新
    bool ProcessSingleUpdate(...);     // 处理单个更新
    void ApplyUnsafeUpdate(...);       // 应用不安全更新
    bool IsUpdateSafe(...);            // 检查更新是否安全
    
    // 向后兼容的别名
    void BatchUpdates(...);           // 原始函数名
};
```

### 3. 主要改进

#### 3.1 代码组织
- **模块化设计**：将复杂的逻辑分解为多个职责明确的方法
- **类封装**：将相关的数据和方法封装在一个类中
- **清晰的接口**：提供明确的公共方法接口

#### 3.2 可维护性
- **单一职责**：每个方法只负责一个特定的功能
- **易于测试**：可以单独测试每个方法
- **易于扩展**：可以轻松添加新的功能或修改现有功能

#### 3.3 可读性
- **有意义的方法名**：方法名清楚地表达了其功能
- **详细的文档**：每个方法都有完整的文档说明
- **逻辑分离**：将复杂的逻辑分解为更小的、可理解的部分

## 使用方法

### 1. 基本使用

```cpp
#include "inter_executor.h"

// 创建 InterExecutor 实例
InterExecutor executor(data_graph, matching_instance);

// 初始化计数器
size_t num_v_updates = 0, num_e_updates = 0, unsafe_updates = 0;
size_t count = 0, positive_num_results_last = 0, negative_num_results_last = 0;
std::atomic_bool reach_time_limit(false);

// 处理批量更新
executor.ProcessBatchUpdates(
    num_v_updates, num_e_updates, unsafe_updates, count,
    positive_num_results_last, negative_num_results_last,
    reach_time_limit
);
```

### 2. 自定义窗口大小

```cpp
// 使用自定义窗口大小
executor.ProcessBatchUpdates(
    num_v_updates, num_e_updates, unsafe_updates, count,
    positive_num_results_last, negative_num_results_last,
    reach_time_limit,
    32  // 自定义窗口大小
);
```

### 3. 处理单个更新

```cpp
InsertUnit update('e', true, 1, 2, 1);  // 添加边 (1,2) 标签为 1
bool is_safe = executor.ProcessSingleUpdate(update);
```

### 4. 向后兼容

```cpp
// 仍然可以使用原始的 BatchUpdates3 函数名
executor.BatchUpdates3(
    num_v_updates, num_e_updates, unsafe_updates, count,
    positive_num_results_last, negative_num_results_last,
    reach_time_limit
);
```

## 文件结构

```
ParaCOSM/CSM/core/
├── inter_executor.h              # 头文件
├── inter_executor.cpp            # 实现文件
├── inter_executor_example.cpp    # 使用示例
└── INTER_EXECUTOR_README.md      # 本文档
```

## 向后兼容性

为了保持向后兼容性，重构后的代码提供了：

1. **相同的函数签名**：`BatchUpdates3` 方法保持与原始函数完全相同的参数
2. **相同的功能**：所有原始功能都得到保留
3. **相同的性能**：重构不会影响性能，实际上可能略有提升

## 扩展建议

### 1. 配置选项
- 可以添加配置文件支持，允许动态调整窗口大小
- 可以添加性能监控和日志记录

### 2. 并行处理
- 可以进一步优化并行处理逻辑
- 可以添加线程池支持

### 3. 错误处理
- 可以添加更完善的错误处理机制
- 可以添加异常处理

### 4. 性能优化
- 可以添加缓存机制
- 可以优化内存使用

## 总结

这次重构将原始的 `BatchUpdates3` 函数转换为一个结构良好的 `InterExecutor` 类，提供了：

- 更好的代码组织和可维护性
- 清晰的接口和文档
- 向后兼容性
- 易于扩展和测试的架构

重构后的代码保持了所有原始功能，同时为未来的改进奠定了良好的基础。
