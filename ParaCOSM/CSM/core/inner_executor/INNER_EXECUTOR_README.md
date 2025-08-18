## InnerExecutor 与内核函数文档

### 概览

`core/inner_executor` 模块提供了基于并行化的子图匹配执行器 `InnerExecutor`，以及若干可直接调用的匹配内核（SymBi / GraphFlow / TurboFlux / ParaCOSM Kernel）。该模块既可作为对外的统一接口（继承自 `matching` 基类），也可按需调用内核级别的函数以实现定制化的并行枚举。

### 主要特性

- **继承式设计**：`InnerExecutor` 直接继承自 `matching` 基类，复用原有执行统计、图接口与配置
- **多内核支持**：支持 SymBi、GraphFlow、TurboFlux 及 ParaCOSM Kernel 的并行版本
- **线程安全与并行化**：OpenMP 进行任务并行，`tbb::concurrent_queue` 支持任务拆分与动态负载均衡
- **可配置**：线程数、是否同态匹配（homomorphism）、自动调参、顶点顺序等均可配置
- **工程友好**：数据结构初始化、资源管理与异常处理开箱即用

### 目录结构

- `inner_executor.h/.cpp`：`InnerExecutor` 的对外接口与实现
- `FindMatchesKernel.h/.cpp`：多种并行匹配内核及其辅助过程
- `FindMatchesKernel_simplified.cpp`：简化/示例化实现，便于阅读与验证

### 类与函数

#### 类：`InnerExecutor`（位于 `core/inner_executor/inner_executor.h`）

- 构造参数：
  - `Graph& query_graph, Graph& data_graph`
  - `uint max_num_results`
  - `bool print_prep, print_enum, homomorphism`
  - `std::vector<std::vector<uint>> orders`
  - `size_t num_threads, auto_tuning`

- 核心方法（覆盖基类或新增）：
  - `void Preprocessing()`：构建索引与候选结构，初始化线程本地数据
  - `void InitialMatching()`：初始化匹配向量与统计
  - `size_t Execute()`：主执行入口（可按需扩展真实求解逻辑）
  - `size_t ExecuteParaCOSMKernel2(uint depth, std::vector<uint>& m, std::vector<uint>& extendable)`：直接调用 ParaCOSM Kernel2 的示例入口
  - 图修改接口（覆盖自基类）：`AddEdge/RemoveEdge/AddVertex/RemoveVertex`
  - 统计/工具：`GetMemoryCost(...)`、`Classify(...)`

- 关键数据结构（节选）：
  - `eidx_`：查询图边索引映射
  - `DCS_ / d1 / d2 / n1 / np1 / n2 / nc2`：候选与约束相关索引
  - `local_vec_m_ / local_vec_visited_local_`：线程本地的匹配/访问状态
  - `vertex_vector_ / job_queue_`：用于任务拆分与并行处理

#### 内核函数（位于 `core/inner_executor/FindMatchesKernel.h/.cpp`）

- SymBi 内核
  - `void Parallel_Symbi_FindMatches(...)`
  - `void Parallel_Symbi_FindMatches_Simplified(...)`
  - `void Parallel_Symbi_FindMatches_ParaCOSM_Kernel(uint depth, std::vector<uint>& m, std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results)`

- GraphFlow 内核
  - `void Parallel_Graphflow_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)`

- TurboFlux 内核
  - `void Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results)`

- 典型辅助过程（节选）
  - `Parallel_FindMatches_local3_simplified(...)`
  - `process_vertex_layer1_simplified(...)`

以上内核流程均遵循相似的“候选检索 → 约束过滤（索引/连通/访问）→ 递归扩展/任务拆分 → 回溯”的结构，并在关键分支处支持任务队列拆分以提升并行吞吐。

### 使用示例

#### 1) 作为统一执行器使用 `InnerExecutor`

```cpp
#include "core/inner_executor/inner_executor.h"
#include "graph_storage/graph.h"

int main() {
    Graph query_graph, data_graph;

    InnerExecutor executor(
        query_graph,
        data_graph,
        /* max_num_results = */ 10000000,
        /* print_prep = */ false,
        /* print_enum = */ false,
        /* homomorphism = */ false,
        /* orders = */ {},
        /* num_threads = */ 8,
        /* auto_tuning = */ 0
    );

    executor.Preprocessing();
    executor.InitialMatching();
    size_t results = executor.Execute();

    return 0;
}
```

#### 2) 直接调用简化内核（SymBi Simplified）

```cpp
#include "core/inner_executor/FindMatchesKernel.h"
#include "graph_storage/graph.h"

void run_kernel(const Graph& query, const Graph& data) {
    std::vector<uint> m(query.NumVertices(), UNMATCHED);
    std::vector<Parrllel_SymBi::ExtendableVertex> extendable(query.NumVertices());

    size_t num_results = 0;
    size_t num_threads = 8;
    bool homomorphism = false;

    // 下列结构需按实际预处理结果进行构造与填充
    std::vector<std::vector<uint>> eidx(query.NumVertices(), std::vector<uint>(query.NumVertices()));
    std::vector<ska::flat_hash_map<uint, std::vector<uint>>> DCS(query.NumEdges() * 2);
    std::vector<ska::flat_hash_map<uint, bool>> d2(query.NumVertices());
    std::vector<ska::flat_hash_map<uint, uint>> n2(query.NumEdges() * 2);
    std::vector<std::vector<bool>> local_visited(num_threads, std::vector<bool>(data.NumVertices(), false));
    std::vector<std::vector<uint>> treeNode_neighbors(query.NumVertices());

    Parallel_Symbi_FindMatches_Simplified(
        /* depth */ 0,
        m,
        extendable,
        num_results,
        query,
        data,
        eidx,
        DCS,
        d2,
        n2,
        local_visited,
        homomorphism,
        num_threads,
        treeNode_neighbors
    );
}
```

### 线程与任务调度

- 使用 OpenMP 进行候选并行与层间并行：`#pragma omp parallel for ...`
- 使用 `tbb::concurrent_queue` 实现任务拆分（深度阈值后将分支入队，尾阶段进行抢占式消费）
- 线程本地状态：`local_vec_m_`、`local_vec_visited_local_`，减少锁争用
- 建议按数据规模设置 `num_threads`，必要时开启任务拆分以提升尾部利用率

### 配置项与建议

- **线程数**：`InnerExecutor(num_threads=...)` 或传入内核的 `NUMTHREAD`
- **同态匹配**：`homomorphism = true/false`（false 时启用访问判重）
- **顶点顺序**：通过 `orders_` 或内核的 `order_index` 控制扩展序
- **预处理**：保证 `eidx / DCS / d2 / n2 / treeNode_neighbors` 与图一致且已排序，以支撑 `lower_bound`/索引过滤

### 依赖

- OpenMP（并行编程）
- Intel TBB：`tbb::concurrent_queue`
- `ska::flat_hash_map`（见 `graph_storage/storage_hash_map.hpp`）
- 图与匹配基类：`graph_storage/graph.h`、`matching_executor/matching.h`

### 构建

- 在顶层或 `CSM/` 的 `CMakeLists.txt` 中启用该模块（确保 OpenMP/TBB 可用）
- 典型选项：`-fopenmp`、链接 TBB；不同平台的具体链接方式请参考项目 CMake 配置

### 实现说明

- 候选检索优先使用已构建的 `DCS` 与 `eidx`，再对邻接做 `lower_bound` 交验证
- 访问判重使用线程本地布尔数组，深度回溯时成对恢复
- 任务拆分仅在设定深度后进行，避免过早产生过多细粒度任务

### 后续工作

- 完善 `Execute()` 主流程的实际求解与调参逻辑
- 丰富/优化候选索引构建（`BuildDCS()`）与边索引（`BuildEdgeIndex()`）
- 增加性能统计与剖析工具，完善自适应线程/切分策略


