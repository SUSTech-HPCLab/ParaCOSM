# GPU Kernel 优化技术参考 (来自 GraphMiner/G²Miner)

> 本文档总结了 GraphMiner (OSDI 2022) 中 GPU 图模式挖掘 kernel 的核心优化技术，供 ParaCOSM GPU kernel 优化参考。

---

## 一、并行粒度策略

### 1.1 Warp-centric Edge-parallel（推荐默认策略）

```cuda
// 每个 warp (32线程) 协作处理 edgelist 中的一条边
for (eidType eid = warp_id; eid < ne; eid += num_warps) {
    auto v = g.get_src(eid);
    auto u = g.get_dst(eid);
    count += intersect_num(g.N(v), v_size, g.N(u), u_size);
}
```

**优势**：
- 避免 warp divergence：同一 warp 内 32 线程执行相同控制流
- 负载均衡：边数远大于顶点数，任务粒度更均匀
- 内存合并：warp 内线程顺序访问邻居列表

**适用场景**：中等度数图的大部分 pattern matching kernel。

### 1.2 CTA-centric Edge-parallel（高度数场景）

```cuda
// 整个 thread block (256线程) 协作处理一条边
for (eidType eid = blockIdx.x; eid < ne; eid += gridDim.x) {
    if (threadIdx.x == 0) { v = g.get_src(eid); u = g.get_dst(eid); }
    __syncthreads();
    count += g.cta_intersect_cache(v, u);  // 全 block 协作交集
}
```

**适用场景**：超高度数节点（邻居列表 >1000），一个 warp 不足以并行处理。

### 1.3 Vertex-parallel（负载不均，慎用）

```cuda
for (auto v = warp_id; v < nv; v += num_warps) {
    for (auto e = 0; e < v_size; e++) {
        count += intersect_num(v_ptr, v_size, g.N(u), u_size);
    }
}
```

**问题**：幂律图中高度数顶点成瓶颈，导致严重负载不均。仅在顶点度数分布均匀时考虑。

---

## 二、集合交集优化（最关键的原语）

### 2.1 2-Phase Binary Search with Shared Memory Caching

这是 GraphMiner 中**性能最好**的集合交集实现。

```cuda
// Phase 1: 将 search 列表均匀采样 32 个 pivot 到 shared memory
cache[warp_lane * WARP_SIZE + thread_lane] = search[thread_lane * search_size / WARP_SIZE];
__syncwarp();

// Phase 2: 每个线程取一个 key，先在 cache 中粗定位，再在 global memory 小区间精确查找
for (auto i = thread_lane; i < lookup_size; i += WARP_SIZE) {
    auto key = lookup[i];
    if (binary_search_2phase(search, cache, key, search_size))
        num += 1;
}
```

**`binary_search_2phase` 实现**：

```cuda
__device__ bool binary_search_2phase(T *list, T *cache, T key, T size) {
    int p = (threadIdx.x / WARP_SIZE) * WARP_SIZE;
    // Phase 1: shared memory 中 O(log 32) = 5 次比较定位区间
    int bottom = 0, top = WARP_SIZE;
    while (top > bottom + 1) {
        int mid = (top + bottom) / 2;
        T y = cache[p + mid];           // ← shared memory，~1 cycle
        if (key == y) return true;
        if (key < y) top = mid;
        if (key > y) bottom = mid;
    }
    // Phase 2: global memory 中 O(log(N/32)) 次比较
    bottom = bottom * size / WARP_SIZE;
    top = top * size / WARP_SIZE - 1;
    while (top >= bottom) {
        int mid = (top + bottom) / 2;
        T y = list[mid];                // ← global memory，~100+ cycles
        if (key == y) return true;
        if (key < y) top = mid - 1;
        else bottom = mid + 1;
    }
    return false;
}
```

**性能分析**：
- 普通 binary search：`O(log N)` 次 global memory 访问
- 2-phase：`O(log 32) + O(log(N/32))` = `5 + O(log N - 5)` 次，但前 5 次在 shared memory
- 实测：对 degree=1000 的列表，global memory 访问从 ~10 次降到 ~5 次（减少 50%）

### 2.2 小表遍历 + 大表搜索

```cuda
T* lookup = a;  // 遍历的表
T* search = b;  // 搜索的表
if (size_a > size_b) { lookup = b; search = a; ... }
```

**总复杂度**：$O(\min(|A|,|B|) \cdot \log(\max(|A|,|B|)))$，比 merge 在非对称大小时更优。

### 2.3 Warp 协作并行查找

```cuda
// warp 内 32 个线程同时取 32 个不同的 key 并行查找
for (auto i = thread_lane; i < lookup_size; i += WARP_SIZE) {
    auto key = lookup[i];  // 每个线程取自己的 key
    if (binary_search_2phase(search, cache, key, search_size))
        num += 1;
}
// 隐式：32 路查找并行，吞吐量是串行的 32 倍
```

### 2.4 带上界的 Early Termination

```cuda
// 利用邻居列表已排序性质，超过 upper_bound 时提前退出
for (auto i = thread_lane; i < lookup_size; i += WARP_SIZE) {
    auto key = lookup[i];
    int is_smaller = key < upper_bound ? 1 : 0;
    int found = 0;
    if (is_smaller && binary_search_2phase(...)) found = 1;
    num += found;
    unsigned mask = __ballot_sync(active, is_smaller);
    if (mask != FULL_MASK) break;  // 所有线程都超界，终止
}
```

**`__ballot_sync` 技巧**：一次 warp vote 判断是否所有线程的 key 都超界，如果是则立刻退出。

### 2.5 Compact Output with Ballot (交集结果写入)

```cuda
// 当需要实际写出交集结果（非仅计数）时
unsigned mask = __ballot_sync(active, found);
auto idx = __popc(mask << (WARP_SIZE - thread_lane - 1));  // 计算写入位置
if (found) c[count[warp_lane] + idx - 1] = key;           // 无冲突写入
if (thread_lane == 0) count[warp_lane] += __popc(mask);    // 更新总数
```

利用 `__ballot_sync` + `__popc` 实现 warp 内无冲突的 stream compaction，避免 atomic。

---

## 三、计算优化

### 3.1 Counting-Only Formula

对于某些 pattern，可以用组合公式避免显式枚举深层循环：

```cuda
// Diamond: 两个三角形共享一条边
// 对边 (v0,v1)，公共邻居数 n → diamond 数 = C(n,2) = n*(n-1)/2
AccType count = intersect_num(g.N(v0), v0_size, g.N(v1), v1_size);
auto n = warp_reduce<AccType>(count);
if (thread_lane == 0) counter += n * (n-1) / 2;
```

**节省**：从 $O(n^2)$ 的嵌套循环降为 $O(n)$ 的交集 + $O(1)$ 的公式。

### 3.2 DAG 对称性消除

```cuda
// 预处理：将无向图转 DAG，边只保留 src < dst 方向
ntasks = gg.init_edgelist(g, /*sym_break=*/true);  // edgelist 大小减半
```

**效果**：
- 任务量减半
- 无需运行时 `if (v > u) continue` 对称性检查
- 减少 warp divergence

### 3.3 Warp Reduce

```cuda
template <typename T>
__device__ T warp_reduce(T val) {
    val += SHFL_DOWN(val, 16);
    val += SHFL_DOWN(val, 8);
    val += SHFL_DOWN(val, 4);
    val += SHFL_DOWN(val, 2);
    val += SHFL_DOWN(val, 1);
    val  = SHFL(val, 0);       // broadcast 结果到所有线程
    return val;
}
```

5 次 shuffle 操作完成 warp 内归约，无需 shared memory。

---

## 四、归约与写回

### 4.1 CUB BlockReduce + atomicAdd

```cuda
typedef cub::BlockReduce<AccType, BLOCK_SIZE> BlockReduce;
__shared__ typename BlockReduce::TempStorage temp_storage;

// 每个线程累积局部 counter，最后一次性归约
AccType block_num = BlockReduce(temp_storage).Sum(counter);
if (threadIdx.x == 0) atomicAdd(total, block_num);
```

**优势**：每个 block 只发一次 atomicAdd（而非每个线程发一次），减少 atomic 竞争 256 倍。

---

## 五、内存管理

### 5.1 Frontier List 动态预算

```cuda
// 每个 warp 需要 (k-3) 层中间结果空间，每层 max_degree 个 vidType
size_t per_block_vlist_size = nwarps * (k-3) * max_deg * sizeof(vidType);

// 根据剩余显存动态决定可启动的 block 数
auto nb = (memsize - mem_graph) / per_block_vlist_size;
if (nb < nblocks) nblocks = nb;  // 避免 OOM

size_t list_size = nblocks * per_block_vlist_size;
cudaMalloc(&frontier_list, list_size);

// 每个 warp 通过 offset 访问自己的独占空间
vidType *vlist = &vlists[int64_t(warp_id) * int64_t(max_deg)];
```

### 5.2 Unified Memory Fallback

```cuda
if (mem_all > mem_gpu) {
    cudaMallocManaged(&d_colidx, nnz * sizeof(vidType));
    std::copy(h_colidx, h_colidx + nnz, d_colidx);
    // 可选 prefetch: cudaMemPrefetchAsync(d_colidx, ..., 0, NULL);
}
```

### 5.3 Shared Memory 布局

```cuda
__shared__ vidType list_size[WARPS_PER_BLOCK];            // 每 warp 的交集结果大小
__shared__ T cache[BLOCK_SIZE];                           // 2-phase search 的采样 pivot
__shared__ typename BlockReduce::TempStorage temp_storage; // CUB 归约
```

---

## 六、Occupancy 与启动配置

### 6.1 `__launch_bounds__` 编译期优化

```cuda
__global__ void __launch_bounds__(256, 8)  // max 256 threads/block, min 8 blocks/SM
warp_edge(...) { ... }
```

**效果**：编译器据此优化寄存器分配：
- 256 threads × 8 blocks = 2048 threads/SM（满占用率）
- 编译器会限制每线程寄存器数 ≤ 65536/(256×8) = 32 个

### 6.2 运行时 Occupancy 计算

```cuda
int max_blocks_per_SM = maximum_residency(warp_edge, nthreads, 0);
size_t max_blocks = max_blocks_per_SM * deviceProp.multiProcessorCount;
nblocks = std::min(6 * max_blocks, nblocks);  // 6x 超额订阅隐藏延迟
```

**6x 超额订阅**：比 SM 实际容量多 6 倍的 block，确保有足够 warp 隐藏内存延迟。

---

## 七、Iterative DFS（通用 k-pattern 支持）

```cuda
// 用显式栈模拟递归 DFS，一个 kernel 支持任意 k
int depth = 1;
Status state = Idle;
vidType stack[MAX_PATTERN_SIZE];
vidType idx[MAX_PATTERN_SIZE - 2];
__shared__ vidType vlist_sizes[WARPS_PER_BLOCK][MAX_PATTERN_SIZE - 2];

while (1) {
    if (depth == k-2) {
        // 叶子层：直接计数
        counter += intersect_num(vlist, v_size, g.N(u), u_size);
        depth--;  // 回溯
    } else if (state == Extending) {
        // 中间层：计算交集，存入 frontier_list
        auto count = intersect(v_ptr, v_size, g.N(u), u_size, out_list);
        vlist_sizes[warp_lane][depth-1] = count;
        idx[depth-1] = 0;
    }
    // 迭代或回溯
    if (idx[depth-1] == vlist_sizes[...]) { depth--; }
    else { stack[++depth] = vlist[idx[depth-1]++]; state = Extending; }
}
```

**vs 嵌套循环**：
- 嵌套循环：每个 k 值需要单独硬编码 kernel，但无额外控制开销
- 迭代 DFS：一个 kernel 支持所有 k，但有栈管理开销（~5-10% 性能损失）

---

## 八、Multi-GPU 策略

```cuda
// Edgelist 均匀分片
eidType n_tasks_per_gpu = (ne - 1) / ndevices + 1;

// 每个 GPU 持有完整 CSR 图，但只处理分配的 edgelist 段
for (int i = 0; i < ndevices; i++) {
    threads.push_back(std::thread([&,i]() {
        cudaSetDevice(i);
        d_graphs[i].init(g, i, ndevices);          // 全图 CSR
        d_graphs[i].copy_edgelist_to_device(...);   // 分片 edgelist
        warp_edge<<<nblocks, nthreads>>>(num_tasks[i], d_graphs[i], d_count[i]);
    }));
}
// 汇总
for (int i = 0; i < ndevices; i++) total += h_counts[i];
```

---

## 九、ParaCOSM 适用性总结

| 优化技术 | 适用于 ParaCOSM 的场景 | 优先级 |
|---|---|---|
| 2-phase binary search cache | `gpu_bfs_search.cu` 中邻居交集 | ⭐⭐⭐ 高 |
| Warp-centric edge-parallel | BFS expand kernel 的并行策略 | ⭐⭐⭐ 高 |
| `__launch_bounds__` + occupancy | 所有 GPU kernel | ⭐⭐⭐ 高 |
| Ballot compact output | 候选集过滤后写回 | ⭐⭐ 中 |
| Counting-only formula | 对称 pattern 的剪枝 | ⭐⭐ 中 |
| DAG 对称性消除 | 减少搜索空间 | ⭐⭐ 中 |
| Frontier list 内存预算 | 避免大图 OOM | ⭐⭐ 中 |
| CUB BlockReduce | 替代 naive atomicAdd | ⭐ 低（已在用） |
| Iterative DFS | 支持可变深度 pattern | ⭐ 低 |

---

## 参考来源

- 代码: `/home/v-haibinlai/haibin/GraphMiner/`
- 论文: Xuhao Chen, Arvind. "Efficient and Scalable Graph Pattern Mining on GPUs", OSDI 2022
- 关键文件:
  - `include/set_intersect.cuh` — 集合交集所有变体
  - `include/search.cuh` — binary_search_2phase 实现
  - `src/triangle/gpu_kernels/bs_warp_edge.cuh` — 最简经典 kernel
  - `src/clique/gpu_kernels/clique4_warp_edge.cuh` — 多层嵌套 pattern
  - `src/clique/gpu_kernels/edge_warp_iterative.cuh` — 迭代 DFS
  - `src/sgl/gpu_kernels/diamond_count.cuh` — counting-only 公式优化
