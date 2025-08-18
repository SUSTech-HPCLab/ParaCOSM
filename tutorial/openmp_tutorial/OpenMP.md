
# OpenMP: Hello World! Hello World! Hello World! 

Good Learning repo and blogs:

https://github.com/muatik/openmp-examples.git

[OpenMP中几个容易混淆的函数（线程数量/线程ID/线程最大数）以及并行区域线程数量的确定_openmp设置线程数目-CSDN博客]
(https://blog.csdn.net/gengshenghong/article/details/7003110)

[OPENMP学习笔记（5）——运行时库函数_ompgetnumthreads-CSDN博客](https://blog.csdn.net/qq_23858785/article/details/97105972)

[Front Page :: OpenMP Validation & Verification](https://crpl.cis.udel.edu/ompvv/)



## Hello World! Hello World! Hello World! 

```CPP
int main() {

    #pragma omp parallel for
    for (int i = 0; i < 3; ++i) {
        cout << "Hello World!";
    }
    cout << endl;
    
    return 0;
}
```

Compile:
```CPP
g++ main.cpp -o hello_world -fopenmp
```

Output same as title

## OpenMP Configuration API

`omp_get_max_threads()` 
Get Maximum threads
```CPP
int max_threads = omp_get_max_threads();
```

`omp_set_num_threads()` 
set threads num. Use outside of parallel region, otherwise: "User Error 1001: omp_set_num_threads should only be called in serial regions"
```
omp_set_num_threads(5);
```


`omp_get_thread_num`
Returns the thread number (ID). The ID here refers to the thread's ID within the OpenMP team. In OpenMP, the threads within a team are assigned IDs in sequential order: 0, 1, 2, ...
Outside of a parallel region, it returns the ID of the master thread, which is 0. Inside a parallel region, each time this function is called, it returns the ID of the currently executing thread.


## How did gcc build OpenMP? 

https://www.haibinlaiblog.top/index.php/gcc%e6%98%af%e6%80%8e%e4%b9%88%e5%ae%9e%e7%8e%b0openmp%e7%9a%84%ef%bc%9f/


## Reference API


| Function                          | Description                                                                 |
|-----------------------------------|-----------------------------------------------------------------------------|
| OMP_SET_NUM_THREADS               | Sets the number of threads to be used in the next parallel region            |
| OMP_GET_NUM_THREADS               | Returns the number of threads in the parallel region executing the call      |
| OMP_GET_MAX_THREADS               | Returns the maximum value that can be returned by OMP_GET_NUM_THREADS        |
| OMP_GET_THREAD_NUM                | Returns the thread number within the team (not to be confused with total threads) |
| OMP_GET_THREAD_LIMIT              | Returns the maximum number of OpenMP threads available to the program        |
| OMP_GET_NUM_PROCS                 | Returns the number of processors available to the program                    |
| OMP_IN_PARALLEL                   | Determines whether the code being executed is in a parallel region           |
| OMP_SET_DYNAMIC                   | Enables or disables dynamic adjustment of the number of threads in parallel regions |
| OMP_GET_DYNAMIC                   | Determines whether dynamic thread adjustment is enabled                      |
| OMP_SET_NESTED                    | Enables or disables nested parallelism                                      |
| OMP_GET_NESTED                    | Determines whether nested parallelism is disabled                           |
| OMP_SET_SCHEDULE                  | Sets the loop scheduling policy when "runtime" is used as the schedule type in OpenMP directives |
| OMP_GET_SCHEDULE                  | Returns the loop scheduling policy when "runtime" is used as the schedule type |
| OMP_SET_MAX_ACTIVE_LEVELS         | Sets the maximum number of nested parallel regions                          |
| OMP_GET_MAX_ACTIVE_LEVELS         | Returns the maximum number of nested parallel regions                       |
| OMP_GET_LEVEL                     | Returns the nesting level of the current parallel region                    |
| OMP_GET_ANCESTOR_THREAD_NUM       | Returns the thread number of the ancestor thread at a given nesting level    |
| OMP_GET_TEAM_SIZE                 | Returns the size of the thread team at a given nesting level                |
| OMP_GET_ACTIVE_LEVEL              | Returns the number of nested active parallel regions containing the calling task |
| OMP_IN_FINAL                      | Returns true if the routine is executed in a final task region, false otherwise |
| OMP_INIT_LOCK                     | Initializes a lock associated with a lock variable                          |
| OMP_DESTROY_LOCK                  | Disassociates a given lock variable from all locks                          |
| OMP_SET_LOCK                      | Acquires ownership of a lock                                                |
| OMP_UNSET_LOCK                    | Releases a lock                                                            |
| OMP_TEST_LOCK                     | Attempts to set a lock but does not block if the lock is unavailable         |
| OMP_INIT_NEST_LOCK                | Initializes a nested lock associated with a lock variable                   |
| OMP_DESTROY_NEST_LOCK             | Disassociates a given nested lock variable from all locks                   |
| OMP_SET_NEST_LOCK                 | Acquires ownership of a nested lock                                         |
| OMP_UNSET_NEST_LOCK               | Releases a nested lock                                                     |
| OMP_TEST_NEST_LOCK                | Attempts to set a nested lock but does not block if the lock is unavailable  |
| OMP_GET_WTIME                     | Provides a portable wall-clock timing routine                               |
| OMP_GET_WTICK                     | Returns the number of seconds between successive clock ticks (double-precision) |


---


| libfunction                          | target                                                       |
|----------------------------------|------------------------------------------------------------|
| OMP_SET_NUM_THREADS              | 设置在下一个并行区域中使用的线程数                         |
| OMP_GET_NUM_THREADS              | 返回当前处于执行调用的并行区域中的线程数                   |
| OMP_GET_MAX_THREADS              | 返回调用OMP_GET_NUM_THREADS函数可以返回的最大值            |
| OMP_GET_THREAD_NUM               | 返回组内线程的线程号（译注：不要和线程总数搞混）            |
| OMP_GET_THREAD_LIMIT             | 返回可用于程序的最大OpenMP线程数                           |
| OMP_GET_NUM_PROCS                | 返回程序可用的处理器数                                     |
| OMP_IN_PARALLEL                  | 用于确定正在执行的代码是否是并行的                         |
| OMP_SET_DYNAMIC                  | 启动或者禁用可执行并行区域的线程数（由运行时系统）的动态调整 |
| OMP_GET_DYNAMIC                  | 用于确定是否启动了动态线程调整                             |
| OMP_SET_NESTED                   | 用于启用或者禁用嵌套并行                                   |
| OMP_GET_NESTED                   | 用于确定嵌套并行是否被弃用                                 |
| OMP_SET_SCHEDULE                 | 当“运行时”被用作OpenMP指令中的调度类型时，设置循环调度策略 |
| OMP_GET_SCHEDULE                 | 当“运行时”被用作OpenMP指令中的调度类型时，返回循环调度策略 |
| OMP_SET_MAX_ACTIVE_LEVELS        | 设置嵌套并行区域的最大数量                                 |
| OMP_GET_MAX_ACTIVE_LEVELS        | 返回嵌套并行区域的最大数量                                 |
| OMP_GET_LEVEL                    | 返回嵌套并行区域的最大数量                                 |
| OMP_GET_ANCESTOR_THREAD_NUM      | 给定当前线程的嵌套级别，返回其祖先线程的线程号             |
| OMP_GET_TEAM_SIZE                | 给定当前线程的嵌套级别，返回其线程组的大小                 |
| OMP_GET_ACTIVE_LEVEL             | 返回包含调用任务的的嵌套活动并行区域的数量                 |
| OMP_IN_FINAL                     | 如果在最终任务区域中执行该例程，则返回true；否则返回false   |
| OMP_INIT_LOCK                    | 初始化与锁变量相关联的锁                                   |
| OMP_DESTROY_LOCK                 | 解除给定的锁变量与所有锁的关联                             |
| OMP_SET_LOCK                     | 获取锁的所有权                                             |
| OMP_UNSET_LOCK                   | 释放锁                                                     |
| OMP_TEST_LOCK                    | 尝试设置锁，但是如果锁不可用，则不会阻止                   |
| OMP_INIT_NEST_LOCK               | 初始化与锁定变量关联的嵌套锁                               |
| OMP_DESTROY_NEST_LOCK            | 将给定的嵌套锁变量与所有锁解除关联                         |
| OMP_SET_NEST_LOCK                | 获取嵌套锁的所有权                                         |
| OMP_UNSET_NEST_LOCK              | 释放嵌套锁                                                 |
| OMP_TEST_NEST_LOCK               | 尝试设置嵌套锁，但如果锁不可用，则不会阻止                 |
| OMP_GET_WTIME                    | 提供便携式挂钟计时程序                                      |
| OMP_GET_WTICK                    | 返回连续时钟之间的秒数（双精度浮点值）                      |



----


## How did LLVM build OpenMP?


这个问题是我在做研究时，为了提升程序性能而做的探索。没想到这个知识后来也用上了。


## 基础知识：编译器的结构

在开启我们的探索旅程前，我们要简单介绍一下编译器的结构。

![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211026888.png)


编译器的“前端”和“后端” 负责编译过程的不同阶段。

- **前端**：负责将源代码解析成抽象语法树并生成中间代码，确保代码的语法和语义正确。
  - 词法分析、语法分析、语义分析和中间代码生成。
  
- **后端**：负责优化中间代码并生成目标机器代码，最终输出可执行文件。
  - 优化、目标代码生成、汇编、链接。

### 前端（Frontend）
编译器的前端负责从源代码到中间表示的转换，通常包括以下几个步骤：

1. **词法分析（Lexical Analysis）**：
   - 前端的第一步是将源代码分解成一系列的标记（tokens）。每个标记代表一个语法单元（如变量、运算符、关键字等）。
   - 例如，源代码 `int x = 10;` 会被分解为 `int`、`x` 、$=$ `10` 和 `;`。

2. **语法分析（Syntax Analysis）**：
   - 语法分析器将词法分析得到的标记流转换为抽象语法树（AST）。这棵树反映了源代码的结构和逻辑关系。
   - 语法分析器根据语法规则检测源代码是否符合语法要求。

3. **语义分析（Semantic Analysis）**：
   - 检查源代码的逻辑错误（如类型错误、变量未定义等）。
   - 生成符号表，记录变量、函数、类等的属性（如类型、作用域等）。

4. **中间代码生成（Intermediate Code Generation）**：
   - 在前端的最后，编译器会将抽象语法树转化为一种或多种中间表示（IR）。中间代码通常是比源代码更接近机器语言，但仍然是平台无关的代码，便于后端处理。（但是这部分LLVM的IR是最可读的）

AI总结：前端的主要任务是确保源代码的语法和语义正确，并生成适合优化和目标代码生成的中间表示（IR）。

### 后端（Backend）
编译器的后端负责将前端生成的中间表示转换为目标机器代码。后端主要包括以下几个步骤：

1. **优化（Optimization）**：
   - 首先会对中间代码进行各种优化。优化可以在多个层次上进行，包括循环优化、数据流分析、指令选择等。
   - 优化可以是局部的（如优化某个函数的性能）或全局的（如对整个程序的内存使用进行优化）。

2. **目标代码生成（Code Generation）**：
   - 将优化后的中间代码转换为目标机器代码。这是编译器的最后一步，生成的代码可以直接在特定的硬件平台上运行。
   - 目标代码生成通常包括将中间代码转换为汇编语言，随后汇编器将其转化为最终的机器码。

3. **汇编和链接（Assembly and Linking）**：
   - 汇编器将汇编语言代码转换为二进制机器码。
   - 链接器将多个目标文件、库文件链接成一个完整的可执行文件。它还负责处理符号引用和重定位。


## 编译器对支持OpenMP的工作

OpenMP（Open Multi-Processing）是一种用于共享内存并行编程的应用程序接口（API）。它允许开发人员利用多核处理器的优势，通过编写并行代码来加速计算任务。在编译器中，OpenMP的支持通常通过解析OpenMP指令并将其转换为适合并行执行的代码来实现。

![1745249338823.jpg](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211029647.jpg)


编译器对OpenMP的支持主要在：
- **识别并解析OpenMP指令（如pragma）**。
- **进行并行化分析和决策**，包括循环并行化、任务分配等。
- **管理线程创建和同步**，包括共享数据和私有数据的管理。
- **生成并行代码**，并利用OpenMP运行时库来管理线程和同步。
- **执行并行优化**，以提高程序的性能。

### 1. **编译器指令支持：识别、转换**
OpenMP使用编译指令（pragma）来指示编译器在特定代码块上应用并行化。编译器需要识别这些指令，并在代码的适当位置插入并行执行的控制结构。例如：

```c
#pragma omp parallel
{
    // 并行执行的代码块
}
```

编译器需要解析这种形式的`#pragma`指令，理解它的含义并对其后续代码进行并行化优化。`#pragma omp parallel`指示编译器将紧随其后的代码块并行执行。

### 2. **并行化决策分析**
编译器分析OpenMP指令，并根据程序结构决定哪些部分可以并行化。这个过程称为并行化分析。编译器通常会进行以下操作：
- **循环并行化**：如`#pragma omp parallel for`指令指示编译器将`for`循环的迭代分配到多个线程并行执行。
  
```CPP
  #pragma omp parallel for
  for (int i = 0; i < N; i++) {
      // 循环体
  }
```

  编译器会识别循环，并对循环进行拆分，生成适合多线程执行的代码。
  
- **任务分配和同步**：当程序中有多个任务（如计算任务）需要并行执行时，编译器会根据OpenMP指令将任务划分并分配给多个线程，同时管理线程之间的同步。

### 3. **线程管理和同步**
OpenMP程序通常涉及多个线程的创建与管理，编译器需要插入适当的代码来管理线程的生命周期、同步机制、数据共享等。这包括：
- **线程创建**：根据指令，编译器会为并行区域创建多个线程。
- **同步机制**：如临界区、屏障（`#pragma omp barrier`）、锁等，编译器会确保并行区域中的数据访问是安全的。
- **数据共享与私有化**：编译器通过指令（如`shared`、`private`）来决定哪些数据是共享的，哪些数据是每个线程私有的。

例如：

```c
#pragma omp parallel for shared(a) private(i)
for (int i = 0; i < N; i++) {
    // 每个线程执行的代码
}
```

### 4. **生成并行代码**

一旦编译器解析完OpenMP指令并进行相应的并行化优化，它会生成底层的并行代码。通常，编译器会生成线程创建、数据划分、同步等操作的代码。这些操作通常使用操作系统提供的线程库（如POSIX线程、Windows线程）或者专门的线程库（如Intel Threading Building Blocks）来实现。

- **OpenMP运行时库**：编译器在生成并行代码时，通常会依赖OpenMP运行时库，这个库负责线程的管理、调度、同步等操作。编译器生成的代码会调用这个库中的函数，以便管理线程的生命周期和同步。

### 5. **优化与性能提升**
编译器会通过不同的优化技术来提高并行程序的性能：
- **循环优化**：如循环拆分、循环合并等，减少冗余计算和同步。
- **负载均衡**：编译器会尽可能将工作均匀地分配给各个线程，以避免某些线程过载或空闲。
- **数据访问优化**：编译器可能会优化数据访问模式，减少内存访问冲突，提升内存访问效率。


### 编译器支持的OpenMP版本
现代编译器（如GCC、Clang、Intel编译器等）通常支持OpenMP的多个版本，并且不断改进对OpenMP指令的支持。具体的编译器版本支持的OpenMP特性可以通过编译器的文档进行查询。

比如GCC对OpenMP的支持可在这里查到：
[OpenMP (Using the GNU Compiler Collection (GCC))](https://gcc.gnu.org/onlinedocs/gcc/OpenMP.html)
![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211035011.png)


### 编译器启用OpenMP支持
要启用编译器对OpenMP的支持，通常需要通过编译选项来开启。例如，在GCC中，可以通过添加`-fopenmp`选项来启用OpenMP支持：

```bash
gcc -fopenmp my_program.c -o my_program
```

在`intel icpx`编译器中，需要启用 `-qopenmp` 开启支持


## LLVM对OpenMP的支持情况


[LLVM OpenMP: Parallel (fork/join)](https://openmp.llvm.org/doxygen/group__PARALLEL.html)

在 LLVM 编译器生态中，**OpenMP 的支持是通过 Clang 前端 + LLVM OpenMP 运行时库（libomp）实现的**。



### LLVM 对 OpenMP 的完整支持链
```bash
Clang 前端 → LLVM IR 生成 → OpenMP 运行时库（libomp） → 系统线程管理
           │
           └─ 自动插入并行运行时调用
```

#### 1. **编译器前端支持**
- **Clang** 负责解析 OpenMP 的编译制导语句（如 `#pragma omp parallel`）
- 将 OpenMP 指令转换为 LLVM 中间表示（IR），插入对 `libomp` 的运行时调用

  ```cpp
  // 原始代码
  #pragma omp parallel for
  for(int i=0; i<100; i++){ ... }

  // 转换后的 IR 伪代码
  __kmpc_fork_call(&loc, 1, microtask, &loop_bound);
  ```

#### 2. **运行时库 (libomp)**
- LLVM 自主开发的 OpenMP 运行时实现
- 核心功能：
  ```c
  // 线程池管理
  void __kmpc_fork_call(ident_t *loc, int argc, kmpc_micro microtask, ...);

  // 任务调度
  void __kmpc_for_static_init_4(ident_t *loc, int gtid, int schedule, ...);

  // 同步原语
  void __kmpc_barrier(ident_t *loc, int gtid);
  ```

#### 3. **架构支持**
| 特性             | LLVM 15 支持状态           |
| -------------- | ---------------------- |
| OpenMP 标准版本    | 5.2 (完整支持 4.5+)        |
| GPU Offloading | ✅ AMDGPU/NVPTX         |
| Target 指令      | ✅ `#pragma omp target` |
| 任务依赖           | ✅                      |
| SIMD 指令        | ✅ AVX2/AVX-512         |




## 关键函数解析


#### 1. **核心函数功能**
- **`__kmpc_fork_call`**  
  OpenMP **并行区域入口函数**，负责将串行代码分支为并行执行。  
  当编译器遇到 `#pragma omp parallel` 时，会生成对此函数的调用，完成以下操作：
  - 激活线程池中的工作线程
  - 将并行区域代码封装为微任务（microtask）
  - 分配共享/私有变量内存空间

- **`__kmpc_barrier`**  
  实现 **线程同步屏障**，确保所有线程到达同步点后才继续执行。  
  对应 OpenMP 中的隐式屏障（如 `parallel` 区域结束时）或显式 `#pragma omp barrier`。

#### 2. **命名规范 `kmpc`**
- **`kmp`**：源自 **Kuck & Associates Multi-Processing**（KAI 公司的并行技术遗产，后被 Intel 收购）
- **`c`**：表示属于 **C 语言接口**（另有 Fortran 接口使用 `kmpf` 前缀）
- 完整前缀含义：**KMP C Interface**

---

### 函数API定义与参数解析

#### 1. **`__kmpc_fork_call` 函数原型**
```c
void __kmpc_fork_call(
    ident_t    *loc,      // 源代码位置标识符
    int         argc,     // 共享参数数量
    microtask_t microtask,// 并行区域函数指针
    ...                   // 可变参数（共享变量地址）
);
```

- **参数详解**：

| 参数          | 类型            | 作用                         |
| ----------- | ------------- | -------------------------- |
| `loc`       | `ident_t*`    | 记录源代码位置（文件名、行号等），用于调试和性能分析 |
| `argc`      | `int`         | 传递给微任务的共享变量数量              |
| `microtask` | `microtask_t` | 实际并行执行的函数指针（编译器生成的包装函数）    |
| `...`       | `void*`       | 共享变量的地址列表（按参数顺序传递）         |

- **执行流程**：
  1. 主线程检查线程池状态
  2. 唤醒或创建工作线程
  3. 分配任务参数到共享内存区
  4. 各线程执行 `microtask` 函数
  5. 回收线程资源（可选）

源代码我们稍后会分析。
[LLVM OpenMP: Parallel (fork/join)](https://openmp.llvm.org/doxygen/group__PARALLEL.html#details)
[LLVM OpenMP: runtime/src/kmp_csupport.cpp Source File](https://openmp.llvm.org/doxygen/kmp__csupport_8cpp_source.html)


#### 2. **`__kmpc_barrier` 函数原型**
```c
void __kmpc_barrier(
    ident_t *loc,       // 源代码位置标识符
    int      gtid       // Global Thread ID（全局线程ID）
);
```

- **参数详解**：

| 参数     | 类型         | 作用                         |
| ------ | ---------- | -------------------------- |
| `loc`  | `ident_t*` | 记录源代码位置（文件名、行号等），用于调试信息    |
| `gtid` | `int`      | 当前线程的全局唯一ID（0为主线程，其他为工作线程） |

- **同步机制**：
  - 采用 **原子计数器 + 条件变量** 实现
  - 每个线程到达屏障时递减计数器
  - 最后一个线程触发条件变量广播唤醒所有线程

---



#### 规范定义来源
- **OpenMP Runtime API**  
  这些函数属于 OpenMP 标准的 **内部实现细节**，未在官方标准文档中公开规范，但各编译器实现遵循通用模式：

| 编译器实现      | 对应运行时库   | 头文件位置                        |
| ---------- | -------- | ---------------------------- |
| LLVM/Clang | libomp   | `omp.h` → `kmp.h`            |
| GCC        | libgomp  | `omp.h` → `gomp-constants.h` |
| Intel ICC  | libiomp5 | `omp.h` → `kmp_os.h`         |


### 调优

#### 1. **调试与性能分析**
- **查看运行时调用**：
  ```bash
  # 使用 GDB 跟踪调用链
  (gdb) break __kmpc_fork_call
  (gdb) backtrace
  # 查看并行区域参数
  (gdb) p *loc@4
  ```

- **性能计数器监控**：
  ```bash
  # Linux perf 统计屏障开销
  perf stat -e 'omp:barrier_wait' ./program
  ```

- VTune Profiler：

![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211040395.png)



### 问题排查

#### 1. **线程泄漏**
- **现象**：程序退出时卡死
- **诊断**：
  ```bash
  # 检查未回收的线程
  pstack <pid> | grep __kmpc_barrier
  ```
- **解决**：确保所有并行区域正确关闭

#### 2. **死锁**
- **常见原因**：
  - 不平衡的条件分支（如某些线程未到达屏障）
  - 在任务（task）区域内错误使用屏障
- **调试方法**：
  ```bash
  export OMP_DEBUG=1    # 启用运行时调试日志
  export KMP_BLOCKTIME=0
  ```

---

### 附：OpenMP 运行时调用层次

```c
// 用户代码
#pragma omp parallel
{
  // 并行区域
}

// 编译器生成代码
void .omp_outlined.(void *data) {
  // 用户代码逻辑
}

void main() {
  __kmpc_fork_call(&loc, 0, .omp_outlined.);
}

// 运行时调用链
__kmpc_fork_call
  ├─ __kmpc_begin                   // 初始化并行环境
  ├─ __kmp_launch_thread            // 启动工作线程
  └─ __kmpc_end                     // 清理资源
```


这是Pytorch CPU下使用GNU OpenMP + Intel MKL 下的程序调用截图。欢迎贡献LLVM的截图。

![d3ea2b5726da1fd0b48c39ab06147cf.png](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211049937.png)



这些底层函数构成了 OpenMP 并行化的核心机制，理解它们的工作原理对于调试性能关键型代码和深度优化至关重要。普通开发者通常无需直接调用这些接口，但掌握其原理有助于更高效地使用 OpenMP。


## 核心代码 `__kmpc_fork_call` 解析

源代码：
```C
/*!
@ingroup PARALLEL
@param loc  source location information
@param argc  total number of arguments in the ellipsis
@param microtask  pointer to callback routine consisting of outlined parallel
construct
@param ...  pointers to shared variables that aren't global
 
Do the actual fork and call the microtask in the relevant number of threads.
*/
void __kmpc_fork_call(ident_t *loc, kmp_int32 argc, kmpc_micro microtask, ...) {
  int gtid = __kmp_entry_gtid();
 
#if (KMP_STATS_ENABLED)
  // If we were in a serial region, then stop the serial timer, record
  // the event, and start parallel region timer
  stats_state_e previous_state = KMP_GET_THREAD_STATE();
  if (previous_state == stats_state_e::SERIAL_REGION) {
    KMP_EXCHANGE_PARTITIONED_TIMER(OMP_parallel_overhead);
  } else {
    KMP_PUSH_PARTITIONED_TIMER(OMP_parallel_overhead);
  }
  int inParallel = __kmpc_in_parallel(loc);
  if (inParallel) {
    KMP_COUNT_BLOCK(OMP_NESTED_PARALLEL);
  } else {
    KMP_COUNT_BLOCK(OMP_PARALLEL);
  }
#endif
 
  // maybe to save thr_state is enough here
  {
    va_list ap;
    va_start(ap, microtask);
 
#if OMPT_SUPPORT
    ompt_frame_t *ompt_frame;
    if (ompt_enabled.enabled) {
      kmp_info_t *master_th = __kmp_threads[gtid];
      ompt_frame = &master_th->th.th_current_task->ompt_task_info.frame;
      ompt_frame->enter_frame.ptr = OMPT_GET_FRAME_ADDRESS(0);
    }
    OMPT_STORE_RETURN_ADDRESS(gtid);
#endif
 
#if INCLUDE_SSC_MARKS
    SSC_MARK_FORKING();
#endif
    __kmp_fork_call(loc, gtid, fork_context_intel, argc,
                    VOLATILE_CAST(microtask_t) microtask, // "wrapped" task
                    VOLATILE_CAST(launch_t) __kmp_invoke_task_func,
                    kmp_va_addr_of(ap));
#if INCLUDE_SSC_MARKS
    SSC_MARK_JOINING();
#endif
    __kmp_join_call(loc, gtid
#if OMPT_SUPPORT
                    ,
                    fork_context_intel
#endif
    );
 
    va_end(ap);
 
#if OMPT_SUPPORT
    if (ompt_enabled.enabled) {
      ompt_frame->enter_frame = ompt_data_none;
    }
#endif
  }
 
#if KMP_STATS_ENABLED
  if (previous_state == stats_state_e::SERIAL_REGION) {
    KMP_EXCHANGE_PARTITIONED_TIMER(OMP_serial);
    KMP_SET_THREAD_STATE(previous_state);
  } else {
    KMP_POP_PARTITIONED_TIMER();
  }
#endif // KMP_STATS_ENABLED
}
 
```


关于`__kmpc_fork_call`实现机制

`__kmpc_fork_call`函数用于在OpenMP中创建一个并行区域（fork），在该区域内并行地执行给定的微任务（`microtask`）。函数通过多个线程来处理任务，并提供了一些计时和性能监控的功能。

它的函数基本思想就是，把我们要并行的函数识别为一个个微任务，然后线程池里的线程就可以每次从任务池里拿任务+做任务。对此我们需要微任务的函数指针，微任务的输入参数，以及环境变量或者共享变量。

执行图解：
![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211104007.png)

知乎里这个大佬讲的非常清楚：
[OpenMP task construct 实现原理以及源码分析 - 知乎](https://zhuanlan.zhihu.com/p/611558324)
![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211056423.png)


### 1. 函数参数

```c
void __kmpc_fork_call(ident_t *loc, kmp_int32 argc, kmpc_micro microtask, ...) {
  int gtid = __kmp_entry_gtid(); // 获取当前线程的全局线程ID (Global Thread ID)
  // ...
}
```
- `loc`: 该参数传递源代码的位置，用于错误报告和调试。
- `argc`: 传递给微任务的参数数量。用于后续处理变量的传递。
- `microtask`: 一个指向回调函数的指针，回调函数是并行执行的微任务。
- `...`: 表示一个可变参数列表，包含并行区域中要共享的变量。

### 2. 获取线程ID

```c
int gtid = __kmp_entry_gtid();
```
此行获取当前线程的线程ID（GTID, Global Thread ID）。GTID用于标识在OpenMP运行时系统中每个线程。

### 3. 性能统计和计时功能 `KMP_STATS_ENABLED`

在启用了统计功能的情况下，代码会记录当前并行区域的性能数据：
```c
#if (KMP_STATS_ENABLED)
  // If we were in a serial region, then stop the serial timer, record
  // the event, and start parallel region timer
  stats_state_e previous_state = KMP_GET_THREAD_STATE();
  if (previous_state == stats_state_e::SERIAL_REGION) {
    KMP_EXCHANGE_PARTITIONED_TIMER(OMP_parallel_overhead);
  } else {
    KMP_PUSH_PARTITIONED_TIMER(OMP_parallel_overhead);
  }
  int inParallel = __kmpc_in_parallel(loc);
  if (inParallel) {
    KMP_COUNT_BLOCK(OMP_NESTED_PARALLEL);
  } else {
    KMP_COUNT_BLOCK(OMP_PARALLEL);
  }
#endif
```


- **功能**：
  1. **状态记录**：`KMP_GET_THREAD_STATE()` 获取线程当前状态（串行或并行）。
  2. **计时器切换**：若从串行区域进入并行，交换计时器以统计并行开销。
  3. **嵌套并行检测**：通过 `__kmpc_in_parallel` 判断是否处于嵌套并行区域。
  4. **性能计数器**：`KMP_COUNT_BLOCK` 记录并行区域的启动次数（普通或嵌套）。



### 4. 可变参数处理 `OMPT_SUPPORT`
```c
va_list ap;
va_start(ap, microtask);
```
  - **可变参数**：通过 `va_start` 和后续的 `kmp_va_addr_of(ap)` 传递共享变量地址。

### 5. OMPT支持
```c
#if OMPT_SUPPORT
    ompt_frame_t *ompt_frame;
    if (ompt_enabled.enabled) {
      kmp_info_t *master_th = __kmp_threads[gtid];
      ompt_frame = &master_th->th.th_current_task->ompt_task_info.frame;
      ompt_frame->enter_frame.ptr = OMPT_GET_FRAME_ADDRESS(0);
    }
    OMPT_STORE_RETURN_ADDRESS(gtid);
#endif
```
- 如果启用了OMPT（OpenMP工具接口），会保存调用堆栈和返回地址等调试信息
	-   **帧信息记录**：`ompt_frame` 跟踪任务执行栈帧。
    - **返回地址存储**：`OMPT_STORE_RETURN_ADDRESS` 用于调试工具定位并行区域入口。




### 6. SSC Markers（性能分析）
```c
#if INCLUDE_SSC_MARKS
    SSC_MARK_FORKING();
#endif
```
- 如果启用了SSC（可能是某种性能分析框架），会在并行区域开始时插入标记，以便性能工具进行跟踪。

### 7. 核心调用
```c
__kmp_fork_call(loc, gtid, fork_context_intel, argc,
                VOLATILE_CAST(microtask_t) microtask, 
                VOLATILE_CAST(launch_t) __kmp_invoke_task_func,
                kmp_va_addr_of(ap));
```

- **这是实际的“fork”操作，创建并行区域，并启动每个线程来执行给定的回调函数（即微任务）**。
- `__kmp_fork_call`实际上会分配线程并执行每个线程的任务。

   **`__kmp_fork_call`**：
     - **线程团队创建**：激活或创建工作线程。
     - **任务分发**：将 `microtask` 和参数传递给线程池。
     - **参数传递**：`kmp_va_addr_of(ap)` 将可变参数转换为统一格式。

### 8. SSC Joining Markers（性能分析）
```c
#if INCLUDE_SSC_MARKS
    SSC_MARK_JOINING(); // 标记并行区域开始（用于Intel® VTune™等工具）
#endif
```

- 在并行区域结束时，插入一个结束标记。


### 9. 等待任务完成
```c
__kmp_join_call(loc, gtid); // 同步所有线程
```
- 等待所有线程完成任务，这意味着调用`__kmp_join_call`来确保所有分支的工作都执行完毕。

- **关键函数**：
   **`__kmp_join_call`**：
     - **屏障同步**：确保所有线程完成微任务执行。
     - **资源回收**：释放临时分配的内存和线程状态。

### 10. 结束可变参数处理
```c
va_end(ap); // 清理可变参数
```
- 结束对可变参数列表的处理。

### 11. 恢复OMPT状态
```c
#if OMPT_SUPPORT
    if (ompt_enabled.enabled) {
      ompt_frame->enter_frame = ompt_data_none;  // 重置帧信息
    }
#endif
```
- 如果启用了OMPT，恢复OMPT框架的状态。OMPT 是 OpenMP 提供的一个标准化接口，帮助开发者创建性能分析、调试或监控工具。
[4 OMPT Interface](https://www.openmp.org/spec-html/5.0/openmpch4.html)

### 12. 性能计时恢复
```c
#if KMP_STATS_ENABLED
  if (previous_state == stats_state_e::SERIAL_REGION) {
    KMP_EXCHANGE_PARTITIONED_TIMER(OMP_serial);
    KMP_SET_THREAD_STATE(previous_state);
  } else {
    KMP_POP_PARTITIONED_TIMER();
  }
#endif // KMP_STATS_ENABLED
```
- 最后，恢复之前的计时状态，结束并行区域的性能统计。


`__kmpc_fork_call` 函数主要用于在OpenMP并行计算中创建并行区域并执行任务。它不仅执行了多线程任务，还包括对性能的监控（通过计时和统计功能），并且支持调试和分析工具（如OMPT和SSC）。


##  核心调用 `__kmp_fork_call` 解析

[LLVM OpenMP: runtime/src/kmp.h File Reference](https://openmp.llvm.org/doxygen/kmp_8h.html#ab07305eec212868cbaea04361bbfce72)

![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211058473.png)

源码过长，放在本文最后。

`__kmp_fork_call` 是 OpenMP 运行时库的一部分，负责在遇到并行区域时处理 fork 操作。

![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211059524.png)





### `__kmp_fork_call` 的核心作用

`__kmp_fork_call` 是 LLVM OpenMP 运行时库中用于处理 OpenMP 并行区域（如 `#pragma omp parallel`）的核心函数。它负责创建和管理线程团队（team），以执行用户指定的并行任务（microtask）。其主要作用包括：

1. **初始化并行环境**：
   - 检查并初始化 OpenMP 运行时环境（如全局线程池、锁等）。
   - 设置主线程（primary thread）的信息，包括线程 ID（`gtid`）、当前任务和团队。

2. **确定线程数量**：
   - 根据用户指定的线程数（`OMP_NUM_THREADS` 或 `num_threads` 指令）、运行时配置（如 `OMP_DYNAMIC`）和任务限制（`task_thread_limit`），决定新并行区域的线程数量（`nthreads`）。
   - 如果线程数为 1，则进入序列化执行（调用 `__kmp_serial_fork_call`）。

3. **分配线程团队**：
   - 调用 `__kmp_allocate_team` 创建新的线程团队，分配线程资源。
   - 配置团队的属性，如线程绑定策略（`proc_bind`）、调度策略（`sched`）和内部控制变量（ICVs，如 `nproc`）。

4. **设置任务参数**：
   - 处理传递给微任务（`microtask`）的参数（`argc` 和 `argv`），这些参数通常是共享变量的地址。
   - 为团队中的线程设置调用栈、任务函数（`microtask`）和其他上下文信息。

5. **启动并行执行**：
   - 调用 `__kmp_fork_team_threads` 启动线程池中的线程，分配任务。
   - 主线程执行微任务（通过 `team->t.t_invoke(gtid)`），其他线程在后台并行执行。

6. **支持 OMPT 和调试**：
   - 如果启用了 OMPT（OpenMP Tools Interface），触发回调（如 `ompt_callback_parallel_begin` 和 `ompt_callback_parallel_end`），通知工具并行区域的开始和结束。
   - 支持性能分析工具（如 Intel VTune）通过 ITT（Instrumentation and Tracing Technology）接口记录并行区域信息。

7. **管理嵌套并行和 teams 构造**：
   - 处理 OpenMP 的嵌套并行（如 `#pragma omp parallel` 内的并行区域）和 `teams` 构造（如 `#pragma omp teams`）。
   - 维护嵌套级别（`level` 和 `active_level`）和线程限制，确保符合 OpenMP 规范。

8. **释放资源**：
   - 在并行区域结束后，释放锁（如 `forkjoin_lock`）并更新线程状态。
   - 如果启用了热团队（hot team），缓存团队以便后续重用。

---

### 代码结构和关键逻辑

以下是代码中关键部分的简要分析：

#### 1. **初始化和环境检查**
 检查运行时是否已初始化，并从软暂停中恢复。
 获取当前线程的信息和父级团队结构。
```c
if (!TCR_4(__kmp_init_parallel))
    __kmp_parallel_initialize();
__kmp_resume_if_soft_paused();
```
- 检查是否已初始化 OpenMP 运行时（`__kmp_init_parallel`），若未初始化则调用 `__kmp_parallel_initialize`。
- 恢复可能的暂停状态（`__kmp_resume_if_soft_paused`），确保运行时处于活跃状态。



#### 2. **获取主线程信息**
通过 `__kmp_fork_in_teams` 特殊处理嵌套在 `teams` 构造中的并行区域。
```c
master_th = __kmp_threads[gtid];
parent_team = master_th->th.th_team;
master_tid = master_th->th.th_info.ds.ds_tid;
root = master_th->th.th_root;
```
- 从全局线程数组 `__kmp_threads` 获取主线程（`master_th`）信息。
- 提取当前团队（`parent_team`）、线程 ID（`master_tid`）和根线程（`root`）。

#### 3. **确定线程数量**

使用 `master_set_numthreads` 或 ICV（内部控制变量）计算线程数（`nthreads`）。
考虑 `task_thread_limit` 对目标区域的影响，并检查嵌套的活跃级别。

```c
if ((!enter_teams &&
     (parent_team->t.t_active_level >=
      master_th->th.th_current_task->td_icvs.max_active_levels)) ||
    (__kmp_library == library_serial)) {
    nthreads = 1;
} else {
    nthreads = master_set_numthreads
                   ? master_set_numthreads
                   : get__nproc_2(parent_team, master_tid);
    nthreads = task_thread_limit > 0 && task_thread_limit < nthreads
                   ? task_thread_limit
                   : nthreads;
}
```
- 如果处于序列化模式（`__kmp_library == library_serial`）或嵌套级别超限（`t_active_level >= max_active_levels`），设置 `nthreads = 1`，进入序列化执行。
- 否则，根据用户设置（`master_set_numthreads`）或默认线程数（`get__nproc_2`）确定线程数，并应用任务线程限制（`task_thread_limit`）。

#### 4. **分配锁和线程资源**

如果 `nthreads == 1`，调用 `__kmp_serial_fork_call` 串行执行，而不创建团队。
**锁获取：** 使用 `__kmp_forkjoin_lock` 安全地保留线程。
- **ICV 传播：** 将内部控制变量（`nproc`、`proc-bind`）复制到新团队中。
- **团队分配：** 使用 `__kmp_allocate_team` 分配并初始化新的 `kmp_team_t` 结构。
（什么是ICV？）
在 OpenMP 中，**ICV** 是 **Internal Control Variables（内部控制变量）** 的缩写。它们是由 OpenMP 运行时用于控制并行执行行为的变量。每个线程和团队通常都会有一组自己的 ICV，用于决定如何分配资源、调度任务、处理线程绑定等。

常见的 ICV 变量包括：

1. **nproc（number of processors）：** 表示可用于并行执行的处理器数量。它通常决定了并行区域的线程数或任务划分的策略。
2. **proc-bind：** 指定线程绑定策略，控制线程如何绑定到物理处理器或核心上。常见的绑定策略包括将线程绑定到特定的核心上，或者将线程分散到多个核心。

**ICV 传播** 是指将这些内部控制变量从父团队（parent team）传递到新创建的子团队，以确保子团队在并行执行时能保持一致的执行策略和资源分配。这样做能够确保新的团队在执行时继承父团队的设置（如处理器数量和线程绑定策略），从而保证 OpenMP 程序的正确性和性能。

例如，在创建一个新团队时，可能需要将父团队的 ICV（如 `nproc` 和 `proc-bind`）复制到新团队中，以便新团队能够使用相同的资源管理和线程调度策略。

```c
if (nthreads > 1) {
    __kmp_acquire_bootstrap_lock(&__kmp_forkjoin_lock);
    nthreads = __kmp_reserve_threads(root, parent_team, master_tid, nthreads, enter_teams);
    if (nthreads == 1) {
        __kmp_release_bootstrap_lock(&__kmp_forkjoin_lock);
    }
}
```
- 如果需要多个线程，获取全局锁（`__kmp_forkjoin_lock`）以同步线程分配。
- 调用 `__kmp_reserve_threads` 预留线程资源，动态调整线程数。
- 如果最终线程数为 1，释放锁并进入序列化模式。



#### 5. **创建线程团队**
```c
team = __kmp_allocate_team(root, nthreads, nthreads, ...);
```
- 调用 `__kmp_allocate_team` 创建新的线程团队，设置团队的线程数、绑定策略（`proc_bind`）和 ICVs。
- 团队结构（`kmp_team_t`）包含线程信息、任务函数和上下文。

#### 6. **设置任务参数**
将来自父团队或 `va_list` 的参数复制到新团队的 `argv` 中。
```c
argv = (void **)team->t.t_argv;
if (ap) {
    for (i = argc - 1; i >= 0; --i) {
        void *new_argv = va.arg(kmp_va_deref(ap), void *);
        KMP_CHECK_UPDATE(*argv, new_argv);
        argv++;
    }
} else {
    for (i = 0; i < argc; ++i) {
        KMP_CHECK_UPDATE(argv[i], team->t.t_parent->t.t_argv[i]);
    }
}
```

- 参数通常是共享变量的地址，按顺序存储在 `team->t.t_argv` 中。

#### 7. **启动线程**

- 调用 `__kmp_fork_team_threads` 创建并唤醒工作线程。
- 为新团队的线程设置 ICV 复制。
```c
__kmp_fork_team_threads(root, team, master_th, gtid, !ap);
```
- 调用 `__kmp_fork_team_threads` 启动团队中的线程，分配任务给工作线程。
- 主线程负责执行微任务（`team->t.t_invoke(gtid)`）。

#### 8. **OMPT 支持**

主线程直接执行用户的并行区域（微任务）。
```c
if (ompt_enabled.enabled) {
    ompt_callbacks.ompt_callback(ompt_callback_parallel_begin)(...);
}
```
- 如果启用了 OMPT，触发并行区域开始的回调（`ompt_callback_parallel_begin`），通知工具并行区域的启动。
- 类似地，`__kmp_join_ompt` 在并行区域结束时触发 `ompt_callback_parallel_end`。

#### 9. **释放资源**
释放 fork-join 锁并处理 ITT/OMPT 跟踪。
```c
__kmp_release_bootstrap_lock(&__kmp_forkjoin_lock);
```
- 在线程启动后释放全局锁，允许其他并行区域继续执行。

- `__kmp_forkjoin_lock` 在线程保留时获取，并在团队设置后释放。这个锁对于防止线程过度订阅至关重要。
- 使用 `KMP_CHECK_UPDATE` 表明对共享团队/线程结构进行原子或有序写入，以确保线程之间的可见性。

---

### 线程池管理

线程池管理的核心逻辑在 `__kmp_fork_call` 和相关函数（如 `__kmp_fork_team_threads` 和 `__kmp_allocate_team`）中实现：

1. **线程池结构**：
   - 线程池通过全局数组 `__kmp_threads` 管理，每个线程有 `kmp_info_t` 结构，记录其状态、ID 和任务。
   - 团队（`kmp_team_t`）包含一组线程，缓存为“热团队”（hot team）以减少创建开销。
![](https://blog-1327458544.cos.ap-guangzhou.myqcloud.com/N2025/202504211115853.png)


2. **热团队（Hot Team）**：
   - 代码中的 `KMP_NESTED_HOT_TEAMS` 宏启用热团队支持，缓存团队在 `master_th->th.th_hot_teams` 中。
   - 热团队在并行区域结束后保持活跃，供后续并行区域重用。

3. **线程分配**：
   - `__kmp_reserve_threads` 根据可用线程和限制动态分配线程。
   - 如果线程不足，可能减少 `nthreads` 或序列化执行。

4. **线程重用**：
   - 线程在并行区域结束后进入等待状态（通过 `__kmp_suspend_initialize_thread`），等待下一次任务分配。
   - 环境变量 `OMP_THREAD_LIMIT` 和 `KMP_ALL_THREADS` 控制线程池的最大规模。

---

### `__kmp_join_ompt` 的作用

`__kmp_join_ompt` 是 OMPT 相关的辅助函数，负责在并行区域结束时：

1. **触发 OMPT 回调**：
   - 调用 `ompt_callback_parallel_end`，通知工具并行区域的结束。
   - 传递并行区域数据（`parallel_data`）、任务数据（`task_info->task_data`）和调用上下文（`codeptr`）。

2. **恢复线程状态**：
   - 调用 `__kmp_join_restore_state`，将线程状态恢复为 `ompt_state_work_serial`（序列化区域）或 `ompt_state_work_parallel`（并行区域）。
   - 清除任务框架（`task_info->frame.enter_frame = ompt_data_none`）。

---

### AI总结

`__kmp_fork_call` 是 LLVM OpenMP 运行时库中处理并行区域的核心函数，负责初始化环境、分配线程团队、设置任务参数和启动并行执行。它通过线程池管理（热团队、线程重用）和 OMPT 支持实现高效的并行计算。线程池管理的关键代码在 `kmp_runtime.cpp` 中，尤其集中在 `__kmp_allocate_team` 和 `__kmp_fork_team_threads` 函数中。


---


### 三、调试与优化技巧

#### 1. 运行时环境控制
```bash
# 设置线程数
export OMP_NUM_THREADS=4

# 绑定线程到物理核心
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# 调试运行时
export OMP_DISPLAY_ENV=VERBOSE
export OMP_DEBUG=1
```

#### 2. 性能分析工具链
```bash
# 使用 LLVM 的 OpenMP 性能计数器
clang -fopenmp -fopenmp-targets=nvptx64 -fopenmp-version=51 \
     -Rpass=openmp-opt -Rpass-analysis=openmp-opt

# 生成优化报告
clang -fsave-optimization-record -fopenmp
```

#### 3. 常见问题诊断
**问题现象**：未找到 OpenMP 符号
```bash
error: undefined reference to `__kmpc_fork_call'
```
**解决方案**：
```bash
# 显式链接 libomp
clang -fopenmp -lomp main.c
# 或指定静态链接
clang -fopenmp -Wl,-Bstatic -lomp -Wl,-Bdynamic
```

**问题现象**：SIMD 向量化失败
```bash
warning: loop not vectorized: failed explicitly specified vectorization
```
**优化策略**：
```cpp
#pragma omp simd collapse(2) safelen(16)
for(int i=0; i<M; i++){
  for(int j=0; j<N; j++){
    // 确保内存连续访问模式
    arr[i][j] = ... 
  }
}
```



## 核心fork源码

https://openmp.llvm.org/doxygen/group__PARALLEL.html

第1927行
