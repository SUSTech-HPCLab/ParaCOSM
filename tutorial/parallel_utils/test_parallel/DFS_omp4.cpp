#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <functional>
#include <atomic>
#include <unordered_set>
#include <iostream>

// 节点结构体
struct Node {
    std::vector<Node*> children;
    int data; // 节点数据
    Node(int val) : data(val) {}
};

// 初始化图结构
Node* initialize_graph() {
    Node* root = new Node(1);
    Node* node2 = new Node(2);
    Node* node3 = new Node(3);
    Node* node4 = new Node(4);
    Node* node5 = new Node(5);

    root->children.push_back(node2);
    root->children.push_back(node3);
    node2->children.push_back(node4);
    node2->children.push_back(node5);

    return root;
}

// 线程池类
class ThreadPool {
public:
    ThreadPool(size_t num_threads) : stop(false) {
        for(size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back([this] {
                while(true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        condition.wait(lock, [this] {
                            return stop || !tasks.empty();
                        });
                        if(stop && tasks.empty()) return;
                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                    task_count--;
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace([this, f]() {
                f();
            });
            task_count++;
        }
        condition.notify_one();
    }

    void wait_all() {
        while (task_count > 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
    std::atomic<int> task_count{0};
};

// 并发访问控制类
class ConcurrentVisited {
    static constexpr int NUM_BUCKETS = 16;
    std::vector<std::mutex> mutexes;
    std::vector<std::unordered_set<Node*>> buckets;

public:
    ConcurrentVisited() : mutexes(NUM_BUCKETS), buckets(NUM_BUCKETS) {}

    bool try_visit(Node* node) {
        size_t hash = std::hash<Node*>()(node);
        size_t bucket = hash % NUM_BUCKETS;
        std::lock_guard<std::mutex> lock(mutexes[bucket]);
        if (buckets[bucket].count(node)) return false;
        buckets[bucket].insert(node);
        return true;
    }
};

// 处理节点的具体逻辑
void process_node(Node* node) {
    // 实现具体的节点处理逻辑
    std::cout << "Processing node with data: " << node->data << std::endl;
}

// 并行DFS函数
void parallel_dfs(Node* node, int depth, int max_depth,
                  ThreadPool& pool, ConcurrentVisited& visited) {
    if (!visited.try_visit(node)) return;

    process_node(node);

    if (depth < max_depth) {
        // 并行处理子节点
        for (Node* child : node->children) {
            pool.enqueue([child, depth, max_depth, &pool, &visited]() {
                parallel_dfs(child, depth + 1, max_depth, pool, visited);
            });
        }
    } else {
        // 串行递归处理
        for (Node* child : node->children) {
            parallel_dfs(child, depth + 1, max_depth, pool, visited);
        }
    }
}

int main() {
    Node* root = initialize_graph(); // 初始化图结构
    const int max_depth = 3;         // 根据实际情况调整
    ThreadPool pool(std::thread::hardware_concurrency());
    ConcurrentVisited visited;

    pool.enqueue([&] {
        parallel_dfs(root, 0, max_depth, pool, visited);
    });

    pool.wait_all(); // 等待所有任务完成
    return 0;
}