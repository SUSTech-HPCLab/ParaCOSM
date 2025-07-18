#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <functional>
#include <atomic>
#include <unordered_set>
#include <iostream>
#include "ThreadPool.h"

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
    Node* node6 = new Node(6);
    Node* node7 = new Node(7);
    Node* node8 = new Node(8);
    Node* node9 = new Node(9);
    Node* node10 = new Node(10);
    Node* node11 = new Node(11);
    Node* node12 = new Node(12);
    Node* node13 = new Node(13);

    root->children.push_back(node2);
    root->children.push_back(node3);
    root->children.push_back(node4);
    node2->children.push_back(node5);
    node2->children.push_back(node6);
    // node3->children.push_back(node6);
    node3->children.push_back(node7);
    // node4->children.push_back(node7);
    // node5->children.push_back(node6);
    // node5->children.push_back(node7);
    // node6->children.push_back(node7);
    node3->children.push_back(node8);
    // node7->children.push_back(node8);
    node8->children.push_back(node9);
    node4->children.push_back(node10);
    // node9->children.push_back(node10);
    node4->children.push_back(node11);
    // node10->children.push_back(node11);
    node10->children.push_back(node13);
    node11->children.push_back(node12);
    // node11->children.push_back(node13);


    return root;
}


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
        // std::lock_guard<std::mutex> lock(mutexes[bucket]);
        if (buckets[bucket].count(node)) return false;
        buckets[bucket].insert(node);
        return true;
    }
};

// 处理节点的具体逻辑
void process_node(Node* node) {
    // 实现具体的节点处理逻辑
    std::cout << "Processing node with data: " << node->data << std::endl;
    // print thread
    std::cout << "Thread " << std::this_thread::get_id() << '\n';
    std::this_thread::sleep_for(std::chrono::milliseconds(800)); // 成功！
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
    const int max_depth = 2;         // 根据实际情况调整
    ThreadPool pool(std::thread::hardware_concurrency());
    ConcurrentVisited visited;

    pool.enqueue([&] {
        parallel_dfs(root, 1, max_depth, pool, visited);
    });

    pool.wait_all(); // 等待所有任务完成
    return 0;
}