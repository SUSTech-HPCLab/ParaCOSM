#include <iostream>
#include <vector>
#include <stack>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <memory>
#include <functional>

constexpr int DEPTH_THRESHOLD = 3;

struct Node {
    std::vector<std::unique_ptr<Node>> children;
    
    const std::vector<std::unique_ptr<Node>>& get_children() const {
        return children;
    }
};

class ThreadPool {
public:
    explicit ThreadPool(size_t threads) : stop(false) {
        for(size_t i = 0; i < threads; ++i) {
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
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker : workers)
            if(worker.joinable())
                worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

class CountDownLatch {
public:
    explicit CountDownLatch(int count) : count(count) {}
    
    void wait() {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return count == 0; });
    }
    
    void count_down() {
        std::unique_lock<std::mutex> lock(mtx);
        if(--count == 0) {
            cv.notify_all();
        }
    }

private:
    std::mutex mtx;
    std::condition_variable cv;
    int count;
};

class ParallelDFS {
public:
    explicit ParallelDFS(Node& root) 
        : pool(std::thread::hardware_concurrency()),
          latch(1) {
        pool.enqueue([this, &root] {
            dfs(root, 0, true);
        });
    }

    void run() {
        latch.wait();
    }

private:
    ThreadPool pool;
    CountDownLatch latch;
    std::mutex print_mutex;

    void process_node(const Node& node) {
        std::lock_guard<std::mutex> lock(print_mutex);
        std::cout << "Processing node on thread: " 
                  << std::this_thread::get_id() << '\n';
    }

    void dfs(Node& node, int depth, bool is_root) {
        process_node(node);
        
        const auto& children = node.get_children();
        if(children.empty()) {
            if(is_root) latch.count_down();
            return;
        }

        if(depth < DEPTH_THRESHOLD) {
            for(auto& child : children) {
                dfs(*child, depth + 1, false);
            }
            if(is_root) latch.count_down();
        } else {
            auto child_latch = std::make_shared<CountDownLatch>(children.size());
            
            for(auto& child : children) {
                Node* raw_child = child.get();  // C++17兼容捕获方式
                pool.enqueue([this, raw_child, depth, child_latch] {
                    iterative_dfs(raw_child, depth + 1);
                    child_latch->count_down();
                });
            }
            
            child_latch->wait();
            if(is_root) latch.count_down();
        }
    }

    void iterative_dfs(Node* start_node, int start_depth) {
        struct StackFrame {
            Node* node;
            int depth;
        };
        
        std::stack<StackFrame> stack;
        stack.push({start_node, start_depth});

        while(!stack.empty()) {
            auto current = stack.top().node;  // C++17兼容解包方式
            auto depth = stack.top().depth;
            stack.pop();
            
            process_node(*current);
            
            for(auto& child : current->get_children()) {
                if(depth + 1 < DEPTH_THRESHOLD) {
                    stack.push({child.get(), depth + 1});
                } else {
                    Node* raw_child = child.get();
                    pool.enqueue([this, raw_child, depth] {
                        iterative_dfs(raw_child, depth + 1);
                    });
                }
            }
        }
    }
};

int main() {
    Node root;
    // 构建测试树结构
    auto& children = root.children;
    for(int i = 0; i < 3; ++i) {
        children.emplace_back(std::make_unique<Node>());
        for(int j = 0; j < 3; ++j) {
            children[i]->children.emplace_back(std::make_unique<Node>());
        }
    }

    ParallelDFS dfs(root);
    dfs.run();
}
