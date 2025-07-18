#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <functional>
#include <atomic>
#include <unordered_set>
#include <iostream>
#include <future>
#include <map>
#include <chrono>

// 任务优先级
enum class TaskPriority {
    HIGH,
    MEDIUM,
    LOW
};

// 任务结构体
struct Task {
    std::function<void()> func;
    TaskPriority priority;

    Task(std::function<void()> f, TaskPriority p) : func(f), priority(p) {}
};

// 比较任务优先级
struct CompareTask {
    bool operator()(const Task& lhs, const Task& rhs) const {
        return static_cast<int>(lhs.priority) > static_cast<int>(rhs.priority);
    }
};

// 线程池类
class ThreadPool {
public:
    ThreadPool(size_t num_threads) : stop(false) {
        for(size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back([this] {
                while(true) {
                    Task task([]{}, TaskPriority::LOW);
                    {
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        condition.wait(lock, [this] {
                            return stop || !tasks.empty();
                        });
                        if(stop && tasks.empty()) return;
                        task = std::move(tasks.top());
                        tasks.pop();
                    }
                    task.func();
                    task_count--;
                }
            });
        }
    }

    template<class F>
    auto enqueue(F&& f, TaskPriority priority = TaskPriority::MEDIUM) -> std::future<decltype(f())> {
        auto task_ptr = std::make_shared<std::packaged_task<decltype(f())()>>(std::forward<F>(f));
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace([task_ptr]() { (*task_ptr)(); }, priority);
            task_count++;
        }
        condition.notify_one();
        return task_ptr->get_future();
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
    std::priority_queue<Task, std::vector<Task>, CompareTask> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
    std::atomic<int> task_count{0};
};


