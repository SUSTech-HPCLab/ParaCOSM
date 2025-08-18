#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <string>

// 使用无序映射存储图的邻接表
typedef std::unordered_map<unsigned int, std::unordered_set<unsigned int>> Graph;

void bfs(const Graph& graph, unsigned int start_node) {
    // 用于保存已访问的节点
    std::unordered_set<unsigned int> visited;
    // BFS 队列
    std::queue<unsigned int> q;

    // 从 start_node 开始搜索
    visited.insert(start_node);
    q.push(start_node);

    // 记录最终可以到达的邻居节点数量
    unsigned int reachable_count = 0;

    while (!q.empty()) {
        unsigned int node = q.front();
        q.pop();

        // 遍历该节点的邻居
        for (const long long& neighbor : graph.at(node)) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                q.push(neighbor);
                reachable_count++;
            }
        }
    }

    // 输出最终可以访问的邻居节点数量
    std::cout << "Total reachable neighbors: " << reachable_count << std::endl;
}

Graph read_graph_from_file(const std::string& filename) {
    Graph graph;
    std::ifstream file(filename);
    std::string line;
    std::unordered_map<unsigned int, std::string> vertices;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        char type;
        ss >> type;

        if (type == 'v') {
            unsigned int vertex_id;
            std::string label;
            ss >> vertex_id ;
            // 读取顶点时，不需要做处理，只需记住顶点即可
            // vertices[vertex_id] = label;
        } else if (type == 'e') {
            unsigned int v1, v2;
            // std::string edge_label;
            ss >> v1 >> v2;
            // 添加边到邻接表
            graph[v1].insert(v2);
            graph[v2].insert(v1); // 因为是无向图，所以需要将v2也添加到v1的邻接表中
        }
    }

    return graph;
}

int main() {
    // 读取图文件
    std::string filename = "test_cache/random_graph_correct_01.graph";
    // std::string filename = "/home/haibin/CSM-Benchmark/ContinuousSubgraphMatching/parallel_tuils/test_cache/random_graph.graph";
    Graph graph = read_graph_from_file(filename);

    std::cout << "Graph loaded with " << graph.size() << " vertices." << std::endl;

    // 从 idx = 1 的节点开始进行 BFS
    unsigned int start_node = 1;
    bfs(graph, start_node);

    return 0;
}

/*
没有排列
(base)  sudo perf stat --repeat 5 -e cache-misses,cache-references,instructions,cycles ./bfs 
Graph loaded with 100000 vertices.
Total reachable neighbors: 99999
Graph loaded with 100000 vertices.
Total reachable neighbors: 99999
Graph loaded with 100000 vertices.
Total reachable neighbors: 99999
Graph loaded with 100000 vertices.
Total reachable neighbors: 99999
Graph loaded with 100000 vertices.
Total reachable neighbors: 99999

 Performance counter stats for './bfs' (5 runs):

        14,958,145      cache-misses              #   34.570 % of all cache refs      ( +-  4.43% )
        43,126,634      cache-references                                              ( +-  0.21% )
     8,040,874,252      instructions              #    1.55  insn per cycle           ( +-  0.01% )
     5,007,094,356      cycles                                                        ( +-  1.76% )

            1.8561 +- 0.0319 seconds time elapsed  ( +-  1.72% )

 Performance counter stats for './bfs' (50 runs):

        14,719,356      cache-misses              #   34.023 % of all cache refs      ( +-  0.68% )
        43,272,895      cache-references                                              ( +-  0.21% )
     8,040,573,755      instructions              #    1.58  insn per cycle           ( +-  0.06% )
     4,961,384,663      cycles                                                        ( +-  0.48% )

           1.83158 +- 0.00852 seconds time elapsed  ( +-  0.47% )

认真排列后？
50 次
 Performance counter stats for './bfs' (50 runs):

        16,805,688      cache-misses              #   38.029 % of all cache refs      ( +-  0.74% )
        43,967,379      cache-references                                              ( +-  0.29% )
     8,066,497,158      instructions              #    1.51  insn per cycle           ( +-  0.08% )
     5,319,082,472      cycles                                                        ( +-  0.61% )

            1.9207 +- 0.0122 seconds time elapsed  ( +-  0.64% )

 Performance counter stats for './bfs' (50 runs):

        15,875,974      cache-misses              #   36.809 % of all cache refs      ( +-  0.49% )
        43,152,126      cache-references                                              ( +-  0.02% )
     8,041,213,039      instructions              #    1.59  insn per cycle           ( +-  0.00% )
     5,173,903,087      cycles                                                        ( +-  0.19% )

           1.81201 +- 0.00354 seconds time elapsed  ( +-  0.20% )



*/ 