#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

void bfs(const SparseMatrix<int>& graph, int start_vertex) {
    vector<bool> visited(graph.rows(), false);
    queue<int> q;

    visited[start_vertex] = true;
    q.push(start_vertex);

    while (!q.empty()) {
        int vertex = q.front();
        q.pop();
        cout << "Visited: " << vertex << endl;

        for (SparseMatrix<int>::InnerIterator it(graph, vertex); it; ++it) {
            int neighbor = it.col();
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }
}

int main() {
    // 输入文件名
    string input_file = "data.graph";

    // 读取文件内容
    ifstream file(input_file);
    if (!file.is_open()) {
        cerr << "Error opening file: " << input_file << endl;
        return 1;
    }

    // 提取顶点和边信息
    vector<int> vertices;
    vector<pair<int, int>> edges;
    string line;
    while (getline(file, line)) {
        if (line[0] == 'v') {
            int vertex;
            sscanf(line.c_str(), "v %d", &vertex);
            vertices.push_back(vertex);
        } else if (line[0] == 'e') {
            int v1, v2;
            sscanf(line.c_str(), "e %d %d", &v1, &v2);
            edges.emplace_back(v1, v2);
        }
    }
    file.close();

    // 创建稀疏矩阵
    int num_vertices = *max_element(vertices.begin(), vertices.end()) + 1;
    SparseMatrix<int> sparse_matrix(num_vertices, num_vertices);

    for (const auto& edge : edges) {
        sparse_matrix.insert(edge.first, edge.second) = 1;
        sparse_matrix.insert(edge.second, edge.first) = 1; // 如果是无向图
    }

    // 打印稀疏矩阵
    cout << "Sparse Matrix:" << endl;
    cout << sparse_matrix << endl;

    // 执行 BFS
    int start_vertex = 0; // 从顶点 0 开始 BFS
    cout << "BFS starting from vertex " << start_vertex << ":" << endl;
    bfs(sparse_matrix, start_vertex);

    return 0;
}