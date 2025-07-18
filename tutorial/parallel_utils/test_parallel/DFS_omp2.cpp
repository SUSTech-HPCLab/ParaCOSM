#include<bits/stdc++.h>
#include<omp.h>
using namespace std;

int main(){
    
    //Input Parameters
    int n;
    cout << "Enter number of nodes: ";
    cin >> n;
    vector<vector<int>> G(n + 1);
    int edges;
    cout << "Enter number of edges: ";
    cin >> edges;
    cout << "---Input for Acyclic Directed Graph---\n";
    while(edges > 0){
        int parent, child;
        cout << "Enter parent and child: ";
        cin >> parent >> child;
        G[parent].push_back(child);
        edges -= 1;
    }
    int src, dst;
    cout << "Enter source and destination: ";
    cin >> src >> dst;
    if(src == dst){
        cout << src << "\n";
        return 0;
    }
    
    //Defining threads for the branching factor from source node
    int num_threads = G[src].size();
    if(num_threads == 0){
        cout << "No Solution";
        return 0;
    }
    vector<pair<int,string>> branch_paths;
    for(auto x : G[src]) branch_paths.push_back({x, to_string(src) + "->" + to_string(x)});
    
    //Starting Parallel DFS
    string solution = "";
    bool found = false;
    double start_time = omp_get_wtime();
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        //Letting each branch take a stack
        int thread_id = omp_get_thread_num();
        stack<pair<int,pair<int,string>>> dfs_parallel_stack;
        dfs_parallel_stack.push({0,{branch_paths[thread_id].first,branch_paths[thread_id].second}});
        
        //Running till stack is empty
        while(!dfs_parallel_stack.empty()){
            
            if(found) break;
            pair<int,pair<int,string>> cur = dfs_parallel_stack.top();
            dfs_parallel_stack.pop();
            int node = cur.second.first, next_node_idx = cur.first;
            string current_path = cur.second.second;
            
            //Critical Section to access solution variable
            #pragma omp critical
            {
                if(node == dst){
                    solution = current_path;
                    found = true;
                }
            }
            if(found) break;

            
            //next node condition
            if(next_node_idx < G[node].size()){
                dfs_parallel_stack.push({next_node_idx + 1,{node, current_path}});
                int next_node = G[node][next_node_idx];
                dfs_parallel_stack.push({0,{next_node,current_path + "->" + to_string(next_node)}});
            }
        }
    }
    double end_time = omp_get_wtime();
    cout << "-----Parallel DFS-----\n";
    cout << "Number of Threads: " << num_threads << "\n";
    if(solution.size() == 0) cout << "No solution available\n";
    else cout << "Path is: " << solution << "\n";
    cout << "Computed in " << end_time - start_time << " units of time\n";
    
    //Starting Normal DFS
    solution = "";
    start_time = omp_get_wtime();
    stack<pair<int,pair<int,string>>> dfs_normal_stack;
    dfs_normal_stack.push({0,{src, to_string(src)}});
    
    //Running till stack is empty
    while(!dfs_normal_stack.empty()){
        pair<int,pair<int,string>> cur = dfs_normal_stack.top();
        dfs_normal_stack.pop();
        int node = cur.second.first, next_node_idx = cur.first;
        string current_path = cur.second.second;
        
        //Breaking the while loop if destination found
        if(solution.size() != 0) break;
        if(node == dst) {
            solution = current_path;
            break;
        }
        
        //next node condition
        if(next_node_idx < G[node].size()){
            dfs_normal_stack.push({next_node_idx + 1,{node, current_path}});
            int next_node = G[node][next_node_idx];
            dfs_normal_stack.push({0,{next_node,current_path + "->" + to_string(next_node)}});
        }
    }
    end_time = omp_get_wtime();
    cout << "-----Normal DFS-----\n";
    if(solution.size() == 0) cout << "No solution available\n";
    else cout << "Path is: " << solution << "\n";
    cout << "Computed in " << end_time - start_time << " units of time\n";
    return 0;
}