#include <bits/stdc++.h>
#include "iostream"
#include "omp.h"
using namespace std;

//this is class of graph contains a map for visited
//nodes and adj map for list of edges per node.
class Graph {
public:
    map<int, bool> visited;
    map<int, vector<int> > adj;

    // function to add an edge to graph
    void addEdge(int v, int w);

    // DFS traversal of the vertices
    // reachable from v
    void DFS(int v);
};

//This method is for adding edges
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
}

// DFS traversal of the vertices
// reachable from v
void Graph::DFS(int v)
{
    //in this section we are checking that if a node is not visited
    //we will turn it as visited and will print it
    //this section is critical since we are writing
    //to the visited map
#pragma omp critical
    {
        if(!visited[v]){
            visited[v] = true;
            //printing the vertex number and threads that are executing
            //at the same time
            cout<<v<<"----"<<omp_get_num_threads()<<endl;
        }
    }

    //in this section we start parallelize the algorithm
    //and after that run each one in specific task
    //the benefit af task is that it don't need to run exactly after running
    for (int i = 0; i < adj[v].size();i++) {
        if (!visited[adj[v][i]]) {
#pragma omp task
            DFS(adj[v][i]);
        }
    }
}

// Driver code
int main()
{
    Graph g;

    ifstream fin;

    // in this part i'm opening a file contain a graph
    //each line is edges of that line number vertex
    //consider that vertex numbers start  from 0
    fin.open("/home/hossein/CLionProjects/untitled/graph.txt");
    string line ;
    int counter = 0;
    while (fin) {
        // Read a Line from File
        getline(fin, line);

        size_t pos = 0;
        std::string token;
        string delimiter = " ";
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            int number = stoi(token);
            g.addEdge(counter,number);
            line.erase(0, pos + delimiter.length());
        }
        counter++;
    }

    // Close the file
    fin.close();
    //start clock to calculate the duration of execution
    clock_t start = clock();

    //start point is node that the algorithm start from
    int startPoint = 297;

    //First we print the start node and make it visited
    g.visited[startPoint] = true;
    cout<<startPoint<<endl;
    //here we start the algorithm
#pragma omp parallel
#pragma omp single
    {
        g.DFS(startPoint);
    };

    clock_t stop = clock();
    //printing the duration of execution of the algorithm
    double e = (double)(stop-start)/CLOCKS_PER_SEC;
    cout<<e;
    return 0;
}