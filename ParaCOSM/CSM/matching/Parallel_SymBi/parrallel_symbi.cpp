#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
#include <stack>
#include <cstring> // memcpy

#include <thread>
#include <mutex>
#include <tbb/task_group.h>
#include <tbb/task_arena.h>

#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/Parallel_SymBi/parallel.h"

#define Print_Time_Nano(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() << " ns" << std::endl


/**
 * @brief Constructor for the Parrllel_SymBi class.
 * 
 * This constructor initializes the data structures and parameters required for the 
 * parallel subgraph matching algorithm. It sets up the query and data graphs, 
 * initializes thread-related variables, and prepares the internal state for processing.
 * 
 * @param query_graph The query graph to be matched.
 * @param data_graph The data graph in which subgraph matches are to be found.
 * @param max_num_results The maximum number of matches to find.
 * @param print_prep Flag to enable or disable printing of preprocessing results.
 * @param print_enum Flag to enable or disable printing of enumeration results.
 * @param homo Flag indicating whether the matching is homomorphic or isomorphic.
 * @param orders Predefined orders for processing query vertices.
 * @param NUMTHREAD The number of threads to use for parallel processing.
 * @param auto_tuning Flag to enable or disable automatic tuning of thread count.
 * 
 * @details
 * 1. Initializes various data structures, including edge indices (`eidx_`), tree nodes (`treeNode_`),
 *    serialized tree representation (`serialized_tree_`), and predefined orders (`pre_defined_order_`).
 * 2. Sets up thread-related variables, such as `NUMTHREAD` and `BIG_THREAD`, based on the query graph size.
 * 3. Allocates memory for thread-local data structures, including `local_vec_m`, `local_vec_extendable`, 
 *    and `local_vec_visited_local`.
 * 4. Prepares the edge index mapping (`eidx_`) for the query graph.
 * 5. Configures predefined orders and backward neighbors if provided.
 */
Parrllel_SymBi::Parrllel_SymBi(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,
        std::vector<std::vector<uint>> orders, size_t NUMTHREAD, size_t auto_tuning)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo
    )
, eidx_(query_.NumVertices())
, treeNode_(query_.NumVertices())
, q_root_(0u)
, serialized_tree_(query_.NumVertices())
, pre_defined_order_()
, pre_defined_backward_nbr_()

, DCS_(query_.NumEdges() * 2)
, d1(query_.NumVertices())
, d2(query_.NumVertices())
, n1(query_.NumEdges() * 2)
, np1(query_.NumVertices())
, n2(query_.NumEdges() * 2)
, nc2(query_.NumVertices()) // vector
, Q1{}
, Q2{}
, thread_pool(query_.NumVertices())
, NUMTHREAD(NUMTHREAD)
, auto_tuning(auto_tuning)
{

    // this->NUM_THREAD = NUMTHREAD;

    NUM_THREAD = 32;

    vertex_vector.resize(0);

    if(auto_tuning == 1){
        std::cout << "Auto-tuning is enabled." << std::endl;
        if(query_.NumVertices() < 7){
            if(NUMTHREAD > query_.NumVertices()){
                NUMTHREAD = query_.NumVertices();
                std::cout << "NUMTHREAD is set to " << NUMTHREAD << "for auto-tuning" << std::endl;
            }
        }
    }else{
        std::cout << "Auto-tuning is disabled." << std::endl;
    }


    if(query_.NumVertices() > 32){
        BIG_THREAD = query_.NumVertices() + 2;
    }else{
        BIG_THREAD = NUMTHREAD + 2;
    }

    local_vec_extendable = std::vector<std::vector<Parrllel_SymBi::ExtendableVertex>>(BIG_THREAD, std::vector<Parrllel_SymBi::ExtendableVertex>(query_.NumVertices()));
    local_vec_m = std::vector<std::vector<uint>>(BIG_THREAD, std::vector<uint>(query_.NumVertices())); 
    local_vec_visited_local = std::vector<std::vector<bool>>(BIG_THREAD, std::vector<bool>(data_.NumVertices(), false));


    num_counter = 0;

    uint edge_pos = 0;

#ifdef DEBUG
    std::chrono::high_resolution_clock::time_point start2, lstart2;
#endif

    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        eidx_[i].resize(query_.NumVertices());
        auto& q_nbrs = query_.GetNeighbors(i);

        for (uint j = 0; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            eidx_[i][u_other] = edge_pos++;
        }
    }

#ifdef DEBUG
    lstart2 = My_Get_Time();
#endif

    if (orders.size() == query_.NumEdges())
    {   
        pre_defined_order_.resize(query_.NumEdges() * 2);

        for (uint i = 0; i < orders.size(); i++)
        {
            pre_defined_order_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ] = orders[i];
        }
    }
    else if (orders.size() == query_.NumEdges() + 1)
    {
        
        std::vector<std::vector<bool>> is_tree_edge (query_.NumVertices());
        for (uint i = 0; i < query_.NumVertices(); i++)
        {
            is_tree_edge[i].resize(query_.NumVertices(), false);
        }
        
        for (uint i = 0; i < orders[0].size(); i += 2)
        {
            is_tree_edge[orders[0][i]][orders[0][i + 1]] = true;
            is_tree_edge[orders[0][i + 1]][orders[0][i]] = true;
        }
        
        pre_defined_order_.resize(query_.NumEdges() * 2);
        pre_defined_backward_nbr_.resize(query_.NumEdges() * 2);


        for (uint i = 1; i < orders.size(); i++)
        {
            pre_defined_order_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ] = orders[i];

            auto& backwards = pre_defined_backward_nbr_[eidx_
                [std::min(orders[i][0], orders[i][1])]
                [std::max(orders[i][0], orders[i][1])]
            ];
            backwards.resize(query_.NumVertices());

            for (uint j = 1; j < orders[i].size(); j++)
            {
                for (uint k = 0; k < j; k++)
                {
                    if (
                        std::get<2>(query_.GetEdgeLabel(orders[i][j], orders[i][k])) != NOT_EXIST &&
                        is_tree_edge[orders[i][j]][orders[i][k]]
                    ) {
                        backwards[orders[i][j]] = orders[i][k];
                    }
                }
            }
        }
    }
#ifdef DEBUG
    Print_Time("Test time for init: ", lstart2); // Test time for 5: 0ms
#endif
}

void Parrllel_SymBi::Preprocessing()
{
    double start_time, end_time;

    // BuildDAG()
    start_time = omp_get_wtime();  // start
    BuildDAG();  
    end_time = omp_get_wtime();  // end
    std::cout << "Time taken to build DAG: " << end_time - start_time << " seconds." << std::endl; 

    // BuildDCS() 
    start_time = omp_get_wtime();  // start
    BuildDCS();  
    end_time = omp_get_wtime();  // end
    std::cout << "Time taken to build DCS: " << end_time - start_time << " seconds." << std::endl; 

}


/**
 * @brief Builds a Directed Acyclic Graph (DAG) for the query graph.
 * 
 * This function constructs a DAG representation of the query graph by analyzing the frequency
 * of edges and vertices. It determines the root vertex of the query graph and organizes the
 * vertices in a Breadth-First Search (BFS) order to create a spanning tree.
 * 
 * @details
 * 1. Calculates the frequency of each edge in the query graph by comparing it with the data graph.
 * 2. Identifies the root vertex of the query graph based on edge and vertex frequencies.
 * 3. Constructs a spanning tree of the query graph using BFS, storing forward and backward edges.
 * 
 * @note The DAG is used to optimize the subgraph matching process by ensuring a topological order.
 */
void Parrllel_SymBi::BuildDAG()
{
    struct EdgeFreq{
        uint v1; uint v2; uint freq;
        EdgeFreq()
        : v1(NOT_EXIST), v2(NOT_EXIST), freq(NOT_EXIST)
        {}
        EdgeFreq(uint v1_arg, uint v2_arg, uint freq_arg)
        : v1(v1_arg), v2(v2_arg), freq(freq_arg)
        {}
    };

    // get the frequency of each query edge
    std::vector<std::vector<size_t>> edge_freqs(query_.NumVertices());
    for (auto& vec: edge_freqs) {
        vec.resize(query_.NumVertices(), 0ul);
    }
    size_t min_freq = ULONG_MAX;
    std::array<uint, 2> min_vs {};

    for (uint u = 0; u < query_.NumVertices(); u++) // query
    {
        const auto& q_nbrs = query_.GetNeighbors(u);
        const auto& q_nbr_labels = query_.GetNeighborLabels(u);

        for (uint i = 0; i < q_nbrs.size(); i++)
        {
            const auto& u_other = q_nbrs[i];
            const auto& u_other_label = q_nbr_labels[i];
            
            if (u > u_other)
            {
                edge_freqs[u][u_other] = edge_freqs[u_other][u];
                continue;
            }

            for (size_t v = 0; v < data_.NumVertices(); v++) // data
            {
                if (data_.GetVertexLabel(v) != NOT_EXIST)
                {
                    const auto& d_nbrs = data_.GetNeighbors(v);
                    const auto& d_nbr_labels = data_.GetNeighborLabels(v);

                    for (uint j = 0; j < d_nbrs.size(); j++)
                    {
                        const auto& v_other = d_nbrs[j];
                        const auto& v_other_label = d_nbr_labels[j];
                        
                        if (
                            data_.GetVertexLabel(v) == query_.GetVertexLabel(u) &&
                            data_.GetVertexLabel(v_other) == query_.GetVertexLabel(u_other) &&
                            v_other_label == u_other_label
                        ) {
                            edge_freqs[u][u_other]++;
                        }
                    }
                }
            }
            if (edge_freqs[u][u_other] < min_freq)
            {
                    min_freq = edge_freqs[u][u_other];
                    min_vs[0] = u;
                    min_vs[1] = u_other;
            }
        }
    }

    // find vertex with low frequency
    std::array<uint, 2> u_freq {};
    for (size_t v = 0; v < data_.NumVertices(); v++)
    if (data_.GetVertexLabel(v) != NOT_EXIST)
    {
        if (data_.GetVertexLabel(v) == query_.GetVertexLabel(min_vs[0]))
            u_freq[0] ++;
        if (data_.GetVertexLabel(v) == query_.GetVertexLabel(min_vs[1]))
            u_freq[1] ++;
    }
    if (u_freq[0] > u_freq[1])
        q_root_ = min_vs[1];
    else if (u_freq[0] < u_freq[1])
        q_root_ = min_vs[0];
    else if (query_.GetDegree(min_vs[0]) > query_.GetDegree(min_vs[1]))
        q_root_ = min_vs[0];
    else
        q_root_ = min_vs[1];

    // build a spaning tree with the BFS order
    std::vector<bool> visited(query_.NumVertices(), false);
    std::vector<bool> exist_on_tree(query_.NumVertices(), false);
    std::queue<uint> bfs_queue;

    bfs_queue.push(q_root_);
    exist_on_tree[q_root_] = true;
    uint order_pos = 0u;
    serialized_tree_[order_pos++] = q_root_;

    while (!bfs_queue.empty())
    {
        const uint u = bfs_queue.front();
        bfs_queue.pop();
        visited[u] = true;

        const auto& q_nbrs = query_.GetNeighbors(u);
        const auto& q_nbr_labels = query_.GetNeighborLabels(u);

        for (uint j = 0; j < q_nbrs.size(); j++)
        {
            const uint u_other = q_nbrs[j];
            const uint u_other_label = q_nbr_labels[j];
            
            if (!exist_on_tree[u_other])
            {
                bfs_queue.push(u_other);
                exist_on_tree[u_other] = true;
                serialized_tree_[order_pos++] = u_other;
            }
            if (!visited[u_other])
            {
                treeNode_[u].forwards_.push_back(u_other);
                treeNode_[u].forward_labels_.push_back(u_other_label);
            }
            else
            {
                treeNode_[u].backwards_.push_back(u_other);
                treeNode_[u].backward_labels_.push_back(u_other_label);
            }
            treeNode_[u].neighbors_.push_back(u_other);
        }
    }  
    if (print_preprocessing_results_)
    {
        std::cout << "DAG: " << std::endl;
        for (uint i = 0; i < query_.NumVertices(); ++i)
        {
            std::cout << serialized_tree_[i] << ": (backwards: ";
            for (auto j: treeNode_[serialized_tree_[i]].backwards_)
                std::cout << j << " ";

            std::cout << ") (forwards: ";
            for (auto j: treeNode_[serialized_tree_[i]].forwards_)
                std::cout << j << " ";

            std::cout << ")" << std::endl;
        }
        std::cout << std::endl;
    }
}


/**
 * @brief Builds the Dynamic Candidate Sets (DCS) for subgraph matching.
 * 
 * This function constructs a hierarchical index structure that maintains candidate vertex sets
 * for each query vertex. The construction process follows a two-phase approach:
 * 1. Top-down phase: Identifies candidate vertices that have connections to their parent vertices.
 * 2. Bottom-up phase: Refines candidate sets by ensuring connectivity with child vertices.
 * 
 * @details
 * The function first processes vertices in a top-down manner (following the serialized spanning tree)
 * to establish parent-child relationships and populate initial candidate sets. During this phase,
 * it identifies vertices that satisfy the "downward" connectivity constraints (d1).
 * 
 * Then, in the bottom-up phase, it processes vertices from leaves to root to enforce "upward" 
 * connectivity constraints (d2) and further refines the candidate sets.
 * 
 * The resulting DCS index efficiently supports incremental subgraph matching by:
 * - Storing candidate vertices for each query vertex
 * - Maintaining connectivity information between candidates
 * - Tracking validity states (d1, d2) for each candidate
 * 
 * This index enables fast filtering during the enumeration phase and efficient updates
 * when the data graph changes.
 */
void Parrllel_SymBi::BuildDCS()
{

    // top-down
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        const uint query_label = query_.GetVertexLabel(u);
        const auto& q_nbrs = query_.GetNeighbors(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++){
            if (data_.GetVertexLabel(j) == query_label){

                const auto& d_nbrs = data_.GetNeighbors(j);
                const auto& d_elabels = data_.GetNeighborLabels(j);

                // add j to hash maps
                for (uint k = 0; k < q_nbrs.size(); k++)
                {
                    {
                        DCS_[eidx_[u][q_nbrs[k]]][j] = {};
                    }
                }
                
                for (uint k = 0; k < treeNode_[u].backwards_.size(); k++)
                {
                    const uint u_other = treeNode_[u].backwards_[k];
                    const uint u_other_elabel = treeNode_[u].backward_labels_[k];

                    for (uint m = 0; m < d_nbrs.size(); m++)
                    {
                        const uint v_other = d_nbrs[m];
                        const uint v_other_elabel = d_elabels[m];
                        if (
                            query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                            u_other_elabel == v_other_elabel
                        ) {
                            {
                                DCS_[eidx_[u][u_other]][j].push_back(v_other);
                                if (d1[u_other][v_other] > 0)
                                    n1[eidx_[u][u_other]][j] += 1;
                            }
                        }
                    }
                    if (n1[eidx_[u][u_other]][j])
                    {
                        np1[u][j] += 1;
                    }
                }
                if (np1[u][j] == treeNode_[u].backwards_.size())
                {
                    d1[u][j] = 1;
                }
            }
        }
    }

    // bottom-up
    
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[query_.NumVertices() - 1 - i];
        // print number of serial tree
        std::cout << "u: " << u << std::endl;
        const uint query_label = query_.GetVertexLabel(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++)
            if (data_.GetVertexLabel(j) == query_label){
                for (uint k = 0; k < treeNode_[u].forwards_.size(); k++){
                    const uint u_other = treeNode_[u].forwards_[k];
                    const uint u_other_elabel = treeNode_[u].forward_labels_[k];

                    const auto& d_nbrs = data_.GetNeighbors(j);
                    const auto& d_elabels = data_.GetNeighborLabels(j);

                    {
                    for (uint m = 0; m < d_nbrs.size(); m++)
                    {
                        const uint v_other = d_nbrs[m];
                        const uint v_other_elabel = d_elabels[m];
                        if (query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                            u_other_elabel == v_other_elabel) 
                        {
                            DCS_[eidx_[u][u_other]][j].push_back(v_other); // DCS Update
                            if (d2[u_other][v_other] > 0)
                                n2[eidx_[u][u_other]][j] += 1;
                        }
                    }
                    if (n2[eidx_[u][u_other]][j]){
                        nc2[u][j] += 1;
                    }
                    }
                }

                {
                if (nc2[u][j] == treeNode_[u].forwards_.size() && d1[u][j] == 1)
                {
                    d2[u][j] = 1;
                }
                }
            }
    }


    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        for (uint k = 0; k < treeNode_[u].backwards_.size(); k++)
        {
            const uint u_other = treeNode_[u].backwards_[k];
            
            for (auto & [v, d_nbrs]: DCS_[eidx_[u][u_other]])
            {
                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];

                    if (d2[u_other][v_other] > 0) {
                        n2[eidx_[u][u_other]][v] += 1;
                    }
                }
            }
        }
    }

}


/**
 * @brief Performs the initial matching phase of the subgraph matching algorithm.
 * 
 * This function initializes the matching process by iterating over all candidate vertices
 * for the root query vertex (`q_root_`). For each candidate, it attempts to find matches
 * by recursively exploring the search space using the `FindMatches` function.
 * 
 * @details
 * 1. Iterates over all candidate vertices for the root query vertex (`q_root_`).
 * 2. For each candidate vertex:
 *    - Initializes the mapping vector `m` and marks the candidate as visited.
 *    - Initializes the `extendable` vector to track the state of extendable query vertices.
 *    - Updates the `extendable` vector based on the neighbors of the root query vertex.
 *    - Calls the `FindMatches` function to recursively explore matches starting from the root.
 *    - Restores the state for backtracking by unmarking the candidate as visited and resetting the mapping.
 * 
 * @note This function is single-threaded and serves as the entry point for the matching process.
 */
void Parrllel_SymBi::InitialMatching() 
{
    
    for (auto & [v, d_nbrs]: DCS_[eidx_[q_root_][query_.GetNeighbors(q_root_).front()]])
    {
        std::vector<uint> m (query_.NumVertices(), UNMATCHED);
        if (d2[q_root_][v] == false) continue;

        m[q_root_] = v;
        visited_[v] = true;

        std::vector<ExtendableVertex> extendable(query_.NumVertices());
        for (auto& u_other: treeNode_[q_root_].neighbors_)
        {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[q_root_][u_other]][v] < extendable[u_other].E)
            {
                extendable[u_other].E = n2[eidx_[q_root_][u_other]][v];
                extendable[u_other].u_min = q_root_;
            }
            extendable[u_other].matched_nbrs ++;
        }
        //         vertex m, extendable, Result counting
        FindMatches(1, m, extendable, num_initial_results_);

        visited_[v] = false;
#ifdef DEBUG
        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            for(size_t j = 0; j < local_vec_visited_local[i].size(); j++){
                std::cout << local_vec_visited_local[i][j] << " ";
            }
        }
#endif
        m[q_root_] = UNMATCHED;
    }
}

// #################### RUNTIME METHODS #################### // 

/**
 * @brief Handles the top-down insertion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the addition of an edge 
 * between two vertices in the data graph in a top-down manner. It increments the edge count 
 * and updates the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_c The child query vertex of `u` in the spanning tree.
 * @param v The data vertex mapped to `u`.
 * @param v_c The data vertex mapped to `u_c`.
 * 
 * @details
 * 1. Checks if the edge count `n1` for the edge between `u_c` and `u` mapped to `v_c` is zero.
 *    - If zero, increments the backward edge count `np1` for `u_c` and `v_c`.
 * 2. If `np1[u_c][v_c]` equals the number of backward edges from `u_c` to its parents:
 *    - Activates `d1[u_c][v_c]` (indicating no data conflict for backward edges).
 *    - Pushes the pair `(u_c, v_c)` into the queue `Q1`.
 * 3. If `nc2[u_c][v_c]` equals the number of forward edges from `u_c` to its children:
 *    - Activates `d2[u_c][v_c]` (indicating no data conflict for forward edges).
 *    - Pushes the pair `(u_c, v_c)` into the queue `Q2`.
 * 4. Increments the edge count `n1` for the edge between `u_c` and `u` mapped to `v_c`.
 */
void Parrllel_SymBi::InsertionTopDown(uint u, uint u_c, uint v, uint v_c)
{
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        np1[u_c][v_c] += 1;
        if (np1[u_c][v_c] == treeNode_[u_c].backwards_.size())
        {
            d1[u_c][v_c] = 1;
            //  It constructs the element in place, meaning it directly constructs the element within the container 
            // without the need for a temporary object, which can be more efficient than using functions like push_back or push
            Q1.emplace(u_c, v_c);
            if (nc2[u_c][v_c] == treeNode_[u_c].forwards_.size())
            {
                d2[u_c][v_c] = 1;
                Q2.emplace(u_c, v_c);
            }
        }
    }
    n1[eidx_[u_c][u]][v_c] += 1;
}


/**
 * @brief Handles the parallel top-down insertion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the addition of an edge 
 * between two vertices in the data graph in a top-down manner. It increments the edge count 
 * and updates the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_c The child query vertex of `u` in the spanning tree.
 * @param v_c The data vertex mapped to `u_c`.
 * @param Q1_para A reference to the queue used for processing updates in the top-down direction.
 * @param Q2_para A reference to the queue used for processing updates in the bottom-up direction.
 * 
 * @details
 * 1. Checks if the edge count `n1` for the edge between `u_c` and `u` mapped to `v_c` is zero.
 *    - If zero, increments the backward edge count `np1` for `u_c` and `v_c`.
 * 2. If `np1[u_c][v_c]` equals the number of backward edges from `u_c` to its parents:
 *    - Activates `d1[u_c][v_c]` (indicating no data conflict for backward edges).
 *    - Pushes the pair `(u_c, v_c)` into the queue `Q1_para`.
 * 3. If `nc2[u_c][v_c]` equals the number of forward edges from `u_c` to its children:
 *    - Activates `d2[u_c][v_c]` (indicating no data conflict for forward edges).
 *    - Pushes the pair `(u_c, v_c)` into the queue `Q2_para`.
 * 4. Increments the edge count `n1` for the edge between `u_c` and `u` mapped to `v_c`.
 */
void Parrllel_SymBi::InsertionTopDown_para(uint u, uint u_c,  uint v_c, 
std::queue<std::pair<uint, uint>> &Q1_para, std::queue<std::pair<uint, uint>> &Q2_para)
{
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        // conflict
        np1[u_c][v_c] += 1;
        if (np1[u_c][v_c] == treeNode_[u_c].backwards_.size())
        {
            // no data conflict
            d1[u_c][v_c] = 1;
            Q1_para.emplace(u_c, v_c);
            if (nc2[u_c][v_c] == treeNode_[u_c].forwards_.size())
            {
                //  no data conflict
                d2[u_c][v_c] = 1;
                Q2_para.emplace(u_c, v_c);
            }
        }
    }
    // conflict
    n1[eidx_[u_c][u]][v_c] += 1;
}

/**
 * @brief Handles the bottom-up insertion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the addition of an edge 
 * between two vertices in the data graph in a bottom-up manner. It increments the edge count 
 * and updates the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_p The parent query vertex of `u` in the spanning tree.
 * @param v The data vertex mapped to `u`.
 * @param v_p The data vertex mapped to `u_p`.
 * 
 * @details
 * 1. Checks if the edge count `n2` for the edge between `u_p` and `u` mapped to `v_p` is zero.
 *    - If zero, increments the forward edge count `nc2` for `u_p` and `v_p`.
 * 2. If `d1[u_p][v_p]` is active and `nc2[u_p][v_p]` equals the number of forward edges from `u_p` to its children:
 *    - Activates `d2[u_p][v_p]`.
 *    - Pushes the pair `(u_p, v_p)` into the queue `Q2`.
 * 3. Increments the edge count `n2` for the edge between `u_p` and `u` mapped to `v_p`.
 */
void Parrllel_SymBi::InsertionBottomUp(uint u, uint u_p, uint v, uint v_p)
{
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        nc2[u_p][v_p] += 1;
        if (d1[u_p][v_p] && nc2[u_p][v_p] == treeNode_[u_p].forwards_.size())
        {
            d2[u_p][v_p] = 1;
            Q2.emplace(u_p, v_p);
        }
    }
    n2[eidx_[u_p][u]][v_p] += 1;
}


/**
 * @brief Handles the parallel bottom-up insertion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the addition of an edge 
 * between two vertices in the data graph in a bottom-up manner. It increments the edge count 
 * and updates the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_p The parent query vertex of `u` in the spanning tree.
 * @param v_p The data vertex mapped to `u_p`.
 * @param Q2_para A reference to the queue used for processing updates in the bottom-up direction.
 * 
 * @details
 * 1. Increments the edge count `nc2` for the edge between `u_p` and `u` mapped to `v_p`.
 * 2. If the edge count becomes non-zero, checks if `d1[u_p][v_p]` is active and if `nc2[u_p][v_p]` 
 *    equals the number of forward edges from `u_p` to its children.
 *    - If both conditions are met, activates `d2[u_p][v_p]` and pushes the pair `(u_p, v_p)` into the queue `Q2_para`.
 */
void Parrllel_SymBi::InsertionBottomUp_para(uint u, uint u_p, uint v_p, std::queue<std::pair<uint, uint>> &Q2_para)
{
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        nc2[u_p][v_p] += 1;
#ifdef DEBUG
        std::cout  << nc2[u_p][v_p] << std::endl;
#endif
        
        if (d1[u_p][v_p] && nc2[u_p][v_p] == treeNode_[u_p].forwards_.size())
        {
            #ifdef DEBUG
            if (d2[u_p][v_p]){
                std::cout << "d2: write again"  << std::endl;
            }
            #endif

            d2[u_p][v_p] = 1; // 没有数据冲突
            Q2_para.emplace(u_p, v_p);

            #ifdef DEBUG
            std::cout << "Q2_para: (" << u_p << ", " << v_p << ")" << std::endl;
            #endif
        }
    }
    n2[eidx_[u_p][u]][v_p] += 1;

    #ifdef DEBUG
    std::cout <<  n2[eidx_[u_p][u]][v_p] << std::endl; 
    #endif
}

/**
 * @brief Handles the top-down deletion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the removal of an edge 
 * between two vertices in the data graph. It decrements the edge count and updates 
 * the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_c The child query vertex of `u` in the spanning tree.
 * @param v The data vertex mapped to `u`.
 * @param v_c The data vertex mapped to `u_c`.
 * 
 * @details
 * 1. Decrements the edge count `n1` for the edge between `u_c` and `u` mapped to `v_c`.
 * 2. If the edge count becomes zero:
 *    - Checks if `d1[u_c][v_c]` is active (value is 1).
 *      - If active, pushes the pair `(u_c, v_c)` into the queue `Q1` and deactivates `d1[u_c][v_c]`.
 *      - If `d2[u_c][v_c]` is also active, pushes the pair `(u_c, v_c)` into the queue `Q2` and deactivates `d2[u_c][v_c]`.
 * 3. Decrements the count `np1[u_c][v_c]` for the number of backward edges from `u_c` to its parents.
 */
void Parrllel_SymBi::DeletionTopDown(uint u, uint u_c, uint v, uint v_c)
{
    n1[eidx_[u_c][u]][v_c] -= 1;
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        if (d1[u_c][v_c] == 1)
        {
            if (d2[u_c][v_c] == 1)
            {
                Q2.emplace(u_c, v_c);
                d2[u_c][v_c] = 0;
            }
            Q1.emplace(u_c, v_c);
            d1[u_c][v_c] = 0;
        }
        np1[u_c][v_c] -= 1;
    }
}


/**
 * @brief Handles the bottom-up deletion of an edge in the data graph.
 * 
 * This function updates the internal data structures to reflect the removal of an edge 
 * between two vertices in the data graph. It decrements the edge count and updates 
 * the corresponding data structures if certain conditions are met.
 * 
 * @param u The current query vertex being processed.
 * @param u_p The parent query vertex of `u` in the spanning tree.
 * @param v The data vertex mapped to `u`.
 * @param v_p The data vertex mapped to `u_p`.
 * 
 * @details
 * 1. Decrements the edge count `n2` for the edge between `u_p` and `u` mapped to `v_p`.
 * 2. If the edge count becomes zero, checks if `d2[u_p][v_p]` is active (value is 1).
 *    - If active, pushes the pair `(u_p, v_p)` into the queue `Q2` and deactivates `d2[u_p][v_p]`.
 * 3. Decrements the count `nc2[u_p][v_p]` for the number of forward edges from `u_p` to its children.
 */
void Parrllel_SymBi::DeletionBottomUp(uint u, uint u_p, uint v, uint v_p)
{
    n2[eidx_[u_p][u]][v_p] -= 1;
    if (n2[eidx_[u_p][u]][v_p] == 0)
    {
        if (d2[u_p][v_p] == 1)
        {
            Q2.emplace(u_p, v_p);
            d2[u_p][v_p] = 0;
        }
        nc2[u_p][v_p] -= 1;
    }
}


//  Main work function for single thread
/**
 * @brief Recursively finds subgraph matches between the query graph and the data graph.
 * 
 * This function implements a depth-first search algorithm to find subgraph matches. It iteratively
 * extends the current partial match by mapping query vertices to data vertices, while ensuring
 * that all constraints are satisfied.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. The function selects a query vertex to extend the match. It prioritizes non-isolated vertices
 *    (vertices with matched neighbors). If no such vertex exists, it selects an isolated vertex.
 * 2. For the selected vertex, it iterates over its candidate data vertices and performs the following checks:
 *    - Index constraints: Ensures the candidate satisfies the index constraints.
 *    - Joinability: Ensures the candidate can connect to already matched neighbors.
 *    - Visitation: Ensures the candidate has not been visited before (if homomorphism is disabled).
 * 3. If all checks pass, the candidate is added to the match, and the function recurses to the next depth.
 * 4. If a complete match is found, the result counter is incremented.
 * 5. The function backtracks by removing the candidate from the match and restoring the previous state.
 * 
 * @note
 * - This is a single-threaded implementation. For large graphs, consider using the parallel version.
 * - The function uses backtracking to explore all possible matches.
 */
void Parrllel_SymBi::FindMatches(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results)
{
    
    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // check if a extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    // if no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local


    #ifdef DEBUG
    bool candidate_empty = true;
    #endif

    // enumerate each neighbor of m[u_min]
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
#ifdef DEBUG
        num_intermediate_results_before_index_check_++;
#endif
        if (d2[u][v] == 0) continue;

#ifdef DEBUG
        num_intermediate_results_after_index_check_++;
#endif

        // 2. check if joinable
        bool joinable = true;

        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }

#ifdef DEBUG
        printf("i_cnt: %d\n", i_cnt);
#endif

        if (!joinable) continue;

#ifdef DEBUG
        num_intermediate_results_after_joinability_check_++;
        candidate_empty = false;
#endif

        // 

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        // num_intermediate_results_after_visit_check_++;

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true; // imp:
        std::vector<ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
            {
                temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs ++;
        }

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;

        #ifdef DEBUG
            std::cout << num_results << std::endl;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;

                std::cout << "num_results by haibin: " << num_results << std::endl;
            }
        #endif

        }
        else
        {
            #ifdef DEBUG
            size_t num_results_before_recursion = num_results;
            #endif

            FindMatches(depth + 1, m, temp_extendable, num_results);
            #ifdef DEBUG            
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
            #endif
        }

        visited_[v] = false;
        m[u] = UNMATCHED;

        #ifdef LIMIT
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
        #endif
    }
    #ifdef DEBUG
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
    #endif
}

/**
 * @brief Parallel version of the FindMatches function.
 * 
 * This function is a parallelized version of the FindMatches function, designed to find subgraph matches
 * between the query graph and the data graph. It uses multiple threads to explore different branches of
 * the search space simultaneously.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 */
void Parrllel_SymBi::Parallel_FindMatches(uint depth, std::vector<uint>& m, 
        std::vector<ExtendableVertex>& extendable, size_t &num_results)
{
    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

#ifdef DEBUG
    std::cout << depth << std::endl;
#endif

    uint u;

// struct ExtendableVertex
// This structure is used to track the extension state of each query vertex 
// during the subgraph matching process.
// It contains information about the number of extendable edges, the minimum 
// extendable vertex, and the number of matched neighbors.

// struct ExtendableVertex {
//     uint E; // Number of extendable edges: Represents the number of vertices in 
// the data graph that the current query vertex can extend to.
//     uint u_min; // Minimum extendable vertex: Indicates the query vertex with the 
// fewest extendable edges that the current vertex can extend to.
//     uint matched_nbrs; // Number of matched neighbors: Represents the number of 
// neighbors of the current query vertex that have already been matched to vertices in the data graph.

//     // Constructor: Initializes the structure with default values.
//     // E is set to UINT_MAX to indicate no extendable edges initially.
//     // u_min is set to NOT_EXIST to indicate no minimum extendable vertex initially.
//     // matched_nbrs is set to 0 as no neighbors are matched initially.
//     ExtendableVertex() : E(UINT_MAX), u_min(NOT_EXIST), matched_nbrs(0) {}
// };



    // init: extendable.size() == query_.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    // If no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]
    // bool candidate_empty = true;

#ifdef DEBUG
    std::cout << DCS_[eidx_[u_min][u]][m[u_min]].size() << std::endl;
#endif

    size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();


    // change the thread number
#ifdef DEBUG
    auto start_nano = std::chrono::high_resolution_clock::now();
#endif

    if (total_size < 3){
        
        for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
        {
            process_vertex(u, u_min, v_idx, m, extendable, num_results, depth);
        }
    }
    else if(total_size < NUM_THREAD){
        

        std::vector<size_t> local_num_results(total_size, 0);

            // Use OpenMP to parallelize the loop over all candidate vertices.
// The number of threads is set to `total_size`, which corresponds to the number of candidate vertices
            #pragma omp parallel for num_threads(total_size)
            for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {
                //  local_num_result for each thread
                size_t thread_id = omp_get_thread_num();
                
                // std::vector<uint> local_m(m);
                // std::vector<ExtendableVertex> local_extendable(extendable);
                // std::vector<bool> local_visited_local(visited_local);

                local_vec_m[thread_id] = m;
                local_vec_extendable[thread_id] = extendable;
                // local_vec_visited_local[thread_id] = visited_local;
                
                // process
                process_vertex(
                    u, u_min, v_idx, 
                    local_vec_m[thread_id], 
                    local_vec_extendable[thread_id], 
                    local_num_results[thread_id], 
                    // local_vec_visited_local[thread_id], 
                    depth);

            }

            for (size_t i = 0; i < local_num_results.size(); ++i) {
                num_results += local_num_results[i];
            }
        
    } else { 

            // OpenMP
            
            std::vector<size_t> local_num_results(NUM_THREAD, 0);
            
            #pragma omp parallel for num_threads(NUM_THREAD)  
            for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {
                // local_num_result
                size_t thread_id = omp_get_thread_num();

                // copy m, extendable, visited_local 
                
                // std::vector<uint> local_m(m);
                // std::vector<ExtendableVertex> local_extendable(extendable);
                // std::vector<bool> local_visited_local(visited_local);

                local_vec_m[thread_id] = m;
                local_vec_extendable[thread_id] = extendable;
                // local_vec_visited_local[thread_id] = visited_local;
                
                // process
                process_vertex(
                    u, u_min, v_idx, 
                    local_vec_m[thread_id], 
                    local_vec_extendable[thread_id], 
                    local_num_results[thread_id], 
                    // local_vec_visited_local[thread_id], 
                    depth);
                // process_vertex(u, u_min, v_idx, local_m, 
                //     local_extendable, local_num_results[thread_id], 
                //     local_visited_local, depth);

                
            }
           // Combine the results from all threads
            for (size_t i = 0; i < local_num_results.size(); ++i) {
                num_results += local_num_results[i];
            }

            // num_counter++;

    }

    #ifdef DEBUG
    Print_Time_Nano("My_Duration on OpenMP: ", start_nano);
    #endif

    

}
/**
 * @brief Parallel version of the FindMatches function for a specific vertex.
 * 
 * This function is a parallelized version of the FindMatches function, designed to find subgraph matches
 * between the query graph and the data graph for a specific vertex. It uses multiple threads to explore
 * different branches of the search space simultaneously.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 */
void Parrllel_SymBi::Parallel_FindMatches3(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{
    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    // extendable.size() == query_.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    // If no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();


    // opt 1: single thread
    if (total_size < 3){
        
        for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
        {
            size_t thread_id = omp_get_thread_num();
            process_vertex_layer1(u, u_min, v_idx, m, extendable, num_results, depth,thread_id);
        }
    }
    // opt 2: small size small thread
    else if(total_size < NUM_THREAD){

        std::vector<size_t> local_num_results(total_size, 0);
        
        #pragma omp parallel for num_threads(total_size)
        for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {

            size_t thread_id = omp_get_thread_num();

            local_vec_m[thread_id] = m;
            local_vec_extendable[thread_id] = extendable;
            // local_vec_visited_local[thread_id] = visited_local;
            
            process_vertex_layer1(
                u, u_min, v_idx, 
                local_vec_m[thread_id], 
                local_vec_extendable[thread_id], 
                local_num_results[thread_id], 
                depth,
                thread_id
            );

        }

        for (size_t i = 0; i < local_num_results.size(); ++i) {
            num_results += local_num_results[i];
        }
        // num_results += std::accumulate(local_num_results.begin(), local_num_results.end(), 0);

    } 
    // opt 3: NUM_THREAD
    else if(total_size < BIG_THREAD){

        

        std::vector<size_t> local_num_results(NUM_THREAD, 0);
        
        #pragma omp parallel for num_threads(NUM_THREAD)  // Incremental Matching: 18576.3ms
        for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {
            //  local_num_result
            size_t thread_id = omp_get_thread_num();

            local_vec_m[thread_id] = m;
            local_vec_extendable[thread_id] = extendable;
            
            // process
            process_vertex_layer1(
                u, u_min, v_idx, 
                local_vec_m[thread_id], 
                local_vec_extendable[thread_id], 
                local_num_results[thread_id], 
                depth,
                thread_id
            );

        }

        for (size_t i = 0; i < local_num_results.size(); ++i) {
            num_results += local_num_results[i];
        }
    }
    // Opt 4: Large Thread
    else
    {

        
        // save only the data

        size_t NUM_BIG_THREAD = 16;

        std::vector<size_t> local_num_results(NUM_BIG_THREAD, 0);
        
        #pragma omp parallel for num_threads(NUM_THREAD)  
        for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {

            size_t thread_id = omp_get_thread_num();

            local_vec_m[thread_id] = m;
            local_vec_extendable[thread_id] = extendable;
            
            auto v = DCS_[eidx_[u_min][u]][local_vec_m[thread_id][u_min]][v_idx];

            // 1. Check index
            if (d2[u][v] == 0) continue;

            // 2. Check if joinable
            bool joinable = true;
            for (auto& u_other: treeNode_[u].neighbors_) {
                if (local_vec_m[thread_id][u_other] == UNMATCHED || u_other == u_min) continue;

                auto it = std::lower_bound(DCS_[eidx_[u_other][u]][local_vec_m[thread_id][u_other]].begin(), DCS_[eidx_[u_other][u]][local_vec_m[thread_id][u_other]].end(), v);
                if (it == DCS_[eidx_[u_other][u]][local_vec_m[thread_id][u_other]].end() || *it != v) {
                    joinable = false;
                    break;
                }
            }

            if (!joinable) continue;

            // 3. Check if visited
            if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

            // 4. Add a vertex mapping
            local_vec_m[thread_id][u] = v;
            local_vec_visited_local[thread_id][v] = true; // imp:

            std::vector<ExtendableVertex> temp_extendable(local_vec_extendable[thread_id]); // local
            for (auto& u_other: treeNode_[u].neighbors_) {
                if (local_vec_m[thread_id][u_other] != UNMATCHED) continue;

                if (n2[eidx_[u][u_other]][local_vec_m[thread_id][u]] < temp_extendable[u_other].E) {
                    temp_extendable[u_other].E = n2[eidx_[u][u_other]][local_vec_m[thread_id][u]];
                    temp_extendable[u_other].u_min = u;
                }
                    temp_extendable[u_other].matched_nbrs++;
            }

            if (depth == query_.NumVertices() - 1) { // match complete
                local_num_results[thread_id]++;
            } else {
                
                uint u2;

                auto depth2 = depth + 1;
                // depth++;

                uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
                uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;    
            
                for (uint i = 0; i < temp_extendable.size(); i++)
                {
                    
                    if (local_vec_m[thread_id][i] != UNMATCHED) continue;
            
                    // __builtin_prefetch(&extendable[i], 1, 0);
                    // Check if an extendable query vertex is isolated or not
                    if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                    {
                        if (temp_extendable[i].E < isolate_minE)
                        {
                            isolate_minE = temp_extendable[i].E; // local
                            isolate_u = i; // local
                        }
                    }
                    else
                    {
                        if (temp_extendable[i].E < non_isolate_minE)
                        {
                            non_isolate_minE = temp_extendable[i].E; // local
                            non_isolate_u = i; // local
                        }
                    }
                }
            
                // If no non-isolated vertex exists, then choose an isolated one
                if (non_isolate_minE == NOT_EXIST)
                    u2 = isolate_u; // local
                else
                    u2 = non_isolate_u; // local
                    
                auto u_min2 = temp_extendable[u2].u_min;
                temp_extendable[u2] = {}; // local
            
            
                // Main parallel for loop
                for (size_t v_idx2 = 0; v_idx2 < DCS_[eidx_[u_min2][u2]][local_vec_m[thread_id][u_min2]].size(); ++v_idx2)
                {
                    process_vertex_layer1(
                        u2, u_min2, v_idx2, 
                        local_vec_m[thread_id], 
                        temp_extendable, 
                        local_num_results[thread_id], 
                        depth2,
                        thread_id
                    );
                    
                }
            }

            // for backtrack
            local_vec_visited_local[thread_id][v] = false;
            local_vec_m[thread_id][u] = UNMATCHED;

        }

        for (size_t i = 0; i < local_num_results.size(); ++i) {
            num_results += local_num_results[i];
        }
    }

    // Print_Time_Nano("My_Duration on OpenMP: ", start_nano);

}



/**
 * @brief Processes a single vertex in the subgraph matching algorithm for layer 1.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * @param thread_id The ID of the thread processing this vertex.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited by the current thread.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Parrllel_SymBi::process_vertex_layer1(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth, size_t thread_id) {

    
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);

        // #ifdef TASK_SPILT
        if(query_.NumVertices() < 7){
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }else{

        if(depth < 5){
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }else{

            // the place for task spilt 




            // if(kongxian){spilt, 逻辑如下:}

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                {
                    if (temp_extendable[i].E < isolate_minE2)
                    {
                        isolate_minE2 = temp_extendable[i].E; // local
                        isolate_u2 = i; // local
                    }
                }
                else
                {
                    if (temp_extendable[i].E < non_isolate_minE2)
                    {
                        non_isolate_minE2 = temp_extendable[i].E; // local
                        non_isolate_u2 = i; // local
                    }
                }
            }

            if (non_isolate_minE2 == NOT_EXIST)
                u2 = isolate_u2; // local
            else
                u2 = non_isolate_u2; // local
                
            uint u_min2 = temp_extendable[u2].u_min;
            temp_extendable[u2] = {}; // local
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            // std::cout << "total_size2: " << total_size2 << std::endl;
            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                // process_vertex(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                auto m3 = m;
                // queue version
                // vertex_queue.emplace(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                
                // vertex_vector.emplace_back(u2, u_min2, v_idx2, m3, temp_extendable,  depth2, v);
                // job
                job_queue.push(std::make_tuple(u2, u_min2, v_idx2, m3, temp_extendable, depth2, v));

            
            }

                    // this is for backtrack
            // local_vec_visited_local[thread_id][v] = false;
            // for(int i=0; i<local_vec_visited_local.size(); i++){
            //     local_vec_visited_local[i][v] = false;
            // }
            // m[u] = UNMATCHED;
        }
    }



    }
// #else
// Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);}

    // #endif
   

    // }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}






inline bool Parrllel_SymBi::process_vertex_layer_local(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth, size_t thread_id) {

    
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);

        #ifdef TASK_SPILT

        if(thread_id != 0){
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }else{

            // the place for task spilt 




            // if(kongxian){spilt, 逻辑如下:}

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                {
                    if (temp_extendable[i].E < isolate_minE2)
                    {
                        isolate_minE2 = temp_extendable[i].E; // local
                        isolate_u2 = i; // local
                    }
                }
                else
                {
                    if (temp_extendable[i].E < non_isolate_minE2)
                    {
                        non_isolate_minE2 = temp_extendable[i].E; // local
                        non_isolate_u2 = i; // local
                    }
                }
            }

            if (non_isolate_minE2 == NOT_EXIST)
                u2 = isolate_u2; // local
            else
                u2 = non_isolate_u2; // local
                
            uint u_min2 = temp_extendable[u2].u_min;
            temp_extendable[u2] = {}; // local
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            // std::cout << "total_size2: " << total_size2 << std::endl;
            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                // process_vertex(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                auto m3 = m;
                // queue version
                // vertex_queue.emplace(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                
                // vertex_vector.emplace_back(u2, u_min2, v_idx2, m3, temp_extendable,  depth2, v);
                // job
                job_queue.push(std::make_tuple(u2, u_min2, v_idx2, m3, temp_extendable, depth2, v));

            
            }
        }



    }
#else
Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);}

    #endif
        // Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);

    // }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}





/**
 * @brief Processes a single vertex in the subgraph matching algorithm.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Parrllel_SymBi::process_vertex(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    size_t thread_id = omp_get_thread_num();
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        Parallel_FindMatches_local2(depth + 1, m, temp_extendable, num_results);
    }

    // 这部分是准备backtrack而生的
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}



/**
 * @brief Processes a single vertex in the subgraph matching algorithm.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Parrllel_SymBi::process_vertex_visit(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth,  std::vector<bool> &visited_local) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    // size_t thread_id = omp_get_thread_num();
    if (!homomorphism_ && visited_local[v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    visited_local[v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        Parallel_FindMatches_local_MMM(depth + 1, m, temp_extendable, num_results, visited_local);
    }

    // for backtrack
    visited_local[v] = false;
    m[u] = UNMATCHED;

    return true;
}



/**
 * @brief Processes a single vertex in the subgraph matching algorithm.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 */
inline bool Parrllel_SymBi::process_vertex2(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth) {

    
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    size_t thread_id = omp_get_thread_num();
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    visited_[v] = true; // imp:
    // std::cout << "u: " << u << " v: " << v << std::endl;

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        FindMatches(depth + 1, m, temp_extendable, num_results);
    }

    // 这部分是准备backtrack而生的
    visited_[v] = false;
    m[u] = UNMATCHED;

    return true;
}



/**
 * @brief Parallel version of the FindMatches function with local visitation tracking.
 * 
 * This function is a parallelized version of the FindMatches function, designed to find subgraph matches
 * between the query graph and the data graph. It uses multiple threads to explore different branches of
 * the search space simultaneously while maintaining a local visitation state for each thread.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param visited_local A vector of booleans tracking the visitation state of data vertices for the current thread.
 * 
 * @details
 * 1. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 2. Iterates over candidate data vertices for the selected query vertex.
 * 3. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 4. If a candidate vertex is valid, updates the match and recursively explores further matches.
 * 5. Performs backtracking to restore the state for the next candidate.
 */
inline void Parrllel_SymBi::Parallel_FindMatches_local2(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{

    size_t thread_id = omp_get_thread_num();

uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

uint u;

for (uint i = 0; i < extendable.size(); i++)
{
    
    if (m[i] != UNMATCHED) continue;

    // __builtin_prefetch(&extendable[i], 1, 0);
    // Check if an extendable query vertex is isolated or not
    if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
    {
        if (extendable[i].E < isolate_minE)
        {
            isolate_minE = extendable[i].E; // local
            isolate_u = i; // local
        }
    }
    else
    {
        if (extendable[i].E < non_isolate_minE)
        {
            non_isolate_minE = extendable[i].E; // local
            non_isolate_u = i; // local
        }
    }
}

// If no non-isolated vertex exists, then choose an isolated one
if (non_isolate_minE == NOT_EXIST)
    u = isolate_u; // local
else
    u = non_isolate_u; // local
    
uint u_min = extendable[u].u_min;
extendable[u] = {}; // local

// Enumerate each neighbor of m[u_min]

// Main parallel for loop
for (auto& v_idx: DCS_[eidx_[u_min][u]][m[u_min]])

{
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) continue;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (
            it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
            *it != v
        ) {
            joinable = false;
            break;
        }
    }

    if (!joinable) continue;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

    // 4. Add a vertex mapping
    {
        m[u] = v;
        local_vec_visited_local[thread_id][v] = true; // imp:
    }

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
        {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
        temp_extendable[u_other].matched_nbrs ++;
    }

    if (depth == query_.NumVertices() - 1) // match complete
    {
        num_results++;
    }
    else
    {
        Parallel_FindMatches_local2(depth + 1, m, temp_extendable, num_results);
    }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;
    }
}


/**
 * @brief Parallel version of the FindMatches function with local visitation tracking.
 * 
 * This function is a parallelized version of the FindMatches function, designed to find subgraph matches
 * between the query graph and the data graph. It uses multiple threads to explore different branches of
 * the search space simultaneously while maintaining a local visitation state for each thread.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param visited_local A vector of booleans tracking the visitation state of data vertices for the current thread.
 * 
 * @details
 * 1. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 2. Iterates over candidate data vertices for the selected query vertex.
 * 3. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 4. If a candidate vertex is valid, updates the match and recursively explores further matches.
 * 5. Performs backtracking to restore the state for the next candidate.
 */
inline void Parrllel_SymBi::Parallel_FindMatches_local_MMM(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local)
{

    // size_t thread_id = omp_get_thread_num();

uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

uint u;

for (uint i = 0; i < extendable.size(); i++)
{
    
    if (m[i] != UNMATCHED) continue;

    // __builtin_prefetch(&extendable[i], 1, 0);
    // Check if an extendable query vertex is isolated or not
    if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
    {
        if (extendable[i].E < isolate_minE)
        {
            isolate_minE = extendable[i].E; // local
            isolate_u = i; // local
        }
    }
    else
    {
        if (extendable[i].E < non_isolate_minE)
        {
            non_isolate_minE = extendable[i].E; // local
            non_isolate_u = i; // local
        }
    }
}

// If no non-isolated vertex exists, then choose an isolated one
if (non_isolate_minE == NOT_EXIST)
    u = isolate_u; // local
else
    u = non_isolate_u; // local
    
uint u_min = extendable[u].u_min;
extendable[u] = {}; // local

// Enumerate each neighbor of m[u_min]

// Main parallel for loop
for (auto& v_idx: DCS_[eidx_[u_min][u]][m[u_min]])

{
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) continue;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (
            it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
            *it != v
        ) {
            joinable = false;
            break;
        }
    }

    if (!joinable) continue;

    // 3. Check if visited
    if (!homomorphism_ && visited_local[v]) continue;

    // 4. Add a vertex mapping
    {
        m[u] = v;
        visited_local[v] = true; // imp:
    }

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
        {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
        temp_extendable[u_other].matched_nbrs ++;
    }

    if (depth == query_.NumVertices() - 1) // match complete
    {
        num_results++;
    }
    else
    {
        Parallel_FindMatches_local_MMM(depth + 1, m, temp_extendable, num_results, visited_local);
    }

    // for backtrack
    visited_local[v] = false;
    m[u] = UNMATCHED;
    }
}


inline void Parrllel_SymBi::Parallel_FindMatches_local3(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results, size_t thread_id)
{

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    for (uint i = 0; i < extendable.size(); i++)
    {
        
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    // If no non-isolated vertex exists, then choose an isolated one
    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]

    // Main parallel for loop
    for (size_t v_idx = 0; v_idx < DCS_[eidx_[u_min][u]][m[u_min]].size(); ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

        // 1. Check index
        if (d2[u][v] == 0) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        {
            m[u] = v;
            local_vec_visited_local[thread_id][v] = true; // imp:
        }

        std::vector<ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
            {
                temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs ++;
        }

        if (depth == query_.NumVertices() - 1) // match complete
        {
            num_results++;
        }
        else
        {

                Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results, thread_id);
        }

        // for backtrack
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
    }
}



/**
 * @brief Processes a single vertex in the subgraph matching algorithm.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param visited_local A vector of booleans tracking the visitation state of data vertices.
 * @param depth The current recursion depth.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Parrllel_SymBi::process_vertex(uint u, uint u_min, size_t v_idx, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results, 
    std::vector<bool> &visited_local, uint depth) {
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    if (!homomorphism_ && visited_local[v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    visited_local[v] = true; // imp:
    // std::cout << "u: " << u << " v: " << v << std::endl;

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {

        Parallel_FindMatches_local2(depth + 1, m, temp_extendable, num_results, visited_local);
    }

    visited_local[v] = false;
    m[u] = UNMATCHED;

    return true;
}




/**
 * @brief Parallel version of the FindMatches function using an explicit stack for recursion.
 * 
 * This function implements a parallelized subgraph matching algorithm using an explicit stack
 * to simulate recursion. It avoids the overhead of recursive function calls and allows for
 * better control over the matching process.
 * 
 * @param initial_depth The initial recursion depth, representing the number of matched query vertices.
 * @param initial_m A vector representing the initial mapping of query vertices to data vertices.
 *                  `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param initial_extendable A vector containing the initial state of extendable query vertices,
 *                           including the number of extendable edges, the minimum extendable vertex,
 *                           and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param initial_visited A vector of booleans tracking the visitation state of data vertices.
 * 
 * @details
 * 1. Uses a `StackFrame` structure to store the state of each recursive call, including the depth,
 *    current mapping, extendable vertices, visited vertices, and loop variables.
 * 2. Simulates recursion by pushing and popping frames from a stack.
 * 3. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 4. Iterates over candidate data vertices for the selected query vertex.
 * 5. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 6. If a candidate vertex is valid, updates the match and pushes a new frame onto the stack.
 * 7. Performs backtracking by restoring the state from the stack.
 */
inline void Parrllel_SymBi::Parallel_FindMatches_local2(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results, std::vector<bool> &visited_local)
{
uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

uint u;

for (uint i = 0; i < extendable.size(); i++)
{
    
    if (m[i] != UNMATCHED) continue;

    __builtin_prefetch(&extendable[i], 1, 0);
    // Check if an extendable query vertex is isolated or not
    if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
    {
        if (extendable[i].E < isolate_minE)
        {
            isolate_minE = extendable[i].E; // local
            isolate_u = i; // local
        }
    }
    else
    {
        if (extendable[i].E < non_isolate_minE)
        {
            non_isolate_minE = extendable[i].E; // local
            non_isolate_u = i; // local
        }
    }
}

// If no non-isolated vertex exists, then choose an isolated one
if (non_isolate_minE == NOT_EXIST)
    u = isolate_u; // local
else
    u = non_isolate_u; // local
    
uint u_min = extendable[u].u_min;
extendable[u] = {}; // local


// Enumerate each neighbor of m[u_min]
// Main parallel for loop
for (size_t v_idx = 0; v_idx < DCS_[eidx_[u_min][u]][m[u_min]].size(); ++v_idx)
{
    __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) continue;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (
            it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
            *it != v
        ) {
            joinable = false;
            break;
        }
    }

    if (!joinable) continue;

    // 3. Check if visited
    if (!homomorphism_ && visited_local[v]) continue;

    // 4. Add a vertex mapping
    {
        m[u] = v;
        visited_local[v] = true; // imp:
    }

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
        {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
        temp_extendable[u_other].matched_nbrs ++;
    }

    if (depth == query_.NumVertices() - 1) // match complete
    {
        num_results++;
    }
    else
    {
        Parallel_FindMatches_local2(depth + 1, m, temp_extendable, num_results, visited_local);
    }

    visited_local[v] = false;
    m[u] = UNMATCHED;
    }
}


/**
 * @brief Parallel version of the FindMatches function using an explicit stack for recursion.
 * 
 * This function implements a parallelized subgraph matching algorithm using an explicit stack
 * to simulate recursion. It avoids the overhead of recursive function calls and allows for
 * better control over the matching process.
 * 
 * @param initial_depth The initial recursion depth, representing the number of matched query vertices.
 * @param initial_m A vector representing the initial mapping of query vertices to data vertices.
 *                  `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param initial_extendable A vector containing the initial state of extendable query vertices,
 *                           including the number of extendable edges, the minimum extendable vertex,
 *                           and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param initial_visited A vector of booleans tracking the visitation state of data vertices.
 * 
 * @details
 * 1. Uses a `StackFrame` structure to store the state of each recursive call, including the depth,
 *    current mapping, extendable vertices, visited vertices, and loop variables.
 * 2. Simulates recursion by pushing and popping frames from a stack.
 * 3. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 4. Iterates over candidate data vertices for the selected query vertex.
 * 5. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 6. If a candidate vertex is valid, updates the match and pushes a new frame onto the stack.
 * 7. Performs backtracking by restoring the state from the stack.
 */
void Parrllel_SymBi::FindMatches2(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{
if (reach_time_limit) return;

uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

uint u;
for (uint i = 0; i < extendable.size(); i++)
{
    if (m[i] != UNMATCHED) continue;

    // check if a extendable query vertex is isolated or not
    if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
    {
        if (extendable[i].E < isolate_minE)
        {
            isolate_minE = extendable[i].E;
            isolate_u = i;
        }
    }
    else
    {
        if (extendable[i].E < non_isolate_minE)
        {
            non_isolate_minE = extendable[i].E;
            non_isolate_u = i;
        }
    }
}

// if no non-isolated vertex exists, then choose an isolated one
if (non_isolate_minE == NOT_EXIST)
    u = isolate_u;
else
    u = non_isolate_u;
uint u_min = extendable[u].u_min;
extendable[u] = {};

// enumerate each neighbor of m[u_min]
bool candidate_empty = true;

for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
{
    // 1. check index
    num_intermediate_results_before_index_check_++;
    if (d2[u][v] == 0) continue;
    num_intermediate_results_after_index_check_++;

    // 2. check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (
            it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
            *it != v
        ) {
            joinable = false;
            break;
        }
    }
    if (!joinable) continue;
    num_intermediate_results_after_joinability_check_++;

    candidate_empty = false;

    // 3. check if visited
    if (!homomorphism_ && visited_[v]) continue;
    num_intermediate_results_after_visit_check_++;

    // 4. add a vertex mapping
    m[u] = v;
    visited_[v] = true;
    std::vector<ExtendableVertex> temp_extendable(extendable);
    for (auto& u_other: treeNode_[u].neighbors_)
    {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E)
        {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
        temp_extendable[u_other].matched_nbrs ++;
    }

    if (depth == query_.NumVertices() - 1)
    {
        num_results++;
        if (print_enumeration_results_)
        {
            for (auto j: m)
                std::cout << j << " ";
            std::cout << std::endl;
        }
    }
    else
    {
        size_t num_results_before_recursion = num_results;
        FindMatches2(depth + 1, m, temp_extendable, num_results);
        if (num_results == num_results_before_recursion)
        {
            num_intermediate_results_without_results_++;
        }
    }

    visited_[v] = false;
    m[u] = UNMATCHED;
    if (num_results >= max_num_results_) return;
    if (reach_time_limit) return;
}
if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


void Parrllel_SymBi::FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{
    // if (reach_time_limit) return;

    uint u = pre_defined_order_[order_index][depth];

    uint u_min = NOT_EXIST;
    if (pre_defined_backward_nbr_.empty())
    {

        int u_min_size = std::numeric_limits<int>::max();  

        // #pragma omp parallel for reduction(min: u_min_size) reduction(min: u_min)
        for (size_t i = 0; i < treeNode_[u].neighbors_.size(); ++i)
        {
            auto& u_other = treeNode_[u].neighbors_[i]; // var
            if (m[u_other] == UNMATCHED) continue;

            int size = DCS_[eidx_[u_other][u]][m[u_other]].size();
            if (size < u_min_size)
            {
                u_min_size = size;
                u_min = u_other;
            }
        }
        // #pragma omp barrier // . wait
    }
    else
    {
        u_min = pre_defined_backward_nbr_[order_index][u];
    }

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;

    // DCS: std::vector<std::unordered_map<uint, std::vector<uint>>> 
    //                      ~var  ~var  
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;
        
        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;
            }
        }
        else
        {
            size_t num_results_before_recursion = num_results;
            FindMatches(order_index, depth + 1, m, num_results);
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


/**
 * @brief Parallel version of the FindMatches function using an explicit stack for recursion.
 * 
 * This function implements a parallelized subgraph matching algorithm using an explicit stack
 * to simulate recursion. It avoids the overhead of recursive function calls and allows for
 * better control over the matching process.
 * 
 * @param initial_depth The initial recursion depth, representing the number of matched query vertices.
 * @param initial_m A vector representing the initial mapping of query vertices to data vertices.
 *                  `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param initial_extendable A vector containing the initial state of extendable query vertices,
 *                           including the number of extendable edges, the minimum extendable vertex,
 *                           and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param initial_visited A vector of booleans tracking the visitation state of data vertices.
 * 
 * @details
 * 1. Uses a `StackFrame` structure to store the state of each recursive call, including the depth,
 *    current mapping, extendable vertices, visited vertices, and loop variables.
 * 2. Simulates recursion by pushing and popping frames from a stack.
 * 3. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 4. Iterates over candidate data vertices for the selected query vertex.
 * 5. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 6. If a candidate vertex is valid, updates the match and pushes a new frame onto the stack.
 * 7. Performs backtracking by restoring the state from the stack.
 */
void Parrllel_SymBi::Parallel_FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = pre_defined_order_[order_index][depth];

    uint u_min = NOT_EXIST;
    if (pre_defined_backward_nbr_.empty())
    {

        int u_min_size = std::numeric_limits<int>::max();  //

        // #pragma omp parallel for reduction(min: u_min_size) reduction(min: u_min)
        for (size_t i = 0; i < treeNode_[u].neighbors_.size(); ++i)
        {
            auto& u_other = treeNode_[u].neighbors_[i]; // var
            if (m[u_other] == UNMATCHED) continue;

            int size = DCS_[eidx_[u_other][u]][m[u_other]].size();
            if (size < u_min_size)
            {
                u_min_size = size;
                u_min = u_other;
            }
        }
        // #pragma omp barrier // . wait
    }
    else
    {
        u_min = pre_defined_backward_nbr_[order_index][u];
    }

    // enumerate each neighbor of m[u_min]
    bool candidate_empty = true;

    // DCS: std::vector<std::unordered_map<uint, std::vector<uint>>> 
    //                      ~var  ~var  
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        num_intermediate_results_before_index_check_++;
        if (d2[u][v] == 0) continue;
        num_intermediate_results_after_index_check_++;

        // 2. check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_)
        {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (
                it == DCS_[eidx_[u_other][u]][m[u_other]].end() ||
                *it != v
            ) {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;
        

        // 4. add a vertex mapping
        m[u] = v;
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
            if (print_enumeration_results_)
            {
                for (auto j: m)
                    std::cout << j << " ";
                std::cout << std::endl;
            }
        }
        else
        {
            size_t num_results_before_recursion = num_results;
            Parallel_FindMatches(order_index, depth + 1, m, num_results);
            if (num_results == num_results_before_recursion)
            {
                num_intermediate_results_without_results_++;
            }
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
    if (candidate_empty) num_intermediate_results_with_empty_candidate_set_++;
}


inline void Parrllel_SymBi::Update_Graph(Graph& data_, int v1, int v2, int label) {

    data_.AddEdge(v1, v2, label);

}



void Parrllel_SymBi::AddEdgeAsync(uint v1, uint v2, uint label)
{
    // Graph Update
    std::thread t1(&Parrllel_SymBi::Update_Graph, this, std::ref(data_), v1, v2, label);
    t1.join();
}


/**
 * @brief Classifies whether an edge in the data graph matches a query edge.
 * 
 * This function checks if the edge (`v1 -> v2`) in the data graph matches any edge in the query graph
 * based on vertex labels and edge labels. It also ensures that the edge follows the directed acyclic
 * graph (DAG) structure of the query graph.
 * 
 * @param v1 The source vertex of the edge in the data graph.
 * @param v2 The destination vertex of the edge in the data graph.
 * @param label The label of the edge in the data graph.
 * 
 * @return `true` if the edge matches a query edge, `false` otherwise.
 * 
 * @details
 * 1. Iterates over all query vertices (`u1`) to find a vertex whose label matches `v1`.
 * 2. For each matching query vertex (`u1`), iterates over all other query vertices (`u2`) to find a vertex
 *    whose label matches `v2`.
 * 3. Checks if the edge label between `u1` and `u2` matches the label of the edge (`v1 -> v2`).
 * 4. Ensures that the edge follows the DAG structure of the query graph by checking the direction of the edge.
 * 5. Returns `true` if all conditions are satisfied; otherwise, returns `false`.
 */
bool Parrllel_SymBi::Classify(uint v1, uint v2, uint label){
    auto v1_label = data_.GetVertexLabel(v1);

    for (uint u1 = 0; u1 < query_.NumVertices(); u1++){ // take one u1 from query
        //  Get v1 that label matched to query label
        if (v1_label == query_.GetVertexLabel(u1)){

            // Get v2 that label matched to query vertex
            for (uint u2 = 0; u2 < query_.NumVertices(); u2++){ // take one u2 from query
                if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){


                    if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                    
                    // Prune 1: detect if the edge is in the query graph
                    // if (label == std::get<2>(query_.GetEdgeLabel(u1, u2))){return false;}                    
                    
                    bool reversed = false;

                    if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                    {
                        std::swap(u1, u2); // 
                        std::swap(v1, v2); //
                        reversed = true;
                    }
                    if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1)
                         != treeNode_[u2].backwards_.end())
                    {
                        // Prune 2: detect if the edge is in the DCS graph
                        // for more complex query graph
                        // detect if DCS has update. If so, that means new FM will occur
                        return false;

                        // Prune 3: detect if the edge is in the FM path
                        // for more complex query graph
                        // If so, that means new FM will occur
                        bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

                        if(old_p_d1 || old_p_d2 || old_c_d2){
                            return false;
                        }


                    }
                    if (reversed)
                    {
                        std::swap(u1, u2);
                        std::swap(v1, v2);
                    }
                }
            }
        }
    }

    return true;

}


/**
 * @brief Adds an edge to the data graph and updates the subgraph matching index.
 * 
 * This function handles the addition of an edge to the data graph. It updates the internal
 * data structures to reflect the addition and re-evaluates the subgraph matches affected by
 * the new edge.
 * 
 * @param v1 The source vertex of the edge to be added in the data graph.
 * @param v2 The destination vertex of the edge to be added in the data graph.
 * @param label The label of the edge to be added.
 * 
 * @details
 * 1. Updates the data graph by adding the edge (`v1 -> v2`) with the specified label.
 * 2. Iterates over all query edges that match the added edge (`v1 -> v2`).
 * 3. For each matching query edge:
 *    - Updates the subgraph matching index by inserting the edge into the corresponding data structures.
 *    - Updates the mapping and visitation state for the affected vertices.
 *    - Recomputes the extendable state for the affected query vertices.
 *    - Calls the `Parallel_FindMatches` function to re-evaluate the matches.
 * 4. Handles backtracking to restore the state for further processing.
 */
void Parrllel_SymBi::AddEdge_Single(uint v1, uint v2, uint label)
{
        // enumerate all query edges that matches v1 --> v2
        {
            // auto start_nano = std::chrono::high_resolution_clock::now();
            std::queue<std::pair<uint, uint>> Q1_para, Q2_para;

            auto v1_label = data_.GetVertexLabel(v1);

   
            for (uint u1 = 0; u1 < query_.NumVertices(); u1++){ // take one u1 from query
                //  Get v1 that label matched to query label
                if (v1_label == query_.GetVertexLabel(u1)){
                    // Get v2 that label matched to query vertex
                    for (uint u2 = 0; u2 < query_.NumVertices(); u2++) // take one u2 from query
                        if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){
                            if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                            
                            bool reversed = false;
                            if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                            {
                                std::swap(u1, u2); 
                                std::swap(v1, v2); 
                                reversed = true;
                            }
                            if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
                            {
                               
                                auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
                      
                                DCS_[eidx_[u1][u2]][v1].insert(it, v2);

                                it = std::lower_bound(DCS_[eidx_[u2][u1]][v2].begin(), DCS_[eidx_[u2][u1]][v2].end(), v1);
                                DCS_[eidx_[u2][u1]][v2].insert(it, v1);

                                bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

                                if (old_p_d1)
                                    InsertionTopDown_para(u1, u2,  v2, Q1_para, Q2_para);

                                if (old_c_d2)
                                    InsertionBottomUp_para(u2, u1, v1, Q2_para);                               
                                if (old_p_d2)
                                    n2[eidx_[u2][u1]][v2] += 1;

                                while (!Q1_para.empty())
                                {
                                    auto [u_queue, v_queue] = Q1_para.front();
                                    Q1_para.pop();
                                    // #pragma omp for
                                    for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                                        for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                                        {
                                            InsertionTopDown_para(u_queue, u_c_queue,  v_c_queue, Q1_para, Q2_para);
                                        }
                                }

                                while (!Q2_para.empty())
                                {
                                    auto [u_queue, v_queue] = Q2_para.front();
                                    Q2_para.pop();
                                    auto size = treeNode_[u_queue].backwards_.size(); 

                                    for (size_t i = 0; i < size; ++i) {
                 
                                        auto& u_p_queue = treeNode_[u_queue].backwards_[i];
                                        auto& eidx_val = eidx_[u_queue][u_p_queue];  

                                        for (size_t j = 0; j < DCS_[eidx_val][v_queue].size(); ++j) {
                                            auto v_p_queue = DCS_[eidx_val][v_queue][j];

                                            InsertionBottomUp_para(u_queue, u_p_queue, v_p_queue, Q2_para);
                                        }
                                    }
                                    
                                    for (size_t i = 0; i < treeNode_[u_queue].forwards_.size(); ++i) {
                                        auto& u_c_queue = treeNode_[u_queue].forwards_[i];
                                        for (size_t j = 0; j < DCS_[eidx_[u_queue][u_c_queue]][v_queue].size(); ++j) {
                                            auto& v_c_queue = DCS_[eidx_[u_queue][u_c_queue]][v_queue][j];

                                            n2[eidx_[u_c_queue][u_queue]][v_c_queue] += 1;
                                        }
                                    }
                                }
                            }

                            // thread N end, join 
                            if (reversed)
                            {
                                std::swap(u1, u2);
                                std::swap(v1, v2);
                            }
                        }
                }

                //  Print_Time_Nano("My_Duration on first loop ", start_nano);
            }


            // enumerate all query edges that matches v1 --> v2
            size_t num_results = 0ul;

                {
                    // #pragma omp parallel for
                    for (uint u1 = 0; u1 < query_.NumVertices(); u1++){
                        for (uint u2 = 0; u2 < query_.NumVertices(); u2++){
                            if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2) 
                                && data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
                            {
                                if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                                    
                                std::vector<uint> m(query_.NumVertices(), UNMATCHED);
                        
                                bool reversed = false;
                                if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                                {
                                    std::swap(u1, u2);
                                    std::swap(v1, v2);
                                    reversed = true; // continue;
                                }
                                if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end()
                                    && d2[u1][v1] == 1 && d2[u2][v2] == 1) // not filted
                                {
                                    m[u1] = v1; // begin search
                                    m[u2] = v2;
                                        
                                    // std::vector<bool> visited_local(data_.NumVertices(), false);
                                    visited_[v1] = true;
                                    visited_[v2] = true;
        
                                    std::vector<ExtendableVertex> extendable(query_.NumVertices()); // 6
                                    for (auto u: {u1, u2})  // all the neighbors of u is e(M(u`),v) belongs to E(g) 
                                    {
                                        for (auto& u_other: treeNode_[u].neighbors_)
                                        {
                                            if (m[u_other] != UNMATCHED) continue;
        
                                            if (n2[eidx_[u][u_other]][m[u]] < extendable[u_other].E)
                                            {
                                                extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                                                extendable[u_other].u_min = u;
                                            }
                                            extendable[u_other].matched_nbrs ++;
                                        }
                                    }
        
                                    // size_t local_num_results = 0; 
        
                                    if (pre_defined_order_.empty()){
                                        Parallel_FindMatches4(2, m, extendable, num_results); // 1s
                                    }
                                    else{
                                        // FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, 
                                        // num_results);
                                        Parallel_FindMatches4(2, m, extendable, num_results); // 1s
                                    }
        
                                    visited_[v1] = false;
                                    visited_[v2] = false;

                                    m[u1] = UNMATCHED;
                                    m[u2] = UNMATCHED;

                                }
                                if (reversed)
                                {
                                    std::swap(u1, u2);
                                    std::swap(v1, v2);
                                }
                            }
                        }
                        // }

                    } // openmp stop
        }
        // END_ENUMERATION:
        num_positive_results_ += num_results;
    }
} // end of section




/**
 * @brief Adds an edge to the data graph and updates the subgraph matching index.
 * 
 * This function handles the addition of an edge to the data graph. It updates the internal
 * data structures to reflect the addition and re-evaluates the subgraph matches affected by
 * the new edge.
 * 
 * @param v1 The source vertex of the edge to be added in the data graph.
 * @param v2 The destination vertex of the edge to be added in the data graph.
 * @param label The label of the edge to be added.
 * 
 * @details
 * 1. Updates the data graph by adding the edge (`v1 -> v2`) with the specified label.
 * 2. Iterates over all query edges that match the added edge (`v1 -> v2`).
 * 3. For each matching query edge:
 *    - Updates the subgraph matching index by inserting the edge into the corresponding data structures.
 *    - Updates the mapping and visitation state for the affected vertices.
 *    - Recomputes the extendable state for the affected query vertices.
 *    - Calls the `Parallel_FindMatches` function to re-evaluate the matches.
 * 4. Handles backtracking to restore the state for further processing.
 */
void Parrllel_SymBi::AddEdge(uint v1, uint v2, uint label)
{

    // data_.AddEdge(v1, v2, label); // add edge to data graph
 
        // enumerate all query edges that matches v1 --> v2

        // #pragma omp task
        {
            // auto start_nano = std::chrono::high_resolution_clock::now();

            std::queue<std::pair<uint, uint>> Q1_para, Q2_para;

            auto v1_label = data_.GetVertexLabel(v1);

            for (uint u1 = 0; u1 < query_.NumVertices(); u1++){ // take one u1 from query

            // record time:
            // one loop for 150 time with 300ns

                //  Get v1 that label matched to query label
                if (v1_label == query_.GetVertexLabel(u1)){

                    // Get v2 that label matched to query vertex
                    for (uint u2 = 0; u2 < query_.NumVertices(); u2++) // take one u2 from query
                        if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){

                            #ifdef DEBUG
                            std::cout << "(v1,v2):" << v1 << " " << v2 << " (u1,u2): " << u1 << " "<< u2 << std::endl;
                            #endif

                            if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                            
                            bool reversed = false;

                            if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                            {
                                std::swap(u1, u2); // u1, u2 conflict
                                std::swap(v1, v2); // v1, v2 follows the DAG order
                                reversed = true;
                            }
                            if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
                            {

                                auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);

                                DCS_[eidx_[u1][u2]][v1].insert(it, v2);

                                it = std::lower_bound(DCS_[eidx_[u2][u1]][v2].begin(), DCS_[eidx_[u2][u1]][v2].end(), v1);
                                DCS_[eidx_[u2][u1]][v2].insert(it, v1);

                                bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

                                if (old_p_d1)
                                    InsertionTopDown_para(u1, u2,  v2, Q1_para, Q2_para);

                                if (old_c_d2)
                                    InsertionBottomUp_para(u2, u1, v1, Q2_para);                               
                                if (old_p_d2)
                                    n2[eidx_[u2][u1]][v2] += 1;

                                // std::cout << Q1_para.size() << std::endl;
                                // std::cout << Q2_para.size() << std::endl;

                                while (!Q1_para.empty())
                                {
                                    auto [u_queue, v_queue] = Q1_para.front();
                                    Q1_para.pop();
                                    // #pragma omp for
                                    for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                                        for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                                        {
                                            InsertionTopDown_para(u_queue, u_c_queue,  v_c_queue, Q1_para, Q2_para);
                                            // print Q2 queue
                                            // std::cout << "Q1_para: (" << u_queue << ", " << v_queue << ")" << std::endl; 
                                            // if (reach_time_limit) return;
                                        }
                                }

                                // thread N end, join
                                // int size = Q2.size();

                                while (!Q2_para.empty())
                                {
                                    auto [u_queue, v_queue] = Q2_para.front();
                                    Q2_para.pop();

                                    auto size = treeNode_[u_queue].backwards_.size(); // from 1 to 2

                                    for (size_t i = 0; i < size; ++i) {
                                        // processing the index
                                        auto& u_p_queue = treeNode_[u_queue].backwards_[i];
                                        auto& eidx_val = eidx_[u_queue][u_p_queue];  // eidx only for u_queue  u_p_queue

                                        // inner
                                        for (size_t j = 0; j < DCS_[eidx_val][v_queue].size(); ++j) {
                                            auto v_p_queue = DCS_[eidx_val][v_queue][j];
                                            // insert
                                            InsertionBottomUp_para(u_queue, u_p_queue, v_p_queue, Q2_para);
                                        }
                                    }
                                    
                                    // #pragma omp parallel for 
                                    for (size_t i = 0; i < treeNode_[u_queue].forwards_.size(); ++i) {
                                        auto& u_c_queue = treeNode_[u_queue].forwards_[i];
                                        for (size_t j = 0; j < DCS_[eidx_[u_queue][u_c_queue]][v_queue].size(); ++j) {
                                            auto& v_c_queue = DCS_[eidx_[u_queue][u_c_queue]][v_queue][j];
                                            // After testing, we don't need atomic here 我们没有数据冲突
                                            // #pragma omp atomic
                                            n2[eidx_[u_c_queue][u_queue]][v_c_queue] += 1;
                                        }
                                    }

                                }
                            }

                            // thread N end, join 
                            if (reversed)
                            {
                                std::swap(u1, u2);
                                std::swap(v1, v2);
                            }
                        }
                }

                //  Print_Time_Nano("My_Duration on first loop ", start_nano);
            }

            // enumerate all query edges that matches v1 --> v2
            size_t num_results = 0ul;

            // #pragma omp parallel 
            {
                // #pragma omp single
                {
                    
                    {
                    // #pragma omp parallel for
                    for (uint u1 = 0; u1 < query_.NumVertices(); u1++){

                        // std::cout << omp_get_thread_num() << std::endl;
                        // std::this_thread::sleep_for(std::chrono::milliseconds(400)); 

                        // std::cout <<  u1 << std::endl;
                        // if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
                        {
                            for (uint u2 = 0; u2 < query_.NumVertices(); u2++){
                                if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2) 
                                    && data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
                                {

                                    // std::cout << omp_get_thread_num();
                                    // if(omp_get_thread_num() == 0){
                                    //     std::cout << std::endl;
                                    // }
                                    #ifdef DEBUG
                                    std::this_thread::sleep_for(std::chrono::milliseconds(300));
                                    #endif
        
                                    if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                                    
                                    // auto thread_id = omp_get_thread_num();
                                    
                                    std::vector<uint> m(query_.NumVertices(), UNMATCHED);
                        
                                    bool reversed = false;
                                    if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                                    {
                                        std::swap(u1, u2);
                                        std::swap(v1, v2);
                                        reversed = true; // continue;
                                    }
                                    if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end()
                                        && d2[u1][v1] == 1 && d2[u2][v2] == 1) // not filted
                                    {
                                        m[u1] = v1; // begin search
                                        m[u2] = v2;
                                        
                                        // std::vector<bool> visited_local(data_.NumVertices(), false);
                                        // visited_[v1] = true;
                                        // local_vec_visited_local[0][v1] = true;
                                        // visited_[v2] = true;
                                        // local_vec_visited_local[0][v2] = true;
                                        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                                            local_vec_visited_local[i][v1] = true;
                                            local_vec_visited_local[i][v2] = true; // wcl
                                        }
        
                                        std::vector<ExtendableVertex> extendable(query_.NumVertices()); // 6
                                        for (auto u: {u1, u2})  // all the neighbors of u is e(M(u`),v) belongs to E(g) 
                                        {
                                            for (auto& u_other: treeNode_[u].neighbors_)
                                            {
                                                if (m[u_other] != UNMATCHED) continue;
        
                                                if (n2[eidx_[u][u_other]][m[u]] < extendable[u_other].E)
                                                {
                                                    extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                                                    extendable[u_other].u_min = u;
                                                }
                                                extendable[u_other].matched_nbrs ++;
                                            }
                                        }
    
        
                                        if (pre_defined_order_.empty()){
                                            
                                            #ifdef DEBUG
                                            std::cout << "Parallel" << std::endl;
                                            auto time_start = std::chrono::high_resolution_clock::now();
                                            #endif

                                            // if(auto_tuning){
                                                if((query_.NumVertices() <=7 ) && (query_.NumEdges() > 9) ){
                                                    FindMatches(2, m, extendable, num_results);
                                                }else{
                                                    Parallel_FindMatches4(2, m, extendable, num_results); //
                                                }


                                            #ifdef DEBUG
                                            Print_Time_Nano("My_Duration on Matching: ", time_start);
                                            #endif
                                        }
                                        else{
                                            // #ifdef DEBUG
                                            // std::cout << "Sequential" << std::endl;
                                            FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, 
                                            num_results);
                                            // #endif
                                            // Parallel_FindMatches4(2, m, extendable, num_results); //
                                        }

    
                                        
                                        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                                            local_vec_visited_local[i][v1] = false;
                                            local_vec_visited_local[i][v2] = false; // wcl
                                        }
        

                                        m[u1] = UNMATCHED;
                                        m[u2] = UNMATCHED;
                                        //  Check time
                                        // if (num_results >= max_num_results_) goto END_ENUMERATION;
                                        // if (reach_time_limit) goto END_ENUMERATION;
        
                                        // Print_Time_Nano("My_Duration on Matching: ", start_nano);
                                    }
                                    if (reversed)
                                    {
                                        std::swap(u1, u2);
                                        std::swap(v1, v2);
                                    }
                                }
                            }
                        }

                    } // openmp stop


                }
            }
        }

            // END_ENUMERATION:
            num_positive_results_ += num_results;

            // std::cout << "Parallel_end" << std::endl; // ok

        }

        
        
    } // end of section


/**
 * @brief Removes an edge from the data graph and updates the subgraph matching index.
 * 
 * This function handles the removal of an edge from the data graph. It updates the internal
 * data structures to reflect the removal and re-evaluates the subgraph matches affected by
 * the removed edge.
 * 
 * @param v1 The source vertex of the edge to be removed in the data graph.
 * @param v2 The destination vertex of the edge to be removed in the data graph.
 * 
 * @details
 * 1. Iterates over all query edges that match the removed edge (`v1 -> v2`).
 * 2. For each matching query edge:
 *    - Checks if the edge exists in the subgraph matching index.
 *    - Updates the mapping and visitation state for the affected vertices.
 *    - Recomputes the extendable state for the affected query vertices.
 *    - Calls the `Parallel_FindMatches_delete` function to re-evaluate the matches.
 * 3. Updates the subgraph matching index to remove the edge.
 * 4. Handles backtracking to restore the state for further processing.
 * 
 * @note This function is designed to handle dynamic updates to the data graph.
 */
void Parrllel_SymBi::RemoveEdge(uint v1, uint v2)
{
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    size_t num_results = 0ul;
    if (max_num_results_ == 0) goto END_ENUMERATION;

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
        for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
        if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
        {
            auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
            if (
                it == DCS_[eidx_[u1][u2]][v1].end() ||
                *it != v2
            ) continue;
            
            bool reversed = false;
            if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
            {
                std::swap(u1, u2);
                std::swap(v1, v2);
                reversed = true;
            }
            if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end()
                && d2[u1][v1] == 1 && d2[u2][v2] == 1)
            {
                m[u1] = v1;
                m[u2] = v2;
                visited_[v1] = true;
                visited_[v2] = true;

                std::vector<ExtendableVertex> extendable(query_.NumVertices());
                for (auto u: {u1, u2})
                {
                    for (auto& u_other: treeNode_[u].neighbors_)
                    {
                        if (m[u_other] != UNMATCHED) continue;

                        if (n2[eidx_[u][u_other]][m[u]] < extendable[u_other].E)
                        {
                            extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                            extendable[u_other].u_min = u;
                        }
                        extendable[u_other].matched_nbrs ++;
                    }
                }
                if (pre_defined_order_.empty())
                    // FindMatches(2, m, extendable, num_results);
                    Parallel_FindMatches_delete(2, m, extendable, num_results);
                else
                    // FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, 
                    //         num_results);
                    Parallel_FindMatches_delete(2, m, extendable, num_results);

                visited_[v1] = false;
                visited_[v2] = false;
                m[u1] = UNMATCHED;
                m[u2] = UNMATCHED;
                if (num_results >= max_num_results_) goto END_ENUMERATION;
                if (reach_time_limit) return;
            }
            if (reversed)
            {
                std::swap(u1, u2);
                std::swap(v1, v2);
            }
        }
    }

    END_ENUMERATION:
    num_negative_results_ += num_results;

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
        if (
            it == DCS_[eidx_[u1][u2]][v1].end() ||
            *it != v2
        ) continue;
        
        bool reversed = false;
        // if <u2,v2> is a parent of <u1,v1>, then we need to swap u1 and u2
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
            reversed = true;
        }
        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
        {
            auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
            DCS_[eidx_[u1][u2]][v1].erase(it);
            it = std::lower_bound(DCS_[eidx_[u2][u1]][v2].begin(), DCS_[eidx_[u2][u1]][v2].end(), v1);
            DCS_[eidx_[u2][u1]][v2].erase(it);

            bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

            if (old_c_d2)
                DeletionBottomUp(u2, u1, v2, v1);
            if (old_p_d2)
                n2[eidx_[u2][u1]][v2] -= 1;
            if (old_p_d1)
                DeletionTopDown(u1, u2, v1, v2);

            while (!Q1.empty())
            {
                auto [u_queue, v_queue] = Q1.front();
                Q1.pop();
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                    for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                    {
                        DeletionTopDown(u_queue, u_c_queue, v_queue, v_c_queue);
                        if (reach_time_limit) return;
                    }
            }

            while (!Q2.empty())
            {
                auto [u_queue, v_queue] = Q2.front();
                Q2.pop();
                for (auto& u_p_queue : treeNode_[u_queue].backwards_)
                    for (auto& v_p_queue : DCS_[eidx_[u_queue][u_p_queue]][v_queue])
                    {
                        DeletionBottomUp(u_queue, u_p_queue, v_queue, v_p_queue);
                        if (reach_time_limit) return;
                    }
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                    for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                    {
                        n2[eidx_[u_c_queue][u_queue]][v_c_queue] -= 1;
                        if (reach_time_limit) return;
                    }
            }
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }

    data_.RemoveEdge(v1, v2);
}



// Add and remove vertex don't need to do much work

void Parrllel_SymBi::AddVertex(uint id, uint label)
{
    #pragma omp critical
    {
    data_.AddVertex(id, label);

    visited_.resize(id + 1, false);

    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        if (data_.GetVertexLabel(id) == query_.GetVertexLabel(u))
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]][id] = {};
            }
            if (treeNode_[u].backwards_.empty())
                d1[u][id] = 1;
            if (d1[u][id] == 1 && treeNode_[u].forwards_.empty())
                d2[u][id] = 1;
        }
    }
    }
}

void Parrllel_SymBi::RemoveVertex(uint id)
{
    #pragma omp critical
    {
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        if (data_.GetVertexLabel(id) == query_.GetVertexLabel(u))
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]].erase(id);
                if (n1[eidx_[u][q_nbrs[k]]].find(id) != n1[eidx_[u][q_nbrs[k]]].end())
                    n1[eidx_[u][q_nbrs[k]]].erase(id);
                if (n2[eidx_[u][q_nbrs[k]]].find(id) != n2[eidx_[u][q_nbrs[k]]].end())
                    n2[eidx_[u][q_nbrs[k]]].erase(id);
            }
            if (np1[u].find(id) != np1[u].end())
                np1[u].erase(id);
            if (nc2[u].find(id) != nc2[u].end())
                nc2[u].erase(id);
            if (d1[u].find(id) != d1[u].end())
                d1[u].erase(id);
            if (d2[u].find(id) != d2[u].end())
                d2[u].erase(id);
        }
    }

    data_.RemoveVertex(id);
    }
}


/**
 * @brief Parallel version of the FindMatches function using an explicit stack for recursion.
 * 
 * This function implements a parallelized subgraph matching algorithm using an explicit stack
 * to simulate recursion. It avoids the overhead of recursive function calls and allows for
 * better control over the matching process.
 * 
 * @param initial_depth The initial recursion depth, representing the number of matched query vertices.
 * @param initial_m A vector representing the initial mapping of query vertices to data vertices.
 *                  `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param initial_extendable A vector containing the initial state of extendable query vertices,
 *                           including the number of extendable edges, the minimum extendable vertex,
 *                           and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param initial_visited A vector of booleans tracking the visitation state of data vertices.
 * 
 * @details
 * 1. Uses a `StackFrame` structure to store the state of each recursive call, including the depth,
 *    current mapping, extendable vertices, visited vertices, and loop variables.
 * 2. Simulates recursion by pushing and popping frames from a stack.
 * 3. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 4. Iterates over candidate data vertices for the selected query vertex.
 * 5. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 6. If a candidate vertex is valid, updates the match and pushes a new frame onto the stack.
 * 7. Performs backtracking by restoring the state from the stack.
 */
void Parrllel_SymBi::Parallel_FindMatches_local(uint initial_depth, std::vector<uint> initial_m,
    std::vector<ExtendableVertex> initial_extendable, size_t &num_results, std::vector<bool> initial_visited)
{
    struct StackFrame {
        uint depth;
        std::vector<uint> m;
        std::vector<ExtendableVertex> extendable;
        std::vector<bool> visited_local;
        bool is_processing_loop;
        uint u;
        uint u_min;
        size_t v_idx;
        size_t v_count;

        StackFrame(uint d, std::vector<uint> m_vec, std::vector<ExtendableVertex> ext, std::vector<bool> vis,
                   bool ipl, uint u_val, uint umin_val, size_t vidx, size_t vcnt)
            : depth(d), m((m_vec)), extendable(std::move(ext)), visited_local(std::move(vis)),
              is_processing_loop(ipl), u(u_val), u_min(umin_val), v_idx(vidx), v_count(vcnt) {}
    };

    std::stack<StackFrame> stack;
    stack.emplace(initial_depth, initial_m, initial_extendable, initial_visited,
                  false, NOT_EXIST, NOT_EXIST, 0, 0);

    while (!stack.empty()) {
        auto frame = stack.top();
        stack.pop();

        // uint u = NOT_EXIST;

        if (!frame.is_processing_loop) {
            // Phase 1: Select u and u_min
            uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
            uint non_isolate_minE = NOT_EXIST;
            uint isolate_minE = NOT_EXIST;

            for (uint i = 0; i < frame.extendable.size(); ++i) {
                if (frame.m[i] != UNMATCHED) continue;

                if (frame.extendable[i].matched_nbrs == query_.GetNeighbors(i).size()) {
                    // Isolated vertex
                    if (frame.extendable[i].E < isolate_minE) {
                        isolate_minE = frame.extendable[i].E;
                        isolate_u = i;
                    }
                } else {
                    // Non-isolated vertex
                    if (frame.extendable[i].E < non_isolate_minE) {
                        non_isolate_minE = frame.extendable[i].E;
                        non_isolate_u = i;
                    }
                }
            }

            // Determine u
            frame.u = (non_isolate_minE != NOT_EXIST) ? non_isolate_u : isolate_u;

            if (frame.u == NOT_EXIST) {
                continue; // No extendable vertex
            }

            frame.u_min = frame.extendable[frame.u].u_min;
            frame.v_count = DCS_[eidx_[frame.u_min][frame.u]][frame.m[frame.u_min]].size();
            frame.v_idx = 0;
            frame.is_processing_loop = true;

            stack.push(frame);
        } else {
            // Phase 2: Process v_idx loop
            for (; frame.v_idx < frame.v_count; ++frame.v_idx) {
                auto v = DCS_[eidx_[frame.u_min][frame.u]][frame.m[frame.u_min]][frame.v_idx];

                // Check d2 condition
                if (d2[frame.u][v] == 0) continue;

                // Check joinable condition
                bool joinable = true;
                for (auto& u_other : treeNode_[frame.u].neighbors_) {
                    if (frame.m[u_other] == UNMATCHED || u_other == frame.u_min) continue;

                    auto& dcs_list = DCS_[eidx_[u_other][frame.u]][frame.m[u_other]];
                    auto it = std::lower_bound(dcs_list.begin(), dcs_list.end(), v);
                    if (it == dcs_list.end() || *it != v) {
                        joinable = false;
                        break;
                    }
                }
                if (!joinable) continue;

                // Check visited condition
                if (!homomorphism_ && frame.visited_local[v]) continue;

                // Prepare new state
                auto new_m = frame.m;
                new_m[frame.u] = v;
                auto new_visited = frame.visited_local;
                new_visited[v] = true;

                // Update extendable vertices
                auto temp_extendable = frame.extendable;
                temp_extendable[frame.u] = {}; // Reset current u

                for (auto& u_other : treeNode_[frame.u].neighbors_) {
                    if (new_m[u_other] != UNMATCHED) continue;

                    uint eidx = eidx_[frame.u][u_other];
                    if (n2[eidx][v] < temp_extendable[u_other].E) {
                        temp_extendable[u_other].E = n2[eidx][v];
                        temp_extendable[u_other].u_min = frame.u;
                    }
                    temp_extendable[u_other].matched_nbrs++;
                }

                // Check termination condition
                if (frame.depth == query_.NumVertices() - 1) {
                    num_results++;
                    // std::cout << num_results << std::endl;
                } else {
                    // Push new stack frame for recursion
                    // Parallel_FindMatches_local2(frame.depth + 1,frame.m,temp_extendable,num_results,new_visited);
                    stack.emplace(frame.depth + 1, new_m, temp_extendable, new_visited,
                                  false, NOT_EXIST, NOT_EXIST, 0, 0);
                }

                // Push current frame back with incremented v_idx
                frame.v_idx++;
                stack.push(frame);

                break;
            }
        }
    }
}



/**
 * @brief Parallel version of the FindMatches function using an explicit stack for recursion.
 * 
 * This function implements a parallelized subgraph matching algorithm using an explicit stack
 * to simulate recursion. It avoids the overhead of recursive function calls and allows for
 * better control over the matching process.
 * 
 * @param initial_depth The initial recursion depth, representing the number of matched query vertices.
 * @param initial_m A vector representing the initial mapping of query vertices to data vertices.
 *                  `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param initial_extendable A vector containing the initial state of extendable query vertices,
 *                           including the number of extendable edges, the minimum extendable vertex,
 *                           and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * @param initial_visited A vector of booleans tracking the visitation state of data vertices.
 * 
 * @details
 * 1. Uses a `StackFrame` structure to store the state of each recursive call, including the depth,
 *    current mapping, extendable vertices, visited vertices, and loop variables.
 * 2. Simulates recursion by pushing and popping frames from a stack.
 * 3. Selects the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 4. Iterates over candidate data vertices for the selected query vertex.
 * 5. Checks index constraints, joinability, and visitation status for each candidate vertex.
 * 6. If a candidate vertex is valid, updates the match and pushes a new frame onto the stack.
 * 7. Performs backtracking by restoring the state from the stack.
 */
void Parrllel_SymBi::Parallel_FindMatches_local_stack(uint initial_depth, std::vector<uint> initial_m,
    std::vector<ExtendableVertex> initial_extendable, size_t &num_results, std::vector<bool> initial_visited)
{
    struct StackFrame {
        uint depth;
        std::vector<uint> m;
        std::vector<ExtendableVertex> extendable;
        std::vector<bool> visited_local;
        bool is_processing_loop;
        uint u;
        uint u_min;
        size_t v_idx;
        size_t v_count;

        StackFrame(uint d, std::vector<uint> m_vec, std::vector<ExtendableVertex> ext, std::vector<bool> vis,
                   bool ipl, uint u_val, uint umin_val, size_t vidx, size_t vcnt)
            : depth(d), m((m_vec)), extendable(std::move(ext)), visited_local(std::move(vis)),
              is_processing_loop(ipl), u(u_val), u_min(umin_val), v_idx(vidx), v_count(vcnt) {}
    };

    std::stack<StackFrame> stack;
    stack.emplace(initial_depth, initial_m, initial_extendable, initial_visited,
                  false, NOT_EXIST, NOT_EXIST, 0, 0);

    while (!stack.empty()) {
        auto frame = stack.top();
        stack.pop();

        if (!frame.is_processing_loop) {
            // Phase 1: Select u and u_min
            uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
            uint non_isolate_minE = NOT_EXIST;
            uint isolate_minE = NOT_EXIST;

            for (uint i = 0; i < frame.extendable.size(); ++i) {
                if (frame.m[i] != UNMATCHED) continue;

                if (frame.extendable[i].matched_nbrs == query_.GetNeighbors(i).size()) {
                    // Isolated vertex
                    if (frame.extendable[i].E < isolate_minE) {
                        isolate_minE = frame.extendable[i].E;
                        isolate_u = i;
                    }
                } else {
                    // Non-isolated vertex
                    if (frame.extendable[i].E < non_isolate_minE) {
                        non_isolate_minE = frame.extendable[i].E;
                        non_isolate_u = i;
                    }
                }
            }

            // Determine u
            frame.u = (non_isolate_minE != NOT_EXIST) ? non_isolate_u : isolate_u;

            if (frame.u == NOT_EXIST) {
                continue; // No extendable vertex
            }

            frame.u_min = frame.extendable[frame.u].u_min;
            frame.v_count = DCS_[eidx_[frame.u_min][frame.u]][frame.m[frame.u_min]].size();
            frame.v_idx = 0;
            frame.is_processing_loop = true;

            stack.push(frame);
        } else {
            // Phase 2: Process v_idx loop
            for (; frame.v_idx < frame.v_count; ++frame.v_idx) {
                auto v = DCS_[eidx_[frame.u_min][frame.u]][frame.m[frame.u_min]][frame.v_idx];

                // Check d2 condition
                if (d2[frame.u][v] == 0) continue;

                // Check joinable condition
                bool joinable = true;
                for (auto& u_other : treeNode_[frame.u].neighbors_) {
                    if (frame.m[u_other] == UNMATCHED || u_other == frame.u_min) continue;

                    auto& dcs_list = DCS_[eidx_[u_other][frame.u]][frame.m[u_other]];
                    auto it = std::lower_bound(dcs_list.begin(), dcs_list.end(), v);
                    if (it == dcs_list.end() || *it != v) {
                        joinable = false;
                        break;
                    }
                }
                if (!joinable) continue;

                // Check visited condition
                if (!homomorphism_ && frame.visited_local[v]) continue;

                // Prepare new state
                auto new_m = frame.m;
                new_m[frame.u] = v;
                auto new_visited = frame.visited_local;
                new_visited[v] = true;

                // Update extendable vertices
                auto temp_extendable = frame.extendable;
                temp_extendable[frame.u] = {}; // Reset current u

                for (auto& u_other : treeNode_[frame.u].neighbors_) {
                    if (new_m[u_other] != UNMATCHED) continue;

                    uint eidx = eidx_[frame.u][u_other];
                    if (n2[eidx][v] < temp_extendable[u_other].E) {
                        temp_extendable[u_other].E = n2[eidx][v];
                        temp_extendable[u_other].u_min = frame.u;
                    }
                    temp_extendable[u_other].matched_nbrs++;
                }

                // Check termination condition
                if (frame.depth == query_.NumVertices() - 1) {
                    num_results++;
                    // std::cout << num_results << std::endl;
                } else {
                    // Push new stack frame for recursion
                    stack.emplace(frame.depth + 1, new_m, temp_extendable, new_visited,
                                  false, NOT_EXIST, NOT_EXIST, 0, 0);
                }

                // Push current frame back with incremented v_idx
                frame.v_idx++;
                stack.push(frame);

                break;
            }
        }
    }
}


/**
 * @brief Parallel version of the FindMatches function for handling deletions.
 * 
 * This function is designed to handle subgraph matching in parallel when edges are deleted
 * from the data graph. It clears the vertex and job queues, identifies the next query vertex
 * to process, and uses parallel threads to explore the search space.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. Clears the `vertex_vector` and `job_queue` to reset the state for processing deletions.
 * 2. Identifies the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 3. Iterates over candidate data vertices for the selected query vertex.
 * 4. Uses OpenMP to parallelize the exploration of the search space.
 * 5. Updates the match and recursively explores further matches.
 * 6. Aggregates results from all threads and updates the global match count.
 */
void Parrllel_SymBi::Parallel_FindMatches4(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{

    vertex_vector.clear();
    job_queue.clear();

    size_t NUMT = NUM_THREAD;

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    // extendable.size() == query_.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        // 获取最小邻居的节点
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]
    // bool candidate_empty = true;

    size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();

    // std::vector<std::tuple<uint, uint, size_t, std::vector<uint>,
    //  std::vector<ExtendableVertex>,  uint , uint>> vertex_vector;
    //  std::vector<bool>> 
     
    std::vector<size_t> local_num_result(NUMT,0);

    for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

        // 1. Check index
        if (d2[u][v] == 0) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_) {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), 
                            DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        size_t thread_id = omp_get_thread_num();
        if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        m[u] = v;
        for(int i=0; i<local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = true;
        }
        // local_vec_visited_local[thread_id][v] = true; // imp:

        #ifdef DEBUG
        std::cout << "thread_id: " << thread_id << std::endl;
        std::cout << "u: " << u << " v: " << v << std::endl;
        #endif

        std::vector<ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_[u].neighbors_) {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
                temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                temp_extendable[u_other].u_min = u;
            }
                temp_extendable[u_other].matched_nbrs++;
        }

        if (depth == query_.NumVertices() - 1) { // match complete
            num_results++;
        } else {


            

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                {
                    if (temp_extendable[i].E < isolate_minE2)
                    {
                        isolate_minE2 = temp_extendable[i].E; // local
                        isolate_u2 = i; // local
                    }
                }
                else
                {
                    if (temp_extendable[i].E < non_isolate_minE2)
                    {
                        non_isolate_minE2 = temp_extendable[i].E; // local
                        non_isolate_u2 = i; // local
                    }
                }
            }

            if (non_isolate_minE2 == NOT_EXIST)
                u2 = isolate_u2; // local
            else
                u2 = non_isolate_u2; // local
                
            uint u_min2 = temp_extendable[u2].u_min;
            temp_extendable[u2] = {}; // local
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            // std::cout << "total_size2: " << total_size2 << std::endl;
            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, temp_extendable,  depth2, v);

                // job_queue.push({u2, u_min2, v_idx2, m, temp_extendable,  depth2, v});

            
            }
        }

        // this is for backtrack
        // local_vec_visited_local[thread_id][v] = false;
        for(int i=0; i<local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = false;
        }
        m[u] = UNMATCHED;

    }


    if(vertex_vector.size() < NUMT){
    // if(job_queue.unsafe_size() < NUMT){
        // NUMT = job_queue.unsafe_size();
        // NUMT = vertex_vector.size();
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }

    // std::cout << "NUMT"<< NUMT << std::endl;
// SOLUTION 1
    #pragma omp parallel for num_threads(NUMT) schedule(auto)
    for (size_t i = 0; i < vertex_vector.size() ; ++i) {
        auto [u3, u_min, v_idx, m, extendable,  depth3, v3] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();

        local_vec_visited_local[thread_id][v3] = true;
        m[u3] = v3;
        // process_vertex_layer1(u3, u_min, v_idx, m, extendable, local_num_result[thread_id], depth3, thread_id);
        process_vertex_layer_local(u3, u_min, v_idx, m, extendable, local_num_result[thread_id], depth3, thread_id);

        local_vec_visited_local[thread_id][v3] = false;
        m[u3] = UNMATCHED;

        if(!job_queue.empty() 
        // && (i > vertex_vector.size() - NUMT)
        ){
            std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                 std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto [u, u_min, v_idx, m, extendable,  depth, v] = job;
                local_vec_visited_local[thread_id][v] = true;
                m[u] = v;
                // FindMatches_task(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                process_vertex_layer_local(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                local_vec_visited_local[thread_id][v] = false; 
                m[u] = UNMATCHED;
            }
        }

        // if(){}
    }

    #ifdef TASK_SPILT
    #pragma omp parallel num_threads(NUMT)
    {
        size_t thread_id = omp_get_thread_num();
        std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                   std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
    
        while (true) {
            // take a job from the queue
            if (!job_queue.try_pop(job)) {
                // if empty, go out
                if (job_queue.empty()) {
                    break;
                }
                continue; // find a job, go
            }
    
            // deal with the job
            auto [u, u_min, v_idx, m, extendable, depth, v] = job;
            local_vec_visited_local[thread_id][v] = true;
            m[u] = v;
            process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
            local_vec_visited_local[thread_id][v] = false;
            m[u] = UNMATCHED;
        }
    }
    #endif

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

    // std::cout << "local_num_result.size(): " << num_results << std::endl;

}


/**
 * @brief Parallel version of the FindMatches function for handling deletions.
 * 
 * This function is designed to handle subgraph matching in parallel when edges are deleted
 * from the data graph. It clears the vertex and job queues, identifies the next query vertex
 * to process, and uses parallel threads to explore the search space.
 * 
 * @param depth The current recursion depth, representing the number of matched query vertices.
 * @param m A vector representing the current mapping of query vertices to data vertices. 
 *          `m[i]` indicates the data vertex mapped to query vertex `i`.
 * @param extendable A vector containing information about extendable query vertices, including
 *                   the number of extendable edges, the minimum extendable vertex, and matched neighbors.
 * @param num_results A reference to the counter for the number of matches found.
 * 
 * @details
 * 1. Clears the `vertex_vector` and `job_queue` to reset the state for processing deletions.
 * 2. Identifies the next query vertex to process based on the number of extendable edges and matched neighbors.
 * 3. Iterates over candidate data vertices for the selected query vertex.
 * 4. Uses OpenMP to parallelize the exploration of the search space.
 * 5. Updates the match and recursively explores further matches.
 * 6. Aggregates results from all threads and updates the global match count.
 */
void Parrllel_SymBi::Parallel_FindMatches_delete(uint depth, std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, size_t &num_results)
{

    vertex_vector.clear();
    job_queue.clear();

    size_t NUMT = NUMTHREAD;

    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;

    uint u;

    // extendable.size() == query_.NumVertices()
    for (uint i = 0; i < extendable.size(); i++)
    {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        // 获取最小邻居的节点
        if (extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
        {
            if (extendable[i].E < isolate_minE)
            {
                isolate_minE = extendable[i].E; // local
                isolate_u = i; // local
            }
        }
        else
        {
            if (extendable[i].E < non_isolate_minE)
            {
                non_isolate_minE = extendable[i].E; // local
                non_isolate_u = i; // local
            }
        }
    }

    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u; // local
    else
        u = non_isolate_u; // local
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {}; // local

    // Enumerate each neighbor of m[u_min]
    // bool candidate_empty = true;

    size_t total_size = DCS_[eidx_[u_min][u]][m[u_min]].size();
     
    std::vector<size_t> local_num_result(NUMT,0);


    for (size_t v_idx = 0; v_idx < total_size; ++v_idx)
    {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
        auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

        // 1. Check index
        if (d2[u][v] == 0) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_[u].neighbors_) {
            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), 
                            DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
            if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        size_t thread_id = omp_get_thread_num();
        if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        m[u] = v;
        // for(int i=0; i<local_vec_visited_local.size(); i++){
        //     local_vec_visited_local[i][v] = true;
        // }
        // local_vec_visited_local[thread_id][v] = true; // imp:
        #ifdef DEBUG
        std::cout << "thread_id: " << thread_id << std::endl;
        std::cout << "u: " << u << " v: " << v << std::endl;
        #endif

        std::vector<ExtendableVertex> temp_extendable(extendable); // local
        for (auto& u_other: treeNode_[u].neighbors_) {
            if (m[u_other] != UNMATCHED) continue;

            if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
                temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
                temp_extendable[u_other].u_min = u;
            }
                temp_extendable[u_other].matched_nbrs++;
        }

        if (depth == query_.NumVertices() - 1) { // match complete
            num_results++;
        } else {

            uint depth2 = depth + 1;

            uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
            uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
        
            uint u2;

            for (uint i = 0; i < temp_extendable.size(); i++)
            {
                if (m[i] != UNMATCHED) continue;
        
                // Check if an extendable query vertex is isolated or not
                if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
                {
                    if (temp_extendable[i].E < isolate_minE2)
                    {
                        isolate_minE2 = temp_extendable[i].E; // local
                        isolate_u2 = i; // local
                    }
                }
                else
                {
                    if (temp_extendable[i].E < non_isolate_minE2)
                    {
                        non_isolate_minE2 = temp_extendable[i].E; // local
                        non_isolate_u2 = i; // local
                    }
                }
            }

            if (non_isolate_minE2 == NOT_EXIST)
                u2 = isolate_u2; // local
            else
                u2 = non_isolate_u2; // local
                
            uint u_min2 = temp_extendable[u2].u_min;
            temp_extendable[u2] = {}; // local
        
            // Enumerate each neighbor of m[u_min]
            // bool candidate_empty = true;
        
            size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
            {

                // process_vertex(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);

                // queue version
                // vertex_queue.emplace(u2, u_min2, v_idx2, m, temp_extendable, num_results, depth2);
                
                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, temp_extendable,  depth2, v);

            
            }
        }

        // this is for backtrack
        // local_vec_visited_local[thread_id][v] = false;
        // for(int i=0; i<local_vec_visited_local.size(); i++){
        //     local_vec_visited_local[i][v] = false;
        // }
        m[u] = UNMATCHED;

    }


    if(vertex_vector.size() < NUMT){
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }
    // std::cout << "wtf" << NUMT << std::endl;

    #pragma omp parallel for num_threads(NUMT) schedule(dynamic, 1)
    for (size_t i = 0; i < vertex_vector.size(); ++i) {
        auto [u, u_min, v_idx, m, extendable,  depth, v] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();

        local_vec_visited_local[thread_id][v] = true;
        m[u] = v;
        process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);

        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;

        if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
            std::tuple<unsigned int, unsigned int, unsigned long, std::vector<unsigned int>,
                 std::vector<Parrllel_SymBi::ExtendableVertex>, unsigned int, unsigned int> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto [u, u_min, v_idx, m, extendable,  depth, v] = job;
                local_vec_visited_local[thread_id][v] = true;
                m[u] = v;
                // FindMatches_task(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                process_vertex_layer1(u, u_min, v_idx, m, extendable, local_num_result[thread_id], depth, thread_id);
                local_vec_visited_local[thread_id][v] = false; 
                m[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

}


/**
 * @brief Processes a single vertex in the subgraph matching algorithm for a specific task.
 * 
 * This function checks whether a candidate vertex can be added to the current match
 * by verifying index constraints, joinability, and visitation status. If the vertex
 * is valid, it updates the match and recursively explores further matches.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The query vertex with the fewest extendable edges.
 * @param v_idx The index of the candidate vertex in the data graph.
 * @param m A vector representing the current mapping of query vertices to data vertices.
 * @param extendable A vector containing information about extendable query vertices.
 * @param num_results A reference to the counter for the number of matches found.
 * @param depth The current recursion depth.
 * @param thread_id The ID of the thread processing this vertex.
 * 
 * @return `true` if the vertex was successfully processed, `false` otherwise.
 * 
 * @details
 * 1. Retrieves the candidate vertex `v` from the data graph.
 * 2. Checks if the vertex satisfies the index constraints (`d2[u][v]`).
 * 3. Verifies if the vertex is joinable with already matched neighbors.
 * 4. Ensures the vertex has not been visited by the current thread.
 * 5. If all checks pass, updates the match and recursively explores further matches.
 * 6. Performs backtracking to restore the state for the next candidate.
 */
inline bool Parrllel_SymBi::FindMatches_task(
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth, size_t thread_id) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] == UNMATCHED || u_other == u_min) continue;

        auto it = std::lower_bound(DCS_[eidx_[u_other][u]][m[u_other]].begin(), DCS_[eidx_[u_other][u]][m[u_other]].end(), v);
        if (it == DCS_[eidx_[u_other][u]][m[u_other]].end() || *it != v) {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    std::vector<ExtendableVertex> temp_extendable(extendable); // local
    for (auto& u_other: treeNode_[u].neighbors_) {
        if (m[u_other] != UNMATCHED) continue;

        if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
            temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
            temp_extendable[u_other].u_min = u;
        }
            temp_extendable[u_other].matched_nbrs++;
    }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {

        // vertex_vector.emplace_back(u, u_min, v_idx, m, temp_extendable, depth + 1, v);
        // Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        uint depth2 = depth + 1;

        if(depth2 < (query_.NumVertices() -1) / 3){
        

        uint non_isolate_u2 = NOT_EXIST, isolate_u2 = NOT_EXIST;
        uint non_isolate_minE2 = NOT_EXIST, isolate_minE2 = NOT_EXIST;
    
        uint u2;

        for (uint i = 0; i < temp_extendable.size(); i++)
        {
            if (m[i] != UNMATCHED) continue;
    
            // Check if an extendable query vertex is isolated or not
            if (temp_extendable[i].matched_nbrs == query_.GetNeighbors(i).size())
            {
                if (temp_extendable[i].E < isolate_minE2)
                {
                    isolate_minE2 = temp_extendable[i].E; // local
                    isolate_u2 = i; // local
                }
            }
            else
            {
                if (temp_extendable[i].E < non_isolate_minE2)
                {
                    non_isolate_minE2 = temp_extendable[i].E; // local
                    non_isolate_u2 = i; // local
                }
            }
        }

        if (non_isolate_minE2 == NOT_EXIST)
            u2 = isolate_u2; // local
        else
            u2 = non_isolate_u2; // local
            
        uint u_min2 = temp_extendable[u2].u_min;
        temp_extendable[u2] = {}; // local
    
        // Enumerate each neighbor of m[u_min]
        // bool candidate_empty = true;
        size_t total_size2 = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();

            // #pragma omp critical
            {
                for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
                {            
                    job_queue.push(std::make_tuple(u2, u_min2, v_idx2, m, temp_extendable, depth2, v));
                }
            }
        }
        else{
            Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);
        }

        
#ifdef MORE_TUNE
if(total_size < NUM_THREAD){

    std::vector<size_t> local_num_results(total_size, 0);
        
        #pragma omp parallel for num_threads(total_size)  // Incremental Matching: 18576.3ms
        for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {
            // thread's local_num_result
            size_t thread_id = omp_get_thread_num();

            local_vec_m[thread_id] = m;
            local_vec_extendable[thread_id] = extendable;
            // local_vec_visited_local[thread_id] = visited_local;
            
            process_vertex(
                u, u_min, v_idx, 
                local_vec_m[thread_id], 
                local_vec_extendable[thread_id], 
                local_num_results[thread_id], 
                // local_vec_visited_local[thread_id], 
                depth);
        }

        for (size_t i = 0; i < local_num_results.size(); ++i) {
            num_results += local_num_results[i];
        }
    
    // process_vertex(u, u_min, 0, m, extendable, num_results, visited_local, depth);
} else { 

    // auto start_nano = std::chrono::high_resolution_clock::now();


    // 40266
        
        std::vector<size_t> local_num_results(NUM_THREAD, 0);
        
        #pragma omp parallel for num_threads(NUM_THREAD)  // Incremental Matching: 18576.3ms
        for (size_t v_idx = 0; v_idx < total_size ; ++v_idx) {
            size_t thread_id = omp_get_thread_num();

            local_vec_m[thread_id] = m;
            local_vec_extendable[thread_id] = extendable;
            // local_vec_visited_local[thread_id] = visited_local;
            
            process_vertex(
                u, u_min, v_idx, 
                local_vec_m[thread_id], 
                local_vec_extendable[thread_id], 
                local_num_results[thread_id], 
                // local_vec_visited_local[thread_id], 
                depth);

        }

        for (size_t i = 0; i < local_num_results.size(); ++i) {
            num_results += local_num_results[i];
        }

        // num_counter++;
}

Print_Time_Nano("My_Duration on OpenMP: ", start_nano);
#endif


    }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}


// Evauluation
/**
 * @brief Evaluates the memory cost of the subgraph matching index.
 * 
 * This function calculates the memory cost of the index by counting the number of vertices
 * and edges stored in the data structures used for subgraph matching. It also provides
 * detailed statistics about the number of candidate vertices and edges, as well as the
 * number of valid candidates that satisfy the constraints.
 * 
 * @param num_edges A reference to store the total number of edges in the index.
 * @param num_vertices A reference to store the total number of vertices in the index.
 * 
 * @details
 * 1. Iterates over all vertices in the query graph and their corresponding children or parents.
 * 2. Counts the number of vertices and edges stored in the index for each query vertex.
 * 3. Separately counts the number of vertices and edges that satisfy `d1` and `d2` constraints.
 * 4. Outputs detailed statistics about the index, including:
 *    - Total number of vertices and edges in the index.
 *    - Number of candidate vertices and edges.
 *    - Number of valid candidate vertices and edges.
 * 
 * @note This function is primarily used for debugging and performance evaluation.
 */
void Parrllel_SymBi::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;

    std::vector<std::unordered_set<uint>> all_vertices_in_index(query_.NumVertices());
    std::vector<std::unordered_set<uint>> d1_vertices_in_index(query_.NumVertices());
    std::vector<std::unordered_set<uint>> d2_vertices_in_index(query_.NumVertices());

    std::vector<size_t> num_all_edges(query_.NumEdges() * 2, 0);
    std::vector<size_t> num_d1_edges(query_.NumEdges() * 2, 0);
    std::vector<size_t> num_d2_edges(query_.NumEdges() * 2, 0);

    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        if (!children.empty())
        {
            for (auto& u_other: children)
            {
                const uint eidx = eidx_[i][u_other];
                for (const auto& [v, d_nbrs]: DCS_[eidx])
                {
                    if (d_nbrs.empty()) continue;

                    all_vertices_in_index[i].insert(v);
                    if (d1[i][v])
                    {
                        d1_vertices_in_index[i].insert(v);
                    }
                    if (d2[i][v])
                    {
                        d2_vertices_in_index[i].insert(v);
                    }

                    for (uint j = 0; j < d_nbrs.size(); j++)
                    {
                        num_all_edges[eidx] += 1;
                        if (d1[u_other][d_nbrs[j]])
                        {
                            num_d1_edges[eidx]++;
                        }
                        if (d2[i][v])
                        {
                            num_d2_edges[eidx]++;
                        }
                    }
                }
            }
        }
        else
        {
            const auto& parents = treeNode_[i].backwards_;
            
            for (auto& u_other: parents)
            {
                const uint eidx = eidx_[i][u_other];
                for (const auto& [v, d_nbrs]: DCS_[eidx])
                {
                    if (d_nbrs.empty()) continue;

                    all_vertices_in_index[i].insert(v);
                    if (d1[i][v])
                        d1_vertices_in_index[i].insert(v);
                    if (d2[i][v])
                        d2_vertices_in_index[i].insert(v);
                }
            }
        }
    }

    size_t total_num_candidate_vs = 0ul;
    size_t total_num_valid_candidate_vs = 0ul;
    size_t total_num_candidate_es = 0ul;
    size_t total_num_valid_candidate_es = 0ul;

    std::cout << "# vertices in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << all_vertices_in_index[i].size() << " ";
        num_vertices += all_vertices_in_index[i].size();
    }
    std::cout << "\n# d1 vertex in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << d1_vertices_in_index[i].size() << " ";
        total_num_candidate_vs += d1_vertices_in_index[i].size();
    }
    std::cout << "\n# d2 vertex in index: "; 
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        std::cout << i << ": " << d2_vertices_in_index[i].size() << " ";
        total_num_valid_candidate_vs += d2_vertices_in_index[i].size();
    }

    std::cout << "\n# edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_all_edges[eidx] << " ";
            
            num_edges += num_all_edges[eidx];
        }
    }
    std::cout << "\n# d1 edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_d1_edges[eidx] << " ";

            total_num_candidate_es += num_d1_edges[eidx];
        }
    }
    std::cout << "\n# d2 edges in index: ";
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const auto& children = treeNode_[i].forwards_;

        for (auto& u_other: children)
        {
            const uint eidx = eidx_[i][u_other];
            std::cout << i << "-" << u_other << ": " <<
                num_d2_edges[eidx] << " ";

            total_num_valid_candidate_es += num_d2_edges[eidx];
        }
    }
    std::cout << "\n# candidates vertices: " << total_num_candidate_vs;
    std::cout << "\n# valid candidates vertices: " << total_num_valid_candidate_vs;
    std::cout << "\n# candidates edges: " << total_num_candidate_es;
    std::cout << "\n# valid candidates edges: " << total_num_valid_candidate_es << std::endl;

    std::cout << "num_counter:" << num_counter << std::endl;
}


