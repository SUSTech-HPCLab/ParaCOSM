#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
#include <cstring> // memcpy
#include <stack>

#include <thread>
#include <mutex>

#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/Parallel_TurboFlux/parallel_turboflux.h"

Parallel_TurboFlux::Parallel_TurboFlux(Graph& query_graph, Graph& data_graph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,     size_t NUMTHREAD,
        size_t auto_tuning)
: matching(
    query_graph, data_graph, max_num_results, 
    print_prep, print_enum, homo)
, eidx_(query_.NumVertices())
, treeNode_(query_.NumVertices())
, q_root_(0u)
, serialized_tree_(query_.NumVertices())
, order_vs_((query_.NumEdges() + 1) * 2)
, backward_vs_((query_.NumEdges() + 1) * 2)

, join_check_vs_((query_.NumEdges() + 1) * 2)
, join_check_labels_((query_.NumEdges() + 1) * 2)

, DCS_(query_.NumEdges() * 2)
, d1(query_.NumVertices())
, d2(query_.NumVertices())
, n1(query_.NumEdges() * 2)
, np1(query_.NumVertices())
, n2(query_.NumEdges() * 2)
, nc2(query_.NumVertices())
, NUMTHREAD(NUMTHREAD)
, auto_tuning(auto_tuning)
, Q1{}
, Q2{}
{
    auto BIG_THREAD =  NUMTHREAD + 2;

    vertex_vector.resize(0);

    if(query_.NumVertices() > 32){
        BIG_THREAD = query_.NumVertices() + 2;
    }else{
        BIG_THREAD = NUMTHREAD + 2;
    }

    if(auto_tuning == 1){
        if(query_.NumVertices() < 7){
            if(NUMTHREAD > query_.NumVertices()){
                NUMTHREAD = query_.NumVertices();
                std::cout << "NUMTHREAD is set to " << NUMTHREAD << "for auto-tuning" << std::endl;
            }
        }
    }

    local_vec_m = std::vector<std::vector<uint>>(BIG_THREAD, std::vector<uint>(query_.NumVertices())); 
    local_vec_visited_local = std::vector<std::vector<bool>>(BIG_THREAD, std::vector<bool>(data_.NumVertices(), false));

    uint edge_pos = 0;
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

    for (uint i = 0; i < (query_.NumEdges() + 1) * 2; i++)
    {
        order_vs_[i].resize(query_.NumVertices());
        backward_vs_[i].resize(query_.NumVertices());
        join_check_vs_[i].resize(query_.NumVertices());
        join_check_labels_[i].resize(query_.NumVertices());
    }
}

void Parallel_TurboFlux::Preprocessing()
{
    double start_time, end_time;
    start_time = omp_get_wtime();  // start
    BuildDAG();
    end_time = omp_get_wtime();  // end
    std::cout << "Time taken to build DAG: " << end_time - start_time << " seconds." << std::endl;
    start_time = omp_get_wtime();  // start 
    BuildDCS();
    end_time = omp_get_wtime();  // end
    GenerateMatchingOrder();
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
void Parallel_TurboFlux::BuildDAG()
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
    for (auto& vec: edge_freqs) vec.resize(query_.NumVertices(), 0ul);
    size_t min_freq = ULONG_MAX;
    std::array<uint, 2> min_vs {};

    for (uint u = 0; u < query_.NumVertices(); u++)
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

            for (size_t v = 0; v < data_.NumVertices(); v++)
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

    // build a spaning tree with the greedy approach
    std::vector<bool> visited(query_.NumVertices(), false);
    std::vector<EdgeFreq> edge_queue;

    visited[q_root_] = true;
    {
        const auto& q_root_nbrs = query_.GetNeighbors(q_root_);
        for (uint i = 0; i < q_root_nbrs.size(); i++)
            edge_queue.emplace_back(q_root_, q_root_nbrs[i], 
                edge_freqs[q_root_][q_root_nbrs[i]]);
    }

    while (!edge_queue.empty())
    {
        std::vector<EdgeFreq> next_level;
        while (!edge_queue.empty())
        {
            EdgeFreq local_min_info = *std::min_element(
                edge_queue.begin(), edge_queue.end(),
                [](const EdgeFreq &v1, const EdgeFreq &v2){
                    return v1.freq < v2.freq;
            });
            auto it = std::remove_if(edge_queue.begin(), edge_queue.end(),
            [&local_min_info](const EdgeFreq& v){
                return v.v2 == local_min_info.v2;
            });

            edge_queue.resize(it - edge_queue.begin());

            if (visited[local_min_info.v2]) continue;
            visited[local_min_info.v2] = true;

            auto edge_label = std::get<2>(query_.GetEdgeLabel(local_min_info.v1, local_min_info.v2));
            treeNode_[local_min_info.v1].forwards_.push_back(local_min_info.v2);
            treeNode_[local_min_info.v1].forward_labels_.push_back(edge_label);
            treeNode_[local_min_info.v2].backwards_.push_back(local_min_info.v1);
            treeNode_[local_min_info.v2].backward_labels_.push_back(edge_label);

            const auto& q_root_nbrs = query_.GetNeighbors(local_min_info.v2);
            for (uint i = 0; i < q_root_nbrs.size(); i++)
            {
                const auto& u_other = q_root_nbrs[i];
                if (!visited[u_other])
                    next_level.emplace_back(local_min_info.v2, u_other, 
                        edge_freqs[local_min_info.v2][u_other]);
            }

        }
        std::swap(edge_queue, next_level);
    }
    // compute serialized tree
    uint write_pos = 0u;
    std::queue<uint> bfs_queue;
    bfs_queue.push(q_root_);

    while (!bfs_queue.empty())
    {
        const uint size = bfs_queue.size();
        for (uint i = 0; i < size; i++)
        {
            const uint front = bfs_queue.front();
            bfs_queue.pop();
            serialized_tree_[write_pos++] = front;
            for (const uint& u_other: treeNode_[front].forwards_)
            {
                bfs_queue.push(u_other);
            }
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
 * @brief Builds the Dynamic Candidate Sets (DCS) for parallel subgraph matching.
 * 
 * This function constructs the hierarchical index structure that maintains candidate vertices
 * and their relationships for efficient subgraph matching. The construction follows a three-phase approach:
 * 
 * @details
 * 1. Top-down Phase:
 *    Processes query vertices in a top-down manner following the serialized spanning tree.
 *    For each query vertex u and corresponding data vertex v with matching labels:
 *    - Initializes empty DCS entries for all neighbors
 *    - For each backward edge (parent relationship), finds matching data vertices
 *    - Populates the DCS entries with connected vertices that satisfy label constraints
 *    - Updates validity counters (n1, np1) and marks vertices as valid (d1=1) when appropriate
 * 
 * 2. Bottom-up Phase:
 *    Processes query vertices from leaves to root:
 *    - For each forward edge (child relationship), finds matching data vertices
 *    - Updates the DCS entries and validity counters (n2, nc2)
 *    - Marks vertices as fully valid (d2=1) when they satisfy both upward and downward constraints
 * 
 * 3. Refinement Phase:
 *    Performs a final top-down pass to update n2 counters for consistency
 *    across the entire indexing structure.
 * 
 * The resulting DCS structure efficiently supports parallel subgraph matching by:
 * - Providing quick access to candidate vertices for each query vertex
 * - Maintaining parent-child relationships between candidates
 * - Tracking validity states (d1, d2) to enable efficient filtering
 * - Supporting parallelized enumeration through thread-safe accessors
 */
void Parallel_TurboFlux::BuildDCS()
{
    // top-down
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[i];
        const uint query_label = query_.GetVertexLabel(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++)
        if (data_.GetVertexLabel(j) == query_label)
        {
            // add j to hash maps
            const auto& q_nbrs = query_.GetNeighbors(u);
            for (uint k = 0; k < q_nbrs.size(); k++)
            {
                DCS_[eidx_[u][q_nbrs[k]]][j] = {};
            }
            for (uint k = 0; k < treeNode_[u].backwards_.size(); k++)
            {
                const uint u_other = treeNode_[u].backwards_[k];
                const uint u_other_elabel = treeNode_[u].backward_labels_[k];

                const auto& d_nbrs = data_.GetNeighbors(j);
                const auto& d_elabels = data_.GetNeighborLabels(j);

                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];
                    const uint v_other_elabel = d_elabels[m];
                    if (
                        query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                        u_other_elabel == v_other_elabel
                    ) {
                        DCS_[eidx_[u][u_other]][j].push_back(v_other);
                        if (d1[u_other][v_other] > 0)
                            n1[eidx_[u][u_other]][j] += 1;
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

    // bottom-up
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        const uint u = serialized_tree_[query_.NumVertices() - 1 - i];
        const uint query_label = query_.GetVertexLabel(u);
        
        for (size_t j = 0; j < data_.NumVertices(); j++)
        if (data_.GetVertexLabel(j) == query_label)
        {
            for (uint k = 0; k < treeNode_[u].forwards_.size(); k++)
            {
                const uint u_other = treeNode_[u].forwards_[k];
                const uint u_other_elabel = treeNode_[u].forward_labels_[k];

                const auto& d_nbrs = data_.GetNeighbors(j);
                const auto& d_elabels = data_.GetNeighborLabels(j);

                for (uint m = 0; m < d_nbrs.size(); m++)
                {
                    const uint v_other = d_nbrs[m];
                    const uint v_other_elabel = d_elabels[m];
                    if (
                        query_.GetVertexLabel(u_other) == data_.GetVertexLabel(v_other) &&
                        u_other_elabel == v_other_elabel
                    ) {
                        DCS_[eidx_[u][u_other]][j].push_back(v_other);
                        if (d2[u_other][v_other] > 0)
                            n2[eidx_[u][u_other]][j] += 1;
                    }
                }
                if (n2[eidx_[u][u_other]][j])
                {
                    nc2[u][j] += 1;
                }
            }
            if (nc2[u][j] == treeNode_[u].forwards_.size() && d1[u][j] == 1)
            {
                d2[u][j] = 1;
            }
        }
    }

    // top-down modify n2
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
 * @brief Recursively counts the number of valid paths in the query graph.
 * 
 * This function traverses the query graph in a top-down manner, counting the number of 
 * valid paths that can be constructed from a given vertex pair (u, v). It uses dynamic 
 * programming to avoid redundant calculations for previously visited vertices.
 * 
 * @param u The current query vertex being processed.
 * @param v The corresponding data vertex that matches u.
 * @param num_explicit_pathes Vector storing the count of explicit paths for each query vertex.
 * @param num_dp Dynamic programming table storing counts for previously processed vertex pairs.
 * 
 * @details
 * The function works by:
 * 1. Examining each forward edge (child relationship) from the current query vertex u.
 * 2. For each child u_c of u:
 *    - If (u_c, v) has not been processed before:
 *      - Iterate through all candidate vertices v_c in DCS that match with u_c and are
 *        connected to v.
 *      - For each valid candidate (satisfying d2 constraint), increment the path count
 *        and recursively process the child vertex pair (u_c, v_c).
 *    - If (u_c, v) has been processed before, reuse the previously calculated count.
 * 3. The results are used to determine the optimal matching order for subgraph enumeration.
 * 
 * This approach efficiently counts paths while avoiding exponential explosion by using
 * memoization to store and reuse intermediate results.
 */
void Parallel_TurboFlux::CountDownwards(uint u, uint v, std::vector<uint>& num_explicit_pathes,
        std::vector<std::unordered_map<uint, uint>>& num_dp)
{
    for (uint i = 0; i < treeNode_[u].forwards_.size(); i++)
    {
        const auto& u_c = treeNode_[u].forwards_[i];

        if (num_dp[u_c].find(v) == num_dp[u_c].end())
        {
            for (uint j = 0; j < DCS_[eidx_[u][u_c]][v].size(); j++)
            {
                const auto& v_c = DCS_[eidx_[u][u_c]][v][j];

                if (d2[u_c][v_c])
                {
                    num_explicit_pathes[u_c] ++;
                    num_dp[u_c][v] ++;
                    CountDownwards(u_c, v_c, num_explicit_pathes, num_dp);
                }
            }
        }
        else
        {
            num_explicit_pathes[u_c] += num_dp[u_c][v];
        }
    }
}



/**
 * @brief Generates optimal matching orders for subgraph enumeration.
 * 
 * This function creates multiple matching orders for different edge-triggered scenarios
 * to optimize the subgraph enumeration process. It uses path statistics to determine the
 * most efficient order to explore vertices during matching.
 * 
 * @details
 * The function works in four main phases:
 * 
 * 1. Path Counting:
 *    - Calls CountDownwards() to determine the number of valid paths for each query vertex
 *    - Uses dynamic programming to efficiently count paths without redundant calculations
 * 
 * 2. Initial Order Generation:
 *    - Identifies leaf vertices (those with no children in the spanning tree)
 *    - Selects vertices in increasing order of path counts (fewest paths first)
 *    - Ensures parent vertices are processed only after all their children
 *    - Creates a globally optimal matching order (initial_order)
 * 
 * 3. Edge-Specific Order Creation:
 *    - For each edge in the query graph, creates a specialized matching order
 *    - Considers both tree edges and non-tree edges
 *    - Handles different topological relationships (parent-child, cross edges)
 *    - Follows an upward-then-downward traversal pattern from the affected vertices
 * 
 * 4. Joinability Constraints Setup:
 *    - For each vertex in each order, identifies which previously matched vertices
 *      need to be checked for connectivity (join_check_vs_)
 *    - Records the corresponding edge labels for validation (join_check_labels_)
 * 
 * The resulting matching orders enable efficient edge-triggered subgraph enumeration
 * by prioritizing vertices that are most constrained and likely to prune the search space.
 */
void Parallel_TurboFlux::GenerateMatchingOrder()
{
    std::vector<uint> num_explicit_pathes (query_.NumVertices());
    std::vector<std::unordered_map<uint, uint>> num_dp(query_.NumVertices());
    uint first_root_child = treeNode_[q_root_].forwards_.front();
    for (auto& [v, _]: DCS_[eidx_[q_root_][first_root_child]])
    {
        CountDownwards(q_root_, v, num_explicit_pathes, num_dp);
    }
    num_explicit_pathes[q_root_] = UINT_MAX - 1;
        
    std::unordered_set<uint> leafs, visited;
    for (uint i = 0; i < query_.NumVertices(); i++)
        if (treeNode_[i].forwards_.empty())
            leafs.insert(i);

    std::vector<uint> initial_order(query_.NumVertices());
    uint order_pos = query_.NumVertices() - 1;
    for (uint i = 0; i < query_.NumVertices(); i++)
    {
        uint min_value = UINT_MAX;
        uint min_value_index = NOT_EXIST;
        for (auto& j: leafs)
        {
            if (num_explicit_pathes[j] < min_value)
            {
                min_value = num_explicit_pathes[j];
                min_value_index = j;
            }
        }
        leafs.erase(min_value_index);
        initial_order[order_pos--] = min_value_index;
        visited.insert(min_value_index);

        if (min_value_index != q_root_)
        {
            bool all_children_visited = true;
            auto parent_id = treeNode_[min_value_index].backwards_.front();
            for (auto j: treeNode_[parent_id].forwards_)
                if (visited.find(j) == visited.end())
                {
                    all_children_visited = false;
                    break;
                }
            if (all_children_visited)
                leafs.insert(parent_id);
        }
    }

    // generate incremental matching orders
    for (uint u = 0u; u < query_.NumVertices(); u++)
    {
        auto& q_nbrs = query_.GetNeighbors(u);

        for (uint i = 0; i < q_nbrs.size(); i++)
        {
            uint u_other = q_nbrs[i];
            if (u > u_other) continue;

            bool reversed = false;
            if (!treeNode_[u_other].backwards_.empty() && treeNode_[u_other].backwards_.front() == u)
            {
                std::swap(u, u_other);
                reversed = true;
            }
            if (!treeNode_[u].backwards_.empty() && treeNode_[u].backwards_.front() == u_other)
            {
                std::vector<bool> temp_visited(query_.NumVertices(), false);
                auto& order = order_vs_[eidx_[u][u_other]];
                auto& backwards = backward_vs_[eidx_[u][u_other]];
                order_pos = 0ul;

                uint cur = u_other;
                temp_visited[u_other] = true;
                temp_visited[u] = true;
                order[order_pos++] = u;
                order[order_pos++] = u_other;

                // upwards
                uint pre = cur;
                while (cur != q_root_)
                {
                    cur = treeNode_[cur].backwards_.front();
                    temp_visited[cur] = true;
                    order[order_pos] = cur;
                    backwards[order_pos++] = pre;
                    pre = cur;
                }

                // downwards
                for (uint j = 0; j < query_.NumVertices(); j++)
                {
                    uint local_u = initial_order[j];
                    if (!temp_visited[local_u])
                    {
                        order[order_pos] = local_u;
                        backwards[order_pos++] = treeNode_[local_u].backwards_.front();

                        auto& q_nbrs_local = query_.GetNeighbors(local_u);
                        auto& q_nbr_labels_local = query_.GetNeighborLabels(local_u);

                        for (uint k = 0; k < q_nbrs_local.size(); k++)
                        {
                            const uint local_u_other = q_nbrs_local[k];
                            const uint local_u_other_label = q_nbr_labels_local[k];

                            if (temp_visited[local_u_other] && treeNode_[local_u].backwards_.front() != local_u_other)
                            {
                                join_check_vs_[eidx_[u][u_other]][local_u].push_back(local_u_other);
                                join_check_labels_[eidx_[u][u_other]][local_u].push_back(local_u_other_label);
                            }
                        }
                        temp_visited[local_u] = true;
                    }
                }
            }
            else
            {
                std::vector<bool> temp_visited(query_.NumVertices(), false);
                auto& order = order_vs_[eidx_[u][u_other]];
                auto& backwards = backward_vs_[eidx_[u][u_other]];
                order_pos = 0ul;

                uint cur = u_other;
                temp_visited[u_other] = true;
                temp_visited[u] = true;
                order[order_pos++] = u;
                order[order_pos++] = u_other;

                // upwards
                uint pre = cur;
                while (cur != q_root_)
                {
                    cur = treeNode_[cur].backwards_.front();
                    auto non_tree_edge = query_.GetEdgeLabel(u, cur);
                    if (std::get<2>(non_tree_edge) != NOT_EXIST)
                    {
                        join_check_vs_[eidx_[u][u_other]][cur].push_back(u);
                        join_check_labels_[eidx_[u][u_other]][cur].push_back(std::get<2>(non_tree_edge));
                    }
                    temp_visited[cur] = true;
                    order[order_pos] = cur;
                    backwards[order_pos++] = pre;
                    pre = cur;
                }

                // downwards
                for (uint j = 0; j < query_.NumVertices(); j++)
                {
                    uint local_u = initial_order[j];
                    if (!temp_visited[local_u])
                    {
                        order[order_pos] = local_u;
                        backwards[order_pos++] = treeNode_[local_u].backwards_.front();

                        auto& q_nbrs_local = query_.GetNeighbors(local_u);
                        auto& q_nbr_labels_local = query_.GetNeighborLabels(local_u);

                        for (uint k = 0; k < q_nbrs_local.size(); k++)
                        {
                            const uint local_u_other = q_nbrs_local[k];
                            const uint local_u_other_label = q_nbr_labels_local[k];

                            if (temp_visited[local_u_other] && treeNode_[local_u].backwards_.front() != local_u_other)
                            {
                                join_check_vs_[eidx_[u][u_other]][local_u].push_back(local_u_other);
                                join_check_labels_[eidx_[u][u_other]][local_u].push_back(local_u_other_label);
                            }
                        }
                        temp_visited[local_u] = true;
                    }
                }
            }
            if (reversed)
            {
                std::swap(u, u_other);
            }
        }
    }
    
    // also generate joinablity check vertices for ARTI -> q_root
    std::vector<bool> temp_visited(query_.NumVertices(), false);
    auto& order = order_vs_[eidx_[q_root_][first_root_child]];
    auto& backwards = backward_vs_[eidx_[q_root_][first_root_child]];
    order_pos = 0ul;

    temp_visited[q_root_] = true;
    order[order_pos++] = q_root_;

    // downwards
    for (uint j = 0; j < query_.NumVertices(); j++)
    {
        uint local_u = initial_order[j];
        if (!temp_visited[local_u])
        {
            order[order_pos] = local_u;
            backwards[order_pos++] = treeNode_[local_u].backwards_.front();

            auto& q_nbrs_local = query_.GetNeighbors(local_u);
            auto& q_nbr_labels_local = query_.GetNeighborLabels(local_u);

            for (uint k = 0; k < q_nbrs_local.size(); k++)
            {
                const uint local_u_other = q_nbrs_local[k];
                const uint local_u_other_label = q_nbr_labels_local[k];

                if (temp_visited[local_u_other] && treeNode_[local_u].backwards_.front() != local_u_other)
                {
                    join_check_vs_[eidx_[q_root_][first_root_child]][local_u].push_back(local_u_other);
                    join_check_labels_[eidx_[q_root_][first_root_child]][local_u].push_back(local_u_other_label);
                }
            }
            temp_visited[local_u] = true;
        }
    }

    
    if (print_preprocessing_results_)
    {
        std::cout << "spanning tree: \n";
        for (uint i = 0; i < query_.NumVertices(); ++i)
        {
            if (initial_order[i] == q_root_)
                std::cout << initial_order[i] << " | ";
            else
                std::cout << treeNode_[initial_order[i]].backwards_.front() 
                    << "-" << initial_order[i] << " | ";
        }
        std::cout << "\nmatching order: ";
        std::cout << "\n-vertex(u_min: joinablity check vertices)-\n";
        for (uint k = 0; k < query_.NumVertices(); k++)
        {
            const uint eidx = eidx_[q_root_][first_root_child];
            const uint cur = order_vs_[eidx][k];
            std::cout << cur << "(";
            std::cout << backward_vs_[eidx][k] << ":";
            for (uint l = 0; l < join_check_vs_[eidx][cur].size(); l++)
            {
                std::cout << join_check_vs_[eidx][cur][l] << ",";
            }
            std::cout << ")-";
        }
        std::cout << std::endl;
        for (uint u = 0; u < query_.NumVertices(); u++)
        {
            const auto& q_nbrs = query_.GetNeighbors(u);

            for (uint j = 0; j < q_nbrs.size(); j++)
            {
                uint u_other = q_nbrs[j];
                if (u > u_other) continue;

                bool reversed = false;
                if (!treeNode_[u_other].backwards_.empty() && treeNode_[u_other].backwards_.front() == u)
                {
                    std::swap(u, u_other);
                    reversed = true;
                }
                for (uint k = 0; k < query_.NumVertices(); k++)
                {
                    const uint eidx = eidx_[u][u_other];
                    const uint cur = order_vs_[eidx][k];
                    std::cout << cur << "(";
                    std::cout << backward_vs_[eidx][k] << ":";
                    for (uint l = 0; l < join_check_vs_[eidx][cur].size(); l++)
                    {
                        std::cout << join_check_vs_[eidx][cur][l] << ",";
                    }
                    std::cout << ")-";
                }
                std::cout << std::endl;
                if (reversed)
                {
                    std::swap(u, u_other);
                }
            }
        }
    }
}


/**
 * @brief Finds all initial matches in the data graph.
 * 
 * This function performs the initial matching of the query graph against the data graph
 * by starting from the root vertex of the query graph and recursively exploring all
 * possible matches. It serves as the entry point for the subgraph matching process.
 * 
 * @details
 * The function works as follows:
 * 1. Creates an empty match array to store the current mapping of query vertices to data vertices.
 * 2. Iterates through all candidate vertices for the query root in the Dynamic Candidate Sets (DCS).
 * 3. For each valid candidate (satisfying d2 constraint):
 *    - Maps the root query vertex to this data vertex
 *    - Marks the data vertex as visited
 *    - Recursively explores all possible matches by calling FindMatches()
 *    - Backtracks by unmarking the data vertex and resetting the mapping
 * 
 * The function uses the spanning tree edge from the root to its first child as the 
 * index to determine the matching order. It maintains the recursive state through the
 * match array (m) and tracks visited vertices to prevent duplicate matches when
 * homomorphism is not enabled.
 */
void Parallel_TurboFlux::InitialMatching()
{
    std::vector<uint> m (query_.NumVertices(), UNMATCHED);
    
    for (auto & [v, d_nbrs]: DCS_[eidx_[q_root_][query_.GetNeighbors(q_root_).front()]])
    {
        if (d2[q_root_][v] == false) continue;

        m[q_root_] = v;
        visited_[v] = true;

        FindMatches(eidx_[q_root_][treeNode_[q_root_].forwards_.front()], 1, m, num_initial_results_);

        visited_[v] = false;
        m[q_root_] = UNMATCHED;
    }
}



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
void Parallel_TurboFlux::InsertionTopDown(uint u, uint u_c, uint v, uint v_c)
{
    if (n1[eidx_[u_c][u]][v_c] == 0)
    {
        np1[u_c][v_c] += 1;
        if (np1[u_c][v_c] == treeNode_[u_c].backwards_.size())
        {
            d1[u_c][v_c] = 1;
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
void Parallel_TurboFlux::InsertionBottomUp(uint u, uint u_p, uint v, uint v_p)
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
void Parallel_TurboFlux::DeletionTopDown(uint u, uint u_c, uint v, uint v_c)
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
void Parallel_TurboFlux::DeletionBottomUp(uint u, uint u_p, uint v, uint v_p)
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


/**
 * @brief Recursively finds matches for the query graph in the data graph.
 * 
 * This function explores all possible mappings of query vertices to data vertices
 * using a depth-first search approach. It checks for joinability and visited status
 * to ensure valid mappings.
 * 
 * @param order_index The index of the current matching order being processed.
 * @param depth The current depth in the recursive search.
 * @param m A reference to the current mapping of query vertices to data vertices.
 * @param num_results A reference to the count of valid matches found so far.
 * 
 * @details
 * The function works by:
 * 1. Enumerating each neighbor of `m[u_min]` (the current vertex being processed).
 * 2. For each neighbor:
 *    - Checks if it is a valid candidate (index check).
 *    - Checks if it is joinable with previously matched vertices.
 *    - If valid, adds the mapping and marks the vertex as visited.
 *    - Recursively calls itself to explore further mappings.
 *    - Backtracks by unmarking the vertex and resetting the mapping.
 */
void Parallel_TurboFlux::FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{
    if (reach_time_limit) return;

    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];

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
        for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
        {
            const auto& u_backward = join_check_vs_[order_index][u][i];
            const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
            const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

            if (std::get<2>(d_elabel) != u_backward_elabel)
            {
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
inline void Parallel_TurboFlux::Parallel_FindMatches(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{

    size_t NUMT = NUMTHREAD;

    vertex_vector.clear();
    job_queue.clear();

    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];


    std::vector<size_t> local_num_result(NUMT+2,0);

    // enumerate each neighbor of m[u_min]
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        if (d2[u][v] == 0) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
        {
            const auto& u_backward = join_check_vs_[order_index][u][i];
            const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
            const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

            if (std::get<2>(d_elabel) != u_backward_elabel)
            {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        
        // 4. add a vertex mapping
        m[u] = v;
        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = true;
        }
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            // FindMatches(order_index, depth + 1, m, num_results);
            auto depth2 = depth + 1;
            uint u2 = order_vs_[order_index][depth2];
            uint u_min2 = backward_vs_[order_index][depth2];

            auto siez_DCS = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();


            for (auto v_idx2 = 0; v_idx2 < siez_DCS; ++v_idx2)
            {

                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, 
                     v);

            }
        }

        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = false;
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;

    }

    if(vertex_vector.size() < NUMT){
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }

    #pragma omp parallel for num_threads(NUMT)
    for(size_t i =0; i < vertex_vector.size(); ++i){
        auto& [u3, u_min2, v_idx2, m2, v3] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();
        // visited_[v3] = true;
        // m2[u3] = v3;
        local_vec_visited_local[thread_id][v3] = true;
        m2[u3] = v3;

        ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
            depth + 1, order_index, thread_id);

        // visited_[v3] = false;
        // m2[u3] = UNMATCHED;
        local_vec_visited_local[thread_id][v3] = false;
        m2[u3] = UNMATCHED;

        if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
            std::tuple<uint, uint, size_t, std::vector<uint>,
                 uint , uint> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
                local_vec_visited_local[thread_id][v3] = true;
                m2[u] = v3;
                ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                    depth3, order_index, thread_id);
                local_vec_visited_local[thread_id][v3] = false; 
                m2[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

}


/**
 * @brief Processes a single vertex in the parallel matching algorithm.
 * 
 * This function evaluates whether a candidate vertex from the data graph can be matched to 
 * a query vertex in the current partial matching. It performs filtering based on indexing 
 * constraints, joinability with previously matched vertices, and visitation status.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The parent query vertex of u in the matching order.
 * @param v_idx The index of the candidate vertex in the DCS array.
 * @param m A reference to the current partial matching (query vertex -> data vertex).
 * @param num_results A reference to the counter for valid matches found.
 * @param depth The current depth in the matching process.
 * @param order_index The index of the current matching order being used.
 * @param thread_id The ID of the thread executing this function.
 * 
 * @return Returns true if the candidate vertex was successfully processed, false otherwise.
 * 
 * @details
 * The function follows these steps:
 * 1. Retrieves the candidate vertex from the DCS structure.
 * 2. Checks if it satisfies the index constraint (d2[u][v]).
 * 3. Verifies if the candidate is joinable with all previously matched vertices.
 * 4. Ensures the vertex hasn't been visited in the current search path (for isomorphism).
 * 5. If all checks pass, adds the vertex to the current matching and marks it as visited.
 * 6. Either records a complete match or recursively explores further matches.
 * 7. Backtracks by unmarking the vertex and removing it from the matching.
 * 
 * This thread-safe version uses thread-local data structures to maintain correctness
 * during parallel execution.
 */
inline bool Parallel_TurboFlux::ProcessVertex(uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    size_t &num_results, 
     uint depth, size_t order_index, size_t thread_id) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];     

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable2 = true;
    for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++) {
        const auto& u_backward = join_check_vs_[order_index][u][i];
        const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
        const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

        if (std::get<2>(d_elabel) != u_backward_elabel) {
            joinable2 = false;
            break;
        }
    }
    if (!joinable2) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add vertex mapping and mark visited
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true;

    // Process result or recurse
    if (depth == query_.NumVertices() - 1) {
        num_results++;
    } else {
        FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
}


/**
 * @brief Thread-local version of FindMatches for parallel subgraph matching.
 * 
 * This function performs depth-first search for subgraph matching within a single thread's
 * context. It maintains thread-local data structures to avoid synchronization overhead
 * while ensuring correctness in a multi-threaded environment.
 * 
 * @param order_index The index of the matching order to use.
 * @param depth The current depth in the recursive search process.
 * @param m A vector representing the current partial mapping of query vertices to data vertices.
 * @param num_results A reference to the thread-local counter for matches found.
 * @param thread_id The ID of the thread executing this function instance.
 * 
 * @details
 * The function follows a standard backtracking approach:
 * 
 * 1. Retrieves the next query vertex (u) and its parent (u_min) from the matching order.
 * 
 * 2. Filters candidates through a pipeline of increasingly strict conditions:
 *    a) Index-based filtering: Verifies if the candidate satisfies structural constraints via d2[u][v]
 *    b) Joinability checking: Ensures the candidate connects properly to previously matched vertices
 *    c) Visitation checking: Prevents reusing vertices in the current search path (for isomorphism)
 * 
 * 3. For each valid candidate:
 *    - Adds it to the current mapping and marks it as visited in thread-local structures
 *    - If at maximum depth, increments the match counter
 *    - Otherwise, recursively explores deeper matches
 *    - Backtracks by unmarking the vertex and removing it from the mapping
 * 
 * This thread-local implementation avoids race conditions by using separate visitation arrays
 * for each thread while maintaining the exploration efficiency of the sequential algorithm.
 */
void Parallel_TurboFlux::FindMatches_local(uint order_index, uint depth, std::vector<uint>& m, 
     size_t &num_results, size_t thread_id)
{
    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];

    // enumerate each neighbor of m[u_min]
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        if (d2[u][v] == 0) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
        {
            const auto& u_backward = join_check_vs_[order_index][u][i];
            const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
            const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

            if (std::get<2>(d_elabel) != u_backward_elabel)
            {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && local_vec_visited_local[thread_id][v]) continue;
        
        // 4. add a vertex mapping
        m[u] = v;
        local_vec_visited_local[thread_id][v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
        }
        
        local_vec_visited_local[thread_id][v] = false;
        m[u] = UNMATCHED;
    }
}


/**
 * @brief Processes a single vertex in the parallel matching algorithm using job queue approach.
 * 
 * This function evaluates whether a candidate vertex from the data graph can be matched to
 * a query vertex, and if valid, adds subsequent matching tasks to a concurrent job queue
 * for parallel processing by worker threads.
 * 
 * @param u The current query vertex being processed.
 * @param u_min The parent query vertex of u in the matching order.
 * @param v_idx The index of the candidate vertex in the DCS array.
 * @param m A reference to the current partial matching (query vertex -> data vertex).
 * @param num_results A reference to the thread-local counter for valid matches found.
 * @param depth The current depth in the matching process.
 * @param order_index The index of the current matching order being used.
 * @param thread_id The ID of the thread executing this function.
 * 
 * @return Returns true if the candidate vertex was successfully processed, false otherwise.
 * 
 * @details
 * The function follows these steps:
 * 1. Retrieves the candidate vertex from the DCS structure.
 * 2. Applies filtering based on the index constraint (d2[u][v]).
 * 3. Verifies candidate joinability with previously matched vertices.
 * 4. Ensures the vertex hasn't been visited in the current search path (for isomorphism).
 * 5. If all checks pass, adds the vertex to the current mapping and marks it as visited.
 * 6. For complete matches (at maximum depth), increments the match counter.
 * 7. For partial matches, instead of recursive calls, adds new matching tasks to the job queue.
 * 8. Backtracks by unmarking the vertex and removing it from the mapping.
 * 
 * This function enables work-stealing parallelism by scheduling tasks that can be executed
 * by any available thread, improving load balancing in the parallel matching process.
 */
inline bool Parallel_TurboFlux::ProcessVertex_queue(uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    size_t &num_results, 
     uint depth, size_t order_index, size_t thread_id) {

    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];     

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable2 = true;
    for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++) {
        const auto& u_backward = join_check_vs_[order_index][u][i];
        const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
        const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

        if (std::get<2>(d_elabel) != u_backward_elabel) {
            joinable2 = false;
            break;
        }
    }
    if (!joinable2) return false;

    // 3. Check if visited
    if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add vertex mapping and mark visited
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true;

    // Process result or recurse
    if (depth == query_.NumVertices() - 1) {
        num_results++;
    } else {
        size_t total_size2 = DCS_[eidx_[u_min][u]][m[u_min]].size();
        for (size_t v_idx2 = 0; v_idx2 < total_size2; ++v_idx2)
        {            
            job_queue.push(std::make_tuple(u, u_min, v_idx2, m, depth+1, v));
        }
        // FindMatches_local(order_index, depth + 1, m, num_results, thread_id);
    }

    // Backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
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
inline void Parallel_TurboFlux::Parallel_FindMatches_delete(uint order_index, uint depth, std::vector<uint>& m, size_t &num_results)
{

    size_t NUMT = NUMTHREAD;

    vertex_vector.clear();
    job_queue.clear();

    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];


    std::vector<size_t> local_num_result(NUMT+2,0);

    // enumerate each neighbor of m[u_min]
    for (auto& v: DCS_[eidx_[u_min][u]][m[u_min]])
    {
        // 1. check index
        if (d2[u][v] == 0) continue;

        // 2. check if joinable
        bool joinable = true;
        for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
        {
            const auto& u_backward = join_check_vs_[order_index][u][i];
            const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
            const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

            if (std::get<2>(d_elabel) != u_backward_elabel)
            {
                joinable = false;
                break;
            }
        }
        if (!joinable) continue;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        
        // 4. add a vertex mapping
        m[u] = v;
        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = true;
        }
        visited_[v] = true;

        if (depth == query_.NumVertices() - 1)
        {
            num_results++;
        }
        else
        {
            // FindMatches(order_index, depth + 1, m, num_results);
            auto depth2 = depth + 1;
            uint u2 = order_vs_[order_index][depth2];
            uint u_min2 = backward_vs_[order_index][depth2];

            auto siez_DCS = DCS_[eidx_[u_min2][u2]][m[u_min2]].size();


            for (auto v_idx2 = 0; v_idx2 < siez_DCS; ++v_idx2)
            {
                vertex_vector.emplace_back(u2, u_min2, v_idx2, m, 
                      v);
            }
        }

        for(size_t i = 0; i< local_vec_visited_local.size(); i++){
            local_vec_visited_local[i][v] = false;
        }
        
        visited_[v] = false;
        m[u] = UNMATCHED;

    }

    if(vertex_vector.size() < NUMT){
        if(vertex_vector.size() == 0){
            NUMT = 1;
        }
    }

    #pragma omp parallel for num_threads(NUMT)
    for(size_t i =0; i < vertex_vector.size(); ++i){
        auto& [u3, u_min2, v_idx2, m2,  v3] = vertex_vector[i];
        size_t thread_id = omp_get_thread_num();
        // visited_[v3] = true;
        // m2[u3] = v3;
        local_vec_visited_local[thread_id][v3] = true;
        m2[u3] = v3;

        ProcessVertex(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
            depth + 1, order_index, thread_id);

        // visited_[v3] = false;
        // m2[u3] = UNMATCHED;
        local_vec_visited_local[thread_id][v3] = false;
        m2[u3] = UNMATCHED;

        if(!job_queue.empty() && (i > vertex_vector.size() - NUMT)){
            std::tuple<uint, uint, size_t, std::vector<uint>,
                 uint , uint> job;
            if(job_queue.try_pop(job)){
                size_t thread_id = omp_get_thread_num();
                auto& [u3, u_min2, v_idx2, m2, depth3, v3] = job;
                local_vec_visited_local[thread_id][v3] = true;
                m2[u] = v3;
                ProcessVertex_queue(u3, u_min2, v_idx2, m2, local_num_result[thread_id], 
                    depth3, order_index, thread_id);
                local_vec_visited_local[thread_id][v3] = false; 
                m2[u] = UNMATCHED;
            }
        }
    }

    for (size_t i = 0; i < local_num_result.size(); ++i) {
        num_results += local_num_result[i];
    }

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
inline bool Parallel_TurboFlux::process_vertex_layer1(
    uint order_index,
    uint u, uint u_min, size_t v_idx, 
    std::vector<uint>& m, 
    // std::vector<ExtendableVertex>& extendable, 
    size_t &num_results, 
     uint depth, size_t thread_id) {

    
    // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);
    auto v = DCS_[eidx_[u_min][u]][m[u_min]][v_idx];

    // 1. Check index
    if (d2[u][v] == 0) return false;

    // 2. Check if joinable
    bool joinable = true;
    for (uint i = 0; i < join_check_vs_[order_index][u].size(); i++)
    {
        const auto& u_backward = join_check_vs_[order_index][u][i];
        const auto& u_backward_elabel = join_check_labels_[order_index][u][i];
        const auto& d_elabel = data_.GetEdgeLabel(m[u_backward], v);

        if (std::get<2>(d_elabel) != u_backward_elabel)
        {
            joinable = false;
            break;
        }
    }

    if (!joinable) return false;

    // 3. Check if visited
    // if (!homomorphism_ && local_vec_visited_local[thread_id][v]) return false;

    // 4. Add a vertex mapping
    m[u] = v;
    local_vec_visited_local[thread_id][v] = true; // imp:

    // std::vector<ExtendableVertex> temp_extendable(extendable); // local
    // for (auto& u_other: treeNode_[u].neighbors_) {
    //     if (m[u_other] != UNMATCHED) continue;

    //     if (n2[eidx_[u][u_other]][m[u]] < temp_extendable[u_other].E) {
    //         temp_extendable[u_other].E = n2[eidx_[u][u_other]][m[u]];
    //         temp_extendable[u_other].u_min = u;
    //     }
    //         temp_extendable[u_other].matched_nbrs++;
    // }

    if (depth == query_.NumVertices() - 1) { // match complete
        num_results++;
    } else {
        // __builtin_prefetch(&DCS_[eidx_[u_min][u]][m[u_min]][v_idx], 1, 0);

        // Parallel_FindMatches_local3(depth + 1, m, temp_extendable, num_results,thread_id);

    }

    // for backtrack
    local_vec_visited_local[thread_id][v] = false;
    m[u] = UNMATCHED;

    return true;
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
bool Parallel_TurboFlux::Classify(uint v1, uint v2, uint label){
    auto v1_label = data_.GetVertexLabel(v1);

    for (uint u1 = 0; u1 < query_.NumVertices(); u1++){ // take one u1 from query
        //  Get v1 that label matched to query label
        if (v1_label == query_.GetVertexLabel(u1)){

            // Get v2 that label matched to query vertex
            for (uint u2 = 0; u2 < query_.NumVertices(); u2++){ // take one u2 from query
                if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2)){

                    if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
                   
                    // Prune 1: detect if the edge is in the query graph
                    if (label == std::get<2>(query_.GetEdgeLabel(u1, u2))){return false;} 
                        
                    bool reversed = false;

                    if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) != treeNode_[u1].backwards_.end())
                    {
                        std::swap(u1, u2); // 
                        std::swap(v1, v2); //
                        reversed = true;
                    }
                    if(std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) != treeNode_[u2].backwards_.end())
                    {

                        // Prune 2: detect if the edge is in the DCS graph
                        // for more complex query graph
                        // detect if DCS has update. If so, that means new FM will occur
                        // return false;


                        bool old_p_d1 = d1[u1][v1], old_p_d2 = d2[u1][v1], old_c_d2 = d2[u2][v2];

                        // Prune 3: detect if the edge is in the FM path
                        // for more complex query graph
                        // If so, that means new FM will occur
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
 * @brief Adds an edge to the data graph and updates the internal structures.
 * 
 * This function adds an edge between two vertices in the data graph and updates
 * the internal data structures to reflect this change. It also enumerates all
 * query edges that match the new edge and updates the corresponding data structures.
 * 
 * @param v1 The source vertex of the edge in the data graph.
 * @param v2 The destination vertex of the edge in the data graph.
 * @param label The label of the edge in the data graph.
 * 
 * @details
 * 1. Adds the edge to the data graph using `data_.AddEdge()`.
 * 2. Enumerates all query edges that match `v1 -> v2` based on vertex labels and edge labels.
 * 3. Updates the internal data structures (`DCS_`, `n1`, `n2`, etc.) based on the new edge.
 */
void Parallel_TurboFlux::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);

    // enumerate all query edges that matches v1 --> v2
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
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
                InsertionTopDown(u1, u2, v1, v2);
            if (old_c_d2)
                InsertionBottomUp(u2, u1, v2, v1);
            if (old_p_d2)
                n2[eidx_[u2][u1]][v2] += 1;

            while (!Q1.empty())
            {
                auto [u_queue, v_queue] = Q1.front();
                Q1.pop();
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    InsertionTopDown(u_queue, u_c_queue, v_queue, v_c_queue);
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
                    InsertionBottomUp(u_queue, u_p_queue, v_queue, v_p_queue);
                    if (reach_time_limit) return;
                }
                for (auto& u_c_queue : treeNode_[u_queue].forwards_)
                for (auto& v_c_queue : DCS_[eidx_[u_queue][u_c_queue]][v_queue])
                {
                    n2[eidx_[u_c_queue][u_queue]][v_c_queue] += 1;
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
    if (max_num_results_ == 0) return;
    
    std::vector<uint> m(query_.NumVertices(), UNMATCHED);

    // enumerate all query edges that matches v1 --> v2
    size_t num_results = 0ul;
    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;
        
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

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; // wcl
            }

            // std::cout << "FindMatches_reduceA" << std::endl;
            Parallel_FindMatches(eidx_[u2][u1], 2, m, num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; // wcl
            }

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) == treeNode_[u1].backwards_.end()
            && std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) == treeNode_[u2].backwards_.end())
        {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; // wcl
            }

            // std::cout << "FindMatches_reduceB" << std::endl;
            Parallel_FindMatches(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, num_results);

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; // wcl
            }

            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
        }
        if (reversed)
        {
            std::swap(u1, u2);
            std::swap(v1, v2);
        }
    }
    }
    END_ENUMERATION:
    num_positive_results_ += num_results;
}


/**
 * @brief Removes an edge from the data graph and updates the internal structures.
 * 
 * This function removes an edge between two vertices in the data graph and updates
 * the internal data structures to reflect this change. It also enumerates all
 * query edges that match the removed edge and updates the corresponding data structures.
 * 
 * @param v1 The source vertex of the edge in the data graph.
 * @param v2 The destination vertex of the edge in the data graph.
 * 
 * @details
 * 1. Removes the edge from the data graph using `data_.RemoveEdge()`.
 * 2. Enumerates all query edges that match `v1 -> v2` based on vertex labels and edge labels.
 * 3. Updates the internal data structures (`DCS_`, `n1`, `n2`, etc.) based on the removed edge.
 */
void Parallel_TurboFlux::RemoveEdge(uint v1, uint v2)
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
        /*auto it = std::lower_bound(DCS_[eidx_[u1][u2]][v1].begin(), DCS_[eidx_[u1][u2]][v1].end(), v2);
        if (
            it == DCS_[eidx_[u1][u2]][v1].end() ||
            *it != v2
        ) continue;*/
        if (std::get<2>(query_.GetEdgeLabel(u1, u2)) == UINT32_MAX) continue;
        
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

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; // wcl
            }

            Parallel_FindMatches_delete(eidx_[u2][u1], 2, m, 
                    num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; // wcl
            }

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;
            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) return;
        }
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2) == treeNode_[u1].backwards_.end()
            && std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1) == treeNode_[u2].backwards_.end())
        {
            m[u1] = v1;
            m[u2] = v2;
            visited_[v1] = true;
            visited_[v2] = true;
            
            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = true;
                local_vec_visited_local[i][v2] = true; // wcl
            }

            Parallel_FindMatches_delete(eidx_[std::min(u1, u2)][std::max(u1, u2)], 2, m, num_results);

            for(size_t i = 0; i< local_vec_visited_local.size(); i++){
                local_vec_visited_local[i][v1] = false;
                local_vec_visited_local[i][v2] = false; // wcl
            }

            visited_[v1] = false;
            visited_[v2] = false;
            m[u1] = UNMATCHED;
            m[u2] = UNMATCHED;
            if (num_results >= max_num_results_) goto END_ENUMERATION;
            if (reach_time_limit) goto END_ENUMERATION;
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


/**
 * @brief Adds a vertex to the data graph and updates the internal structures.
 * 
 * This function adds a vertex to the data graph and updates the internal data structures
 * to reflect this change. It also checks for matches with existing vertices in the query graph
 * and updates the corresponding data structures.
 * 
 * @param id The ID of the vertex being added.
 * @param label The label of the vertex being added.
 */
void Parallel_TurboFlux::AddVertex(uint id, uint label)
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


/**
 * @brief Removes a vertex from the data graph and updates the internal structures.
 * 
 * This function removes a vertex from the data graph and updates the internal data structures
 * to reflect this change. It also checks for matches with existing vertices in the query graph
 * and updates the corresponding data structures.
 * 
 * @param id The ID of the vertex being removed.
 */
void Parallel_TurboFlux::RemoveVertex(uint id)
{
    #pragma omp critical // if we want to remove in parallel, we need to change DCS into concurrent hashmap
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
void Parallel_TurboFlux::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
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
}
