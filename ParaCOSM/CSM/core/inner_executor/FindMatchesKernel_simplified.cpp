#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h> // for openmp

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph_storage/graph.h"
#include "matching_executor/Parallel_GraphFlow/parallel_graphflow.h"
#include "matching_executor/Parallel_SymBi/parallel_symbi.h"
#include "core/FindMatchesKernel.h"

// Simplified version of Parallel_Symbi_FindMatches without global variables
void Parallel_Symbi_FindMatches_Simplified(
    uint depth, 
    std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, 
    size_t& num_results,
    // Additional parameters to replace global variables
    const Graph& query_graph,
    const Graph& data_graph,
    const std::vector<std::vector<uint>>& eidx,
    const std::vector<ska::flat_hash_map<uint, std::vector<uint>>>& DCS,
    const std::vector<ska::flat_hash_map<uint, bool>>& d2,
    const std::vector<ska::flat_hash_map<uint, uint>>& n2,
    const std::vector<std::vector<bool>>& local_vec_visited_local,
    bool homomorphism,
    size_t NUMTHREAD,
    const std::vector<std::vector<uint>>& treeNode_neighbors
) {
    uint non_isolate_u = NOT_EXIST, isolate_u = NOT_EXIST;
    uint non_isolate_minE = NOT_EXIST, isolate_minE = NOT_EXIST;
    uint u;

    // extendable.size() == query_graph.NumVertices()
    for (uint i = 0; i < extendable.size(); i++) {
        if (m[i] != UNMATCHED) continue;

        // Check if an extendable query vertex is isolated or not
        if (extendable[i].matched_nbrs == query_graph.GetNeighbors(i).size()) {
            if (extendable[i].E < isolate_minE) {
                isolate_minE = extendable[i].E;
                isolate_u = i;
            }
        } else {
            if (extendable[i].E < non_isolate_minE) {
                non_isolate_minE = extendable[i].E;
                non_isolate_u = i;
            }
        }
    }

    if (non_isolate_minE == NOT_EXIST)
        u = isolate_u;
    else
        u = non_isolate_u;
        
    uint u_min = extendable[u].u_min;
    extendable[u] = {};

    // Enumerate each neighbor of m[u_min]
    // Use the flat_hash_map access pattern
    auto it = DCS[eidx[u_min][u]].find(m[u_min]);
    if (it == DCS[eidx[u_min][u]].end()) {
        return; // No candidates found
    }
    
    const auto& candidates = it->second;
    size_t total_size = candidates.size();
     
    std::vector<size_t> local_num_result(NUMTHREAD, 0);

    for (size_t v_idx = 0; v_idx < total_size; ++v_idx) {
        auto v = candidates[v_idx];

        // 1. Check index
        auto d2_it = d2[u].find(v);
        if (d2_it == d2[u].end() || !d2_it->second) continue;

        // 2. Check if joinable
        bool joinable = true;
        for (auto& u_other: treeNode_neighbors[u]) {
            if (m[u_other] != UNMATCHED || u_other == u_min) continue;

            auto dcs_it = DCS[eidx[u_other][u]].find(m[u_other]);
            if (dcs_it == DCS[eidx[u_other][u]].end()) {
                joinable = false;
                break;
            }
            
            const auto& neighbor_candidates = dcs_it->second;
            auto it = std::lower_bound(neighbor_candidates.begin(), 
                            neighbor_candidates.end(), v);
            if (it == neighbor_candidates.end() || *it != v) {
                joinable = false;
                break;
            }
        }

        if (!joinable) continue;

        // 3. Check if visited
        size_t thread_id = omp_get_thread_num();
        if (!homomorphism && local_vec_visited_local[thread_id][v]) continue;

        // 4. Add a vertex mapping
        m[u] = v;

        std::vector<Parrllel_SymBi::ExtendableVertex> temp_extendable(extendable);
        for (auto& u_other: treeNode_neighbors[u]) {
            if (m[u_other] != UNMATCHED) continue;

            auto n2_it = n2[eidx[u][u_other]].find(m[u]);
            if (n2_it != n2[eidx[u][u_other]].end() && 
                n2_it->second < temp_extendable[u_other].E) {
                temp_extendable[u_other].E = n2_it->second;
                temp_extendable[u_other].u_min = u;
            }
            temp_extendable[u_other].matched_nbrs++;
        }

        if (depth == query_graph.NumVertices() - 1) { // match complete
            num_results++;
        } else {
            // Recursive call for deeper matching
            // You can implement this part based on your needs
            // For now, we'll just increment depth and continue
            // uint depth2 = depth + 1;

            // ParaCOSM_Kernel2
            // Parallel_Symbi_FindMatches_ParaCOSM_Kernel2(depth + 1, m, temp_extendable, num_results);

            // Continue with the next level of matching
            // This is a simplified version - you may need to expand this
        }

        // Backtrack
        m[u] = UNMATCHED;
    }
}
