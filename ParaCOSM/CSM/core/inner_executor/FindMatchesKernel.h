#ifndef FIND_MATCHES_KERNEL_H
#define FIND_MATCHES_KERNEL_H

#include <vector>
#include <cstdint>

// Forward declarations
class Graph;

// Include the header that defines ExtendableVertex
#include "matching_executor/Parallel_SymBi/parallel_symbi.h"

// Function declaration for Parallel_Symbi_FindMatches
void Parallel_Symbi_FindMatches(
    uint depth, 
    std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, 
    size_t& num_results,
    // Additional parameters to replace global variables
    const Graph& query_graph,
    const Graph& data_graph,
    std::vector<std::vector<uint>>& eidx,
    const std::vector<ska::flat_hash_map<uint, std::vector<uint>>>& DCS,
    const std::vector<ska::flat_hash_map<uint, bool>>& d2,
    const std::vector<ska::flat_hash_map<uint, uint>>& n2,
    const std::vector<std::vector<bool>>& local_vec_visited_local,
    bool homomorphism,
    size_t NUMTHREAD,
    std::vector< std::tuple<uint, uint, size_t, std::vector<uint>,
        std::vector<Parrllel_SymBi::ExtendableVertex>,  uint , uint> >& vertex_vector,
    tbb::concurrent_queue<std::tuple<uint, uint, size_t, std::vector<uint>,
    std::vector<Parrllel_SymBi::ExtendableVertex>,  uint , uint> >& job_queue,
    const std::vector<std::vector<uint>>& treeNode_neighbors
);

// Simplified version without global variables
void Parallel_Symbi_FindMatches_Simplified(
    uint depth, 
    std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, 
    size_t& num_results,
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
);


void Parallel_Symbi_FindMatches_ParaCOSM_Kernel(uint depth, std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results);

void Parallel_Symbi_FindMatches_ParaCOSM_Kernel2(uint depth, std::vector<uint>& m, 
    std::vector<Parrllel_SymBi::ExtendableVertex>& extendable, size_t &num_results);

void Parallel_Graphflow_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);

void Parallel_TurboFlux_FindMatches_ParaCOSM_Kernel(uint order_index, uint depth, std::vector<uint> m, size_t &num_results);


#endif // FIND_MATCHES_KERNEL_H
