#ifndef CORE_GPU_CANDIDATE_FILTER_H
#define CORE_GPU_CANDIDATE_FILTER_H

#include <cstddef>
#include <cstdint>

// Maximum number of joinability constraints (= max query vertex degree)
static constexpr int MAX_JOIN_CONSTRAINTS = 16;

// A joinability constraint: candidate v must be adjacent to constraint_vertex
// with edge label == expected_elabel
struct JoinConstraint {
    uint32_t constraint_vertex;   // m[u_other] — data vertex that v must be adjacent to
    uint32_t expected_elabel;     // expected edge label
    uint32_t nbr_offset;         // offset into flattened neighbor array
    uint32_t nbr_size;           // number of neighbors of constraint_vertex
};

// Task description for GPU candidate filtering
struct CandidateFilterTask {
    // Candidates: neighbors of m[u_min] in data graph
    const uint32_t* candidates;          // neighbor IDs (sorted)
    const uint32_t* candidate_elabels;   // edge labels to m[u_min]
    uint32_t num_candidates;

    // Expected labels
    uint32_t expected_vlabel;            // query vertex u's label
    uint32_t expected_elabel;            // edge label between u and u_min in query

    // Joinability constraints
    JoinConstraint constraints[MAX_JOIN_CONSTRAINTS];
    uint32_t num_constraints;

    // Flattened neighbor lists for all constraints
    // constraints[i].nbr_offset indexes into this array
    const uint32_t* flat_nbrs;           // flattened sorted neighbor lists
    const uint32_t* flat_elabels;        // flattened edge labels
    uint32_t flat_nbrs_total;            // total size

    // Visited vertices (small set, typically depth <= 10)
    const uint32_t* visited_vertices;
    uint32_t num_visited;
    bool check_visited;
};

class GPUCandidateFilter {
public:
    GPUCandidateFilter();
    ~GPUCandidateFilter();

    // Non-copyable
    GPUCandidateFilter(const GPUCandidateFilter&) = delete;
    GPUCandidateFilter& operator=(const GPUCandidateFilter&) = delete;

    /**
     * Initialize with vertex labels (reuses from Phase 1 if available).
     */
    void Init(const uint32_t* vlabels, size_t num_vertices);

    /**
     * Update vertex labels on GPU.
     */
    void UpdateVertexLabels(const uint32_t* vlabels, size_t num_vertices);

    /**
     * Filter candidates on GPU.
     * @param task  The filtering task description
     * @param valid_indices  Output: indices into candidates[] that passed all checks
     * @return Number of valid candidates
     */
    uint32_t FilterCandidates(
        const CandidateFilterTask& task,
        uint32_t* valid_indices
    );

    void Destroy();
    bool IsInitialized() const { return initialized_; }

private:
    bool initialized_ = false;

    // Device: vertex labels
    uint32_t* d_vlabels_ = nullptr;
    size_t vlabels_capacity_ = 0;

    // Device: candidate buffers (pre-allocated, grown as needed)
    uint32_t* d_candidates_ = nullptr;
    uint32_t* d_candidate_elabels_ = nullptr;
    uint32_t* d_flat_nbrs_ = nullptr;
    uint32_t* d_flat_elabels_ = nullptr;
    uint32_t* d_visited_ = nullptr;
    uint32_t* d_valid_indices_ = nullptr;
    uint32_t* d_valid_count_ = nullptr;   // single uint32 on device

    // Pinned host buffers
    uint32_t* h_valid_indices_ = nullptr;
    uint32_t* h_valid_count_ = nullptr;

    // Capacities
    size_t candidates_cap_ = 0;
    size_t flat_nbrs_cap_ = 0;
    size_t visited_cap_ = 0;

    void EnsureCapacity(size_t num_candidates, size_t flat_nbrs_total, size_t num_visited);
};

#endif // CORE_GPU_CANDIDATE_FILTER_H
