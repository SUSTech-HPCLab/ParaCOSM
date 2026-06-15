// Adapter: bridges main framework's matching interface to provide
// versioned/GPU BFS search for NewSP queries.
// The single-thread path is a no-op (use standalone parallel_newsp for that).
// Versioned/GPU paths use the same order-based search as SymBi/CaLiG/TurboFlux.

#include "matching_executor/Parallel_NewSP/newsp_adapter.h"

#include <omp.h>
#include <atomic>
#include <algorithm>
#include <iostream>
#include <queue>
#include "utils/globals.h"

struct Parallel_NewSP_Adapter::Impl {
    // Placeholder — no NewSP-specific state needed for versioned/GPU
};

Parallel_NewSP_Adapter::Parallel_NewSP_Adapter(
    Graph& query_graph, Graph& data_graph, uint max_num_results,
    bool print_prep, bool print_enum, bool homo,
    size_t num_threads, size_t auto_tuning)
: matching(query_graph, data_graph, max_num_results, print_prep, print_enum, homo)
, NUMTHREAD_(num_threads)
, auto_tuning_(auto_tuning)
, impl_(nullptr)
{
}

Parallel_NewSP_Adapter::~Parallel_NewSP_Adapter()
{
    delete impl_;
}

void Parallel_NewSP_Adapter::Preprocessing()
{
    BuildDAGForVersioned();
    GenerateVersionedMatchingOrders();
}

void Parallel_NewSP_Adapter::InitialMatching() {}

void Parallel_NewSP_Adapter::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);
    // CaLiG-style: no internal index to maintain, count is done by single-thread path
}

void Parallel_NewSP_Adapter::RemoveEdge(uint v1, uint v2)
{
    data_.RemoveEdge(v1, v2);
}

void Parallel_NewSP_Adapter::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
}

void Parallel_NewSP_Adapter::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void Parallel_NewSP_Adapter::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0; num_vertices = 0;
}

bool Parallel_NewSP_Adapter::Classify(uint v1, uint v2, uint label)
{
    // Conservative: all edges are unsafe (ensures correctness)
    return false;
}

size_t Parallel_NewSP_Adapter::EnumerateNewEdge(uint v1, uint v2, uint label, size_t /*thread_id*/)
{
    // Not implemented for NewSP adapter — use versioned mode
    return 0;
}

void Parallel_NewSP_Adapter::UpdateIndexForEdge(uint v1, uint v2, uint label)
{
    // No-op for versioned mode
}

void Parallel_NewSP_Adapter::PrepareBatchEnumeration(size_t num_threads)
{
    const size_t nv = data_.NumVertices();
    const size_t nq = query_.NumVertices();

    // Slot pool
    size_t desired_pool = std::max<size_t>(num_threads * 8, 16);
    if (slot_pool_size_ < desired_pool) {
        size_t old = slot_pool_size_;
        slot_m_.resize(desired_pool, std::vector<uint>(nq, UNMATCHED));
        slot_visited_.resize(desired_pool, std::vector<bool>(nv, false));
        for (size_t s = old; s < desired_pool; s++) {
            free_slots_.push(s);
        }
        slot_pool_size_ = desired_pool;
    }
    for (size_t s = 0; s < slot_pool_size_; s++) {
        if (slot_visited_[s].size() < nv) slot_visited_[s].resize(nv, false);
        if (slot_m_[s].size() < nq) slot_m_[s].assign(nq, UNMATCHED);
    }

    // Read tunables
    if (const char* env = std::getenv("PARACOSM_HOT_SPLIT_THRESHOLD")) {
        size_t v = std::strtoul(env, nullptr, 10);
        if (v > 0) hot_split_threshold_ = v;
    }
    if (const char* env = std::getenv("PARACOSM_HOT_CHUNK_SIZE")) {
        size_t v = std::strtoul(env, nullptr, 10);
        if (v > 0) hot_chunk_size_ = v;
    }
}

// ---------------------------------------------------------------------------
// Versioned DAG + matching order (same pattern as SymBi/CaLiG)
// ---------------------------------------------------------------------------

void Parallel_NewSP_Adapter::BuildDAGForVersioned()
{
    const uint nq = query_.NumVertices();
    eidx_.resize(nq);
    treeNode_.resize(nq);
    serialized_tree_.resize(nq);
    tree_parent_.resize(nq);

    uint edge_pos = 0;
    for (uint i = 0; i < nq; i++) {
        eidx_[i].resize(nq);
        auto& q_nbrs = query_.GetNeighbors(i);
        for (uint j = 0; j < q_nbrs.size(); j++) {
            eidx_[i][q_nbrs[j]] = edge_pos++;
        }
    }

    // Find root by minimum frequency
    std::vector<size_t> freq(nq, 0);
    for (uint u = 0; u < nq; u++) {
        uint ql = query_.GetVertexLabel(u);
        for (uint v = 0; v < data_.NumVertices(); v++) {
            if (data_.GetVertexLabel(v) == ql) freq[u]++;
        }
    }
    q_root_ = std::min_element(freq.begin(), freq.end()) - freq.begin();
    tree_parent_[q_root_] = q_root_;

    // BFS spanning tree
    std::vector<bool> visited(nq, false);
    std::vector<bool> on_tree(nq, false);
    std::queue<uint> bfs_q;
    bfs_q.push(q_root_);
    on_tree[q_root_] = true;
    uint opos = 0;
    serialized_tree_[opos++] = q_root_;

    while (!bfs_q.empty()) {
        uint u = bfs_q.front(); bfs_q.pop();
        visited[u] = true;
        auto& nbrs = query_.GetNeighbors(u);
        for (uint j = 0; j < nbrs.size(); j++) {
            uint uo = nbrs[j];
            if (!on_tree[uo]) {
                bfs_q.push(uo);
                on_tree[uo] = true;
                serialized_tree_[opos++] = uo;
                tree_parent_[uo] = u;
            }
            if (!visited[uo]) {
                treeNode_[u].forwards_.push_back(uo);
            } else {
                treeNode_[u].backwards_.push_back(uo);
            }
            treeNode_[u].neighbors_.push_back(uo);
        }
    }
}

void Parallel_NewSP_Adapter::GenerateVersionedMatchingOrders()
{
    const uint nq = query_.NumVertices();
    const auto& initial_order = serialized_tree_;
    uint first_root_child = treeNode_[q_root_].forwards_.empty()
        ? 0u : treeNode_[q_root_].forwards_.front();

    const size_t num_edge_slots = (query_.NumEdges() + 1) * 2;
    order_vs_.resize(num_edge_slots);
    backward_vs_.resize(num_edge_slots);
    join_check_vs_.resize(num_edge_slots);
    join_check_labels_.resize(num_edge_slots);
    for (size_t i = 0; i < num_edge_slots; i++) {
        order_vs_[i].resize(nq);
        backward_vs_[i].resize(nq);
        join_check_vs_[i].resize(nq);
        join_check_labels_[i].resize(nq);
    }

    // Generate for every directed query edge
    for (uint u = 0; u < nq; u++) {
        auto& q_nbrs = query_.GetNeighbors(u);
        for (uint i = 0; i < q_nbrs.size(); i++) {
            uint u_other = q_nbrs[i];
            const uint eidx = eidx_[u][u_other];
            std::vector<bool> tv(nq, false);
            auto& order = order_vs_[eidx];
            auto& bw = backward_vs_[eidx];
            uint op = 0;
            tv[u] = true; tv[u_other] = true;
            order[op++] = u; order[op++] = u_other;

            uint cur = u_other;
            while (cur != q_root_) {
                uint parent = tree_parent_[cur];
                if (tv[parent]) { cur = parent; continue; }
                tv[parent] = true;
                order[op] = parent; bw[op] = cur;
                auto& pn = query_.GetNeighbors(parent);
                auto& pl = query_.GetNeighborLabels(parent);
                for (uint k = 0; k < pn.size(); k++) {
                    if (tv[pn[k]] && pn[k] != cur) {
                        join_check_vs_[eidx][parent].push_back(pn[k]);
                        join_check_labels_[eidx][parent].push_back(pl[k]);
                    }
                }
                op++; cur = parent;
            }

            for (uint j = 0; j < nq; j++) {
                uint lu = initial_order[j];
                if (!tv[lu]) {
                    order[op] = lu;
                    bw[op] = tree_parent_[lu];
                    auto& ln = query_.GetNeighbors(lu);
                    auto& ll = query_.GetNeighborLabels(lu);
                    for (uint k = 0; k < ln.size(); k++) {
                        if (tv[ln[k]] && ln[k] != tree_parent_[lu]) {
                            join_check_vs_[eidx][lu].push_back(ln[k]);
                            join_check_labels_[eidx][lu].push_back(ll[k]);
                        }
                    }
                    op++; tv[lu] = true;
                }
            }
        }
    }

    // ARTI -> root
    if (!treeNode_[q_root_].forwards_.empty()) {
        const uint ae = eidx_[q_root_][first_root_child];
        std::vector<bool> tv(nq, false);
        auto& order = order_vs_[ae];
        auto& bw = backward_vs_[ae];
        uint op = 0;
        tv[q_root_] = true;
        order[op++] = q_root_;
        for (uint j = 0; j < nq; j++) {
            uint lu = initial_order[j];
            if (!tv[lu]) {
                order[op] = lu; bw[op] = tree_parent_[lu];
                auto& ln = query_.GetNeighbors(lu);
                auto& ll = query_.GetNeighborLabels(lu);
                for (uint k = 0; k < ln.size(); k++) {
                    if (tv[ln[k]] && ln[k] != tree_parent_[lu]) {
                        join_check_vs_[ae][lu].push_back(ln[k]);
                        join_check_labels_[ae][lu].push_back(ll[k]);
                    }
                }
                op++; tv[lu] = true;
            }
        }
    }

    // GPU orders
    gpu_orders_.clear();
    for (uint u = 0; u < nq; u++) {
        auto& q_nbrs = query_.GetNeighbors(u);
        for (uint i = 0; i < q_nbrs.size(); i++) {
            uint uo = q_nbrs[i];
            if (u > uo) continue;
            gpu_orders_.push_back(order_vs_[eidx_[u][uo]]);
        }
    }
}

// ---------------------------------------------------------------------------
// Slot pool
// ---------------------------------------------------------------------------

size_t Parallel_NewSP_Adapter::AcquireSlot()
{
    size_t s;
    while (!free_slots_.try_pop(s)) {
        #pragma omp taskyield
    }
    return s;
}

void Parallel_NewSP_Adapter::ReleaseSlot(size_t s)
{
    free_slots_.push(s);
}

// ---------------------------------------------------------------------------
// Versioned search (same as SymBi/CaLiG pattern)
// ---------------------------------------------------------------------------

static inline bool NewSPVersionedJoinCheck(const Graph& data, uint u_b_data_v, uint v,
                                            uint expected_elabel, uint max_ts)
{
    const auto& nbrs = data.GetNeighbors(u_b_data_v);
    auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
    if (it == nbrs.end() || *it != v) return false;
    size_t idx = static_cast<size_t>(it - nbrs.begin());
    if (data.GetNeighborLabels(u_b_data_v)[idx] != expected_elabel) return false;
    if (data.GetEdgeTimestamps(u_b_data_v)[idx] > max_ts) return false;
    return true;
}

size_t Parallel_NewSP_Adapter::EnumerateNewEdgeVersioned(uint v1, uint v2, uint label,
                                                          size_t /*thread_id*/, uint max_timestamp)
{
    if (max_num_results_ == 0) return 0;

    size_t slot = AcquireSlot();
    auto& m = slot_m_[slot];
    auto& visited = slot_visited_[slot];
    size_t num_results = 0;

    for (uint u1 = 0; u1 < query_.NumVertices(); u1++)
    if (data_.GetVertexLabel(v1) == query_.GetVertexLabel(u1))
    {
    for (uint u2 = 0; u2 < query_.NumVertices(); u2++)
    if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(u2))
    {
        if (std::get<2>(query_.GetEdgeLabel(u1, u2)) != label) continue;

        bool reversed = false;
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2)
            != treeNode_[u1].backwards_.end())
        {
            std::swap(u1, u2); std::swap(v1, v2); reversed = true;
        }

        auto run_search = [&](uint order_index) {
            uint u_min = backward_vs_[order_index][2];
            uint vmin = m[u_min];
            const auto& nbrs = data_.GetNeighbors(vmin);
            const size_t n_nbrs = nbrs.size();
            if (n_nbrs < hot_split_threshold_ || !omp_in_parallel()) {
                FindMatches_versioned_v2(order_index, 2, m, visited, num_results, max_timestamp);
                return;
            }
            std::atomic<size_t> sub{0};
            const uint cu1 = u1, cu2 = u2, cv1 = v1, cv2 = v2;
            const size_t chunk = hot_chunk_size_;
            #pragma omp taskgroup
            {
                for (size_t lo = 0; lo < n_nbrs; lo += chunk) {
                    size_t hi = std::min(lo + chunk, n_nbrs);
                    #pragma omp task firstprivate(lo, hi, order_index, cu1, cu2, cv1, cv2, max_timestamp) shared(sub)
                    {
                        if (!reach_time_limit && sub.load(std::memory_order_relaxed) < max_num_results_) {
                            size_t s = AcquireSlot();
                            auto& tm = slot_m_[s]; auto& tv = slot_visited_[s];
                            tm[cu1] = cv1; tm[cu2] = cv2;
                            tv[cv1] = true; tv[cv2] = true;
                            size_t r = 0;
                            FindMatches_versioned_chunk(order_index, 2, tm, tv, r, max_timestamp, lo, hi);
                            tm[cu1] = UNMATCHED; tm[cu2] = UNMATCHED;
                            tv[cv1] = false; tv[cv2] = false;
                            sub.fetch_add(r, std::memory_order_relaxed);
                            ReleaseSlot(s);
                        }
                    }
                }
            }
            num_results += sub.load();
        };

        if (std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1)
            != treeNode_[u2].backwards_.end())
        {
            m[u1] = v1; m[u2] = v2;
            visited[v1] = true; visited[v2] = true;
            run_search(eidx_[u2][u1]);
            visited[v1] = false; visited[v2] = false;
            m[u1] = UNMATCHED; m[u2] = UNMATCHED;
            if (num_results >= max_num_results_ || reach_time_limit) goto DONE_V;
        }
        if (std::find(treeNode_[u1].backwards_.begin(), treeNode_[u1].backwards_.end(), u2)
                == treeNode_[u1].backwards_.end()
            && std::find(treeNode_[u2].backwards_.begin(), treeNode_[u2].backwards_.end(), u1)
                == treeNode_[u2].backwards_.end())
        {
            m[u1] = v1; m[u2] = v2;
            visited[v1] = true; visited[v2] = true;
            run_search(eidx_[std::min(u1, u2)][std::max(u1, u2)]);
            visited[v1] = false; visited[v2] = false;
            m[u1] = UNMATCHED; m[u2] = UNMATCHED;
            if (num_results >= max_num_results_ || reach_time_limit) goto DONE_V;
        }
        if (reversed) { std::swap(u1, u2); std::swap(v1, v2); }
    }
    }
    DONE_V:
    ReleaseSlot(slot);
    return num_results;
}

void Parallel_NewSP_Adapter::FindMatches_versioned_v2(uint order_index, uint depth,
    std::vector<uint>& m, std::vector<bool>& visited,
    size_t& num_results, uint max_ts)
{
    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];
    auto q_elabel = std::get<2>(query_.GetEdgeLabel(u_min, u));
    const uint q_vlabel_u = query_.GetVertexLabel(u);
    const auto& nbrs = data_.GetNeighbors(m[u_min]);
    const auto& nbr_labels = data_.GetNeighborLabels(m[u_min]);
    const auto& nbr_ts = data_.GetEdgeTimestamps(m[u_min]);
    const uint nv = query_.NumVertices();
    const auto& jc_vs = join_check_vs_[order_index][u];
    const auto& jc_labels = join_check_labels_[order_index][u];
    const size_t jc_n = jc_vs.size();

    for (uint idx = 0; idx < nbrs.size(); idx++) {
        if (nbr_ts[idx] > max_ts) continue;
        if (nbr_labels[idx] != q_elabel) continue;
        uint v = nbrs[idx];
        if (data_.GetVertexLabel(v) != q_vlabel_u) continue;
        bool joinable = true;
        for (size_t i = 0; i < jc_n; i++) {
            if (!NewSPVersionedJoinCheck(data_, m[jc_vs[i]], v, jc_labels[i], max_ts)) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;
        if (!homomorphism_ && visited[v]) continue;
        m[u] = v; visited[v] = true;
        if (depth == nv - 1) { num_results++; }
        else { FindMatches_versioned_v2(order_index, depth + 1, m, visited, num_results, max_ts); }
        visited[v] = false; m[u] = UNMATCHED;
        if (num_results >= max_num_results_ || reach_time_limit) return;
    }
}

void Parallel_NewSP_Adapter::FindMatches_versioned_chunk(uint order_index, uint depth,
    std::vector<uint>& m, std::vector<bool>& visited,
    size_t& num_results, uint max_ts, size_t lo, size_t hi)
{
    uint u = order_vs_[order_index][depth];
    uint u_min = backward_vs_[order_index][depth];
    auto q_elabel = std::get<2>(query_.GetEdgeLabel(u_min, u));
    const uint q_vlabel_u = query_.GetVertexLabel(u);
    const auto& nbrs = data_.GetNeighbors(m[u_min]);
    const auto& nbr_labels = data_.GetNeighborLabels(m[u_min]);
    const auto& nbr_ts = data_.GetEdgeTimestamps(m[u_min]);
    const uint nv = query_.NumVertices();
    const auto& jc_vs = join_check_vs_[order_index][u];
    const auto& jc_labels = join_check_labels_[order_index][u];
    const size_t jc_n = jc_vs.size();

    size_t end = std::min<size_t>(hi, nbrs.size());
    for (size_t idx = lo; idx < end; idx++) {
        if (nbr_ts[idx] > max_ts) continue;
        if (nbr_labels[idx] != q_elabel) continue;
        uint v = nbrs[idx];
        if (data_.GetVertexLabel(v) != q_vlabel_u) continue;
        bool joinable = true;
        for (size_t i = 0; i < jc_n; i++) {
            if (!NewSPVersionedJoinCheck(data_, m[jc_vs[i]], v, jc_labels[i], max_ts)) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;
        if (!homomorphism_ && visited[v]) continue;
        m[u] = v; visited[v] = true;
        if (depth == nv - 1) { num_results++; }
        else { FindMatches_versioned_v2(order_index, depth + 1, m, visited, num_results, max_ts); }
        visited[v] = false; m[u] = UNMATCHED;
        if (num_results >= max_num_results_ || reach_time_limit) return;
    }
}
