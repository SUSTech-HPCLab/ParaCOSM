#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include <omp.h>
#include "utils/types.h"
#include "utils/utils.h"
#include "graph_storage/graph.h"

Graph::Graph()
: edge_count_(0)
, vlabel_count_(0)
, elabel_count_(0)
, neighbors_{}
, elabels_{}
, hash_adj_{}
, updates_{}
, updates_vec_{}
, vlabels_{}
{}

void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
        hash_adj_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    
    vlabel_count_ = std::max(vlabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();
    hash_adj_[id].clear();
}

void Graph::AddEdge(uint v1, uint v2, uint label)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;
    
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].emplace(lower, v2);
    elabels_[v1].emplace(elabels_[v1].begin() + dis, label);
    if (!edge_timestamps_.empty()) {
        edge_timestamps_[v1].emplace(edge_timestamps_[v1].begin() + dis, 0u);
    }
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].emplace(lower, v1);
    elabels_[v2].emplace(elabels_[v2].begin() + dis, label);
    if (!edge_timestamps_.empty()) {
        edge_timestamps_[v2].emplace(edge_timestamps_[v2].begin() + dis, 0u);
    }

    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);

    // Maintain hash index for high-degree vertices
    auto maybe_hash_insert = [&](uint v, uint nbr, uint lbl) {
        if (neighbors_[v].size() > HASH_DEGREE_THRESHOLD) {
            if (hash_adj_[v].empty()) {
                // First time crossing threshold: bulk-build from vectors
                hash_adj_[v].reserve(neighbors_[v].size() * 2);
                for (size_t i = 0; i < neighbors_[v].size(); i++)
                    hash_adj_[v][neighbors_[v][i]] = elabels_[v][i];
            } else {
                hash_adj_[v][nbr] = lbl;
            }
        }
    };
    maybe_hash_insert(v1, v2, label);
    maybe_hash_insert(v2, v1, label);
}

void Graph::AddEdgeVersioned(uint v1, uint v2, uint label, uint timestamp)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;
    
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].emplace(lower, v2);
    elabels_[v1].emplace(elabels_[v1].begin() + dis, label);
    edge_timestamps_[v1].emplace(edge_timestamps_[v1].begin() + dis, timestamp);
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].emplace(lower, v1);
    elabels_[v2].emplace(elabels_[v2].begin() + dis, label);
    edge_timestamps_[v2].emplace(edge_timestamps_[v2].begin() + dis, timestamp);

    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
}

// Parallel batch insertion of versioned edges.
// Strategy: bucket inserts per source vertex (sequential scatter, O(E)),
// then per-vertex sorted merge in parallel (independent buckets).
// Note: hash_adj_ is NOT maintained (matches AddEdgeVersioned behavior).
size_t Graph::AddEdgesVersionedBatch(const std::vector<VersionedEdgeBatch>& edges,
                                     size_t num_threads)
{
    if (edges.empty()) return 0;
    const size_t N = neighbors_.size();
    const bool has_ts = !edge_timestamps_.empty();

    // ---- Phase 1: per-vertex insertion counts ----
    std::vector<size_t> off(N + 1, 0);
    for (const auto& e : edges) {
        off[e.v1 + 1]++;
        off[e.v2 + 1]++;
    }
    for (size_t v = 0; v < N; v++) off[v + 1] += off[v];

    // ---- Phase 2: scatter into flat per-vertex bucket ----
    struct Ins { uint nbr, label, ts; };
    std::vector<Ins> flat(off[N]);
    std::vector<size_t> cur = off;  // cursors
    for (const auto& e : edges) {
        flat[cur[e.v1]++] = Ins{e.v2, e.label, e.timestamp};
        flat[cur[e.v2]++] = Ins{e.v1, e.label, e.timestamp};
    }

    // ---- Phase 3: per-vertex parallel sorted merge ----
    std::atomic<size_t> total_added_endpoints{0};
    uint emax = elabel_count_;

    #pragma omp parallel for schedule(dynamic, 64) num_threads(num_threads) reduction(max:emax)
    for (size_t v = 0; v < N; v++) {
        size_t k = off[v + 1] - off[v];
        if (k == 0) continue;

        Ins* ins = flat.data() + off[v];
        std::sort(ins, ins + k, [](const Ins& a, const Ins& b) { return a.nbr < b.nbr; });

        // Drop duplicates within `ins` itself (keep earliest timestamp).
        size_t k_uniq = 0;
        for (size_t a = 0; a < k; a++) {
            if (a > 0 && ins[a].nbr == ins[a - 1].nbr) continue;
            ins[k_uniq++] = ins[a];
        }
        k = k_uniq;

        auto& nbrs = neighbors_[v];
        auto& elab = elabels_[v];
        const size_t m = nbrs.size();

        std::vector<uint> new_nbrs; new_nbrs.reserve(m + k);
        std::vector<uint> new_elab; new_elab.reserve(m + k);
        std::vector<uint> new_ts;   if (has_ts) new_ts.reserve(m + k);
        const std::vector<uint>* old_ts = has_ts ? &edge_timestamps_[v] : nullptr;

        size_t i = 0, j = 0, added_v = 0;
        while (i < m && j < k) {
            if (nbrs[i] < ins[j].nbr) {
                new_nbrs.push_back(nbrs[i]);
                new_elab.push_back(elab[i]);
                if (has_ts) new_ts.push_back((*old_ts)[i]);
                i++;
            } else if (nbrs[i] > ins[j].nbr) {
                new_nbrs.push_back(ins[j].nbr);
                new_elab.push_back(ins[j].label);
                if (has_ts) new_ts.push_back(ins[j].ts);
                if (ins[j].label + 1 > emax) emax = ins[j].label + 1;
                added_v++;
                j++;
            } else {
                // duplicate: keep existing
                new_nbrs.push_back(nbrs[i]);
                new_elab.push_back(elab[i]);
                if (has_ts) new_ts.push_back((*old_ts)[i]);
                i++; j++;
            }
        }
        while (i < m) {
            new_nbrs.push_back(nbrs[i]);
            new_elab.push_back(elab[i]);
            if (has_ts) new_ts.push_back((*old_ts)[i]);
            i++;
        }
        while (j < k) {
            new_nbrs.push_back(ins[j].nbr);
            new_elab.push_back(ins[j].label);
            if (has_ts) new_ts.push_back(ins[j].ts);
            if (ins[j].label + 1 > emax) emax = ins[j].label + 1;
            added_v++;
            j++;
        }

        nbrs = std::move(new_nbrs);
        elab = std::move(new_elab);
        if (has_ts) edge_timestamps_[v] = std::move(new_ts);

        total_added_endpoints.fetch_add(added_v, std::memory_order_relaxed);
    }

    // Each undirected edge contributes to 2 endpoints.
    size_t n_added = total_added_endpoints.load() / 2;
    edge_count_ += static_cast<uint>(n_added);
    elabel_count_ = emax;
    return n_added;
}

void Graph::InitTimestamps()
{
    edge_timestamps_.resize(neighbors_.size());
    for (size_t v = 0; v < neighbors_.size(); v++) {
        edge_timestamps_[v].assign(neighbors_[v].size(), 0u);
    }
}

void Graph::ClearTimestamps()
{
    edge_timestamps_.clear();
    edge_timestamps_.shrink_to_fit();
}

const std::vector<uint>& Graph::GetEdgeTimestamps(uint v) const
{
    return edge_timestamps_[v];
}

uint Graph::GetEdgeTimestamp(uint v1, uint v2) const
{
    if (edge_timestamps_.empty()) return 0;
    // Use the smaller-degree side for binary search
    uint src = (GetDegree(v1) < GetDegree(v2)) ? v1 : v2;
    uint tgt = (src == v1) ? v2 : v1;
    auto it = std::lower_bound(neighbors_[src].begin(), neighbors_[src].end(), tgt);
    if (it != neighbors_[src].end() && *it == tgt) {
        size_t idx = std::distance(neighbors_[src].begin(), it);
        return edge_timestamps_[src][idx];
    }
    return 0;
}

void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);
    elabels_[v1].erase(elabels_[v1].begin() + std::distance(neighbors_[v1].begin(), lower));

    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + std::distance(neighbors_[v2].begin(), lower));

    edge_count_--;

    // Maintain hash index
    auto maybe_hash_erase = [&](uint v, uint nbr) {
        if (!hash_adj_[v].empty()) {
            hash_adj_[v].erase(nbr);
            if (neighbors_[v].size() <= HASH_DEGREE_THRESHOLD)
                hash_adj_[v].clear();
        }
    };
    maybe_hash_erase(v1, v2);
    maybe_hash_erase(v2, v1);
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}

const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    const std::vector<uint> *nbrs;
    const std::vector<uint> *elabel;
    uint other;
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }
    
    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return {v1_label, v2_label, e_label};
        }
    }
    return {v1_label, v2_label, -1};
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)
    if (GetVertexLabel(i) != NOT_EXIST)
    {
        std::queue<uint> bfs_queue;
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}

void Graph::LoadFromFile(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            AddEdge(from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    // updates_.clear();
    updates_vec_.clear();
    // updates_.reserve(10000);
    updates_vec_.reserve(10000);


    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label);
            updates_vec_.emplace_back('v', type == "v", vertex_id, 0u, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            updates_.emplace('e', type == "e", from_id, to_id, label);
            updates_vec_.emplace_back('e', type == "e", from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() <<
        "\n# edges = " << NumEdges() << std::endl;
}