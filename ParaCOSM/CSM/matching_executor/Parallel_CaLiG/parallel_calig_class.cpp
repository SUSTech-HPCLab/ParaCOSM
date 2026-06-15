#include <algorithm>
#include <atomic>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>

#include <omp.h>
#include <tbb/task_group.h>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph_storage/graph.h"
#include "matching_executor/Parallel_CaLiG/parallel_calig_class.h"

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

Parallel_CaLiG::Parallel_CaLiG(
    Graph& query_graph, Graph& data_graph, uint max_num_results,
    bool print_prep, bool print_enum, bool homo,
    size_t num_threads, size_t auto_tuning)
: matching(query_graph, data_graph, max_num_results, print_prep, print_enum, homo)
, NUMTHREAD_(num_threads)
, auto_tuning_(auto_tuning)
{
    // Initialize eidx_ for versioned path
    uint edge_pos = 0;
    eidx_.resize(query_.NumVertices());
    for (uint i = 0; i < query_.NumVertices(); i++) {
        eidx_[i].resize(query_.NumVertices());
        auto& q_nbrs = query_.GetNeighbors(i);
        for (uint j = 0; j < q_nbrs.size(); j++) {
            eidx_[i][q_nbrs[j]] = edge_pos++;
        }
    }

    // Initialize versioned order structures
    const size_t num_edge_slots = (query_.NumEdges() + 1) * 2;
    order_vs_.resize(num_edge_slots);
    backward_vs_.resize(num_edge_slots);
    join_check_vs_.resize(num_edge_slots);
    join_check_labels_.resize(num_edge_slots);
    for (size_t i = 0; i < num_edge_slots; i++) {
        order_vs_[i].resize(query_.NumVertices());
        backward_vs_[i].resize(query_.NumVertices());
        join_check_vs_[i].resize(query_.NumVertices());
        join_check_labels_[i].resize(query_.NumVertices());
    }

    treeNode_.resize(query_.NumVertices());
    serialized_tree_.resize(query_.NumVertices());
    tree_parent_.resize(query_.NumVertices());
}

// ---------------------------------------------------------------------------
// Preprocessing
// ---------------------------------------------------------------------------

void Parallel_CaLiG::Preprocessing()
{
    double t0 = omp_get_wtime();
    AnalyzeQuery();
    GenerateMatchingOrders();
    double t1 = omp_get_wtime();
    std::cout << "CaLiG: AnalyzeQuery + GenerateMO: " << (t1 - t0) << " s" << std::endl;

    ConstructCandidates();
    double t2 = omp_get_wtime();
    std::cout << "CaLiG: ConstructCandidates: " << (t2 - t1) << " s" << std::endl;

    StaticFilter();
    double t3 = omp_get_wtime();
    std::cout << "CaLiG: StaticFilter: " << (t3 - t2) << " s" << std::endl;

    // Build DAG + versioned matching orders for versioned/GPU path
    BuildDAGForVersioned();
    GenerateVersionedMatchingOrders();
    double t4 = omp_get_wtime();
    std::cout << "CaLiG: BuildDAG + VersionedOrders: " << (t4 - t3) << " s" << std::endl;
}

// ---------------------------------------------------------------------------
// Query analysis
// ---------------------------------------------------------------------------

void Parallel_CaLiG::AnalyzeQuery()
{
    const uint Q_size = query_.NumVertices();
    labels_.clear();
    rep_nei_.resize(Q_size);

    for (uint ui = 0; ui < Q_size; ui++) {
        uint lb = query_.GetVertexLabel(ui);
        labels_[lb].push_back(ui);

        std::unordered_map<uint, vec> rep_nei;
        const auto& q_nbrs = query_.GetNeighbors(ui);
        for (uint j = 0; j < q_nbrs.size(); j++) {
            uint uj = q_nbrs[j];
            uint uj_label = query_.GetVertexLabel(uj);
            rep_nei[uj_label].push_back(uj);
        }
        for (auto& [label, vs] : rep_nei) {
            if (vs.size() > 1) {
                rep_nei_[ui][label] = vs;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Core-shell matching order generation
// ---------------------------------------------------------------------------

static bool isInVecU(uint a, std::vector<uint>& v) {
    for (size_t i = 0; i < v.size(); i++) {
        if (a == v[i]) return true;
    }
    return false;
}

void Parallel_CaLiG::GenerateMatchingOrders()
{
    const uint Q_size = query_.NumVertices();

    auto insertNei = [&](u_set& temp, uint i, uint ui, uint uj, Split& sp) {
        const auto& nbrs = query_.GetNeighbors(i);
        for (uint j = 0; j < nbrs.size(); j++) {
            uint n = nbrs[j];
            if (n != ui && n != uj && !isInVecU(n, sp.core) && !isInVecU(n, sp.shell)) {
                temp.insert(n);
            }
        }
    };

    auto neiAllCore = [&](uint i, Split& sp, uint ui, uint uj) -> bool {
        const auto& nbrs = query_.GetNeighbors(i);
        for (uint j = 0; j < nbrs.size(); j++) {
            uint n = nbrs[j];
            if (n != ui && n != uj && !isInVecU(n, sp.core)) {
                return false;
            }
        }
        return true;
    };

    for (uint ui = 0; ui < Q_size; ui++) {
        std::unordered_map<uint, Split> ui_first;
        const auto& q_nbrs = query_.GetNeighbors(ui);

        for (uint j = 0; j < q_nbrs.size(); j++) {
            uint uj = q_nbrs[j];
            Split sp;
            u_set temp;
            insertNei(temp, ui, ui, uj, sp);
            insertNei(temp, uj, ui, uj, sp);

            while (!temp.empty()) {
                u_set temp2 = temp;
                // Stage 1: move all-core-neighbors to shell
                for (auto& ti : temp2) {
                    if (neiAllCore(ti, sp, ui, uj)) {
                        u_set nei_info;
                        const auto& tn = query_.GetNeighbors(ti);
                        for (uint k = 0; k < tn.size(); k++) nei_info.insert(tn[k]);

                        uint v_core = 0;
                        for (auto& c : sp.core) {
                            if (nei_info.count(c)) v_core = c;
                        }
                        sp.c_s_nei[v_core].insert(static_cast<uint>(sp.shell.size()));
                        sp.shell.push_back(ti);
                        sp.shell_nei.push_back(nei_info);
                        temp.erase(ti);
                    }
                }
                // Stage 2: move highest-degree to core
                if (!temp.empty()) {
                    uint nei_num = 0;
                    uint core_v = *temp.begin();
                    for (auto& t : temp) {
                        uint deg = 0;
                        const auto& tn = query_.GetNeighbors(t);
                        for (uint k = 0; k < tn.size(); k++) {
                            if (tn[k] != ui && tn[k] != uj) deg++;
                        }
                        if (deg > nei_num) {
                            nei_num = deg;
                            core_v = t;
                        }
                    }
                    sp.core.push_back(core_v);
                    u_set nei_info;
                    const auto& cn = query_.GetNeighbors(core_v);
                    for (uint k = 0; k < cn.size(); k++) {
                        uint n = cn[k];
                        if (n == ui || n == uj || isInVecU(n, sp.core)) {
                            nei_info.insert(n);
                        }
                    }
                    sp.core_nei.push_back(nei_info);
                    insertNei(temp, core_v, ui, uj, sp);
                    temp.erase(core_v);
                }
            }
            ui_first[uj] = sp;
        }
        matching_order_[ui] = ui_first;
    }
}

// ---------------------------------------------------------------------------
// Candidate construction
// ---------------------------------------------------------------------------

void Parallel_CaLiG::ConstructCandidates()
{
    const size_t nv = data_.NumVertices();
    cidx_.resize(nv);

    for (uint vi = 0; vi < nv; vi++) {
        uint lb = data_.GetVertexLabel(vi);
        if (lb == NOT_EXIST) continue;

        auto it = labels_.find(lb);
        if (it == labels_.end()) continue;

        const auto& d_nbrs = data_.GetNeighbors(vi);

        for (uint ui : it->second) {
            std::unordered_map<uint, u_set> ui_cand;

            // Initialize empty sets for ALL query neighbors of ui
            const auto& q_nbrs = query_.GetNeighbors(ui);
            for (uint k = 0; k < q_nbrs.size(); k++) {
                ui_cand[q_nbrs[k]] = u_set();
            }

            // Populate candidate sets
            for (uint j = 0; j < d_nbrs.size(); j++) {
                uint vj = d_nbrs[j];
                uint vj_lb = data_.GetVertexLabel(vj);
                if (vj_lb == NOT_EXIST) continue;

                for (uint k = 0; k < q_nbrs.size(); k++) {
                    uint uj = q_nbrs[k];
                    if (query_.GetVertexLabel(uj) == vj_lb) {
                        ui_cand[uj].insert(vj);
                    }
                }
            }
            cidx_[vi].cand[ui] = ui_cand;
            cidx_[vi].LI[ui] = true;
        }
    }
}

// ---------------------------------------------------------------------------
// Static filter
// ---------------------------------------------------------------------------

void Parallel_CaLiG::StaticFilter()
{
    for (uint vi = 0; vi < data_.NumVertices(); vi++) {
        TurnOff(vi);
    }
}

// ---------------------------------------------------------------------------
// Index maintenance helpers
// ---------------------------------------------------------------------------

bool Parallel_CaLiG::TryNei(uint th, uint vi, uint ui, u_set& used, vec& to_check)
{
    if (th == to_check.size()) return true;
    uint uj = to_check[th];
    for (auto& vj : cidx_[vi].cand[ui][uj]) {
        if (data_.GetVertexLabel(vj) != query_.GetVertexLabel(uj)) continue;
        if (used.find(vj) == used.end()) {
            used.insert(vj);
            if (TryNei(th + 1, vi, ui, used, to_check)) return true;
            used.erase(vj);
        }
    }
    return false;
}

bool Parallel_CaLiG::CheckNei(uint vi, uint ui)
{
    const auto& q_nbrs = query_.GetNeighbors(ui);
    for (uint j = 0; j < q_nbrs.size(); j++) {
        uint uj = q_nbrs[j];
        auto it = cidx_[vi].cand[ui].find(uj);
        if (it == cidx_[vi].cand[ui].end()) return false;
        if (it->second.empty()) return false;
    }
    auto it = rep_nei_.size() > ui ? rep_nei_[ui].begin() : rep_nei_[ui].end();
    for (auto& [label, rep_vec] : rep_nei_[ui]) {
        u_set used;
        vec to_check = rep_vec;
        if (!TryNei(0, vi, ui, used, to_check)) return false;
    }
    return true;
}

void Parallel_CaLiG::DeleteAndCheck(uint ui, uint vi)
{
    for (auto& [uj, vj_set] : cidx_[vi].cand[ui]) {
        for (auto& vj : vj_set) {
            if (data_.GetVertexLabel(vj) != query_.GetVertexLabel(uj)) continue;
            if (cidx_[vj].LI[uj]) {
                cidx_[vj].cand[uj][ui].erase(vi);
                if (cidx_[vj].cand[uj][ui].empty()) {
                    cidx_[vj].LI[uj] = false;
                    DeleteAndCheck(uj, vj);
                } else {
                    uint lb = query_.GetVertexLabel(ui);
                    if (rep_nei_[uj].find(lb) != rep_nei_[uj].end()) {
                        auto& rep = rep_nei_[uj][lb];
                        u_set used;
                        if (!TryNei(0, vj, uj, used, rep)) {
                            cidx_[vj].LI[uj] = false;
                            DeleteAndCheck(uj, vj);
                        }
                    }
                }
            }
        }
    }
}

void Parallel_CaLiG::TurnOff(uint vi)
{
    // Iterate over a copy of keys to avoid mutation during iteration
    std::vector<uint> keys;
    for (auto& [ui, _] : cidx_[vi].cand) keys.push_back(ui);

    for (uint ui : keys) {
        if (data_.GetVertexLabel(vi) != query_.GetVertexLabel(ui)) {
            cidx_[vi].cand.erase(ui);
            cidx_[vi].LI.erase(ui);
            continue;
        }
        if (!cidx_[vi].LI[ui]) continue;
        if (!CheckNei(vi, ui)) {
            cidx_[vi].LI[ui] = false;
            DeleteAndCheck(ui, vi);
        }
    }
}

void Parallel_CaLiG::AddAndCheck(uint ui, uint vi, vec& temp_v, vec& temp_u)
{
    for (auto& [uj, vj_set] : cidx_[vi].cand[ui]) {
        for (auto& vj : vj_set) {
            if (data_.GetVertexLabel(vj) != query_.GetVertexLabel(uj)) continue;
            cidx_[vj].cand[uj][ui].insert(vi);
            if (!cidx_[vj].LI[uj]) {
                if (CheckNei(vj, uj)) {
                    cidx_[vj].LI[uj] = true;
                    AddAndCheck(uj, vj, temp_v, temp_u);
                } else {
                    temp_v.push_back(vj);
                    temp_u.push_back(uj);
                }
            }
        }
    }
}

void Parallel_CaLiG::TurnOnProcess(uint v1, uint v2)
{
    vec temp_v, temp_u;

    // Keys snapshot to iterate safely
    std::vector<uint> keys;
    for (auto& [ui, _] : cidx_[v1].cand) keys.push_back(ui);

    for (uint ui : keys) {
        if (data_.GetVertexLabel(v1) != query_.GetVertexLabel(ui)) {
            cidx_[v1].cand.erase(ui);
            cidx_[v1].LI.erase(ui);
            continue;
        }
        for (auto& [uj, _] : cidx_[v1].cand[ui]) {
            if (data_.GetVertexLabel(v2) == query_.GetVertexLabel(uj)) {
                cidx_[v1].cand[ui][uj].insert(v2);
                if (cidx_[v1].LI[ui]) {
                    cidx_[v2].cand[uj][ui].insert(v1);
                    AddAndCheck(ui, v1, temp_v, temp_u);
                    for (size_t i = 0; i < temp_v.size(); i++) {
                        if (!cidx_[temp_v[i]].LI[temp_u[i]]) {
                            DeleteAndCheck(temp_u[i], temp_v[i]);
                        }
                    }
                    temp_v.clear();
                    temp_u.clear();
                } else if (CheckNei(v1, ui)) {
                    cidx_[v1].LI[ui] = true;
                    AddAndCheck(ui, v1, temp_v, temp_u);
                    for (size_t i = 0; i < temp_v.size(); i++) {
                        if (!cidx_[temp_v[i]].LI[temp_u[i]]) {
                            DeleteAndCheck(temp_u[i], temp_v[i]);
                        }
                    }
                    temp_v.clear();
                    temp_u.clear();
                }
            }
        }
    }
}

void Parallel_CaLiG::TurnOffProcess(uint v1, uint v2)
{
    std::vector<uint> keys;
    for (auto& [ui, _] : cidx_[v1].cand) keys.push_back(ui);

    for (uint ui : keys) {
        if (data_.GetVertexLabel(v1) != query_.GetVertexLabel(ui)) {
            cidx_[v1].cand.erase(ui);
            cidx_[v1].LI.erase(ui);
            continue;
        }
        for (auto& [uj, _] : cidx_[v1].cand[ui]) {
            if (query_.GetVertexLabel(uj) == data_.GetVertexLabel(v2)) {
                cidx_[v1].cand[ui][uj].erase(v2);
                cidx_[v2].cand[uj][ui].erase(v1);
            }
        }
    }
    TurnOff(v1);
    TurnOff(v2);
}

// ---------------------------------------------------------------------------
// CaLiG search (core-shell enumeration)
// ---------------------------------------------------------------------------

bool Parallel_CaLiG::ShellCand(std::vector<u_set>& result,
    std::unordered_map<uint, uint>& m, const vec& s,
    const std::vector<u_set>& s_n, vec& used)
{
    result.resize(s.size());
    for (size_t i = 0; i < s.size(); i++) {
        auto p_ui = s_n[i].begin();
        if (s_n[i].size() > 1) {
            uint ui_l = *p_ui; ++p_ui;
            uint ui2 = *p_ui;
            // intersection
            result[i].clear();
            auto& s1 = cidx_[m[ui_l]].cand[ui_l][s[i]];
            auto& s2 = cidx_[m[ui2]].cand[ui2][s[i]];
            for (auto& x : s1) {
                if (s2.count(x)) result[i].insert(x);
            }
            for (++p_ui; p_ui != s_n[i].end(); ++p_ui) {
                uint ui3 = *p_ui;
                u_set tmp;
                for (auto& x : result[i]) {
                    if (cidx_[m[ui3]].cand[ui3][s[i]].count(x)) tmp.insert(x);
                }
                result[i] = std::move(tmp);
            }
            if (result[i].empty()) return true;
        } else {
            uint ui = *p_ui;
            result[i] = cidx_[m[ui]].cand[ui][s[i]];
        }
        for (auto& vth : used) {
            result[i].erase(vth);
        }
        if (result[i].empty()) return true;
    }
    return false;
}

int Parallel_CaLiG::NumAdd(int th, const std::vector<u_set>& cand, u_set& used)
{
    int result = 0;
    if (th == static_cast<int>(cand.size()) - 1) {
        int del = 0;
        for (auto& vth : used) {
            if (cand[th].count(vth)) del++;
        }
        return static_cast<int>(cand[th].size()) - del;
    }
    for (auto& vth : cand[th]) {
        if (used.find(vth) == used.end()) {
            used.insert(vth);
            result += NumAdd(th + 1, cand, used);
            used.erase(vth);
        }
    }
    return result;
}

bool Parallel_CaLiG::NotExit(uint shell_v, const u_set& nei, std::unordered_map<uint, uint>& m)
{
    if (nei.size() == 1) return false;
    uint n = *nei.begin();
    for (auto& cand_v : cidx_[m[n]].cand[n][shell_v]) {
        uint count = 0;
        for (auto& ni : nei) {
            count++;
            if (count == 1) continue;
            if (cidx_[m[ni]].cand[ni][shell_v].count(cand_v)) {
                if (count == nei.size()) return false;
            } else {
                break;
            }
        }
    }
    return true;
}

int Parallel_CaLiG::SearchCore(int th, std::unordered_map<uint, uint>& m, vec& used,
    const vec& c, const std::vector<u_set>& c_n,
    const vec& s, const std::vector<u_set>& s_n,
    std::unordered_map<uint, u_set>& c2check)
{
    int result = 0;
    if (th == static_cast<int>(c.size())) {
        std::vector<u_set> candidates;
        if (ShellCand(candidates, m, s, s_n, used)) return 0;
        u_set used_v;
        return NumAdd(0, candidates, used_v);
    }
    u_set candidates;
    auto p_ui = c_n[th].begin();
    uint ui = *p_ui;

    if (c_n[th].size() > 1) {
        auto& ref = cidx_[m[ui]].cand[ui][c[th]];
        candidates = ref;
        for (++p_ui; p_ui != c_n[th].end(); ++p_ui) {
            ui = *p_ui;
            u_set tmp;
            for (auto& x : candidates) {
                if (cidx_[m[ui]].cand[ui][c[th]].count(x)) tmp.insert(x);
            }
            candidates = std::move(tmp);
            if (candidates.empty()) return 0;
        }
        for (auto& vi : used) candidates.erase(vi);
        for (auto& vth : candidates) {
            m[c[th]] = vth;
            if (c2check.count(c[th])) {
                bool to_continue = false;
                for (auto& shell_idx : c2check[c[th]]) {
                    if (NotExit(s[shell_idx], s_n[shell_idx], m)) {
                        to_continue = true; break;
                    }
                }
                if (to_continue) continue;
            }
            used.push_back(vth);
            result += SearchCore(th + 1, m, used, c, c_n, s, s_n, c2check);
            used.pop_back();
            if (reach_time_limit) return result;
        }
    } else {
        candidates = cidx_[m[ui]].cand[ui][c[th]];
        for (auto& vth : candidates) {
            if (!isInVecU(vth, used)) {
                m[c[th]] = vth;
                if (c2check.count(c[th])) {
                    bool to_continue = false;
                    for (auto& shell_idx : c2check[c[th]]) {
                        if (NotExit(s[shell_idx], s_n[shell_idx], m)) {
                            to_continue = true; break;
                        }
                    }
                    if (to_continue) continue;
                }
                used.push_back(vth);
                result += SearchCore(th + 1, m, used, c, c_n, s, s_n, c2check);
                used.pop_back();
                if (reach_time_limit) return result;
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// matching interface: AddEdge / RemoveEdge
// ---------------------------------------------------------------------------

void Parallel_CaLiG::AddEdge(uint v1, uint v2, uint label)
{
    // Ensure cidx_ is large enough
    uint mx = std::max(v1, v2);
    if (mx >= cidx_.size()) cidx_.resize(mx + 1);

    data_.AddEdge(v1, v2, label);

    // Update CaLiG index for both directions
    TurnOnProcess(v1, v2);
    TurnOnProcess(v2, v1);

    // Enumerate matches
    size_t num_results = 0;
    for (auto& [ui, li_val] : cidx_[v1].LI) {
        if (!li_val) continue;
        for (auto& [uj, cand_set] : cidx_[v1].cand[ui]) {
            if (cand_set.find(v2) != cand_set.end()) {
                std::unordered_map<uint, uint> m;
                m[ui] = v1; m[uj] = v2;
                vec used_v = {v1, v2};
                num_results += SearchCore(0, m, used_v,
                    matching_order_[ui][uj].core, matching_order_[ui][uj].core_nei,
                    matching_order_[ui][uj].shell, matching_order_[ui][uj].shell_nei,
                    matching_order_[ui][uj].c_s_nei);
            }
        }
    }
    num_positive_results_ += num_results;
}

void Parallel_CaLiG::RemoveEdge(uint v1, uint v2)
{
    if (TurnOffProcessSafe(v1, v2)) {
        data_.RemoveEdge(v1, v2);
        return;
    }

    // Enumerate matches before removal
    size_t num_results = 0;
    for (auto& [ui, li_val] : cidx_[v1].LI) {
        if (!li_val) continue;
        for (auto& [uj, cand_set] : cidx_[v1].cand[ui]) {
            if (cand_set.find(v2) != cand_set.end()) {
                std::unordered_map<uint, uint> m;
                m[ui] = v1; m[uj] = v2;
                vec used_v = {v1, v2};
                num_results += SearchCore(0, m, used_v,
                    matching_order_[ui][uj].core, matching_order_[ui][uj].core_nei,
                    matching_order_[ui][uj].shell, matching_order_[ui][uj].shell_nei,
                    matching_order_[ui][uj].c_s_nei);
            }
        }
    }
    num_negative_results_ += num_results;

    TurnOffProcess(v1, v2);
    data_.RemoveEdge(v1, v2);
}

void Parallel_CaLiG::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
    if (id >= cidx_.size()) cidx_.resize(id + 1);
}

void Parallel_CaLiG::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void Parallel_CaLiG::InitialMatching() {}

void Parallel_CaLiG::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0; num_vertices = 0;
}

// ---------------------------------------------------------------------------
// Safe-edge classification
// ---------------------------------------------------------------------------

bool Parallel_CaLiG::Classify(uint v1, uint v2, uint label)
{
    uint v1_lb = data_.GetVertexLabel(v1);
    uint v2_lb = data_.GetVertexLabel(v2);

    // Check v1's candidates (regardless of LI status)
    for (auto& [ui, cand_map] : cidx_[v1].cand) {
        if (v1_lb != query_.GetVertexLabel(ui)) continue;
        for (auto& [uj, _] : cand_map) {
            if (query_.GetVertexLabel(uj) == v2_lb) {
                return false;  // unsafe
            }
        }
    }
    // Check v2's candidates
    for (auto& [uj, cand_map] : cidx_[v2].cand) {
        if (v2_lb != query_.GetVertexLabel(uj)) continue;
        for (auto& [ui, _] : cand_map) {
            if (query_.GetVertexLabel(ui) == v1_lb) {
                return false;  // unsafe
            }
        }
    }
    return true;  // safe
}

bool Parallel_CaLiG::TurnOnProcessSafe(uint v1, uint v2)
{
    return Classify(v1, v2, 0);
}

bool Parallel_CaLiG::TurnOffProcessSafe(uint v1, uint v2)
{
    return Classify(v1, v2, 0);
}

// ---------------------------------------------------------------------------
// Batch support: EnumerateNewEdge / UpdateIndexForEdge / PrepareBatchEnumeration
// ---------------------------------------------------------------------------

void Parallel_CaLiG::UpdateIndexForEdge(uint v1, uint v2, uint label)
{
    uint mx = std::max(v1, v2);
    if (mx >= cidx_.size()) cidx_.resize(mx + 1);
    TurnOnProcess(v1, v2);
}

void Parallel_CaLiG::PrepareBatchEnumeration(size_t num_threads)
{
    // Slot pool for versioned search
    size_t desired_pool = std::max<size_t>(num_threads * 8, 16);
    if (slot_pool_size_ < desired_pool) {
        size_t old = slot_pool_size_;
        slot_m_.resize(desired_pool);
        slot_used_.resize(desired_pool);
        for (size_t s = old; s < desired_pool; s++) {
            free_slots_.push(s);
        }
        slot_pool_size_ = desired_pool;
    }
}

size_t Parallel_CaLiG::EnumerateNewEdge(uint v1, uint v2, uint label, size_t thread_id)
{
    if (max_num_results_ == 0) return 0;
    size_t num_results = 0;

    for (auto& [ui, li_val] : cidx_[v1].LI) {
        if (!li_val) continue;
        for (auto& [uj, cand_set] : cidx_[v1].cand[ui]) {
            if (cand_set.find(v2) != cand_set.end()) {
                std::unordered_map<uint, uint> m;
                m[ui] = v1; m[uj] = v2;
                vec used_v = {v1, v2};
                num_results += SearchCore(0, m, used_v,
                    matching_order_[ui][uj].core, matching_order_[ui][uj].core_nei,
                    matching_order_[ui][uj].shell, matching_order_[ui][uj].shell_nei,
                    matching_order_[ui][uj].c_s_nei);
                if (num_results >= max_num_results_ || reach_time_limit) return num_results;
            }
        }
    }
    return num_results;
}

// ---------------------------------------------------------------------------
// Versioned search support (DAG-based, same as SymBi/TurboFlux)
// ---------------------------------------------------------------------------

void Parallel_CaLiG::BuildDAGForVersioned()
{
    const uint Q_size = query_.NumVertices();

    // Find root: vertex with smallest label frequency
    uint min_freq = UINT_MAX;
    q_root_ = 0;
    for (uint u = 0; u < Q_size; u++) {
        uint lb = query_.GetVertexLabel(u);
        uint freq = 0;
        for (uint v = 0; v < data_.NumVertices(); v++) {
            if (data_.GetVertexLabel(v) == lb) freq++;
        }
        if (freq < min_freq) {
            min_freq = freq;
            q_root_ = u;
        }
    }

    // BFS spanning tree
    std::vector<bool> visited(Q_size, false);
    std::vector<bool> on_tree(Q_size, false);
    std::queue<uint> bfs_q;
    tree_parent_[q_root_] = q_root_;
    bfs_q.push(q_root_);
    on_tree[q_root_] = true;
    uint pos = 0;
    serialized_tree_[pos++] = q_root_;

    while (!bfs_q.empty()) {
        uint u = bfs_q.front(); bfs_q.pop();
        visited[u] = true;
        const auto& nbrs = query_.GetNeighbors(u);
        for (uint j = 0; j < nbrs.size(); j++) {
            uint u_other = nbrs[j];
            if (!on_tree[u_other]) {
                bfs_q.push(u_other);
                on_tree[u_other] = true;
                serialized_tree_[pos++] = u_other;
                tree_parent_[u_other] = u;
            }
            if (!visited[u_other]) {
                treeNode_[u].forwards_.push_back(u_other);
            } else {
                treeNode_[u].backwards_.push_back(u_other);
            }
            treeNode_[u].neighbors_.push_back(u_other);
        }
    }
}

void Parallel_CaLiG::GenerateVersionedMatchingOrders()
{
    const auto& initial_order = serialized_tree_;
    uint first_root_child = treeNode_[q_root_].forwards_.empty()
        ? 0u : treeNode_[q_root_].forwards_.front();

    // Generate for every directed edge
    for (uint u = 0; u < query_.NumVertices(); u++) {
        const auto& q_nbrs = query_.GetNeighbors(u);
        for (uint i = 0; i < q_nbrs.size(); i++) {
            uint u_other = q_nbrs[i];
            const uint eidx = eidx_[u][u_other];

            std::vector<bool> temp_visited(query_.NumVertices(), false);
            auto& order = order_vs_[eidx];
            auto& backwards = backward_vs_[eidx];
            uint order_pos = 0;

            temp_visited[u] = true;
            temp_visited[u_other] = true;
            order[order_pos++] = u;
            order[order_pos++] = u_other;

            // Walk upward from u_other toward root
            uint cur = u_other;
            while (cur != q_root_) {
                uint parent = tree_parent_[cur];
                if (temp_visited[parent]) { cur = parent; continue; }
                temp_visited[parent] = true;
                order[order_pos] = parent;
                backwards[order_pos] = cur;
                auto& pn = query_.GetNeighbors(parent);
                auto& pl = query_.GetNeighborLabels(parent);
                for (uint k = 0; k < pn.size(); k++) {
                    if (temp_visited[pn[k]] && pn[k] != cur) {
                        join_check_vs_[eidx][parent].push_back(pn[k]);
                        join_check_labels_[eidx][parent].push_back(pl[k]);
                    }
                }
                order_pos++;
                cur = parent;
            }

            // Fill remaining
            for (uint j = 0; j < query_.NumVertices(); j++) {
                uint local_u = initial_order[j];
                if (!temp_visited[local_u]) {
                    order[order_pos] = local_u;
                    backwards[order_pos] = tree_parent_[local_u];
                    auto& q_nbrs_local = query_.GetNeighbors(local_u);
                    auto& q_nbr_labels_local = query_.GetNeighborLabels(local_u);
                    for (uint k = 0; k < q_nbrs_local.size(); k++) {
                        if (temp_visited[q_nbrs_local[k]] && q_nbrs_local[k] != tree_parent_[local_u]) {
                            join_check_vs_[eidx][local_u].push_back(q_nbrs_local[k]);
                            join_check_labels_[eidx][local_u].push_back(q_nbr_labels_local[k]);
                        }
                    }
                    order_pos++;
                    temp_visited[local_u] = true;
                }
            }
        }
    }

    // ARTI→root order
    if (!treeNode_[q_root_].forwards_.empty()) {
        const uint arti_eidx = eidx_[q_root_][first_root_child];
        std::vector<bool> temp_visited(query_.NumVertices(), false);
        auto& order = order_vs_[arti_eidx];
        auto& backwards = backward_vs_[arti_eidx];
        uint order_pos = 0;
        temp_visited[q_root_] = true;
        order[order_pos++] = q_root_;
        for (uint j = 0; j < query_.NumVertices(); j++) {
            uint local_u = initial_order[j];
            if (!temp_visited[local_u]) {
                order[order_pos] = local_u;
                backwards[order_pos] = tree_parent_[local_u];
                auto& q_nbrs_local = query_.GetNeighbors(local_u);
                auto& q_nbr_labels_local = query_.GetNeighborLabels(local_u);
                for (uint k = 0; k < q_nbrs_local.size(); k++) {
                    if (temp_visited[q_nbrs_local[k]] && q_nbrs_local[k] != tree_parent_[local_u]) {
                        join_check_vs_[arti_eidx][local_u].push_back(q_nbrs_local[k]);
                        join_check_labels_[arti_eidx][local_u].push_back(q_nbr_labels_local[k]);
                    }
                }
                order_pos++;
                temp_visited[local_u] = true;
            }
        }
    }

    // Build gpu_orders_
    gpu_orders_.clear();
    for (uint u = 0; u < query_.NumVertices(); u++) {
        const auto& q_nbrs = query_.GetNeighbors(u);
        for (uint i = 0; i < q_nbrs.size(); i++) {
            uint u_other = q_nbrs[i];
            if (u > u_other) continue;
            gpu_orders_.push_back(order_vs_[eidx_[u][u_other]]);
        }
    }
}

size_t Parallel_CaLiG::AcquireSlot()
{
    size_t s;
    while (!free_slots_.try_pop(s)) {
        #pragma omp taskyield
    }
    return s;
}

void Parallel_CaLiG::ReleaseSlot(size_t s)
{
    free_slots_.push(s);
}

static inline bool CaLiGVersionedJoinCheck(const Graph& data, uint u_b_data_v, uint v,
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

size_t Parallel_CaLiG::EnumerateNewEdgeVersioned(uint v1, uint v2, uint label,
                                                   size_t /*thread_id*/, uint max_timestamp)
{
    if (max_num_results_ == 0) return 0;

    // Use versioned DAG-based search (same as SymBi/TurboFlux)
    const uint nv = query_.NumVertices();
    size_t num_results = 0;

    // Allocate temporary state
    std::vector<uint> m(nv, UNMATCHED);
    std::vector<bool> visited(data_.NumVertices(), false);

    for (uint u1 = 0; u1 < nv; u1++) {
        if (data_.GetVertexLabel(v1) != query_.GetVertexLabel(u1)) continue;
        for (uint u2 = 0; u2 < nv; u2++) {
            if (data_.GetVertexLabel(v2) != query_.GetVertexLabel(u2)) continue;
            if (std::get<2>(query_.GetEdgeLabel(u1, u2)) == NOT_EXIST) continue;

            m[u1] = v1; m[u2] = v2;
            visited[v1] = true; visited[v2] = true;

            // Use the order for directed edge u1→u2
            uint eidx = eidx_[u1][u2];
            FindMatches_versioned_v2(eidx, 2, m, visited, num_results, max_timestamp);

            visited[v1] = false; visited[v2] = false;
            m[u1] = UNMATCHED; m[u2] = UNMATCHED;

            if (num_results >= max_num_results_ || reach_time_limit) return num_results;
        }
    }
    return num_results;
}

void Parallel_CaLiG::FindMatches_versioned_v2(uint order_index, uint depth,
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
            if (!CaLiGVersionedJoinCheck(data_, m[jc_vs[i]], v, jc_labels[i], max_ts)) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;
        if (!homomorphism_ && visited[v]) continue;

        m[u] = v;
        visited[v] = true;

        if (depth == nv - 1) {
            num_results++;
        } else {
            FindMatches_versioned_v2(order_index, depth + 1, m, visited, num_results, max_ts);
        }

        visited[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_ || reach_time_limit) return;
    }
}

void Parallel_CaLiG::FindMatches_versioned_chunk(uint order_index, uint depth,
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
            if (!CaLiGVersionedJoinCheck(data_, m[jc_vs[i]], v, jc_labels[i], max_ts)) {
                joinable = false; break;
            }
        }
        if (!joinable) continue;
        if (!homomorphism_ && visited[v]) continue;

        m[u] = v;
        visited[v] = true;

        if (depth == nv - 1) {
            num_results++;
        } else {
            FindMatches_versioned_v2(order_index, depth + 1, m, visited, num_results, max_ts);
        }

        visited[v] = false;
        m[u] = UNMATCHED;
        if (num_results >= max_num_results_ || reach_time_limit) return;
    }
}
