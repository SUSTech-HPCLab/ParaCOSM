#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <algorithm>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "utils/types.h"
#include "utils/utils.h"

class Edge{
    private:
        uint v1;
        uint v2;
        uint v1Label;
        uint v2Label;
        uint eLabel;
        uint index;
        bool exist;
    public:
        Edge(uint v1, uint v2, uint v1Label, uint v2Label, uint eLabel, uint index):v1(v1),v2(v2),v1Label(v1Label),v2Label(v2Label),eLabel(eLabel),index(index){
            this->exist = true;
        }
        Edge(uint v1Label, uint v2Label, uint eLabel):v1Label(v1Label),v2Label(v2Label),eLabel(eLabel){}
        bool GetExist(){return this->exist;}
        void EdgeDelete(){this->exist = false;}
        uint GetV1() const{return this->v1;}
        uint GetV2() const{return this->v2;}
        uint GetV1Label() const{return this->v1Label;}
        uint GetV2Label() const{return this->v2Label;}
        uint GeteLabel() const{return this->eLabel;}
        uint GetIndex(){return this->index;}
        bool operator == (const Edge& edge)const {
            if((edge.v1Label == this->v1Label && edge.v2Label == this->v2Label && edge.eLabel == this->eLabel) ||
            (edge.v2Label == this->v1Label && edge.v1Label == this->v2Label && edge.eLabel == this->eLabel)){
                return true;
            }
            return false;
        }
        bool operator < (const Edge & edge) const{
            if(this->v1Label < edge.v1Label){
                return true;
            } else if (this->v1Label > edge.v1Label) {
                return false;
            }
    
            if(this->v2Label < edge.v2Label){
                return true;
            } else if (this->v2Label > edge.v2Label) {
                return false;
            }
    
    
            if(this->eLabel < edge.eLabel){
                return true;
            } else if (this->eLabel > edge.eLabel) {
                return false;
            }
    
            return false;
        }
    };

class Graph
{
protected:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;
    std::vector<std::vector<uint>> elabels_;
    // Hash index for high-degree vertices: neighbor_id → edge_label
    static constexpr uint HASH_DEGREE_THRESHOLD = 64;
    std::vector<std::unordered_map<uint, uint>> hash_adj_;
    // Per-edge timestamps for versioned search
    std::vector<std::vector<uint>> edge_timestamps_;

public:
    std::queue<InsertUnit> updates_;
    std::vector<InsertUnit> updates_vec_;
    std::vector<uint> vlabels_;

public:
    Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void AddEdgeVersioned(uint v1, uint v2, uint label, uint timestamp);

    // Batch parallel insert with per-edge timestamps.
    // Each edge: {v1, v2, label, timestamp}. Duplicates (already-existing or
    // repeated) are skipped. Returns the number of newly-inserted edges.
    struct VersionedEdgeBatch { uint v1, v2, label, timestamp; };
    size_t AddEdgesVersionedBatch(const std::vector<VersionedEdgeBatch>& edges,
                                  size_t num_threads);

    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    const std::vector<uint>& GetEdgeTimestamps(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void InitTimestamps();
    void ClearTimestamps();

    // Get timestamp of edge (v1, v2). Returns 0 if no timestamps or edge not found.
    uint GetEdgeTimestamp(uint v1, uint v2) const;

    /// Fast joinability check: is dst a neighbor of src with the given edge label?
    inline bool HasNeighborWithLabel(uint src, uint dst, uint label) const {
        if (src < hash_adj_.size() && !hash_adj_[src].empty()) {
            auto it = hash_adj_[src].find(dst);
            return it != hash_adj_[src].end() && it->second == label;
        }
        const auto& nbrs = neighbors_[src];
        auto it = std::lower_bound(nbrs.begin(), nbrs.end(), dst);
        if (it == nbrs.end() || *it != dst) return false;
        return elabels_[src][it - nbrs.begin()] == label;
    }

    void LoadFromFile(const std::string &path);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
};

#endif //GRAPH_GRAPH
