#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
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
        const uint GetV1() const{return this->v1;}
        const uint GetV2() const{return this->v2;}
        const uint GetV1Label() const{return this->v1Label;}
        const uint GetV2Label() const{return this->v2Label;}
        const uint GeteLabel() const{return this->eLabel;}
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
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void LoadFromFile(const std::string &path);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
};

#endif //GRAPH_GRAPH
