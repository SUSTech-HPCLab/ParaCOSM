#include <vector>
#include <map>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/matching.h"

// enum searchType{
//     init, pos, neg
// };
class CSMPP : public matching
{
private:

struct SearchContext {
    std::vector<int>* match;
    std::vector<std::vector<uint>>* matchCandidate;
    std::vector<bool>* visited;
};


    
    // std::vector<std::pair<Edge, std::vector<std::pair<uint,uint>>>> updateEdgeFindQueryEqual;// record the edge update 
    std::map<Edge,std::vector<std::pair<uint,uint>>> updateEdgeFindQuery; // const after init, ok
    std::map<Edge,std::vector<std::pair<uint,uint>>> initEdgeFindQuery; // ok
   
    std::vector<int>match;
    std::vector<std::vector<int>> match_parallel;

    std::vector<std::vector<uint>> matchCandidate;
    std::vector<std::vector<std::vector<uint>>> matchCandidate_parallel;

    // std::vector<uint> indexCheckCNT;
    
    size_t is_safe;
    size_t Thread_MAX;

    size_t NUM_L;

    //combine part
    std::vector<std::vector<std::vector<uint>>> intersectionResult;
    std::vector<std::vector<std::vector<std::vector<uint>>>> intersection_parallelResult;


    std::vector<std::vector<int>> combineStack;
    std::vector<std::vector<std::vector<uint>>> combineStack_parallel;

    std::vector<std::vector<uint>> headRecord;
    std::vector<std::vector<std::vector<uint>>> headRecord_parallel;

    std::vector<uint> stackSize;
    std::vector<std::vector<uint>> stackSize_parallel;

    std::vector<unFreezeStackType> type;
    std::vector<std::vector<unFreezeStackType>> type_parallel;

    std::vector<int> stackHead;
    std::vector<std::vector<int>> stackHead_parallel;


    bool print_init;

public:
    CSMPP(Graph& data_graph, Graph& query_graph, std::vector<Graph> & multiQueryGraph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, bool print_init);
    ~CSMPP() override{};

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void TimePrint();
 
    void localPopVertex(std::vector<int>& match, std::vector<bool>& visited, uint vertex, uint pos);
    void localMatchVertex(std::vector<int>& match, std::vector<bool>& visited, uint vertex, uint pos);

    void localSearchVertex(SearchContext& context, uint queryIndex, uint edgeIndex, searchType type, uint depth);

    void searchVertex_local(uint queryIndex, uint edgeIndex, searchType type, uint depth);
    void searchVertex_local2(uint queryIndex, uint edgeIndex, searchType type, uint depth, size_t thread_id);
    void Parallel_searchVertex(uint queryIndex, uint edgeIndex, searchType type, uint depth, size_t thread_id);

    void Parallel_searchVertex(uint queryIndex, uint edgeIndex, searchType type, 
        uint depth, size_t thread_id, std::vector<int>& match_local, std::vector<bool>& visited_local);

        void Parallel_searchVertex(uint queryIndex, uint edgeIndex, searchType type, 
            uint depth);

private:

    bool indexCheck(uint data_v, uint query_v, uint queryID);
    bool indexCheck(uint data_v, uint query_v, uint queryID, uint pos);
    void EdgeFindQueryInit();
    void initEdgeFindQueryInit();
    std::vector<std::pair<uint, uint>> UpdateEdgeFineQuery(uint v1, uint v2, uint edgeLabel, searchType type);
    void setMatchVertex(const std::vector<uint> & matchingIndex, const std::vector<int> & vertexs);
    void setMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex, const std::vector<int> & vertexs);
    void unsetMatchVertex(const std::vector<uint> & matchingIndex);
    void unsetMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex);
    void addMatchResult(uint queryIndex, uint edgeIndex, searchType type);



    #if defined(ADD_RESULT_DFS)
    void addMatchResult_DFS(uint queryIndex, uint edgeIndex, searchType type);
    
    void DFS_SEARCH_ALL(const std::vector<uint> & needToCombine, const std::vector<uint> & cacheMap, uint idx, searchType type);
    void DFS_SEARCH_ALL(const std::vector<std::vector<int>> & needToCombine, const std::vector<uint> & NoOverLeafWeight, const std::vector<uint> & isolatedIndex, uint idx, searchType type, size_t weight);
    #endif
    void printMatch();
    void searchInit(uint v1, uint v2, uint label, searchType type);
    void searchVertex(uint queryIndex, uint edgeIndex, searchType type, uint depth);
    void matchVertex(uint vertex, uint pos);
    void matchVertexlocal(uint vertex, uint pos, std::vector<int> & match, std::vector<bool>& visited_ );
    void matchVertexAll(uint vertex, uint pos);
    void matchVertex(std::vector<uint> & Candidate, uint pos);

    void popVertexAll(uint vertex, uint pos);
    
    void matchBatchVertices(uint index, std::vector<uint> & candidates);
    void unmatchBatchVertices(uint index);


    void Parallel_searchVertex_local(uint queryIndex, uint edgeIndex, searchType type, 
        uint depth, size_t thread_id)  ;

    void popVertex(uint vertex, uint pos);
    void popVertex(uint pos);
    void popVertex_local(uint vertex, uint pos, std::vector<int> & match, std::vector<bool>& visited_ );
    void popVertex_local(uint pos, std::vector<int> & match, std::vector<std::vector<std::seed_seq::result_type>>& matchCandidate );
    
    void printMatch(const std::vector<uint> & matchOrder);


    bool vertexPushCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel);
    bool getIntersection(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel);

    bool getIntersection_Parallel(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel);

    std::vector<uint> & getItersectionTop(int depth);
    void combineStackPopTail(int depth);
    void combineStackPopTail(int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    void combinePushBack(uint vertex, uint vertexInWhichCandidateNeighbor, int depth);
    void combinePushBack(int vertex, uint vertexInWhichCandidateNeighbor, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    bool headChange(const std::vector<std::vector<uint>> & needToCombine, int depth);
    bool headChange(const std::vector<std::vector<int>> & needToCombine, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);

    void setVisitedPatch(const std::vector<int> & vertex);
    void setUnVisitedPatch(const std::vector<int> & vertex);

    bool safe_detect(uint v1, uint v2, uint label, searchType type) override;

    void Safe_Update(uint v1, uint v2, uint label) override;

    void Safe_Update_remove(uint v1, uint v2) override;

    static bool cmp(const std::pair<int,int> &p1,const std::pair<int,int> &p2){
        return p1.second < p2.second;
    }


    void matchVertexlocal(std::vector<uint> & Candidate, uint pos, std::vector<int> & match,
        std::vector<std::vector<std::seed_seq::result_type>>& matchCandidate);
    
    void localPopVertex(std::vector<int>& match,  std::vector<std::vector<uint>>& matchCandidate, uint pos);


    void matchVertexAll(std::vector<uint> & Candidate, uint pos);
    void popVertexAll(uint pos);

    void matchVertexlocal(std::vector<uint> & Candidate, uint pos, std::vector<int> & match,
        std::vector<std::vector<std::seed_seq::result_type>>& matchCandidate, size_t thread_id);

        void popVertex_local(uint pos, std::vector<int> & match, 
            std::vector<std::vector<std::seed_seq::result_type>>& matchCandidate, size_t thread_id );

};
