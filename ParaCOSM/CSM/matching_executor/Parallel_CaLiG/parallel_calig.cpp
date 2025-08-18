#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <map>

#include <omp.h>  

#include <chrono>
#include <climits>
#include <functional>
#include <stdlib.h>

// #include <ska::flat_hash_map>
#include "storage_hash_map.hpp"

#define CaLiG_Get_Time() std::chrono::high_resolution_clock::now()

#define CaLiG_Duration(start) std::chrono::duration_cast<\
    std::chrono::microseconds>(CaLiG_Get_Time() - start).count()/(float)1000

#define CaLiG_Print_Time(str, start) std::cout << str << CaLiG_Duration(start) << \
    "ms" << std::endl
#define CaLiG_Print_Time_Nano(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() << " ns" << std::endl
#define CaLiG_Print_Time_Micro(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " us" << std::endl


// std::vector<ska::flat_hash_map<uint, bool>> d1;

using namespace std;
using u_set = unordered_set<uint>;
using vec = vector<uint>;


struct vertex_Q{
    uint label;
    u_set nei;
    ska::flat_hash_map<uint, vec> rep_nei;
};

struct vertex_G{
    uint label;
    u_set nei;
    ska::flat_hash_map<uint, ska::flat_hash_map<uint, u_set>> cand;
    ska::flat_hash_map<uint, bool> LI;
};

vector<vertex_Q> Q;
vector<vertex_G> G;
uint Q_size, G_size;
vec update;
ska::flat_hash_map<uint, vec> labels;

bool isInVec(uint a, vec& v){
    for(uint i=0; i<v.size(); i++){
        if(a == v[i]) return 1;
    }
    return 0;
}

u_set intersection(u_set& us1, u_set& us2){
    u_set result;
    if(us1.size()<=us2.size()){
        for(auto& item : us1){
            if(us2.find(item)!=us2.end()) result.insert(item);
        }
    }
    else{
        for(auto& item : us2){
            if(us1.find(item)!=us1.end()) result.insert(item);
        }
    }
    return result;
}

void dumpG(const string& path){
    ofstream outfile(path);
    for(uint vi=0; vi<G_size; vi++){
        for(auto& candi : G[vi].cand){
            uint ui = candi.first;
            if(G[vi].LI[ui]){
                outfile << "v" << vi << "-u" << ui << " nei:";
                for(auto& nei : candi.second){
                    uint uj = nei.first;
                    outfile << " [u" << uj << "]";
                    for(uint vj : nei.second){
                        if(G[vj].LI[uj]) outfile << " v" << vj;
                    }
                }
                outfile << endl;
            }
        }
    }
    outfile.close();
}



void inputQ(string& qp){
    ifstream qg(qp);
    if(!qg) cerr << "Fail to open query file." << endl;

    char c;
    uint id, id1, id2, lb, dvir;
    vertex_Q u;
    while(qg >> c){
        if(c == 'v'){
            qg >> id >> lb;
            if(lb == -1) lb = 0;
            u.label = lb;
            Q.push_back(u);
        }
        else{
            qg >> id1 >> id2;
            Q[id1].nei.insert(id2);
            Q[id2].nei.insert(id1);
            break;
        }
    }
    while(qg >> c){
        qg >> id1 >> id2;
        Q[id1].nei.insert(id2);
        Q[id2].nei.insert(id1);
    }
    qg.close();
    Q_size = Q.size();

    for(uint ui=0; ui<Q_size; ui++){
        labels[Q[ui].label].emplace_back(ui);
        ska::flat_hash_map<uint, vec> rep_nei;
        for(auto& uni : Q[ui].nei){
            uint& uni_label = Q[uni].label;
            rep_nei[uni_label].emplace_back(uni);
        }
        for(auto& repi : rep_nei){
            if(repi.second.size()>1){
                Q[ui].rep_nei[repi.first] = repi.second;
            }
        }
    }
}

struct split{
    vec core;
    vector<u_set> core_nei;
    ska::flat_hash_map<uint, u_set> c_s_nei;
    vec shell;
    vector<u_set> shell_nei;
};

ska::flat_hash_map<uint, ska::flat_hash_map<uint, split>> matching_order;

void insertNei(u_set &temp, uint i, uint ui, uint uj, split& uj_second){
    for(auto& nei : Q[i].nei){
        uint n = nei;
        if(n!=ui && n!=uj && !isInVec(n, uj_second.core) && !isInVec(n, uj_second.shell)){
            temp.insert(n);
        }
    }
}

bool neiAllCore(uint i, split& uj_second, uint ui, uint uj){
    for(auto& nei : Q[i].nei){
        uint n = nei;
        if(n!=ui && n!=uj && !isInVec(n, uj_second.core)){
            return false;
        }
    }
    return true;
}


// Generate the matching order for each vertex in Q
// The matching order is stored in a hash map where the key is the vertex index
void generateMO(){
    for(uint ui=0; ui<Q_size; ui++) {
        ska::flat_hash_map<uint, split> ui_first;
        for (auto &uj : Q[ui].nei) {
            if (uj >= 0) {
                split uj_second;
                u_set temp;
                insertNei(temp, ui, ui, uj, uj_second);
                insertNei(temp, uj, ui, uj, uj_second);
                while(!temp.empty()){
                    u_set temp2 = temp;
                    for(auto& ti : temp2){
                        if(neiAllCore(ti, uj_second, ui, uj)){
                            u_set nei_info = Q[ti].nei;
                            uint v_core;
                            for(uint & i : uj_second.core){
                                if(nei_info.count(i)){
                                    v_core = i;
                                }
                            }
                            uj_second.c_s_nei[v_core].insert(uj_second.shell.size());
                            uj_second.shell.push_back(ti);
                            uj_second.shell_nei.push_back(nei_info);
                            temp.erase(ti);
                        }
                    }
                    if(!temp.empty()){
                        uint nei_num=0;
                        uint core_v;
                        for(auto& t : temp){
                            u_set nei_noij = Q[t].nei;
                            nei_noij.erase(ui);
                            nei_noij.erase(uj);
                            if(nei_noij.size() > nei_num){
                                nei_num = nei_noij.size();
                                core_v = t;
                            }
                        }
                        uj_second.core.push_back(core_v);
                        u_set nei_info;
                        for(auto& nei : Q[core_v].nei){
                            uint n = nei;
                            if(n == ui || n == uj || isInVec(n, uj_second.core)){
                                nei_info.insert(nei);
                            }
                        }
                        uj_second.core_nei.push_back(nei_info);
                        insertNei(temp, core_v, ui, uj, uj_second);
                        temp.erase(core_v);
                    }
                }
                ui_first[uj] = uj_second;
            }
        }
        matching_order[ui] = ui_first;
    }
}

// void inputG(string& gp){
//     ifstream dg(gp);
//     if(!dg) cerr << "Fail to open graph file." << endl;

//     char c;
//     uint id, id1, id2, lb;
//     vertex_G v;
//     while(dg >> c){
//         if(c == 'v'){
//             dg >> id >> lb;
//             if(labels.count(lb))  v.label = lb;
//             else  v.label = -1;
//             G.push_back(v);
//         }
//         else{
//             G_size = G.size();
//             dg >> id1 >> id2;
//             if(id1<G_size && id2<G_size && G[id1].label!=-1 && G[id2].label!=-1){
//                 G[id1].nei.insert(id2);
//                 G[id2].nei.insert(id1);
//             }
//             break;
//         }
//     }
//     while(dg >> c){
//         dg >> id1 >> id2;
//         if(id1<G_size && id2<G_size && G[id1].label!=-1 && G[id2].label!=-1){
//             G[id1].nei.insert(id2);
//             G[id2].nei.insert(id1);
//         }
//     }
//     dg.close();
// }



void inputG(string& gp) {
    ifstream dg(gp);
    if(!dg) {
        cerr << "Fail to open graph file." << endl;
        return;
    }

    char c;
    int id, id1, id2, lb;
    int max_id = -1;  // 记录遇到的最大ID

    // 第一遍：找出最大ID
    while(dg >> c) {
        if(c == 'v') {
            dg >> id >> lb;
            if(id > max_id) max_id = id;
        }
    }

    // 重置文件指针
    dg.clear();
    dg.seekg(0);

    // 预先分配足够大的空间，填充无效顶点
    G.resize(max_id + 1);
    for(auto& v : G) {
        v.label = -1;  // 标记为无效
    }

    // 第二遍：读取实际顶点和边
    while(dg >> c) {
        if(c == 'v') {
            dg >> id >> lb;
            if(labels.count(lb)) {
                G[id].label = lb;
            }
        }
        else if(c == 'e') {
            dg >> id1 >> id2;
            if(id1 <= max_id && id2 <= max_id && 
               G[id1].label != -1 && G[id2].label != -1) {
                G[id1].nei.insert(id2);
                G[id2].nei.insert(id1);
            }
            break;
        }
    }

    // 处理剩余的边
    while(dg >> c) {
        if(c == 'e') {
            dg >> id1 >> id2;
            if(id1 <= max_id && id2 <= max_id && 
               G[id1].label != -1 && G[id2].label != -1) {
                G[id1].nei.insert(id2);
                G[id2].nei.insert(id1);
            }
        }
    }

    dg.close();
    G_size = G.size();  // 更新全局G_size
}



void constructCand(){
    for(uint vi=0; vi<G_size; vi++){
        uint& lb = G[vi].label;
        if(lb!=-1){
            for(auto& ui: labels[lb]){
                ska::flat_hash_map<uint, u_set> ui_cand;
                for(auto& vj : G[vi].nei){
                    uint vj_lb = G[vj].label;
                    for(auto& uj : Q[ui].nei){
                        if(ui_cand.count(uj) == 0) {
                            ui_cand[uj] = u_set();  
                        }
                        if(uj>=0 && vj>=0 && Q[uj].label==vj_lb){
                            ui_cand[uj].insert(vj);
                        }
                    }
                }
                G[vi].cand[ui] = ui_cand;
                G[vi].LI[ui] = 1;
            }
        }
    }
    //dumpG("./dump/G_init");
}


bool tryNei(uint th, uint vi, uint ui, u_set& used, vec& to_check){
    if(th==to_check.size()){
        return 1;
    }
    else{
        auto& uj = to_check[th];
        for(auto& vj : G[vi].cand[ui][uj]){
            if(G[vj].label!=Q[uj].label) {
                G[vj].cand.erase(uj);
                G[vj].LI.erase(uj);
                continue;
            }
            if(used.find(vj)==used.end()){
                used.insert(vj);
                if(tryNei(th+1, vi, ui, used, to_check)) return 1;
                used.erase(vj);
            }
        }
        return 0;
    }
}

bool checkNei(uint vi, uint ui){
    for(auto& uj : Q[ui].nei){
        if(G[vi].cand[ui].find(uj)==G[vi].cand[ui].end()) return 0;
        else if(!G[vi].cand[ui][uj].empty()) continue;
        else return 0;
    }
    for(auto& rep_nei : Q[ui].rep_nei){
        u_set used;
        if(!tryNei(0, vi, ui, used, rep_nei.second)) return 0;
    }
    return 1;
}

void deleteAndCheck(uint ui, uint vi){

    for(auto& nei : G[vi].cand[ui]){
        uint uj = nei.first;
        for(auto& vj : nei.second){
            if(G[vj].label!=Q[uj].label) {
                G[vj].cand.erase(uj);
                G[vj].LI.erase(uj);
                continue;
            }
            if(G[vj].LI[uj]) {
                G[vj].cand[uj][ui].erase(vi);
                if (G[vj].cand[uj][ui].empty()) {
                    G[vj].LI[uj] = 0;
                    deleteAndCheck(uj, vj);
                }
                else{
                    auto lb = Q[ui].label;
                    if (Q[uj].rep_nei.find(lb)!=Q[uj].rep_nei.end()){
                        auto& rep_nei = Q[uj].rep_nei[lb];
                        u_set used;
                        if(!tryNei(0, vj, uj, used, rep_nei)){
                            G[vj].LI[uj] = 0;
                            deleteAndCheck(uj, vj);
                        }
                    }
                }
            }
        }
    }
}


void deleteAndCheck2(uint u, uint v){

    for(auto& nei : G[v].cand[u]){
        uint uj = nei.first;
        for(auto& vj : nei.second){
            if(G[vj].label!=Q[uj].label) {
                G[vj].cand.erase(uj);
                G[vj].LI.erase(uj);
                continue;
            }
            if(G[vj].LI[uj]) {
                G[vj].cand[uj][u].erase(v);
                if (G[vj].cand[uj][u].empty()) {
                    G[vj].LI[uj] = 0;
                    deleteAndCheck(uj, vj);
                }
                else{
                    auto lb = Q[u].label;
                    if (Q[uj].rep_nei.find(lb)!=Q[uj].rep_nei.end()){
                        auto& rep_nei = Q[uj].rep_nei[lb];
                        u_set used;
                        if(!tryNei(0, vj, uj, used, rep_nei)){
                            G[vj].LI[uj] = 0;
                            deleteAndCheck(uj, vj);
                        }
                    }
                }
            }
        }
    }
}


void turnOff(uint& vi){
    for(auto& candi : G[vi].cand){
        auto& ui = candi.first;
        if(G[vi].label!=Q[ui].label) {
            G[vi].cand.erase(ui);
            G[vi].LI.erase(ui);
            continue;
        }
        if(!G[vi].LI[ui]) continue;
        if ( !checkNei(vi, ui) ){
            G[vi].LI[ui] = 0;
            deleteAndCheck(ui, vi);
        }
    }
}


void turnOffProcess(uint v1, uint v2){
    G[v1].nei.erase(v2);
    G[v2].nei.erase(v1);
    for(auto& candi : G[v1].cand){
        uint ui = candi.first;
        if(G[v1].label!=Q[ui].label){
            G[v1].cand.erase(ui);
            G[v1].LI.erase(ui);
            continue;
        }
        for(auto& ui_nei : candi.second){
            uint uj = ui_nei.first;
            if(Q[uj].label==G[v2].label){
                G[v1].cand[ui][uj].erase(v2);
                G[v2].cand[uj][ui].erase(v1);
            }
        }
    }

    turnOff(v1);
    turnOff(v2);
}

void turnOffProcess_parallel(uint v1, uint v2){

    G[v1].nei.erase(v2);
    G[v2].nei.erase(v1);
    
    // v1 
    auto cand_v1 = G[v1].cand;
    
    #pragma omp parallel
    {

        std::vector<std::tuple<uint, uint, uint>> local_updates; // (ui, uj, v2)
        std::vector<uint> local_erase_ui; // ui
        
        #pragma omp for nowait
        for (int i = 0; i < cand_v1.size(); i++) {
            // 使用迭代器访问 map 的第 i 个元素
            auto it = cand_v1.begin();
            std::advance(it, i);
            uint ui = it->first;
            auto& candi = it->second;
            
            if(G[v1].label != Q[ui].label){
                local_erase_ui.push_back(ui);
                continue;
            }
            
            for(auto& ui_nei : candi){
                uint uj = ui_nei.first;
                if(Q[uj].label == G[v2].label){
                    local_updates.push_back(std::make_tuple(ui, uj, v2));
                }
            }
        }

        #pragma omp critical
        {
            for(auto ui : local_erase_ui){
                G[v1].cand.erase(ui);
                G[v1].LI.erase(ui);
            }
            
            for(auto& update : local_updates){
                uint ui = std::get<0>(update);
                uint uj = std::get<1>(update);
                uint v_val = std::get<2>(update);
                
                G[v1].cand[ui][uj].erase(v_val);
                G[v2].cand[uj][ui].erase(v1);
            }
        }
    }
    
    turnOff(v1);
    turnOff(v2);
}

void addAndCheck(uint ui, uint vi, vec& temp_v, vec& temp_u){
    for(auto& nei : G[vi].cand[ui]){
        uint uj = nei.first;
        for(auto& vj : nei.second){
            if(G[vj].label!=Q[uj].label) {
                G[vj].cand.erase(uj);
                G[vj].LI.erase(uj);
                continue;
            }
            G[vj].cand[uj][ui].insert(vi);
            if(!G[vj].LI[uj]){
                if(checkNei(vj, uj)){
                    G[vj].LI[uj] = 1;
                    addAndCheck(uj, vj, temp_v, temp_u);
                }
                else{
                    temp_v.emplace_back(vj);
                    temp_u.emplace_back(uj);
                }
            }
        }
    }
}

/*
Explanation of Comments
Function Purpose: The initial comment describes the function's role in handling edge additions and updating the CaLiG index.
Step-by-Step Breakdown:
    Step 1: Adding the edge to the graph's neighbor sets.
    Step 2: Iterating over v1's candidates and filtering by label compatibility.
    Step 3: Checking ui's neighbors and updating candidate sets based on v2.
State Update Logic:
    Case 1: When (ui, v1) is already "ON", updates are propagated symmetrically.
    Case 2: When (ui, v1) is not "ON", checks if it should be turned "ON" and propagates.
Propagation Handling: Comments explain how addAndCheck and deleteAndCheck manage state changes and their cascading effects.
*/ 

// Function to handle the addition of an edge (v1, v2) in the data graph G
// and update the CaLiG index's candidate sets and lighting states accordingly
void turnOnProcess(uint& v1, uint& v2) {
    // Temporary vectors to store affected vertex pairs during state propagation
    vec temp_v, temp_u;

    // Step 1: Add the edge (v1, v2) to the data graph G by updating neighbor sets
    G[v1].nei.insert(v2);  // Add v2 to v1's neighbors
    G[v2].nei.insert(v1);  // Add v1 to v2's neighbors (undirected graph)

    // Step 2: Iterate through v1's candidate set to update based on the new edge
    for (auto& candi : G[v1].cand) {
        uint ui = candi.first;  // ui is a query vertex from Q that v1 might match

        // Check label compatibility: v1 must have the same label as ui to be a candidate
        if (G[v1].label != Q[ui].label) {
            G[v1].cand.erase(ui);  // Remove ui from v1's candidate set if labels don't match
            G[v1].LI.erase(ui);    // Remove ui's lighting state as it's no longer relevant
            continue;              // Skip to the next candidate
        }

        // Step 3: Check ui's neighbors in the query graph Q
        for (auto& ui_nei : candi.second) {
            uint uj = ui_nei.first;  // uj is a neighbor of ui in Q

            // If v2's label matches uj's label, v2 could be a candidate for uj
            if (G[v2].label == Q[uj].label) {
                // Update v1's candidate set: v2 is a potential match for uj given (ui, v1)
                G[v1].cand[ui][uj].insert(v2);

                // Case 1: If (ui, v1) is already "ON", propagate the update
                if (G[v1].LI[ui]) {
                    // Symmetrically update v2's candidate set: v1 is a candidate for ui
                    G[v2].cand[uj][ui].insert(v1);

                    // Propagate "ON" state changes starting from (ui, v1)
                    addAndCheck(ui, v1, temp_v, temp_u);

                    // Check affected pairs and turn off invalid states if necessary
                    for (uint i = 0; i < temp_v.size(); i++) {
                        uint& v = temp_v[i];  // Data graph vertex
                        uint& u = temp_u[i];  // Query graph vertex
                        if (!G[v].LI[u]) {    // If (u, v) is not "ON"
                            deleteAndCheck(u, v);  // Turn it "OFF" and propagate
                        }
                    }
                    temp_v.clear();  // Clear temporary storage
                    temp_u.clear();
                }
                // Case 2: If (ui, v1) is not "ON", check if it should be turned "ON"
                else if (checkNei(v1, ui)) {
                    G[v1].LI[ui] = true;  // Set (ui, v1) to "ON" state

                    // Propagate "ON" state changes starting from (ui, v1)
                    addAndCheck(ui, v1, temp_v, temp_u);

                    // Check affected pairs and turn off invalid states if necessary
                    for (uint i = 0; i < temp_v.size(); i++) {
                        uint& v = temp_v[i];  // Data graph vertex
                        uint& u = temp_u[i];  // Query graph vertex
                        if (!G[v].LI[u]) {    // If (u, v) is not "ON"
                            deleteAndCheck(u, v);  // Turn it "OFF" and propagate
                        }
                    }
                    temp_v.clear();  // Clear temporary storage
                    temp_u.clear();
                }
            }
        }
    }
}

void turnOnProcess2(uint& v1, uint& v2){
    vec temp_v, temp_u;
    G[v1].nei.insert(v2);
    G[v2].nei.insert(v1);
    for(auto& candi : G[v1].cand){
        uint ui = candi.first;
        if(G[v1].label!=Q[ui].label){
            G[v1].cand.erase(ui);
            G[v1].LI.erase(ui);
            continue;
        }
        for(auto& ui_nei : candi.second){
            uint uj = ui_nei.first;
            if(G[v2].label == Q[uj].label){
                G[v1].cand[ui][uj].insert(v2);
                if(G[v1].LI[ui]){
                    G[v2].cand[uj][ui].insert(v1);
                    addAndCheck(ui, v1, temp_v, temp_u);
                    for(uint i=0; i<temp_v.size(); i++){
                        uint& v = temp_v[i];
                        uint& u = temp_u[i];
                        if(!G[v].LI[u]){
                            // deleteAndCheck(u, v);
                            for(auto& nei : G[v].cand[u]){
                                uint uj = nei.first;
                                for(auto& vj : nei.second){
                                    if(G[vj].label!=Q[uj].label) {
                                        G[vj].cand.erase(uj);
                                        G[vj].LI.erase(uj);
                                        continue;
                                    }
                                    if(G[vj].LI[uj]) {
                                        G[vj].cand[uj][u].erase(v);
                                        if (G[vj].cand[uj][u].empty()) {
                                            G[vj].LI[uj] = 0;
                                            deleteAndCheck(uj, vj);
                                        }
                                        else{
                                            auto lb = Q[u].label;
                                            if (Q[uj].rep_nei.find(lb)!=Q[uj].rep_nei.end()){
                                                auto& rep_nei = Q[uj].rep_nei[lb];
                                                u_set used;
                                                if(!tryNei(0, vj, uj, used, rep_nei)){
                                                    G[vj].LI[uj] = 0;
                                                    deleteAndCheck(uj, vj);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    temp_v.clear();
                    temp_u.clear();
                }
                else if ( checkNei(v1, ui) ){
                    G[v1].LI[ui] = true;
                    addAndCheck(ui, v1, temp_v, temp_u);
                    for(uint i=0; i<temp_v.size(); i++){
                        uint& v = temp_v[i];
                        uint& u = temp_u[i];
                        if(!G[v].LI[u]){
                            // deleteAndCheck(u, v);
                            for(auto& nei : G[v].cand[u]){
                                uint uj = nei.first;
                                for(auto& vj : nei.second){
                                    if(G[vj].label!=Q[uj].label) {
                                        G[vj].cand.erase(uj);
                                        G[vj].LI.erase(uj);
                                        continue;
                                    }
                                    if(G[vj].LI[uj]) {
                                        G[vj].cand[uj][u].erase(v);
                                        if (G[vj].cand[uj][u].empty()) {
                                            G[vj].LI[uj] = 0;
                                            deleteAndCheck(uj, vj);
                                        }
                                        else{
                                            auto lb = Q[u].label;
                                            if (Q[uj].rep_nei.find(lb)!=Q[uj].rep_nei.end()){
                                                auto& rep_nei = Q[uj].rep_nei[lb];
                                                u_set used;
                                                if(!tryNei(0, vj, uj, used, rep_nei)){
                                                    G[vj].LI[uj] = 0;
                                                    deleteAndCheck(uj, vj);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    temp_v.clear();
                    temp_u.clear();
                }
            }
        }
    }


}

bool turnOnProcess_safe(uint v1, uint v2) {
    // Temporary vectors to store affected vertex pairs during state propagation
    // vec temp_v, temp_u;

    // Step 1: Add the edge (v1, v2) to the data graph G by updating neighbor sets
    G[v1].nei.insert(v2);  // Add v2 to v1's neighbors
    G[v2].nei.insert(v1);  // Add v1 to v2's neighbors (undirected graph)

    // Step 2: Iterate through v1's candidate set to update based on the new edge
    for (auto& candi : G[v1].cand) {
        uint ui = candi.first;  // ui is a query vertex from Q that v1 might match

        // // Check label compatibility: v1 must have the same label as ui to be a candidate
        if (G[v1].label != Q[ui].label) {
            // G[v1].cand.erase(ui);  // Remove ui from v1's candidate set if labels don't match
            // G[v1].LI.erase(ui);    // Remove ui's lighting state as it's no longer relevant
            continue;              // Skip to the next candidate
        }

        if(G[v1].label != Q[ui].label){continue; }
        // Step 3: Check ui's neighbors in the query graph Q
        for (auto& ui_nei : candi.second) {
            return false;

        }
    }

    return true;
}

bool turnOffProcess_safe(uint v1, uint v2) {
    // Temporary vectors to store affected vertex pairs during state propagation
    // vec temp_v, temp_u;

    // Step 1: Add the edge (v1, v2) to the data graph G by updating neighbor sets
    G[v1].nei.insert(v2);  // Add v2 to v1's neighbors
    G[v2].nei.insert(v1);  // Add v1 to v2's neighbors (undirected graph)

    // Step 2: Iterate through v1's candidate set to update based on the new edge
    for (auto& candi : G[v1].cand) {
        uint ui = candi.first;  // ui is a query vertex from Q that v1 might match

        // // Check label compatibility: v1 must have the same label as ui to be a candidate
        if (G[v1].label != Q[ui].label) {
            // G[v1].cand.erase(ui);  // Remove ui from v1's candidate set if labels don't match
            // G[v1].LI.erase(ui);    // Remove ui's lighting state as it's no longer relevant
            continue;              // Skip to the next candidate
        }

        if(G[v1].label != Q[ui].label){continue; }
        // Step 3: Check ui's neighbors in the query graph Q
        for (auto& ui_nei : candi.second) {
            uint uj = ui_nei.first;  // uj is a neighbor of ui in Q

            // If v2's label matches uj's label, v2 could be a candidate for uj
            if (G[v2].label == Q[uj].label) {
                // // Update v1's candidate set: v2 is a potential match for uj given (ui, v1)
                // G[v1].cand[ui][uj].insert(v2);

                // // Case 1: If (ui, v1) is already "ON", propagate the update
                // if (G[v1].LI[ui]) {
                //     // Symmetrically update v2's candidate set: v1 is a candidate for ui
                //     G[v2].cand[uj][ui].insert(v1);

                //     // Propagate "ON" state changes starting from (ui, v1)
                //     addAndCheck(ui, v1, temp_v, temp_u);

                //     // Check affected pairs and turn off invalid states if necessary
                //     for (uint i = 0; i < temp_v.size(); i++) {
                //         uint& v = temp_v[i];  // Data graph vertex
                //         uint& u = temp_u[i];  // Query graph vertex
                //         if (!G[v].LI[u]) {    // If (u, v) is not "ON"
                //             deleteAndCheck(u, v);  // Turn it "OFF" and propagate
                //         }
                //     }
                //     temp_v.clear();  // Clear temporary storage
                //     temp_u.clear();
                // }
                // // Case 2: If (ui, v1) is not "ON", check if it should be turned "ON"
                // else if (checkNei(v1, ui)) {
                //     G[v1].LI[ui] = true;  // Set (ui, v1) to "ON" state

                //     // Propagate "ON" state changes starting from (ui, v1)
                //     addAndCheck(ui, v1, temp_v, temp_u);

                //     // Check affected pairs and turn off invalid states if necessary
                //     for (uint i = 0; i < temp_v.size(); i++) {
                //         uint& v = temp_v[i];  // Data graph vertex
                //         uint& u = temp_u[i];  // Query graph vertex
                //         if (!G[v].LI[u]) {    // If (u, v) is not "ON"
                //             deleteAndCheck(u, v);  // Turn it "OFF" and propagate
                //         }
                //     }
                //     temp_v.clear();  // Clear temporary storage
                //     temp_u.clear();
                // }
                return false;
            }
        }
    }

    return true;
}





bool turnOnProcess_safe_parallel(uint& v1, uint& v2) {
    // Step 1
    G[v1].nei.insert(v2);
    G[v2].nei.insert(v1);
    
    auto cand_v1 = G[v1].cand;
    
    #pragma omp parallel num_threads(32)
    {
        vec local_temp_v, local_temp_u;
        std::vector<std::tuple<uint, uint, uint>> updates; 
        
        #pragma omp for nowait
        for (int i = 0; i < cand_v1.size(); i++) {
            auto it = cand_v1.begin();
            std::advance(it, i);
            uint ui = it->first;
            auto& candi = it->second;
            
            // LABEL
            if (G[v1].label != Q[ui].label) {
                #pragma omp critical
                {
                    G[v1].cand.erase(ui);
                    G[v1].LI.erase(ui);
                }
                continue;
            }
            
            for (auto& ui_nei : candi) {
                
                uint uj = ui_nei.first;
                
                if (G[v2].label == Q[uj].label) {

                    updates.push_back(std::make_tuple(ui, v1, uj));
                    
                    // only one thread at one time
                    // #pragma omp critical
                    {

                        G[v1].cand[ui][uj].insert(v2);
                        
                        // Case 1: "ON"
                        if (G[v1].LI[ui]) {
                            G[v2].cand[uj][ui].insert(v1);
                            
                            // ON
                            addAndCheck(ui, v1, local_temp_v, local_temp_u);
                            
                            for (uint i = 0; i < local_temp_v.size(); i++) {
                                uint& v = local_temp_v[i];
                                uint& u = local_temp_u[i];
                                if (!G[v].LI[u]) {
                                    deleteAndCheck(u, v);
                                }
                            }
                            local_temp_v.clear();
                            local_temp_u.clear();
                        }
                        // Case 2: "OFF"
                        else if (checkNei(v1, ui)) {
                            G[v1].LI[ui] = true;
                            
                            // ON
                            addAndCheck(ui, v1, local_temp_v, local_temp_u);
                            
                            for (uint i = 0; i < local_temp_v.size(); i++) {
                                uint& v = local_temp_v[i];
                                uint& u = local_temp_u[i];
                                if (!G[v].LI[u]) {
                                    deleteAndCheck(u, v);
                                }
                            }
                            local_temp_v.clear();
                            local_temp_u.clear();
                        }
                      
                    }
                }
            }

        }
    }

    return true; 
}

void staticFilter(){
    for(uint vi=0; vi<G_size; vi++){
        turnOff(vi);
    }
    //dumpG("./dump/G_static");
}

void dumpMatch(ska::flat_hash_map<uint, uint>& m){
    cerr << "matching results is：" ;
    for(uint i=0; i<m.size(); i++){
        uint v = m[i];
        cerr << " u" << i << "-v" << v;
    }
    cerr << endl;
}

bool shellCand(vector<u_set>& result, ska::flat_hash_map<uint, uint>& m, const vec& s, const vector<u_set>& s_n, vec& used){
    uint s_size = s.size();
    result.resize(s_size);
    for(uint i = 0; i<s_size; i++) {
        auto p_ui = s_n[i].begin();
        if(s_n[i].size()>1){
            auto ui_l = *p_ui;
            ++p_ui;
            auto ui = * p_ui;
            result[i] = intersection(G[m[ui_l]].cand[ui_l][s[i]], G[m[ui]].cand[ui][s[i]]);
            for(++p_ui;p_ui!=s_n[i].end();++p_ui){
                ui = *p_ui;
                result[i] = intersection(result[i], G[m[ui]].cand[ui][s[i]]);
            }
            if (result[i].empty()) return 1;
        }
        else{
            auto ui = *p_ui;
            result[i] = G[m[ui]].cand[ui][s[i]];
        }
        for(auto& vth : used){
            result[i].erase(vth);
        }
        if (result[i].empty()) return 1;
    }
    return 0;
}

uint numAdd(uint th, const vector<u_set>& cand, u_set& used){
    uint result = 0;
    if(th==cand.size()-1){
        uint del = 0;
        for(auto& vth : used){
            if(cand[th].find(vth)!=cand[th].end()){
                del += 1;
            }
        }
        return cand[th].size() - del;
    }
    for(auto& vth : cand[th]){
        if(used.find(vth)==used.end()){
            used.insert(vth);
            result += numAdd(th+1, cand, used);
            used.erase(vth);
        }
    }
    return result;
}

bool notExit(uint shell, u_set nei, ska::flat_hash_map<uint, uint>& m){
    if(nei.size()==1){
        return false;
    }
    uint n = *nei.begin();
    for (auto cand : G[m[n]].cand[n][shell]){
        uint count = 0;
        for(uint ni : nei){
            count += 1;
            if(count==1){
                continue;
            }
            else{
                if(G[m[ni]].cand[ni][shell].count(cand)){
                    if(count == nei.size()){
                        return false;
                    }
                }
                else{
                    break;
                }
            }
        }
    }
    return true;
}

uint searchCore(uint th, ska::flat_hash_map<uint, uint>& m, vec& used, const vec& c, const vector<u_set>& c_n, const vec& s, const vector<u_set>& s_n, ska::flat_hash_map<uint, u_set>& c2check){
    uint result = 0;
    if(th == c.size()){
        vector<u_set> candidates;
        if(shellCand(candidates, m, s, s_n, used)) return 0;
        u_set used_v;
        return numAdd(0, candidates, used_v);
    }
    u_set candidates;
    auto p_ui = c_n[th].begin();
    auto ui = *p_ui;
    if(c_n[th].size()>1){
        u_set* temp;
        temp = &G[m[ui]].cand[ui][c[th]];
        candidates = intersection(*temp,G[m[ui]].cand[ui][c[th]]);
        for(++p_ui;p_ui!=c_n[th].end();++p_ui){
            ui = *p_ui;
            candidates = intersection(candidates,G[m[ui]].cand[ui][c[th]]);
            if(candidates.empty()) return 0;
        }
        for(auto& vi : used){
            candidates.erase(vi);
        }
        for(auto& vth : candidates){
            m[c[th]] = vth;
            if(c2check.count(c[th])){
                bool to_continue = false;
                for (auto& shell : c2check[c[th]]){
                    if(notExit(s[shell], s_n[shell], m)) {
                        to_continue = true;
                    }
                }
                if(to_continue) {
                    continue;
                }
            }
            used.emplace_back(vth);
            result += searchCore(th+1, m, used, c, c_n, s, s_n, c2check);
            used.pop_back();
        }
    }
    else{
        candidates = G[m[ui]].cand[ui][c[th]];
        for(auto& vth : candidates){
            if(!isInVec(vth, used)){
                m[c[th]] = vth;
                if(c2check.count(c[th])){
                    bool to_continue = false;
                    for (auto& shell : c2check[c[th]]){
                        if(notExit(s[shell], s_n[shell], m)) {
                            to_continue = true;
                        }
                    }
                    if(to_continue) {
                        continue;
                    }

                }
                used.emplace_back(vth);
                result += searchCore(th+1, m, used, c, c_n, s, s_n, c2check);
                used.pop_back();
            }
        }
    }
    return result;
}

uint searchMatch(uint v1, uint v2){
    uint update_result = 0;
    for(auto& li : G[v1].LI){
        if(li.second){
            uint u1 = li.first;
            for(auto& candi : G[v1].cand[u1]){
                uint u2 = candi.first;
                if(u2>=0 && candi.second.find(v2)!=candi.second.end()){
                    ska::flat_hash_map<uint, uint> matching;
                    matching[u1] = v1;
                    matching[u2] = v2;
                    vec core_v;
                    core_v.emplace_back(v1);
                    core_v.emplace_back(v2);

                    update_result += searchCore(0, matching, core_v, matching_order[u1][u2].core,
                                                matching_order[u1][u2].core_nei, matching_order[u1][u2].shell,
                                                matching_order[u1][u2].shell_nei, matching_order[u1][u2].c_s_nei);
                }
            }
        }
    }
    return update_result;
}



uint Parallel_Search(uint v1, uint v2){

    uint update_result = 0;
    
    // Map
    std::vector<std::tuple<uint, uint, uint>> tasks; // (u1, u2, v2)
    
    for(auto& li : G[v1].LI) {
        if(li.second) {
            uint u1 = li.first;
            for(auto& candi : G[v1].cand[u1]) {
                uint u2 = candi.first;
                if(u2 >= 0 && candi.second.find(v2) != candi.second.end()) {
                    tasks.push_back(std::make_tuple(u1, u2, v2));
                }
            }
        }
    }
    
    // reduce
    #pragma omp parallel for reduction(+:update_result)
    for(int i = 0; i < tasks.size(); i++) {
        uint u1 = std::get<0>(tasks[i]);
        uint u2 = std::get<1>(tasks[i]);
        uint v2_val = std::get<2>(tasks[i]);
        
        ska::flat_hash_map<uint, uint> matching;
        matching[u1] = v1;
        matching[u2] = v2;
        vec core_v;
        core_v.emplace_back(v1);
        core_v.emplace_back(v2);
        
        uint local_result = searchCore(0, matching, core_v, matching_order[u1][u2].core,
                                      matching_order[u1][u2].core_nei, matching_order[u1][u2].shell,
                                      matching_order[u1][u2].shell_nei, matching_order[u1][u2].c_s_nei);
        
        update_result += local_result;
    }
    
    return update_result;

}


void inputUpdate(string& path, uint max_num){
    update.clear();
    ifstream infile(path);
    char c;
    uint v1, v2, w;
    uint cnt = 0;
    while (infile >> c >> v1 >> v2){
        if(max_num!=0 && ++cnt>max_num) break;
        update.push_back(v1);
        update.push_back(v2);
    }
}

void updateAndMatching() {
    long long add_matches = 0;
    long long del_matches = 0;
    double update_time = 0.0;
    double DCS_update_time = 0.0;
    double search_time = 0.0;
    
    // OpenMP time
    double start_time, end_time;

    for (uint t = 0; t < update.size(); t += 2) {
        uint v1 = update[t];
        uint v2 = update[t + 1];

        if(v1<0){
            if(G[-v1-1].label==-1 || G[-v2-1].label==-1) continue;
            if(!G[-v1-1].nei.count(-v2-1)) continue;
            if(turnOffProcess_safe(-1-v1, -v2-1)){
                continue;
            }
    
            // Auxiliary data structures update
            start_time = omp_get_wtime();
            del_matches += searchMatch(-v1-1, -v2-1);
            end_time = omp_get_wtime();
            search_time += (end_time - start_time) * 1000; // ms

            // FM
            start_time = omp_get_wtime();
            turnOffProcess(-v1-1, -v2-1);
            end_time = omp_get_wtime();
            DCS_update_time += (end_time - start_time) * 1000; // ms
        }
        else{
            if(G[v1].label==-1 || G[v2].label==-1) continue;
            if(G[v1].nei.count(v2)) continue;
            if(turnOnProcess_safe(v1,v2)){
                continue;
            }
            // Auxiliary data structures update
            start_time = omp_get_wtime();
            turnOnProcess_safe_parallel(v1, v2);
            end_time = omp_get_wtime();
            DCS_update_time += (end_time - start_time) * 1000; // ms

            // FM
            start_time = omp_get_wtime();
            add_matches += searchMatch(v1, v2);
            end_time = omp_get_wtime();
            search_time += (end_time - start_time) * 1000; // ms
        }
    }

    cerr << "added matches: " << add_matches << endl;
    cerr << "deleted matches: " << del_matches << endl;
    cerr << "Auxiliary Data Structure update time: " << DCS_update_time << "ms" << endl;
    cerr << "updated matches: " << add_matches + del_matches << endl;
    cerr << "search time: " << search_time << "ms" << endl;
}


void Parallel_batch_updateAndMatchingSingle() {
    long long add_matches = 0;
    long long del_matches = 0;
    double update_time = 0.0;
    double DCS_update_time = 0;
    double search_time = 0.0;
    
    // OpenMP
    double start_time, end_time;

    size_t window_size = 16;

    for (uint t = 0; t < update.size(); t += 2) {
        uint v1 = update[t];
        uint v2 = update[t + 1];

        if(v1<0){
            if(G[-v1-1].label==-1 || G[-v2-1].label==-1) continue;
            if(!G[-v1-1].nei.count(-v2-1)) continue;
            if(turnOffProcess_safe(-1-v1, -v2-1)){
                continue;
            }
    
            // FM
            start_time = omp_get_wtime();
            del_matches += Parallel_Search(-v1-1, -v2-1);
            // del_matches += searchMatch(-v1-1, -v2-1);
            end_time = omp_get_wtime();
            search_time += (end_time - start_time) * 1000; // ms

            // Auxiliary data structures update   
            start_time = omp_get_wtime();
            turnOffProcess(-v1-1, -v2-1);
            end_time = omp_get_wtime();
            DCS_update_time += (end_time - start_time) * 1000; // ms
        }
        else{
            // bool is_safe = false;
            if(G[v1].label==-1 || G[v2].label==-1) continue;
            if(G[v1].nei.count(v2)) continue;
            if(turnOnProcess_safe(v1,v2)){
                continue;
            }
            // Auxiliary data structures update
            start_time = omp_get_wtime();
            // turnOnProcess(v1, v2);
            turnOnProcess_safe_parallel(v1, v2);
            end_time = omp_get_wtime();
            DCS_update_time += (end_time - start_time) * 1000; //ms

            // FM
            start_time = omp_get_wtime();
            add_matches += Parallel_Search(v1, v2);
            // add_matches += searchMatch(v1, v2);
            end_time = omp_get_wtime();
            search_time += (end_time - start_time) * 1000; //ms
        }
    }
    cerr << "added matches: " << add_matches << endl;
    cerr << "deleted matches: " << del_matches << endl;
    cerr << "Auxiliary Data Structure update time: " << DCS_update_time << "ms" << endl;
    cerr << "updated matches: " << add_matches + del_matches << endl;
    // cerr << "update time: " << update_time << "ms" << endl;
    cerr << "search time: " << search_time << "ms" << endl;
}


void Parallel_batch_updateAndMatchingMoreDetect() {
    long long add_matches = 0;
    long long del_matches = 0;
    double update_time = 0.0;

    double DCS_update_time = 0;

    double search_time = 0.0;
    clock_t time;

    size_t window_size = 16;

    for (uint t = 0; t < update.size(); t += 2) {

        if (t + window_size < update.size()) {
            // #pragma omp parallel for
            for (uint i = t; i < t + window_size; i += 2) {
                uint v1 = update[i];
                uint v2 = update[i + 1];
        
                if(v1<0){
                    if(G[-v1-1].label==-1 || G[-v2-1].label==-1) continue;
                    if(!G[-v1-1].nei.count(-v2-1)) continue;
                    if(turnOffProcess_safe(-1-v1, -v2-1)){
                        continue;
                    }
        
                    // Auxiliary data structures update         
                    time = clock();
                    del_matches += searchMatch(-v1-1, -v2-1);
                    search_time += double(clock() - time)*1000 / CLOCKS_PER_SEC;
        
                    // FM
                    time = clock();
                    turnOffProcess(-v1-1, -v2-1);
                    DCS_update_time += double(clock() - time) * 1000 / CLOCKS_PER_SEC;
                }
        
                else{
                    // bool is_safe = false;
                    if(G[v1].label==-1 || G[v2].label==-1) continue;
                    if(G[v1].nei.count(v2)) continue;
                    if(turnOnProcess_safe(v1,v2)){
                        continue;
                    }
                    // Auxiliary data structures update
                    time = clock();
                    turnOnProcess_safe_parallel(v1, v2);
                    DCS_update_time += double(clock() - time)*1000 / CLOCKS_PER_SEC;
        
                    // FM
                    time = clock();
                    add_matches += searchMatch(v1, v2);
                    search_time += double(clock() - time)*1000 / CLOCKS_PER_SEC;
                }
            }
            t += window_size - 2; // Adjust the loop index to skip the processed window
        }

      
    }
    cerr << "added matches: " << add_matches << endl;
    cerr << "deleted matches: " << del_matches << endl;
    cerr << "Auxiliary Data Structure update time: " << DCS_update_time << "ms" << endl;
    cerr << "updated matches: " << add_matches + del_matches << endl;
    // cerr << "update time: " << update_time << "ms" << endl;
    cerr << "search time: " << search_time << "ms" << endl;
}




int main(int argc, char** argv){

    string g_path = "../dataset/dz_3/initial";
    string q_path = "../dataset/dz_3/Q/8/q10";
    string s_path = "../dataset/dz_3/s";
    uint cnum = 0;
    size_t thread_num = 32;

    for (uint i=1; i<argc; i++) {
        if(string(argv[i]) == "-d")
            g_path = argv[i+1];
        else if(string(argv[i]) == "-q")
            q_path = argv[i+1];
        else if(string(argv[i]) == "-c")
            cnum = atoi(argv[i+1]);
        else if(string(argv[i]) == "-s") {
            s_path = argv[i+1];
        }
        else if(string(argv[i]) == "-t") {
            thread_num = atoi(argv[i+1]);
        }
    }

    omp_set_num_threads(thread_num); // Set the number of threads to use

    clock_t time = clock();
    cerr << "========== Start inputting. ==========" << endl;
    inputQ(q_path);
    generateMO();
    inputG(g_path);
    inputUpdate(s_path, cnum);
    cerr << "Inputting cost (ms): " << double(clock() - time)*1000 / CLOCKS_PER_SEC << endl;
    cerr << "The number of vertices in Q: " << Q_size << endl;
    cerr << "The number of vertices in G: " << G_size << endl;
    cerr << "========== End inputting. ==========" << endl;

    time = clock();
    cerr << "========== Start constructing. ==========" << endl;
    constructCand();
    cerr << "Constructing cost (ms): " << double(clock() - time)*1000 / CLOCKS_PER_SEC << endl;
    cerr << "========== End constructing. ==========" << endl;

    time = clock();
    cerr << "========== Start static filtering. ==========" << endl;
    staticFilter();
    cerr << "Static filtering cost (ms): " << double(clock() - time)*1000 / CLOCKS_PER_SEC << endl;
    cerr << "========== End static filtering. ==========" << endl;

    double start_time = omp_get_wtime();
    if(thread_num == 1){
        // baseline
        cerr << "========== Start updating for single thread ==========" << endl;
        updateAndMatching();
        double end_time = omp_get_wtime();
        cerr << "Updating totally cost (ms): " << (end_time - start_time) * 1000 << endl;
        cerr << "========== End single thread updating. ==========" << endl;
    }else{
        cerr << "========== Start updating for "<< thread_num << " thread ==========" << endl;
        Parallel_batch_updateAndMatchingSingle();
        double end_time = omp_get_wtime();
        cerr << "Updating totally cost (ms): " << (end_time - start_time) * 1000 << endl;
        cerr << "========== End Multi tgihread updating. ==========" << endl;
    
        cerr << "========== DONE!! ==========" << endl;
    }


    return 0;
}