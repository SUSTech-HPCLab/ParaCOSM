#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>

#include <omp.h>

#include "utils/CLI11.hpp"
#include "utils/globals.h"
#include "utils/types.h"

#include "graph/graph.h"
#include "matching/matching.h"


#include "matching/SJTree/sj_tree.h"
#include "matching/TurboFlux/turboflux.h"
#include "matching/GraphFlow/graphflow.h"
#include "matching/SymBi/symbi.h"
#include "matching/Iedyn/iedyn.h"

#include "matching/Parallel_SymBi/parallel.h"
#include "matching/Parallel_TurboFlux/parallel_turboflux.h"
#include "matching/Parallel_GraphFlow/parallel_graphflow.h"


#define Print_Time_Nano_Main(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1e6 << " ms" << std::endl
#define Print_Time_Milli_Main(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << " ms" << std::endl


inline void BatchUpdates(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit){

    const size_t batch_size = 80000;
    const size_t sliding_widow_size = 200;
    const size_t classify_threads = 8;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);

    std::vector<bool> is_safes(batch_size, false);
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);

    bool safe = true;
    size_t idx_update = batch_size + 5;

    // size_t inc_size = data_graph.updates_.size();

    size_t main_size = batch_size / sliding_widow_size ;

    // size_t which_80000 = 0;

    while (!data_graph.updates_.empty() && !reach_time_limit) {
        batch_updates.clear();
        is_safes.assign(batch_size, false);
        num_e_update_vec.assign(batch_size, 0ul);

        size_t actual_threads = std::min(classify_threads, sliding_widow_size);

        // 80000
        #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
        {
            #pragma omp single
            // if(omp_get_thread_num() == 0)
            {
                // 从队列中提取一批更新
                for (size_t i = 0; i < batch_size && !data_graph.updates_.empty(); ++i) {
                    batch_updates.push_back(data_graph.updates_.front());
                    data_graph.updates_.pop();
                }
            }

            // 80000/200
            // size_t main_size = batch_size / sliding_widow_size + 1;
            for(size_t local_batch = 0; local_batch < main_size; ++local_batch){

                // 200/1
                #pragma omp for
                for(size_t local_window = 0; local_window < sliding_widow_size; ++local_window){
                    // size_t cnt = which_80000 * batch_size + local_batch * sliding_widow_size + local_window;
                    size_t cnt =  local_batch * sliding_widow_size + local_window;
                    if(cnt >= batch_updates.size()){
                        continue;
                    }

                    const auto& insert = batch_updates[cnt];
                    if (insert.type == 'e' && insert.is_add) {
                        is_safes[cnt] = mm->Classify(insert.id1, insert.id2, insert.label);
                        num_e_update_vec[cnt]++;


                        if (!is_safes[cnt]) {
                            #pragma omp critical
                            {
                                safe = false;
                                if (idx_update >= cnt) {
                                    idx_update = cnt;
                                    // mm->AddEdge(insert.id1, insert.id2, insert.label);
                                }
                            }
                            #pragma omp cancel for
                        }
                    }
                }

                #pragma omp barrier

                #pragma omp single
                {
                    // #pragma omp critical
                    {
                // std::cout << "safe: " << safe << std::endl;
                        if (!safe) { // correct, without update, 330ms
                            // if(omp_get_thread_num() == 1)
                            {
                                // size_t cnt = which_80000 * batch_size + local_batch * sliding_widow_size;
                                // size_t cnt =  local_batch * sliding_widow_size;
                                for (size_t i = 0; i < sliding_widow_size; i++) {
                                    // std::cout<<i<<std::endl;
                                    if(i < batch_updates.size()){
                                        const auto& insert = batch_updates[local_batch* sliding_widow_size+i];
                                        mm->AddEdge(insert.id1, insert.id2, insert.label);
                                    }
                                    
                                }
                            }
                        }
                    }
                }


            }
        }



        for (size_t i = 0; i < num_e_update_vec.size(); i++) {
            num_e_updates += num_e_update_vec[i];
        }

        // which_80000++;
    }
}


inline void BatchUpdates2(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
        size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
        size_t& negative_num_results_last, std::atomic_bool& reach_time_limit){
    
        const size_t batch_size = 400;
        const size_t classify_threads = 16;

        std::vector<InsertUnit> batch_updates;
        batch_updates.reserve(batch_size);
        size_t idx_update = batch_size + 5;
        std::vector<bool> is_safes(batch_size, false);
        bool safe = true;
        std::vector<size_t> num_e_update_vec(batch_size, 0ul);
        // constexpr size_t inner_batch_size = 200;
    
        while (!data_graph.updates_.empty() && !reach_time_limit) {
            batch_updates.clear();
            safe = true;
            is_safes.assign(batch_size, false);        
            num_e_update_vec.assign(batch_size, 0ul);
            
            size_t actual_threads = std::min(classify_threads, batch_updates.size());
    
            #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
            {
                #pragma omp single
                {
                    // get from queue
                    for (size_t i = 0; i < batch_size && !data_graph.updates_.empty(); ++i) {
                        batch_updates.push_back(data_graph.updates_.front());
                        data_graph.updates_.pop();
                        // num_e_update_vec[i]++;
                        num_e_updates++;
                    }
                }
    
                #pragma omp for 
                for(size_t cnt = 0; cnt < batch_updates.size(); cnt++){
                    if(!safe){
                        continue;
                    }

                    const auto& insert = batch_updates[cnt];
                    if (insert.type == 'e' && insert.is_add) {
                        is_safes[cnt] = mm->Classify(insert.id1, insert.id2, insert.label);
                        // num_e_update_vec[cnt]++;
                        if (!is_safes[cnt]) {
                            // #pragma omp atomic write
                            
                            // #pragma omp critical
                            {
                                safe = false;
                                if(idx_update >= cnt){
                                    #pragma omp critical
                                    idx_update = cnt;
                                }   
                            }
                            // #pragma omp cancel for
                        }
                    }
                }
                // #pragma omp barrier

                // #pragma omp single
                // {
                //     if(!safe){ // correct, without update, 330ms
                //         for(size_t i = 0; i < batch_updates.size(); i++){
                //             const auto& insert = batch_updates[i];
                //             mm->AddEdge(insert.id1, insert.id2, insert.label);
                //         }
                //     }
                // }
            }
    
            if(!safe){ // correct, without update, 330ms
                // const auto& insert = batch_updates[idx_update];
                // mm->AddEdge(insert.id1, insert.id2, insert.label);
                for(size_t i = idx_update; i < batch_updates.size(); i++){
                    const auto& insert = batch_updates[i];
                    mm->AddEdge(insert.id1, insert.id2, insert.label);
                    
                    size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
                    mm->GetNumPositiveResults(positive_num_results_cur);
                    mm->GetNumNegativeResults(negative_num_results_cur);
                    if (positive_num_results_cur != positive_num_results_last || negative_num_results_cur != negative_num_results_last) {
                        positive_num_results_last = positive_num_results_cur;
                        negative_num_results_last = negative_num_results_cur;

                        unsafe_updates++;
                        // std::cout << "Unsafe update number: " << unsafe_updates << "At" << count << std::endl;
                    }
                }

            }


            // sum:
            // num_e_updates += std::accumulate(num_e_update_vec.begin(), num_e_update_vec.end(), 0ul);
            // for(size_t i = 0; i < num_e_update_vec.size(); i++){
            //     num_e_updates += num_e_update_vec[i];
            // }
        }
    }


inline void BatchUpdates3(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit){
    
    size_t sliding_window_base = 0;
    const size_t window_size = 16;
    const size_t update_size = data_graph.updates_vec_.size();

    std::vector<bool> update_safes(window_size+1, false);
    size_t min_safe = window_size + 10;
    bool is_safes = true;

    while(sliding_window_base < update_size){
        is_safes = true;
        min_safe = window_size + 10;

        // #pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < window_size; i++){

            if(i + sliding_window_base >= update_size){
                continue;
            }
            const auto& insert = data_graph.updates_vec_[i+sliding_window_base];

            // std::cout << "insert: " << insert.id1 << " " << insert.id2 << " " << insert.label << std::endl;

            if (insert.type == 'e' && insert.is_add) {
                // std::cout << "insert edge: " << insert.id1 << " " << insert.id2 << " " << insert.label << std::endl;
                update_safes[i] = mm->Classify(insert.id1, insert.id2, insert.label);
                // std::cout << "update_safes[i]: " << update_safes[i] << std::endl;

                if(!update_safes[i]){
                    // #pragma omp critical
                    {
                        is_safes = false;
                        if(min_safe > i){
                            min_safe = i;
                        }
                        // std::cout << "unsafe PHD at "<< min_safe+sliding_window_base << std::endl;
                        break;
                    }
                }
            }
            else if (insert.type == 'v' && insert.is_add)
            {
                mm->AddVertex(insert.id1, insert.label);

            }
            else if (insert.type == 'v' && !insert.is_add)
            {
                // #pragma omp critical
                {
                is_safes = false;
                if(min_safe > i){
                    min_safe = i;
                }
                // std::cout << "unsafe PHD at "<< min_safe+sliding_window_base << std::endl;
                break;
                }
                // mm->RemoveVertex(insert.id1);
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                update_safes[i] = mm->Classify(insert.id1, insert.id2, insert.label);
                if(!update_safes[i]){
                    // #pragma omp critical
                    {
                        is_safes = false;
                        if(min_safe > i){
                            min_safe = i;
                        }
                        break;
                    }
                }
            }

            
     

        }

        // for(size_t i =0; i<update_safes.size(); i++ ){
        //     if(!update_safes[i]){
        //         min_safe = i;
        //         is_safes = false;
        //         break;
        //     } 
        // }

        if(!is_safes){
            // std::cout << "unsafe at "<< min_safe+sliding_window_base << std::endl;
            const auto& insert_unsafe = data_graph.updates_vec_[min_safe+sliding_window_base];

            if(insert_unsafe.type == 'e' && insert_unsafe.is_add){
                mm ->AddEdge(insert_unsafe.id1, insert_unsafe.id2, insert_unsafe.label);
            }
            else if(insert_unsafe.type == 'v' && !insert_unsafe.is_add){
                mm->RemoveVertex(insert_unsafe.id1);
            }
            
            // if(insert_unsafe.type == 'e' && insert_unsafe.is_add){
            //     mm->AddEdge(insert_unsafe.id1, insert_unsafe.id2, insert_unsafe.label);
            // }

            // else if(insert_unsafe.type == 'v' && insert_unsafe.is_add){
            //     mm->AddVertex(insert_unsafe.id1, insert_unsafe.label);
            // }else if(insert_unsafe.type == 'v' && !insert_unsafe.is_add){
            //     mm->RemoveVertex(insert_unsafe.id1);
            // }else if(insert_unsafe.type == 'e' && !insert_unsafe.is_add){
            //     mm->RemoveEdge(insert_unsafe.id1, insert_unsafe.id2);
            // }
            
            sliding_window_base += min_safe + 1; // add 1

            num_e_updates += min_safe + 1;
            
            // std::cout << "sliding_window_base: " << sliding_window_base << "detecting:" << std::endl;

            // std::cout << "sliding_window_base: " << sliding_window_base << "detecting:" << std::endl;
            size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
            mm->GetNumPositiveResults(positive_num_results_cur);
            mm->GetNumNegativeResults(negative_num_results_cur);
            if (positive_num_results_cur != positive_num_results_last || negative_num_results_cur != negative_num_results_last) {
                positive_num_results_last = positive_num_results_cur;
                negative_num_results_last = negative_num_results_cur;

                unsafe_updates++;
                // std::cout << "Unsafe update number: " << unsafe_updates << "At" << count << std::endl;
            }

        }else{
            // std::cout << "safe at "<< sliding_window_base << std::endl;
            sliding_window_base += window_size;
            num_e_updates += window_size;
        }

        // update progress
        // if((update_size - sliding_window_base) % (update_size / 10) < (window_size-1) ){
        //     std::cout << "update progress: " << (sliding_window_base) * 100 / update_size << "%" << std::endl;
        // }

        // std::cout << "END" << std::endl;
    }

}




inline void BatchUpdates_OpenMP(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit){
    
    size_t sliding_window_base = 0;
    const size_t window_size = 160;
    const size_t update_size = data_graph.updates_vec_.size();

    std::vector<bool> update_safes(window_size+1, false);
    size_t min_safe = window_size + 10;
    bool is_safes = true;

    while(sliding_window_base < update_size){
        is_safes = true;
        min_safe = window_size + 10;

        #pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < window_size; i++){

            if(i + sliding_window_base >= update_size){
                continue;
            }
            const auto& insert = data_graph.updates_vec_[i+sliding_window_base];

            if (insert.type == 'e' && insert.is_add) {
                update_safes[i] = mm->Classify(insert.id1, insert.id2, insert.label);

                if(!update_safes[i]){
                    {
                        is_safes = false;
                        if(min_safe > i){
                            min_safe = i;
                        }
                    }
                }
            }
            else if (insert.type == 'v' && insert.is_add)
            {
                mm->AddVertex(insert.id1, insert.label);

            }
            else if (insert.type == 'v' && !insert.is_add)
            {
                #pragma omp critical
                {
                    is_safes = false;
                    if(min_safe > i){
                        min_safe = i;
                    }
                }
                
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                update_safes[i] = mm->Classify(insert.id1, insert.id2, insert.label);
                if(!update_safes[i]){
                    #pragma omp critical
                    {
                        is_safes = false;
                        if(min_safe > i){
                            min_safe = i;
                        }
                    }
                }
            }

            
     

        }

        if(!is_safes){
            const auto& insert_unsafe = data_graph.updates_vec_[min_safe+sliding_window_base];

            if(insert_unsafe.type == 'e' && insert_unsafe.is_add){
                mm ->AddEdge(insert_unsafe.id1, insert_unsafe.id2, insert_unsafe.label);
            }
            else if(insert_unsafe.type == 'v' && !insert_unsafe.is_add){
                mm->RemoveVertex(insert_unsafe.id1);
            }
            
            
            sliding_window_base += min_safe + 1; // add 1

            num_e_updates += min_safe + 1;

            size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
            mm->GetNumPositiveResults(positive_num_results_cur);
            mm->GetNumNegativeResults(negative_num_results_cur);
            if (positive_num_results_cur != positive_num_results_last || negative_num_results_cur != negative_num_results_last) {
                positive_num_results_last = positive_num_results_cur;
                negative_num_results_last = negative_num_results_cur;

                unsafe_updates++;
            }

        }else{

            sliding_window_base += window_size;
            num_e_updates += window_size;
        }

        // update progress
        if((update_size - sliding_window_base) % (update_size / 10) < (window_size-1) ){
            std::cout << "update progress: " << (sliding_window_base) * 100 / update_size << "%" << std::endl;
        }

    }

}




inline void ProcessBatchUpdates(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit){

    const size_t batch_size = 10000;
    const size_t classify_threads = 8;

    std::vector<InsertUnit> batch_updates;
    batch_updates.reserve(batch_size);

    std::vector<bool> is_safes(batch_size, false);
    std::vector<size_t> num_e_update_vec(batch_size, 0ul);

    // for (size_t i = 0; i < batch_size && !data_graph.updates_.empty(); ++i) {
    //     batch_updates.push_back(data_graph.updates_.front());
    //     data_graph.updates_.pop();
    // }



    // size_t actual_threads = std::min(classify_threads, batch_updates.size());

    // #pragma omp parallel for num_threads(actual_threads)
    // for(size_t cnt = 0; cnt < batch_updates.size() ; cnt++){
    //     const auto& insert = batch_updates[cnt];
    //     if (insert.type == 'e' && insert.is_add) {
    //         is_safes[cnt] = mm->Classify(insert.id1, insert.id2, insert.label);
    //         num_e_update_vec[cnt]++;
    //     }
    // }

    bool safe = true;
    // idx = max
    size_t idx_update = batch_size+5;

    size_t actual_threads = std::min(classify_threads, batch_updates.size());

    #pragma omp parallel num_threads(actual_threads + 1) shared(safe, idx_update)
    {
        #pragma omp single
        {
            // 从队列中提取一批更新
            for (size_t i = 0; i < batch_size && !data_graph.updates_.empty(); ++i) {
                batch_updates.push_back(data_graph.updates_.front());
                data_graph.updates_.pop();
            }
        }

        #pragma omp for
        for(size_t cnt = 0; cnt < batch_updates.size(); cnt++){
            const auto& insert = batch_updates[cnt];
            if (insert.type == 'e' && insert.is_add) {
                is_safes[cnt] = mm->Classify(insert.id1, insert.id2, insert.label);
                num_e_update_vec[cnt]++;
                if (!is_safes[cnt]) {
                    #pragma omp critical
                    {
                        safe = false;
                        if(idx_update >= cnt){
                            idx_update = cnt;
                        }   
                    }
                    // #pragma omp cancel for
                }
                
            }
            if (insert.type == 'v' && insert.is_add)
            {
                mm->AddVertex(insert.id1, insert.label);
            }
            else if (insert.type == 'v' && !insert.is_add)
            {
                mm->RemoveVertex(insert.id1);
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                is_safes[cnt] = mm->Classify(insert.id1, insert.id2, insert.label);
                num_e_update_vec[cnt]++;
                if (!is_safes[cnt]) {
                    #pragma omp critical
                    {
                        safe = false;
                        if(idx_update >= cnt){
                            idx_update = cnt;
                        }   
                    }
                    // #pragma omp cancel for
                }
            }
        }
    }



    // kkkk
    if(!safe){ // correct, without update, 330ms
        for(size_t i = 0; i < batch_updates.size(); i++){
            const auto& insert = batch_updates[i];
            mm->AddEdge(insert.id1, insert.id2, insert.label);
        }
    }



    for(size_t i =0; i< num_e_update_vec.size();i++){
        num_e_updates += num_e_update_vec[i];
    }

  
}


inline void SingleThreadUpdate(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit)
{
    size_t num_updates = data_graph.updates_.size();
    std::cout << "Start with " << num_updates << " Update" << std::endl;

    for (size_t i = 0; i < num_updates; i++)
    {
        InsertUnit insert = data_graph.updates_.front();
        data_graph.updates_.pop();

        if (insert.type == 'v' && insert.is_add)
        {
            mm->AddVertex(insert.id1, insert.label);
            num_v_updates++;
        }
        else if (insert.type == 'v' && !insert.is_add)
        {
            mm->RemoveVertex(insert.id1);
            num_v_updates++;
        }
        else if (insert.type == 'e' && insert.is_add)
        {
            mm->AddEdge(insert.id1, insert.id2, insert.label);
            num_e_updates++;
        }
        else if (insert.type == 'e' && !insert.is_add)
        {
            mm->RemoveEdge(insert.id1, insert.id2);
            num_e_updates++;
        }
        if (reach_time_limit) break;

        // count unsafe update
        size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
        mm->GetNumPositiveResults(positive_num_results_cur);
        mm->GetNumNegativeResults(negative_num_results_cur);
        if (positive_num_results_cur != positive_num_results_last || negative_num_results_cur != negative_num_results_last)
        {
            positive_num_results_last = positive_num_results_cur;
            negative_num_results_last = negative_num_results_cur;
            unsafe_updates++;
        }
        count += 1;
    }
}


int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string query_path = "", initial_path = "", stream_path = "", algorithm = "none";
    uint max_num_results = UINT_MAX, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
    bool print_prep = true, print_enum = false, homo = false, report_initial = true;
    std::vector<std::vector<uint>> orders;

    size_t thread_num = 6;

    size_t auto_tuning = 1;

    app.add_option("-q,--query", query_path, "query graph path")->required();
    app.add_option("-d,--data", initial_path, "initial data graph path")->required();
    app.add_option("-u,--update", stream_path, "data graph update stream path")->required();
    app.add_option("-a,--algorithm", algorithm, "algorithm");
    app.add_option("--max-results", max_num_results, "max number of results for one edge update");
    app.add_option("--time-limit", time_limit, "time limit for the incremental matching (second)");
    app.add_option("--print-prep", print_prep, "print processing results or not");
    app.add_option("--print-enum", print_enum, "print enumeration results or not");
    app.add_option("--homo", homo, "using graph homomorphism");
    app.add_option("--report-initial", report_initial, "report the result of initial matching or not");
    app.add_option("--initial-time-limit", initial_time_limit, "time limit for the initial matching (second)");
    app.add_option("--orders", orders, "pre-defined matching orders");
    app.add_option("-t,--thread-num", thread_num, "Number of thread that program will use.");
    app.add_option("--auto-tuning", auto_tuning, "Framework will tune the thread number with query vertex");

    
    CLI11_PARSE(app, argc, argv);

    // std::cout << "use " << thread_num << " threads for incresemental matching " << std::endl;
    
    std::chrono::high_resolution_clock::time_point start;

    start = My_Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;
    Graph query_graph {};
    query_graph.LoadFromFile(query_path);
    query_graph.PrintMetaData();

    Graph data_graph {};
    data_graph.LoadFromFile(initial_path);
    data_graph.PrintMetaData();
    Print_Time("Load Graphs: ", start);

    std::cout << "------------ Preprocessing ------------" << std::endl;
    matching *mm = nullptr;

    //IncIsoMatch *incIsoMatch = nullptr;
    SJTree *sjtree = nullptr;
    Graphflow *graphflow = nullptr;
    TurboFlux *turboflux = nullptr;
    SymBi *symbi = nullptr;
    IEDyn *iedyn = nullptr;

    Parrllel_SymBi *parrallel = nullptr;
    Parallel_TurboFlux *parallel_turboflux = nullptr;
    Parallel_Graphflow *parallel_graphflow = nullptr;
    
    start = My_Get_Time();

    // Single Thread Version
    if (algorithm == "sj-tree")
        mm = sjtree         = new SJTree        (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);
    else if (algorithm == "graphflow")
        mm = graphflow      = new Graphflow     (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);
    else if (algorithm == "turboflux")
        mm = turboflux      = new TurboFlux     (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);
    else if (algorithm == "symbi")
        mm = symbi          = new SymBi         (query_graph, data_graph, max_num_results, print_prep, print_enum, homo, orders);
    else if (algorithm == "iedyn")
        mm = iedyn          = new IEDyn         (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);

    // Parallel Version
    else if (algorithm == "parallel_symbi")
        mm = parrallel      = new Parrllel_SymBi(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, orders, thread_num, auto_tuning);
    else if (algorithm == "parallel_turboflux")
        mm = parallel_turboflux = new Parallel_TurboFlux(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
    else if (algorithm == "parallel_graphflow")
        mm = parallel_graphflow = new Parallel_Graphflow (query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
    else if (algorithm == "none")
        mm                  = new matching      (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);
    else
    {
        std::cout << "Unknown algorithm" << std::endl;
        exit(-1);
    }

    
    mm->Preprocessing();
    Print_Time("Preprocessing: ", start);

    if (report_initial || algorithm == "sj-tree")
    {
        std::cout << "----------- Initial Matching ----------" << std::endl;
        
        start = My_Get_Time();
        auto InitialFun = [&mm]()
        {
            mm->InitialMatching();
        };
        execute_with_time_limit(InitialFun, initial_time_limit, reach_time_limit);
        Print_Time("Initial Matching: ", start);
        
        size_t num_results = 0ul;
        mm->GetNumInitialResults(num_results);
        std::cout << num_results << " matches. \n";

        if (reach_time_limit) return 1;
    }
    

    {
        // #pragma omp single
        {

        std::cout << "--------- Incremental Matching --------" << std::endl;

        if(auto_tuning == 1){
            std::cout << "auto_tuning on:" << auto_tuning << std::endl;
        }else if(auto_tuning == 0){
            std::cout << "auto_tuning off:" << auto_tuning << std::endl;
            std::cout << "use " << thread_num << " threads for incresemental matching " << std::endl;
        }
        // 




        data_graph.LoadUpdateStream(stream_path);
            
        size_t num_v_updates = 0ul, num_e_updates = 0ul;

            // count unsafe updates
        size_t positive_num_results_last = 0ul, negative_num_results_last = 0ul;
        // size_t positive_num_results_cur = 0ul, negative_num_results_cur = 0ul;
        size_t unsafe_updates = 0ul;
            
        size_t count = 0ul;


        // auto timer = std::chrono::high_resolution_clock::now();
        start = My_Get_Time();

        BatchUpdates3(data_graph, mm, num_v_updates, num_e_updates, unsafe_updates, count, 
            positive_num_results_last, negative_num_results_last, reach_time_limit);

        // while (!data_graph.updates_.empty() && !reach_time_limit) {
            // ProcessBatchUpdates(data_graph, mm, num_v_updates, num_e_updates, unsafe_updates, count, 
            //     positive_num_results_last, negative_num_results_last, reach_time_limit);
        // }

        // print ms 
        // Print_Time_Nano_Main("Incremental Matching: ", timer);

        // print ms of matching

        // single

        // auto IncrementalFun = [&]()
        // {
            // SingleThreadUpdate(data_graph, mm, num_v_updates, num_e_updates, 
            //     unsafe_updates, count, positive_num_results_last, negative_num_results_last, reach_time_limit);
        // };
            
            

            // Test if out of time
            // execute_with_time_limit(IncrementalFun, time_limit, reach_time_limit);
            
            Print_Time("Incremental Matching: ", start);

            std::cout << num_v_updates << " vertex updates.\n";
            std::cout << num_e_updates << " edge updates.\n";
            std::cout << unsafe_updates << " unsafe updates.\n";

            size_t positive_num_results = 0ul, negative_num_results = 0ul;
            mm->GetNumPositiveResults(positive_num_results);
            mm->GetNumNegativeResults(negative_num_results);
            std::cout << positive_num_results << " positive matches.\n";
            std::cout << negative_num_results << " negative matches.\n";

            mm->PrintCounter();

            size_t num_edges = 0u, num_vertices = 0ul;
            mm->GetMemoryCost(num_edges, num_vertices);
            std::cout << "\n# edges in index in total: " << num_edges;
            std::cout << "\n# vertices in index in total: " << num_vertices;

            std::cout << "\nPeak Virtual Memory: " << mem::getValue() << " KB";
            std::cout << "\n\n----------------- End -----------------" << std::endl;

            //delete ;
            delete mm;

        }
    }
}