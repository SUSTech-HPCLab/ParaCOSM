#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <filesystem>
#include <algorithm>

#include <omp.h>

#include "utils/CLI11.hpp"
#include "utils/globals.h"
#include "utils/types.h"

#include "graph_storage/graph.h"
#include "matching_executor/matching.h"


#include "matching_executor/SJTree/sj_tree.h"
#include "matching_executor/TurboFlux/turboflux.h"
#include "matching_executor/GraphFlow/graphflow.h"
#include "matching_executor/SymBi/symbi.h"
#include "matching_executor/Iedyn/iedyn.h"

#include "matching_executor/Parallel_SymBi/parallel_symbi.h"
#include "matching_executor/Parallel_TurboFlux/parallel_turboflux.h"
#include "matching_executor/Parallel_GraphFlow/parallel_graphflow.h"
#include "matching_executor/Parallel_CaLiG/parallel_calig_class.h"
#include "matching_executor/Parallel_NewSP/newsp_adapter.h"

#include "core/inter_executor/inter_executor.h"

#define Print_Time_Nano_Main(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1e6 << " ms" << std::endl
#define Print_Time_Milli_Main(str, start) std::cout << str << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << " ms" << std::endl





// // Legacy function wrapper for backward compatibility
// // This maintains the exact same interface as the original BatchUpdates3 function
// inline void BatchUpdates3_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates, 
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last, 
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {
    
//     InterExecutor executor(data_graph, mm);
//     executor.BatchUpdates3(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// // Wrapper for InterExecutor::BatchUpdates (port of free BatchUpdates)
// inline void BatchUpdates_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {

//     InterExecutor executor(data_graph, mm);
//     executor.BatchUpdates(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// // Wrapper for InterExecutor::BatchUpdates2 (port of free BatchUpdates2)
// inline void BatchUpdates2_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {

//     InterExecutor executor(data_graph, mm);
//     executor.BatchUpdates2(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// // Wrapper for InterExecutor::BatchUpdates_OpenMP (port of free BatchUpdates_OpenMP)
// inline void BatchUpdates_OpenMP_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {

//     InterExecutor executor(data_graph, mm);
//     executor.BatchUpdates_OpenMP(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// // Wrapper for InterExecutor::ProcessBatchUpdatesQueue (port of free ProcessBatchUpdates)
// inline void ProcessBatchUpdates_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {

//     InterExecutor executor(data_graph, mm);
//     executor.ProcessBatchUpdatesQueue(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// // Wrapper for InterExecutor::SingleThreadUpdate (port of free SingleThreadUpdate)
// inline void SingleThreadUpdate_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
//     size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
//     size_t& negative_num_results_last, std::atomic_bool& reach_time_limit) {

//     InterExecutor executor(data_graph, mm);
//     executor.SingleThreadUpdate(
//         num_v_updates, num_e_updates, unsafe_updates, count,
//         positive_num_results_last, negative_num_results_last,
//         reach_time_limit
//     );
// }

// Unified wrapper: dispatch by update_mode string and call the corresponding InterExecutor method
inline void RunUpdates_InterExecutor(Graph& data_graph, matching* mm, size_t& num_v_updates, size_t& num_e_updates,
    size_t& unsafe_updates, size_t& count, size_t& positive_num_results_last,
    size_t& negative_num_results_last, std::atomic_bool& reach_time_limit, const std::string& update_mode,
    size_t num_threads = 8) {

    InterExecutor executor(data_graph, mm);

    if (update_mode == "batch") {
        executor.BatchUpdates(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "batch2") {
        executor.BatchUpdates2(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "batch3") {
        executor.BatchUpdates3(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "batch4") {
        executor.BatchUpdates4(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "openmp") {
        executor.BatchUpdates_OpenMP(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "queue") {
        executor.ProcessBatchUpdatesQueue(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "single") {
        executor.SingleThreadUpdate(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "persistent") {
        executor.BatchUpdates_Persistent(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit, num_threads);
    } else if (update_mode == "batch_all") {
        executor.BatchUpdates_AllAtOnce(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit, num_threads);
    } else if (update_mode == "gpu") {
        executor.BatchUpdates_GPU(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "gpu_all") {
        executor.BatchUpdates_GPU_AllAtOnce(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "gpu_bfs") {
        executor.BatchUpdates_GPU_BFS(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "gpu_bfs_single") {
        executor.BatchUpdates_GPU_BFS_Single(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else if (update_mode == "versioned") {
        executor.BatchUpdates_Versioned(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit, num_threads);
    } else if (update_mode == "gpu_bfs_versioned") {
        executor.BatchUpdates_GPU_BFS_Versioned(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    } else {
        // default fallback
        executor.BatchUpdates3(
            num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit);
    }
}


// ---------------------------------------------------------------------------
// Batch helper: run ONE query against a fresh copy of the (already-parsed) data
// graph, applying the update stream, and print a compact one-line result.
// Used by --query-dir so the 50s data-graph parse happens only once.
// `pristine_data` is copied (deep) per call so each query starts from the same
// graph state; the update stream is provided pre-parsed in `stream_updates`.
// Returns positive match count (UINT64_MAX on timeout).
// ---------------------------------------------------------------------------
static uint64_t run_one_query(
    const std::string& query_path,
    const Graph& pristine_data,
    const std::vector<InsertUnit>& stream_updates,
    const std::string& algorithm,
    const std::string& update_mode,
    uint max_num_results, bool print_prep, bool print_enum, bool homo,
    size_t thread_num, size_t auto_tuning, uint time_limit)
{
    Graph query_graph {};
    query_graph.LoadFromFile(query_path);

    // Fresh copy of the data graph for this query (deep copy of vectors — fast,
    // ~1s for LJ vs ~50s to re-parse from disk).
    Graph data_graph = pristine_data;

    matching* mm = nullptr;
    Parallel_Graphflow* pg = nullptr;
    Parrllel_SymBi* ps = nullptr;
    Parallel_TurboFlux* pt = nullptr;
    if (algorithm == "parallel_graphflow")
        mm = pg = new Parallel_Graphflow(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
    else if (algorithm == "parallel_symbi")
        mm = ps = new Parrllel_SymBi(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, {}, thread_num, auto_tuning);
    else if (algorithm == "parallel_turboflux")
        mm = pt = new Parallel_TurboFlux(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
    else {
        std::cout << "[batch] unsupported algorithm for --query-dir: " << algorithm << std::endl;
        return 0;
    }
    (void)pg; (void)ps; (void)pt;

    mm->Preprocessing();

    // Load the (pre-parsed) update stream into this copy.
    data_graph.updates_vec_ = stream_updates;

    size_t num_v = 0, num_e = 0, unsafe = 0, count = 0, pos_last = 0, neg_last = 0;
    std::atomic_bool reach_time_limit{false};

    auto run = [&]() {
        RunUpdates_InterExecutor(data_graph, mm, num_v, num_e, unsafe, count,
            pos_last, neg_last, reach_time_limit, update_mode, thread_num);
    };
    execute_with_time_limit(run, time_limit, reach_time_limit);

    uint64_t pos = 0;
    if (reach_time_limit) { delete mm; return UINT64_MAX; }
    size_t p = 0; mm->GetNumPositiveResults(p); pos = p;
    delete mm;
    return pos;
}


int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string query_path = "", initial_path = "", stream_path = "", algorithm = "none";
    std::string query_dir = "", csv_out = "";
    std::string update_mode = "batch3"; // batch, batch2, batch3, batch4, openmp, queue, single
    uint max_num_results = UINT_MAX, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
    bool print_prep = true, print_enum = false, homo = false, report_initial = true;
    std::vector<std::vector<uint>> orders;

    size_t thread_num = 6;

    size_t auto_tuning = 1;

    app.add_option("-q,--query", query_path, "query graph path");
    app.add_option("--query-dir", query_dir, "directory of query graphs (batch: load data once, run all)");
    app.add_option("--csv-out", csv_out, "append batch results to this CSV file");
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
    app.add_option("-m,--update-mode", update_mode, "Update strategy: batch | batch2 | batch3 | batch4 | openmp | queue | single | persistent | batch_all");

    
    CLI11_PARSE(app, argc, argv);

    // ----- Batch mode: --query-dir loads data + stream ONCE, runs all queries -----
    if (!query_dir.empty()) {
        auto t_load = My_Get_Time();
        std::cout << "----------- Batch: loading data once ------------" << std::endl;
        Graph pristine_data {};
        pristine_data.LoadFromFile(initial_path);
        pristine_data.PrintMetaData();
        pristine_data.LoadUpdateStream(stream_path);
        std::vector<InsertUnit> stream_updates = pristine_data.updates_vec_;
        Print_Time("Load data + stream: ", t_load);

        // collect & sort query files
        std::vector<std::string> qfiles;
        for (auto& e : std::filesystem::directory_iterator(query_dir))
            if (e.is_regular_file()) qfiles.push_back(e.path().string());
        std::sort(qfiles.begin(), qfiles.end());
        std::cout << "[batch] " << qfiles.size() << " queries in " << query_dir << std::endl;

        FILE* csv = nullptr;
        if (!csv_out.empty()) { csv = fopen(csv_out.c_str(), "a"); }

        for (size_t i = 0; i < qfiles.size(); i++) {
            auto t_q = My_Get_Time();
            uint64_t pos = run_one_query(qfiles[i], pristine_data, stream_updates,
                algorithm, update_mode, max_num_results, false, false, homo,
                thread_num, auto_tuning, time_limit);
            double ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                My_Get_Time() - t_q).count();
            std::string qname = std::filesystem::path(qfiles[i]).filename().string();
            if (pos == UINT64_MAX)
                printf("[batch] %4zu/%zu  %-12s  TIMEOUT      (%.0fms)\n", i+1, qfiles.size(), qname.c_str(), ms);
            else
                printf("[batch] %4zu/%zu  %-12s  %15llu  (%.0fms)\n", i+1, qfiles.size(), qname.c_str(),
                       (unsigned long long)pos, ms);
            fflush(stdout);
            if (csv) {
                fprintf(csv, "%s,%s,%s,%s\n", query_dir.c_str(), qname.c_str(),
                        (pos==UINT64_MAX? "TIMEOUT" : std::to_string(pos).c_str()),
                        update_mode.c_str());
                fflush(csv);
            }
        }
        if (csv) fclose(csv);
        std::cout << "[batch] done." << std::endl;
        return 0;
    }

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
    Parallel_CaLiG *parallel_calig = nullptr;
    
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
    else if (algorithm == "parallel_calig")
        mm = parallel_calig = new Parallel_CaLiG(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
    else if (algorithm == "parallel_newsp")
        mm = new Parallel_NewSP_Adapter(query_graph, data_graph, max_num_results, print_prep, print_enum, homo, thread_num, auto_tuning);
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

        // Select update strategy (mode variable removed as it's not used)

        RunUpdates_InterExecutor(
            data_graph, mm, num_v_updates, num_e_updates, unsafe_updates, count,
            positive_num_results_last, negative_num_results_last, reach_time_limit, update_mode,
            thread_num);

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