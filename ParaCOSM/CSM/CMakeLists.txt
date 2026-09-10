cmake_minimum_required(VERSION 3.16)

# -----------------------------
# 1. Project declaration (CXX + CUDA)
# -----------------------------
project(CSM LANGUAGES CXX CUDA)

# -----------------------------
# 2. Compilation options
# -----------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA settings
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # Build a portable binary for the machines used by the project:
    # Volta V100 (sm_70) and Ampere A100 (sm_80). Users can override this with
    # -DCMAKE_CUDA_ARCHITECTURES=<arch-list>.
    set(CMAKE_CUDA_ARCHITECTURES "70;80" CACHE STRING
        "CUDA architectures to build for (70=V100, 80=A100)")
endif()

# -----------------------------
# 3. Dependencies
# -----------------------------
find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(CUDAToolkit REQUIRED)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CUDA Toolkit: ${CUDAToolkit_VERSION}")
message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")

# -----------------------------
# 4. Source files
# -----------------------------
set(MATCH matching_executor)
set(UTILS utils)
set(GRAPH graph_storage)

set(IEDYN ${MATCH}/Iedyn)
set(GRAPHFLOW ${MATCH}/GraphFlow)
set(SYMBI ${MATCH}/SymBi)
set(SJ_TREE ${MATCH}/SJTree)
set(TURBOFLUX ${MATCH}/TurboFlux)

set(PARALLEL_SYMBI ${MATCH}/Parallel_SymBi)
set(PARALLEL_TURBOFLUX ${MATCH}/Parallel_TurboFlux)
set(PARALLEL_GRAPHFLOW ${MATCH}/Parallel_GraphFlow)
set(PARALLEL_NEWSP ${MATCH}/Parallel_NewSP)
set(PARALLEL_CALIG ${MATCH}/Parallel_CaLiG)

set(SRC_FILES
    ${MATCH}/main.cpp
    ${MATCH}/matching.cpp

    ${IEDYN}/iedyn.cpp
    ${SYMBI}/symbi.cpp
    ${TURBOFLUX}/turboflux.cpp
    ${GRAPHFLOW}/graphflow.cpp
    ${SJ_TREE}/sj_tree.cpp

    ${PARALLEL_SYMBI}/parrallel_symbi.cpp
    ${PARALLEL_TURBOFLUX}/parallel_turboflux.cpp
    ${PARALLEL_GRAPHFLOW}/parallel_graphflow.cpp
    ${PARALLEL_CALIG}/parallel_calig_class.cpp
    ${PARALLEL_NEWSP}/newsp_adapter.cpp
    # ${PARALLEL_NEWSP}/parallel_newsp.cpp

    ${GRAPH}/graph.cpp
    ${GRAPH}/induced_graph.cpp

    ${UTILS}/globals.cpp

    core/inner_executor/FindMatchesKernel.cpp
    # core/FindMatchesKernel_simplified.cpp
    core/inner_executor/inner_executor.cpp
    core/inter_executor/inter_executor.cpp

    core/gpu/gpu_classifier.cu
    core/gpu/gpu_candidate_filter.cu
    core/gpu/gpu_search.cu
    core/gpu/gpu_bfs_search.cu
)

# -----------------------------
# 7. Generate executable
# -----------------------------
add_executable(csm ${SRC_FILES})

target_include_directories(csm PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)

target_compile_options(csm PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wno-unused-command-line-argument>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3>
)

# -----------------------------
# 8. Link libraries
# -----------------------------
target_link_libraries(csm PRIVATE 
    TBB::tbb
    OpenMP::OpenMP_CXX
    pthread
    CUDA::cudart
)

# -----------------------------
# 9. Standalone CaLiG executable
# -----------------------------
add_executable(parallel_calig
    ${MATCH}/Parallel_CaLiG/parallel_calig.cpp
)

target_link_libraries(parallel_calig PRIVATE 
    pthread
    OpenMP::OpenMP_CXX
)


add_executable(calig
    ${MATCH}/CaLiG/calig.cpp
)

target_link_libraries(calig PRIVATE 
    pthread
    OpenMP::OpenMP_CXX
)



### -----------------------------
### 10. Standalone Parallel_NewSP executable
### -----------------------------
add_executable(parallel_newsp
    ${PARALLEL_NEWSP}/matching/CSMPP.cpp
    ${PARALLEL_NEWSP}/matching/main.cpp
    ${PARALLEL_NEWSP}/matching/matching.cpp
    ${PARALLEL_NEWSP}/graph/graph.cpp
    ${PARALLEL_NEWSP}/graph/induced_graph.cpp
    ${PARALLEL_NEWSP}/DecisionMakingSystem/DecisionMakingSystem.cpp
    ${PARALLEL_NEWSP}/utils/globals.cpp
    
)

target_include_directories(parallel_newsp PRIVATE
    ${PARALLEL_NEWSP}
    ${PARALLEL_NEWSP}/matching
    ${PARALLEL_NEWSP}/graph
    ${PARALLEL_NEWSP}/DecisionMakingSystem
    ${PARALLEL_NEWSP}/configuration
    ${PARALLEL_NEWSP}/utils
)

target_link_libraries(parallel_newsp PRIVATE
    pthread
    OpenMP::OpenMP_CXX
)


### -----------------------------
### 11. Output directories (unify top-level outputs)
### -----------------------------
set(MAIN_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin)

set_target_properties(csm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${MAIN_OUTPUT_DIR}
)

set_target_properties(parallel_calig PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${MAIN_OUTPUT_DIR}
)

set_target_properties(calig PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${MAIN_OUTPUT_DIR}
)

set_target_properties(parallel_newsp PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${MAIN_OUTPUT_DIR}
)


### -----------------------------
### 12. (Removed: Parallel_NewSP ExternalProject — now built directly as parallel_newsp target above)
### -----------------------------
