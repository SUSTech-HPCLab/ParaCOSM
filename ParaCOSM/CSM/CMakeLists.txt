cmake_minimum_required(VERSION 3.16)

# -----------------------------
# 1. Set Intel oneAPI compiler
# -----------------------------
set(ONEAPI_ROOT $ENV{HOME}/intel/oneapi)
set(INTEL_COMPILER_ROOT ${ONEAPI_ROOT}/compiler/2025.3)
set(CMAKE_CXX_COMPILER ${INTEL_COMPILER_ROOT}/bin/icpx)
set(CMAKE_LINKER ${INTEL_COMPILER_ROOT}/bin/icpx)

# -----------------------------
# 2. Project declaration
# -----------------------------
project(CSM LANGUAGES CXX)

# -----------------------------
# 3. Compilation options
# -----------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Common optimizations and warnings
add_compile_options(
    -g
    -O3
    -Wall
    -Wno-unused-command-line-argument
    -qopenmp
)

# -----------------------------
# 4. Intel TBB paths
# -----------------------------
set(TBB_ROOT ${ONEAPI_ROOT}/tbb/latest)
set(TBB_INC ${TBB_ROOT}/include)

# Auto-detect TBB lib directory (different oneAPI versions use different layouts)
if(EXISTS "${TBB_ROOT}/lib/intel64/gcc4.8")
    set(TBB_LIB_DIR ${TBB_ROOT}/lib/intel64/gcc4.8)
else()
    set(TBB_LIB_DIR ${TBB_ROOT}/lib)
endif()

include_directories(
    ${TBB_INC}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
 
link_directories(${TBB_LIB_DIR})

# -----------------------------
# 5. OpenMP
# -----------------------------
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
endif()

# -----------------------------
# 6. Source files
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
    # ${PARALLEL_NEWSP}/parallel_newsp.cpp

    ${GRAPH}/graph.cpp
    ${GRAPH}/induced_graph.cpp

    ${UTILS}/globals.cpp

    core/inner_executor/FindMatchesKernel.cpp
    # core/FindMatchesKernel_simplified.cpp
    core/inner_executor/inner_executor.cpp
    core/inter_executor/inter_executor.cpp
)

# -----------------------------
# 7. Generate executable
# -----------------------------
add_executable(csm ${SRC_FILES})

# -----------------------------
# 8. Link libraries
# -----------------------------
target_link_libraries(csm PRIVATE 
    tbb
    pthread
    iomp5   # Intel OpenMP runtime
)

# Set rpath so runtime can find Intel OpenMP libraries
set(INTEL_OPENMP_LIB_DIR ${INTEL_COMPILER_ROOT}/lib)
target_link_options(csm PRIVATE "-Wl,-rpath,${INTEL_OPENMP_LIB_DIR}")

# -----------------------------
# 9. Standalone CaLiG executable
# -----------------------------
add_executable(parallel_calig
    ${MATCH}/Parallel_CaLiG/parallel_calig.cpp
)

target_link_libraries(parallel_calig PRIVATE 
    pthread
    iomp5
)

target_link_options(parallel_calig PRIVATE "-Wl,-rpath,${INTEL_OPENMP_LIB_DIR}")


add_executable(calig
    ${MATCH}/CaLiG/calig.cpp
)

target_link_libraries(calig PRIVATE 
    pthread
    iomp5
)

target_link_options(calig PRIVATE "-Wl,-rpath,${INTEL_OPENMP_LIB_DIR}")



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
    iomp5
)

target_link_options(parallel_newsp PRIVATE "-Wl,-rpath,${INTEL_OPENMP_LIB_DIR}")


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
