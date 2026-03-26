cmake_minimum_required(VERSION 3.15)
project(CSM LANGUAGES CXX)

# ========= Compiler 设置 =========
set(CMAKE_CXX_COMPILER icpx)

# 默认编译选项
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Debug/Release 配置
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG   "-g")

# ========= OneAPI & TBB 路径 =========
if(DEFINED ENV{ONEAPI_ROOT})
    set(ONEAPI_ROOT $ENV{ONEAPI_ROOT})
elseif(DEFINED ONEAPI_ROOT)
    # keep as-is (passed via -D)
else()
    set(ONEAPI_ROOT "$ENV{HOME}/intel/oneapi")
endif()

set(TBB_ROOT "${ONEAPI_ROOT}/tbb/latest")
set(TBB_INC "${TBB_ROOT}/include")

if(EXISTS "${TBB_ROOT}/lib/intel64/gcc4.8")
    set(TBB_LIB_DIR "${TBB_ROOT}/lib/intel64/gcc4.8")
else()
    set(TBB_LIB_DIR "${TBB_ROOT}/lib")
endif()

include_directories(${TBB_INC})
link_directories(${TBB_LIB_DIR})

# ========= 源码路径 =========
set(GRAPH_DIR graph)
set(MATCH_DIR matching)
set(UTILS_DIR utils)
set(DMS_DIR DecisionMakingSystem)

include_directories(
    ${CMAKE_SOURCE_DIR}      # <== 加这个，等价于 Makefile 的 -I.
    ${GRAPH_DIR}
    ${MATCH_DIR}
    ${UTILS_DIR}
    ${DMS_DIR}
    configuration
)

# ========= 源文件 =========
set(SOURCES
    ${MATCH_DIR}/main.cpp
    ${MATCH_DIR}/matching.cpp
    ${MATCH_DIR}/CSMPP.cpp
    ${GRAPH_DIR}/graph.cpp
    ${GRAPH_DIR}/induced_graph.cpp
    ${UTILS_DIR}/globals.cpp
    ${DMS_DIR}/DecisionMakingSystem.cpp
)

# ========= 可执行文件 =========
add_executable(csm ${SOURCES})

# ========= 链接库 =========
target_link_libraries(csm
    -qopenmp                 # Intel OpenMP
    -qopenmp-link=static
    tbb
    pthread
)

# ========= 运行时路径 =========
set_target_properties(csm PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build
    LINK_FLAGS "-Wl,-rpath,${TBB_LIB_DIR}"
)
