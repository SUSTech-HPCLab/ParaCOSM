(base) haibin@hpclab03:~/CSM-Benchmark/ContinuousSubgraphMatching$ cat run_parrallel.sh 
# Test Workflow

# Data directory
DIR=/data/haibin/CSM_dataset
# Dataset
DATA_SET=livejournal/30


# Testing Algo
# sj-tree, graphflow, turboflux, symbi, iedyn, parrallel_symbi
# ALGORITHM=symbi 
ALGORITHM=parallel_symbi 
# symbi 40478.3ms
# 6551.23ms
# parallel: 12-16s
# 在深度很浅的情况下，搜索的非常慢
# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph
QUERY_GRAPH=${DIR}/${DATA_SET}/query_graph/sparse_6/Q_32

# low freg 0 27
LOWFREEDATA=/data/sicheng/lsbench/low_freq/data_graph/data.graph
LOWFREEINC=/data/sicheng/lsbench/low_freq/data_graph/insertion.graph
LOWQUERY=/data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32

Hey=/home/haibin/Q_10
# Q32
# Time
TIME_LIMIT=3600

echo "Start Testing ${ALGORITHM} ${LOWFREEDATA} on ${LOWFREEINC} with time limit ${QUERY_GRAPH}"

# Run
sudo  ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} 

# sudo  ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} 
# strace -e trace=futex
# strace -c
# valgrind --tool=memcheck --leak-check=full ./build/csm -a ${ALGORITHM} -d ${LOWFREEDATA} -u ${LOWFREEINC} -q ${LOWQUERY} --time-limit ${TIME_LIMIT}

# sudo sysctl -w kernel.yama.ptrace_scope=0
# vtune -collect hotspots ./build/csm -a parallel_symbi -d /data/sicheng/lsbench/low_freq/data_graph/data.graph -u /data/sicheng/lsbench/low_freq/data_graph/insertion.graph -q /data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32 --time-limit 3600


# sudo gdb ./build/csm 
# run -a parrallel_symbi -d /data/sicheng/lsbench/low_freq/data_graph/data.graph -u /data/sicheng/lsbench/low_freq/data_graph/insertion.graph -q /data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32 --time-limit 3600
# sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

# 数据那么大，发生大更新的问题是什么？PageFault?

# sudo ../CaLiG/CaLiG/calig -d ${DATA_GRAPH} -q ${QUERY_GRAPH} -s ${INSERT_GRAPH} 
# ./calig -d ./github10/initial -q github10/Q/10/q1 -s ./github10/ss

# sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}


# gdb debug:
# # 编译代码时启用调试信息
# g++ -g -o csm your_source_files.cpp

# # 使用gdb运行程序
# gdb ./build/csm

# # 在gdb中运行程序
# (gdb) run -a ${ALGORITHM} -d ${LOWFREEDATA} -u ${LOWFREEINC} -q ${LOWQUERY} --time-limit ${TIME_LIMIT}

# # 当程序遇到段错误时，使用backtrace命令查看调用堆栈
# (gdb) backtrace
