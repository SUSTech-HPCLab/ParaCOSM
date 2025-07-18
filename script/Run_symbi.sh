# Test Workflow

# Data directory
DIR=/home/cc/haibin2/
# Dataset
DATA_SET=livejournal/30

# export KMP_AFFINITY=compact,1,0,granularity=fine
# Testing Algo
# sj-tree, graphflow, turboflux, symbi, iedyn, parrallel_symbi
ALGORITHM=parrallel_symbi 

# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph
QUERY_GRAPH=${DIR}/${DATA_SET}/query_graph/sparse_6/Q_9
QUERY_GRAPH_DIR=${DIR}/${DATA_SET}/query_graph/sparse_6
QUERY_GRAPH_big=${DIR}/${DATA_SET}/query_graph/sparse_6/Q_big

FINAL=${QUERY_GRAPH}

# low freg 0 27
# LOWFREEDATA=/data/sicheng/lsbench/low_freq/data_graph/data.graph
# LOWFREEINC=/data/sicheng/lsbench/low_freq/data_graph/insertion.graph
# LOWQUERY=/data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32

# Time
TIME_LIMIT=180

# echo "Start Testing ${DATA_GRAPH} on ${INSERT_GRAPH} with time limit ${FINAL}"

# # Run
# sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${FINAL} --time-limit ${TIME_LIMIT} 

# sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${FINAL} --time-limit ${TIME_LIMIT} 


for i in $(seq 1 100); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}
    
    # 检查 QUERY_GRAPH 是否存在
    if [ ! -f "$QUERY_GRAPH" ]; then
        echo "Query graph $QUERY_GRAPH does not exist. Skipping..."
        continue
    fi

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} >> logs_txt/symbi_good.txt
    # sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

    # 可选：记录日志
    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done



# valgrind --tool=memcheck --leak-check=full ./build/csm -a ${ALGORITHM} -d ${LOWFREEDATA} -u ${LOWFREEINC} -q ${LOWQUERY} --time-limit ${TIME_LIMIT}

# sudo gdb ./build/csm 
# run -a parrallel_symbi -d /data/sicheng/lsbench/low_freq/data_graph/data.graph -u /data/sicheng/lsbench/low_freq/data_graph/insertion.graph -q /data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32 --time-limit 3600
# sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

# 数据那么大，发生大更新的问题是什么？PageFault?
# ./build/csm -a parrallel_symbi -d /home/cc/haibin/CSM-Benchmark/ContinuousSubgraphMatching/livejournal/30/data_graph/data.graph -u /home/cc/haibin/CSM-Benchmark/ContinuousSubgraphMatching/livejournal/30/data_graph/insertion.graph -q /home/cc/haibin/CSM-Benchmark/ContinuousSubgraphMatching/livejournal/30/query_graph/sparse_6/Q_big9 --time-limit 7200


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


