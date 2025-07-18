# Test Workflow

# Data directory
DIR=/data/haibin/CSM_dataset
# Dataset
DATA_SET=livejournal/30


# Testing Algo
# sj-tree, graphflow, turboflux, symbi, iedyn, parrallel_symbi
ALGORITHM=parrallel_symbi 

# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph
QUERY_GRAPH=${DIR}/${DATA_SET}/query_graph/sparse_6/Q_6

# low freg 0 27
LOWFREEDATA=/data/sicheng/lsbench/low_freq/data_graph/data.graph
LOWFREEINC=/data/sicheng/lsbench/low_freq/data_graph/insertion.graph
LOWQUERY=/data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32

# Time
TIME_LIMIT=3600

echo "Start Testing ${LOWFREEDATA} on ${LOWFREEINC} with time limit ${LOWQUERY}"


export ALGORITHM=parrallel_symbi
export LOWFREEDATA=/data/sicheng/lsbench/low_freq/data_graph/insertion.graph
export LOWFREEINC=/data/sicheng/lsbench/low_freq/data_graph/insertion.graph
export LOWQUERY=/data/sicheng/lsbench/low_freq/query_graph/sparse_6/Q_32
export TIME_LIMIT="120"


# Run
sudo gdb ./build/csm 

# # gdb ./build/csm
# (gdb) set args -a ${ALGORITHM} -d ${LOWFREEDATA} -u ${LOWFREEINC} -q ${LOWQUERY} --time-limit ${TIME_LIMIT}
# (gdb) run


# sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

# 数据那么大，发生大更新的问题是什么？PageFault?

# sudo ../CaLiG/CaLiG/calig -d ${DATA_GRAPH} -q ${QUERY_GRAPH} -s ${INSERT_GRAPH} 
# ./calig -d ./github10/initial -q github10/Q/10/q1 -s ./github10/ss

# sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}