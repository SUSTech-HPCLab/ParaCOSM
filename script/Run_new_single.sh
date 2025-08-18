# Test Workflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/
# Dataset
DATA_SET=amazon/6

export OMP_PROC_BIND=close
export OMP_PLACES=threads
export OMP_DYNAMIC=true

# export KMP_AFFINITY=compact,1,0,granularity=fine
# Testing Algo
# sj-tree, graphflow, turboflux, symbi, iedyn, parallel_symbi, parallel_turboflux, parallel_graphflow
ALGORITHM=graphflow 

# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph
QUERY_GRAPH=${DIR}/${DATA_SET}/query_graph/sparse_6/Q_3
# QUERY_GRAPH_DIR=${DIR}/${DATA_SET}/query_graph/sparse_6

# change 1
QUERY_GRAPH_DIR_ROOT=/home/cc/haibin2/amazon/8_self
FINAL=${QUERY_GRAPH}

OUTPUT_DIR=logs_txt/amazon/Graphflow
OUTPUT=${ALGORITHM}_Amazon_1_self2.txt


# Time
TIME_LIMIT=180

# echo "Start Testing ${DATA_GRAPH} on ${INSERT_GRAPH} with time limit ${FINAL}"

# # Run Test
sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} \
        -q ${FINAL} --time-limit ${TIME_LIMIT} --report-initial off -t 6 --auto-tuning 0

# sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${FINAL} --time-limit ${TIME_LIMIT} 


# 创建所有必要的输出目录
for dir in "sparse_8" "dense_8" "tree_8"; do
    if [ ! -d "${OUTPUT_DIR}/${dir}" ]; then
        echo "Directory ${OUTPUT_DIR}/${dir} does not exist. Creating it..."
        mkdir -p ${OUTPUT_DIR}/${dir}
    fi
done


QUERY_GRAPH_DIR=${QUERY_GRAPH_DIR_ROOT}/sparse

for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}
    
    # 检查 QUERY_GRAPH 是否存在
    if [ ! -f "$QUERY_GRAPH" ]; then
        echo "Query graph $QUERY_GRAPH does not exist. Skipping..."
        continue
    fi

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} --report-initial off -t 6 \
            >> ${OUTPUT_DIR}/sparse_8/${OUTPUT}
    # sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

    # 可选：记录日志
    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done

QUERY_GRAPH_DIR2=${QUERY_GRAPH_DIR_ROOT}/dense


for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR2}/Q_${i}
    
    # 检查 QUERY_GRAPH 是否存在
    if [ ! -f "$QUERY_GRAPH" ]; then
        echo "Query graph $QUERY_GRAPH does not exist. Skipping..."
        continue
    fi

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} --report-initial off -t 16 \
            >> ${OUTPUT_DIR}/dense_8/${OUTPUT}
    # sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

    # 可选：记录日志
    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done

QUERY_GRAPH_DIR3=${QUERY_GRAPH_DIR_ROOT}/tree
# QUERY_GRAPH_DIR3=${DIR}/${DATA_SET}/query_graph/tree_6


for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR3}/Q_${i}
    
    # 检查 QUERY_GRAPH 是否存在
    if [ ! -f "$QUERY_GRAPH" ]; then
        echo "Query graph $QUERY_GRAPH does not exist. Skipping..."
        continue
    fi

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} \
    -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} --report-initial off -t 32 \
            >> ${OUTPUT_DIR}/tree_8/${OUTPUT}
    # sudo ./build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

    # 可选：记录日志
    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done
