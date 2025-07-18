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
ALGORITHM=parallel_graphflow 

# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph

OUTPUT_DIR=logs_txt/amazon/Parallel_graphflow
OUTPUT=${ALGORITHM}_Amazon_9.txt

# Time
TIME_LIMIT=180

# 创建所有必要的输出目录
for dir in "sparse_9" "dense_9" "tree_9"; do
    if [ ! -d "${OUTPUT_DIR}/${dir}" ]; then
        echo "Directory ${OUTPUT_DIR}/${dir} does not exist. Creating it..."
        mkdir -p ${OUTPUT_DIR}/${dir}
    fi
done


QUERY_GRAPH_DIR=/home/cc/haibin2/amazon/9_self/sparse

for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} --report-initial off -t 6 \
            >> ${OUTPUT_DIR}/sparse_9/${OUTPUT}

    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done

QUERY_GRAPH_DIR2=/home/cc/haibin2/amazon/9_self/dense


for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR2}/Q_${i}

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} --report-initial off -t 16 \
            >> ${OUTPUT_DIR}/dense_9/${OUTPUT}

    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done


QUERY_GRAPH_DIR3=/home/cc/haibin2/amazon/9_self/tree


for i in $(seq 1 99); do
    QUERY_GRAPH=${QUERY_GRAPH_DIR3}/Q_${i}

    echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

    # 执行测试
    timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
    --time-limit ${TIME_LIMIT} --report-initial off -t 32 --auto-tuning 0 \
            >> ${OUTPUT_DIR}/tree_9/${OUTPUT}

    echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH" 
    echo "----------------------------------------"
    echo " "
done
