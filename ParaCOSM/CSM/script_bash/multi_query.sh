# Test Workflow
# This file was intended to run lsbench_x25, but it was too big for vtune.
# Data directory
DIR=/data/sicheng/lsbench
# Dataset
DATA_SET=high_freq
QUERY_TYPE=sparse_6

# Testing Algo
# sj-tree, graphflow, turboflux, symbi, iedyn, parrallel_symbi
ALGORITHM=parrallel_symbi 

# Time
TIME_LIMIT=3600

  
# Graph
DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph
QUERY_GRAPH_DIR=${DIR}/${DATA_SET}/query_graph/${QUERY_TYPE}
LOG_DIR=/data/sicheng/log

# write a test_dir array
test_dir=(/data/sicheng/high_freq/query_graph/sparse_6/Q_0)

# Get the file names from QUERY_GRAPH_DIR and set them to QUERY_GRAPH
for file in $QUERY_GRAPH_DIR/*
# for file in ${test_dir[@]}
do
    QUERY_GRAPH=$file
    QUERY_GRAPH_NAME=$(basename $QUERY_GRAPH)
    echo $QUERY_GRAPH_NAME

    LOG_PATH=${LOG_DIR}/${QUERY_TYPE}/high_freq_${QUERY_GRAPH_NAME}.log
    echo $LOG_PATH

    # Run
    # sudo touch ${LOG_PATH}
    # sudo chmod 777 ${LOG_PATH}
    # sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT} > ${LOG_PATH} 2>&1
done

