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

# Time
TIME_LIMIT=3600

# Run
# vtune -collect hotspots $(pwd)/build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}

vtune -collect hotspots $(pwd)/build/csm -a symbi -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} --time-limit ${TIME_LIMIT}