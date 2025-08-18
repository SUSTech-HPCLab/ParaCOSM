# Test Workflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/
# Dataset
DATA_SET=amazon/6
export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH
# Time
TIME_LIMIT=3600

# echo "Start Testing ${DATA_GRAPH} on ${INSERT_GRAPH} with time limit ${FINAL}"

# # Run
./build/csm  -d /home/cc/haibin2/livejournal/data_graph/data.graph -u /home/cc/haibin2/livejournal/data_graph/insertion.graph \
        -q /home/cc/haibin2/livejournal/random_walk/6_self/dense/Q_16 --time-limit ${TIME_LIMIT} --report-initial off

# ./build/csm  -d /home/cc/haibin2/livejournal/30/data_graph/data.graph -u /home/cc/haibin2/livejournal/30/data_graph/insertion.graph \
#         -q /home/cc/haibin2/livejournal/9_self/dense/Q_3 --time-limit ${TIME_LIMIT} --report-initial off
