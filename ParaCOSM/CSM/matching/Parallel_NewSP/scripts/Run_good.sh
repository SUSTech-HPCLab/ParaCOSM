# Test Workflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/
# Dataset
DATA_SET=amazon/6

# Time
TIME_LIMIT=3600

# echo "Start Testing ${DATA_GRAPH} on ${INSERT_GRAPH} with time limit ${FINAL}"

# Run for each folder under 6_self
for folder in /home/cc/haibin2/amazon/10_self/{sparse,dense,tree}; do
    echo "Running CSM for folder: ${folder}"
    sudo ./build/csm -d /home/cc/haibin2/amazon/6/data_graph/data.graph \
        -u /home/cc/haibin2/amazon/6/data_graph/insertion.graph \
        -q "${folder}" --time-limit ${TIME_LIMIT} --report-initial off
done

# /home/cc/haibin2/amazon/9_self