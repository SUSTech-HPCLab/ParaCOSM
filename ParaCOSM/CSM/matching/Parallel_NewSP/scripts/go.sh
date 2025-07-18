DIR=/home/cc/haibin2/
# Dataset
DATA_SET=amazon/6

# Time
TIME_LIMIT=180

# Run for each folder under 6_self
for folder in /home/cc/haibin2/amazon/10_self/{sparse,dense,tree}; do
    echo "Processing folder: ${folder}"
    # Run for each graph file Q_1 to Q_100
    for i in {1..100}; do
        graph_file="${folder}/Q_${i}"
        echo "Running CSM for graph: ${graph_file}"
        timeout 900 ./build/csm -d /home/cc/haibin2/amazon/6/data_graph/data.graph \
            -u /home/cc/haibin2/amazon/6/data_graph/insertion.graph \
            -q "${graph_file}" --time-limit ${TIME_LIMIT} --report-initial off
    done
done