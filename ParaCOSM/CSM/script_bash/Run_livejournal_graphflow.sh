#!/bin/bash

# Test Workflow for Symbi

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/

# export OMP_PROC_BIND=close
# export OMP_PLACES=threads
# export OMP_DYNAMIC=true

# Testing Algo
ALGORITHM=graphflow

# Output base directory
OUTPUT_DIR=logs_txt/orkut/graphflow

# Time limit (in seconds)
TIME_LIMIT=1800

# List of suffixes to test
SUFFIXES=(8)

# Graph files
DATA_GRAPH=/home/cc/haibin2/orkut/data_graph/data_.graph
INSERT_GRAPH=/home/cc/haibin2/orkut/data_graph/insert.graph

# Loop over each suffix
for SUFFIX in "${SUFFIXES[@]}"; do
    # Dataset
    DATA_SET=Orkut/${SUFFIX}

    # Output file
    OUTPUT=${ALGORITHM}_Orkut_${SUFFIX}_P2.txt

    # Define configurations for sparse, dense, tree
    CONFIGS=(
        "dense /home/cc/haibin2/livejournal/random_walk/${SUFFIX}_self/dense"
        "sparse /home/cc/haibin2/livejournal/random_walk/${SUFFIX}_self/sparse"
        "tree ${DIR}/amazon/${SUFFIX}_self/tree"
    )

    # Process each configuration
    for CONFIG in "${CONFIGS[@]}"; do
        # Extract type and query directory
        read TYPE QUERY_GRAPH_DIR <<< "$CONFIG"

        # Set threads
        THREADS=32

        # Create output directory if it doesn't exist
        TARGET_DIR=${OUTPUT_DIR}/${TYPE}_${SUFFIX}
        if [ ! -d "${TARGET_DIR}" ]; then
            echo "Directory ${TARGET_DIR} does not exist. Creating it..."
            mkdir -p ${TARGET_DIR}
        fi

        echo "Starting tests for ${TYPE}_${SUFFIX} with dataset ${DATA_SET}"

        # Run tests for queries 1 to 99
        for i in $(seq 1 50); do
            QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}

            echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

            # Run the test with appropriate parameters
            if [ "$TYPE" = "tree" ]; then
                # Tree uses timeout and auto-tuning off
                timeout 1600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
                    >> ${TARGET_DIR}/${OUTPUT}
            else
                # Sparse and dense use standard command
                timeout 1600 ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
                    >> ${TARGET_DIR}/${OUTPUT}
            fi

            echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH"
            echo "----------------------------------------"
            echo " "
        done

        echo "Completed tests for ${TYPE}_${SUFFIX}"
        echo "========================================="
    done
done

echo "All tests completed for ${ALGORITHM} with suffixes ${SUFFIXES[*]}"