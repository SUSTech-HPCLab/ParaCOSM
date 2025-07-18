#!/bin/bash

# Test Workflow for Graphflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/

export OMP_PROC_BIND=close
export OMP_PLACES=threads
# export OMP_DYNAMIC=true

export OMP_WAIT_POLICY=ACTIVE
# export OMP_DYNAMIC=false

# export TBB_NUM_THREADS=4

# Testing Algo
ALGORITHM=parallel_turboflux

# Output base directory
# OUTPUT_DIR=logs_txt/amazon/Graphflow

# Time limit (in seconds)
TIME_LIMIT=1800

# Graph files
DATA_GRAPH=/home/cc/haibin2/lsbench_x1/data_graph/data.graph
INSERT_GRAPH=/home/cc/haibin2/lsbench_x1/data_graph/insertion.graph

# List of suffixes to test
SUFFIXES=(8)

# Loop over each suffix
for SUFFIX in "${SUFFIXES[@]}"; do
    # Dataset
    DATA_SET=lsbench/${SUFFIX}



    # Output file
    OUTPUT=${ALGORITHM}_lsbench_x1_${SUFFIX}_P2.txt

    # Define configurations for sparse, dense, tree
    CONFIGS=(
        "sparse ${DIR}/lsbench_x1/random_walk/${SUFFIX}_self/sparse"
        "dense ${DIR}/lsbench_x1/random_walk/${SUFFIX}_self/dense"
        "tree ${DIR}/lsbench_x1/random_walk/${SUFFIX}_self/tree"
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
        for i in $(seq 1 100); do
            QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}

            echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

            # Run the test with appropriate parameters
            # if [ "$TYPE" = "tree" ]; then
            #     # Tree uses timeout and auto-tuning off
            #     timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
            #         --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
            #         >> ${TARGET_DIR}/${OUTPUT}
            # else
                # Sparse and dense use standard command
                timeout 3600 ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 >> ${OUTPUT}
                    
            # fi

            echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH"
            # echo "----------------------------------------"
            # echo " "
        done

        echo "Completed tests for ${TYPE}_${SUFFIX}"
        # echo "========================================="
    done
done

# echo "All tests completed for ${ALGORITHM} with suffixes ${SUFFIXES[*]}"

