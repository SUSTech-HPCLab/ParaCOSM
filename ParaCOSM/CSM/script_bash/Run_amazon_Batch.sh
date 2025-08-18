#!/bin/bash

# Test Workflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/

export OMP_PROC_BIND=close
export OMP_PLACES=threads
export OMP_DYNAMIC=true

# Testing Algo
ALGORITHM=turboflux 

# Output base directory
OUTPUT_DIR=logs_txt/amazon/Turboflux

# Time limit (in seconds)
TIME_LIMIT=180

# List of suffixes to test
SUFFIXES=(6 7 8 9 10)

# Loop over each suffix
for SUFFIX in "${SUFFIXES[@]}"; do
    # Dataset
    DATA_SET=amazon/${SUFFIX}

    # Graph files
    DATA_GRAPH=${DIR}/${DATA_SET}/data_graph/data.graph
    INSERT_GRAPH=${DIR}/${DATA_SET}/data_graph/insertion.graph

    # Output file
    OUTPUT=${ALGORITHM}_Amazon_${SUFFIX}_A.txt

    # Define configurations for sparse, dense, tree
    CONFIGS=(
        "sparse ${DIR}/amazon/${SUFFIX}_self/sparse 6"
        "dense ${DIR}/amazon/${SUFFIX}_self/dense 16"
        "tree ${DIR}/amazon/${SUFFIX}_self/tree 32"
    )

    # Process each configuration
    for CONFIG in "${CONFIGS[@]}"; do
        # Extract type, query directory, and threads
        read TYPE QUERY_GRAPH_DIR THREADS <<< "$CONFIG"

        THREADS=32

        # Create output directory if it doesn't exist
        TARGET_DIR=${OUTPUT_DIR}/${TYPE}_${SUFFIX}
        if [ ! -d "${TARGET_DIR}" ]; then
            echo "Directory ${TARGET_DIR} does not exist. Creating it..."
            mkdir -p ${TARGET_DIR}
        fi

        echo "Starting tests for ${TYPE}_${SUFFIX} with dataset ${DATA_SET}"

        # Run tests for queries 1 to 99
        for i in $(seq 1 99); do
            QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}

            echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

            # Run the test with appropriate parameters
            if [ "$TYPE" = "tree" ]; then
                # Tree uses timeout and auto-tuning off
                timeout 600 sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
                    >> ${TARGET_DIR}/${OUTPUT}
            else
                # Sparse and dense use standard command
                sudo ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} \
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

echo "All tests completed for suffixes ${SUFFIXES[*]}"