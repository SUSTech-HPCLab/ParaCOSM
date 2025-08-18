# Parallel Algorithm Running 


## Graphflow, Symbi, Turboflux

### Compile

ensure to load intel oneapi

then

```bash
make
```

### Dataset 

Your dataset file should be:

    CONFIGS=(
        "sparse ${DIR}/${DATASET}/${SUFFIX}_self/sparse"
        "dense ${DIR}/${DATASET}/${SUFFIX}_self/dense"
        "tree ${DIR}/${DATASET}/${SUFFIX}_self/tree"
    )


其中:
- `sparse`、`dense`和`tree`表示不同类型的图结构
- `${DIR}`是数据根目录的路径
- `${DATASET}`是数据集名称（如示例中的"amazon"）
- `${SUFFIX}`表示查询图的大小（如示例中的6-10）
- 目录结构应为`${DIR}/${DATASET}/${SUFFIX}_self/[sparse|dense|tree]`


运行脚本：


```bash

#!/bin/bash

# Test Workflow for Graphflow

export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:$LD_LIBRARY_PATH

# Data directory
DIR=/home/cc/haibin2/

DATASET=amazon

# Set threads
THREADS=32

export OMP_PROC_BIND=close
export OMP_PLACES=threads
export OMP_DYNAMIC=true

# Testing Algo
ALGORITHM=parallel_graphflow

# Output base directory
OUTPUT_DIR=logs_txt/amazon/Parallel_graphflow

# Graph files
DATA_GRAPH=/home/cc/haibin2/amazon/6/data_graph/data.graph
INSERT_GRAPH=/home/cc/haibin2/amazon/6/data_graph/insertion.graph

# List of suffixes to test. 6-10 means test query graph size from 6-10
SUFFIXES=(6 7 8 9 10)

# Loop over each suffix
for SUFFIX in "${SUFFIXES[@]}"; do
    # Dataset
    DATA_SET=${DATASET}/${SUFFIX}

    # Output file
    OUTPUT=${ALGORITHM}_Amazon_${SUFFIX}_P2.txt

    # Define configurations for sparse, dense, tree
    CONFIGS=(
        "sparse ${DIR}/${DATASET}/${SUFFIX}_self/sparse"
        "dense ${DIR}/${DATASET}/${SUFFIX}_self/dense"
        "tree ${DIR}/${DATASET}/${SUFFIX}_self/tree"
    )

    # Process each configuration
    for CONFIG in "${CONFIGS[@]}"; do
        # Extract type and query directory
        read TYPE QUERY_GRAPH_DIR <<< "$CONFIG"

        # Create output directory if it doesn't exist
        TARGET_DIR=${OUTPUT_DIR}/${TYPE}_${SUFFIX}
        if [ ! -d "${TARGET_DIR}" ]; then
            echo "Directory ${TARGET_DIR} does not exist. Creating it..."
            mkdir -p ${TARGET_DIR}
        fi

        # echo "Starting tests for ${TYPE}_${SUFFIX} with dataset ${DATA_SET}"

        # Run tests for queries 1 to 99
        for i in $(seq 1 99); do
            QUERY_GRAPH=${QUERY_GRAPH_DIR}/Q_${i}

            echo "Testing with QUERY_GRAPH: $QUERY_GRAPH"

            # Run the test with appropriate parameters
            if [ "$TYPE" = "tree" ]; then
                # Tree uses timeout and auto-tuning off
                timeout 3600 ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
                    >> ${TARGET_DIR}/${OUTPUT}
            else
                # Sparse and dense use standard command
                timeout 3600 ./build/csm -a ${ALGORITHM} -d ${DATA_GRAPH} -u ${INSERT_GRAPH} -q ${QUERY_GRAPH} \
                    --time-limit ${TIME_LIMIT} --report-initial off -t ${THREADS} --auto-tuning 0 \
                    >> ${TARGET_DIR}/${OUTPUT}
            fi

            echo "Finished testing QUERY_GRAPH: $QUERY_GRAPH"
            # echo "----------------------------------------"
            # echo " "
        done

        echo "Completed tests for ${TYPE}_${SUFFIX}"
        echo "========================================="
    done
done

echo "All tests completed for ${ALGORITHM} with suffixes ${SUFFIXES[*]}"


```


## CaLiG

Ensure that you use CaLiG dataset!

### Compile

1. go into ContinuousSubgraphMatching/matching/CaLiG file

2. see Run.sh, you can compile with intel ICPX


WIP


## NewSP

### Compile 

1. go into newSP_GOOD file

2. `make` with intel icpx

### Run

see newSP_GOOD/run_parallel.sh

WIP
