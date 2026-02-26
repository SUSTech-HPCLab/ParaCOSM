#!/bin/bash
# =============================================================================
# Graphflow parallel test script (clean version)
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Environment and library path
# -----------------------------------------------------------------------------
export LD_LIBRARY_PATH=/home/cc/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH}

# -----------------------------------------------------------------------------
# Global paths and dataset
# -----------------------------------------------------------------------------
BASE_DIR=/your/path/to/ParaCOSM
# Optional dataset: livejournal, amazon, orkut, lsbench
DATA_SET_NAME=livejournal
DATA_GRAPH=${BASE_DIR}/paracosm/${DATA_SET_NAME}/data_graph/data.graph
INSERT_GRAPH=${BASE_DIR}/paracosm/${DATA_SET_NAME}/data_graph/insertion.graph

# -----------------------------------------------------------------------------
# Algorithm and run parameters
# -----------------------------------------------------------------------------
# Algorithm: parallel_symbi, parallel_turboflux, parallel_graphflow
ALGORITHM=parallel_symbi
TIME_LIMIT=1800
THREADS=32
RUN_TIMEOUT=3600

# Optional: output directory (leave empty to print to terminal only)
OUTPUT_DIR=""
# OUTPUT_DIR=logs_txt/livejournal/Graphflow

# -----------------------------------------------------------------------------
# Configurations to test: (type, query graph directory)
# -----------------------------------------------------------------------------
SUFFIXES=(6)

get_configs() {
    local suffix=$1
    # Uncomment or add more configs as needed
    echo "sparse ${BASE_DIR}/livejournal/6/query_graph/sparse_6"
    # echo "sparse ${BASE_DIR}/livejournal/random_walk/6_self/sparse"
    # echo "dense  ${BASE_DIR}/livejournal/6/random_walk/${suffix}_self/dense"
    # echo "tree   ${BASE_DIR}/livejournal/6/random_walk/${suffix}_self/tree"
}

# Query IDs to run (currently only Q_3)
QUERY_IDS=(3)
# For a range use: QUERY_IDS=($(seq 3 5))

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
for SUFFIX in "${SUFFIXES[@]}"; do
    DATA_SET=${DATA_SET_NAME}/${SUFFIX}
    OUTPUT_FILE="${ALGORITHM}_${DATA_SET_NAME}_${SUFFIX}.txt"

    while IFS= read -r CONFIG_LINE; do
        [ -z "$CONFIG_LINE" ] && continue
        read -r TYPE QUERY_GRAPH_DIR <<< "$CONFIG_LINE"

        TARGET_DIR=""
        if [ -n "${OUTPUT_DIR}" ]; then
            TARGET_DIR="${OUTPUT_DIR}/${TYPE}_${SUFFIX}"
            if [ ! -d "${TARGET_DIR}" ]; then
                echo "Creating directory: ${TARGET_DIR}"
                mkdir -p "${TARGET_DIR}"
            fi
        fi

        echo "========== ${TYPE}_${SUFFIX} | dataset: ${DATA_SET} =========="

        for i in "${QUERY_IDS[@]}"; do
            QUERY_GRAPH="${QUERY_GRAPH_DIR}/Q_${i}"
            echo "Running QUERY_GRAPH: ${QUERY_GRAPH}"

            # timeout: run up to RUN_TIMEOUT seconds; after SIGTERM wait 5s then SIGKILL
            CMD=(
                timeout --kill-after=5 "${RUN_TIMEOUT}"
                ./build/bin/csm
                -a "${ALGORITHM}"
                -d "${DATA_GRAPH}"
                -u "${INSERT_GRAPH}"
                -q "${QUERY_GRAPH}"
                --time-limit "${TIME_LIMIT}"
                --report-initial off
                -t "${THREADS}"
                --auto-tuning 0
            )

            if [ -n "${TARGET_DIR}" ]; then
                "${CMD[@]}" >> "${TARGET_DIR}/${OUTPUT_FILE}"
            else
                "${CMD[@]}"
            fi

            echo "Done: ${QUERY_GRAPH}"
        done

        echo "Completed config: ${TYPE}_${SUFFIX}"
    done < <(get_configs "$SUFFIX")
done

echo "All tests completed: ${ALGORITHM}, suffixes ${SUFFIXES[*]}"
