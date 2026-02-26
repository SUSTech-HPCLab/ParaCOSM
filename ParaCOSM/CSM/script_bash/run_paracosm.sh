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
BASE_DIR=
# Optional dataset: livejournal, amazon, orkut, lsbench
DATA_SET_NAME=amazon_dataset

# -----------------------------------------------------------------------------
# Algorithm and run parameters
# -----------------------------------------------------------------------------
# Algorithm: parallel_symbi, parallel_turboflux, parallel_graphflow
ALGORITHM=parallel_graphflow
TIME_LIMIT=1800
THREADS=8
RUN_TIMEOUT=3600

# Optional: output directory (leave empty to print to terminal only)
OUTPUT_DIR="logs_txt/amazon/parallel_graphflow"
# OUTPUT_DIR=logs_txt/livejournal/Graphflow

# -----------------------------------------------------------------------------
# Configurations to test: (type, query graph directory)
# -----------------------------------------------------------------------------
SUFFIXES=(6)

get_configs() {
    local suffix=$1
    # Uncomment or add more configs as needed
    echo "sparse ${BASE_DIR}/paracosm/${DATA_SET_NAME}/${suffix}_self/sparse"
    # echo "dense  ${BASE_DIR}/paracosm/${DATA_SET_NAME}/${suffix}_self/dense"
    # echo "tree   ${BASE_DIR}/paracosm/${DATA_SET_NAME}/${suffix}_self/tree"
    # echo "sparse ${BASE_DIR}/livejournal/random_walk/6_self/sparse"
    # echo "dense  ${BASE_DIR}/livejournal/6/random_walk/${suffix}_self/dense"
    # echo "tree   ${BASE_DIR}/livejournal/6/random_walk/${suffix}_self/tree"
}

# Query IDs to run (currently only Q_3)
QUERY_IDS=($(seq 3 20))
# For a range use: QUERY_IDS=($(seq 3 5))

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
TOTAL_RUNS=0
FAILED_RUNS=0
SKIPPED_RUNS=0

for SUFFIX in "${SUFFIXES[@]}"; do
    DATA_GRAPH=${BASE_DIR}/paracosm/${DATA_SET_NAME}/${SUFFIX}/data_graph/data.graph
    INSERT_GRAPH=${BASE_DIR}/paracosm/${DATA_SET_NAME}/${SUFFIX}/data_graph/insertion.graph
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
            if [ ! -e "${QUERY_GRAPH}" ]; then
                TOTAL_RUNS=$((TOTAL_RUNS + 1))
                SKIPPED_RUNS=$((SKIPPED_RUNS + 1))
                continue
            fi

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

            TOTAL_RUNS=$((TOTAL_RUNS + 1))
            if [ -n "${TARGET_DIR}" ]; then
                if "${CMD[@]}" >> "${TARGET_DIR}/${OUTPUT_FILE}" 2>&1; then
                    echo "Done: ${QUERY_GRAPH}"
                else
                    EXIT_CODE=$?
                    FAILED_RUNS=$((FAILED_RUNS + 1))
                    echo "[ERROR] Failed (${EXIT_CODE}): ${QUERY_GRAPH}" | tee -a "${TARGET_DIR}/${OUTPUT_FILE}"
                fi
            else
                if "${CMD[@]}"; then
                    echo "Done: ${QUERY_GRAPH}"
                else
                    EXIT_CODE=$?
                    FAILED_RUNS=$((FAILED_RUNS + 1))
                    echo "[ERROR] Failed (${EXIT_CODE}): ${QUERY_GRAPH}"
                fi
            fi
        done

        echo "Completed config: ${TYPE}_${SUFFIX}"
    done < <(get_configs "$SUFFIX")
done

echo "All tests completed: ${ALGORITHM}, suffixes ${SUFFIXES[*]}"
echo "Summary: total=${TOTAL_RUNS}, failed=${FAILED_RUNS}, skipped=${SKIPPED_RUNS}, success=$((TOTAL_RUNS - FAILED_RUNS - SKIPPED_RUNS))"
