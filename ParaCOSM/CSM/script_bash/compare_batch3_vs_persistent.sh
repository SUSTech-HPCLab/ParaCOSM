#!/bin/bash
# =============================================================================
# Compare batch3 (default, per-update fork/join) vs persistent (single OMP region)
# for parallel_graphflow on LiveJournal 7_self queries.
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"
QUERY_BASE="${HOME}/haibin/paracosm_data/lj/livejournal/random_walk/7_self"

ALGORITHM="parallel_graphflow"
THREADS=8
TIME_LIMIT=180
INITIAL_TIME_LIMIT=180
RUN_TIMEOUT=600
AUTO_TUNING=0
REPORT_INITIAL="off"

# Which query types to test
QUERY_TYPES=("sparse" "tree" "dense")
# Max queries per type (to keep run time manageable)
MAX_QUERIES=${MAX_QUERIES:-5}

OUTPUT_DIR="${1:-logs_txt/batch3_vs_persistent_lj_7self}"
mkdir -p "${OUTPUT_DIR}"

SUMMARY="${OUTPUT_DIR}/summary.csv"
echo "query,type,mode,initial,positive,negative,updates,time_ms,exit_code" > "${SUMMARY}"

run_one() {
    local query="$1" mode="$2" logfile="$3"
    timeout --kill-after=5 "${RUN_TIMEOUT}" \
        "${CSM_BIN}" \
        -a "${ALGORITHM}" \
        -d "${DATA_GRAPH}" \
        -u "${INSERT_GRAPH}" \
        -q "${query}" \
        --time-limit "${TIME_LIMIT}" \
        --initial-time-limit "${INITIAL_TIME_LIMIT}" \
        --report-initial "${REPORT_INITIAL}" \
        -t "${THREADS}" \
        --auto-tuning "${AUTO_TUNING}" \
        -m "${mode}" \
        > "${logfile}" 2>&1
    return $?
}

parse_log() {
    local logfile="$1"
    local initial positive negative updates time_ms
    initial=$(awk '/^[0-9]+ matches\.\s*$/ {print $1; exit}' "${logfile}")
    positive=$(awk '$2=="positive" && $3=="matches." {print $1; exit}' "${logfile}")
    negative=$(awk '$2=="negative" && $3=="matches." {print $1; exit}' "${logfile}")
    updates=$(awk '$2=="edge" && $3=="updates." {print $1; exit}' "${logfile}")
    time_ms=$(awk '/^Incremental Matching:/ {gsub(/ms/,"",$3); print $3; exit}' "${logfile}")
    echo "${initial:-?},${positive:-?},${negative:-?},${updates:-?},${time_ms:-?}"
}

echo ""
echo "========================================================"
echo " batch3 vs persistent   |  LJ 7_self  |  ${THREADS} threads"
echo "========================================================"

for QTYPE in "${QUERY_TYPES[@]}"; do
    QDIR="${QUERY_BASE}/${QTYPE}"
    [ -d "${QDIR}" ] || { echo "SKIP: ${QDIR} not found"; continue; }

    queries=( $(ls -d "${QDIR}"/Q_* 2>/dev/null | sort -t_ -k2 -n | head -n "${MAX_QUERIES}") )
    echo ""
    echo "--- type: ${QTYPE}  (${#queries[@]} queries) ---"
    printf "%-16s  %-8s  %10s  %10s  %10s  %8s  %10s  %10s  %s\n" \
        "query" "status" "pos_b3" "pos_per" "neg_b3" "neg_per" "ms_b3" "ms_per" "speedup"

    for q in "${queries[@]}"; do
        name=$(basename "$q")
        log_b3="${OUTPUT_DIR}/${name}.${QTYPE}.batch3.log"
        log_per="${OUTPUT_DIR}/${name}.${QTYPE}.persistent.log"

        # Run batch3
        ec_b3=0
        run_one "${q}" "batch3" "${log_b3}" || ec_b3=$?
        parsed_b3=$(parse_log "${log_b3}")

        # Run persistent
        ec_per=0
        run_one "${q}" "persistent" "${log_per}" || ec_per=$?
        parsed_per=$(parse_log "${log_per}")

        IFS=',' read -r init_b3 pos_b3 neg_b3 upd_b3 ms_b3 <<< "${parsed_b3}"
        IFS=',' read -r init_per pos_per neg_per upd_per ms_per <<< "${parsed_per}"

        # Correctness check
        status="PASS"
        if [[ "${pos_b3}" != "${pos_per}" || "${neg_b3}" != "${neg_per}" ]]; then
            status="MISMATCH"
        fi
        if [[ "${ec_b3}" -ne 0 || "${ec_per}" -ne 0 ]]; then
            status="ERROR(b3=${ec_b3},per=${ec_per})"
        fi

        # Speedup
        speedup="N/A"
        if [[ "${ms_b3}" != "?" && "${ms_per}" != "?" && "${ms_per}" != "0" ]]; then
            speedup=$(awk "BEGIN{printf \"%.2fx\", ${ms_b3}/${ms_per}}")
        fi

        printf "%-16s  %-8s  %10s  %10s  %10s  %8s  %10s  %10s  %s\n" \
            "${name}" "${status}" "${pos_b3}" "${pos_per}" "${neg_b3}" "${neg_per}" \
            "${ms_b3}" "${ms_per}" "${speedup}"

        echo "${name},${QTYPE},batch3,${init_b3},${pos_b3},${neg_b3},${upd_b3},${ms_b3},${ec_b3}" >> "${SUMMARY}"
        echo "${name},${QTYPE},persistent,${init_per},${pos_per},${neg_per},${upd_per},${ms_per},${ec_per}" >> "${SUMMARY}"
    done
done

echo ""
echo "Logs: ${OUTPUT_DIR}/"
echo "Summary CSV: ${SUMMARY}"
echo "Done."
