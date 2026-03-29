#!/bin/bash
# =============================================================================
# Screen 200 generated 8-vertex queries on LJ: serial graphflow vs parallel
# Uses -m single (no classify). Short time limit for fast triage.
# =============================================================================
set -euo pipefail

CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"

QUERY_DIR="${1:?Usage: $0 <query_dir> [output_dir]}"
OUTPUT_DIR="${2:-logs_txt/screen_8v_200}"
mkdir -p "${OUTPUT_DIR}"

THREADS=8
TIME_LIMIT=${TIME_LIMIT:-120}
RUN_TIMEOUT=${RUN_TIMEOUT:-300}

SUMMARY="${OUTPUT_DIR}/summary.csv"
# Only write header if file doesn't exist or is empty
if [ ! -s "${SUMMARY}" ]; then
    echo "query,edges,serial_ms,parallel_ms,speedup,pos_serial,pos_parallel,neg_serial,neg_parallel,status" > "${SUMMARY}"
fi

run_one() {
    local algo="$1" threads="$2" query="$3" logfile="$4"
    timeout --kill-after=5 "${RUN_TIMEOUT}" \
        "${CSM_BIN}" \
        -a "${algo}" \
        -d "${DATA_GRAPH}" \
        -u "${INSERT_GRAPH}" \
        -q "${query}" \
        --time-limit "${TIME_LIMIT}" \
        --report-initial off \
        -t "${threads}" \
        --auto-tuning 0 \
        -m single \
        > "${logfile}" 2>&1
    return $?
}

parse_time() { awk '/^Incremental Matching:/ {gsub(/ms/,"",$3); print $3; exit}' "$1"; }
parse_pos()  { awk '$2=="positive" && $3=="matches." {print $1; exit}' "$1"; }
parse_neg()  { awk '$2=="negative" && $3=="matches." {print $1; exit}' "$1"; }

QUERIES=( $(ls "${QUERY_DIR}"/Q_gen_* 2>/dev/null | sort) )
TOTAL=${#QUERIES[@]}

echo ""
echo "================================================================"
echo " Screen 8v queries: serial graphflow(1T) vs parallel(${THREADS}T)"
echo " Queries: ${TOTAL}  TIME_LIMIT=${TIME_LIMIT}s"
echo "================================================================"
printf "%-14s %3s %10s %10s %7s  %s\n" "query" "E" "serial" "parallel" "speedup" "status"
echo "--------------------------------------------------------------"

for idx in "${!QUERIES[@]}"; do
    q="${QUERIES[$idx]}"
    name=$(basename "$q")
    edges=$(grep -c '^e ' "$q")

    # Skip if already in summary
    if grep -q "^${name}," "${SUMMARY}" 2>/dev/null; then
        # Print cached result
        line=$(grep "^${name}," "${SUMMARY}")
        sp=$(echo "$line" | cut -d, -f5)
        st=$(echo "$line" | cut -d, -f10)
        printf "[%d/%d] %-14s %3d %43s  %s\n" "$((idx+1))" "$TOTAL" "$name" "$edges" "cached speedup=${sp}" "$st"
        continue
    fi

    log_s="${OUTPUT_DIR}/${name}.serial.log"
    log_p="${OUTPUT_DIR}/${name}.parallel.log"

    # Serial
    ec_s=0; run_one "graphflow" 1 "$q" "$log_s" || ec_s=$?
    ms_s=$(parse_time "$log_s")
    pos_s=$(parse_pos "$log_s")
    neg_s=$(parse_neg "$log_s")

    # Parallel
    ec_p=0; run_one "parallel_graphflow" "${THREADS}" "$q" "$log_p" || ec_p=$?
    ms_p=$(parse_time "$log_p")
    pos_p=$(parse_pos "$log_p")
    neg_p=$(parse_neg "$log_p")

    status="PASS"
    if [[ $ec_s -ne 0 || $ec_p -ne 0 ]]; then
        status="TIMEOUT(s=${ec_s},p=${ec_p})"
    elif [[ "${pos_s:-x}" != "${pos_p:-y}" || "${neg_s:-x}" != "${neg_p:-y}" ]]; then
        status="MISMATCH"
    fi

    speedup="N/A"
    if [[ -n "${ms_s}" && -n "${ms_p}" && "${ms_p}" != "0" ]]; then
        speedup=$(awk "BEGIN{printf \"%.2f\", ${ms_s}/${ms_p}}")
    fi

    printf "[%d/%d] %-14s %3d %10s %10s %7s  %s\n" "$((idx+1))" "$TOTAL" "$name" "$edges" "${ms_s:-?}" "${ms_p:-?}" "${speedup}" "${status}"
    echo "${name},${edges},${ms_s:-?},${ms_p:-?},${speedup},${pos_s:-?},${pos_p:-?},${neg_s:-?},${neg_p:-?},${status}" >> "${SUMMARY}"
done

echo ""
echo "================================================================"
echo "Queries with speedup > 3.0x:"
awk -F, 'NR>1 && $5+0 > 3.0 {printf "  %-14s edges=%-2s  serial=%-10s  parallel=%-10s  speedup=%sx\n", $1, $2, $3, $4, $5}' "${SUMMARY}"
echo ""
echo "Top 10 by speedup:"
tail -n +2 "${SUMMARY}" | sort -t, -k5 -rn | head -10 | \
    awk -F, '{printf "  %-14s edges=%-2s  serial=%-10s  parallel=%-10s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}'
echo ""
echo "Summary: ${SUMMARY}"
