#!/bin/bash
# =============================================================================
# Quick screen: serial graphflow vs parallel_graphflow on LJ 7_self
# Phase 1: fast scan with short time limit to identify heavy queries
# Phase 2: re-run heavy ones with longer limit
# =============================================================================
set -euo pipefail

CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"
QUERY_BASE="${HOME}/haibin/paracosm_data/lj/livejournal/random_walk/7_self"

THREADS=8
TIME_LIMIT=${TIME_LIMIT:-180}
RUN_TIMEOUT=${RUN_TIMEOUT:-400}

OUTPUT_DIR="${1:-logs_txt/screen_lj_7self}"
mkdir -p "${OUTPUT_DIR}"

SUMMARY="${OUTPUT_DIR}/summary.csv"
echo "query,type,serial_ms,parallel_ms,speedup,pos_serial,pos_parallel,neg_serial,neg_parallel,status" > "${SUMMARY}"

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

echo ""
echo "================================================================"
echo " Serial graphflow (1T) vs Parallel_graphflow (${THREADS}T)"
echo " Dataset: LiveJournal 7_self  |  update-mode: single (no classify)"
echo " TIME_LIMIT=${TIME_LIMIT}s  RUN_TIMEOUT=${RUN_TIMEOUT}s"
echo "================================================================"

# Collect query list from all types; randomize order for diversity
QUERIES=()
for QTYPE in sparse tree dense; do
    QDIR="${QUERY_BASE}/${QTYPE}"
    [ -d "${QDIR}" ] || continue
    for q in $(ls -d "${QDIR}"/Q_* 2>/dev/null | sort -t_ -k2 -n); do
        QUERIES+=("${QTYPE}|${q}")
    done
done

TOTAL=${#QUERIES[@]}
echo "Total queries to screen: ${TOTAL}"
echo ""

printf "%-12s %-6s %10s %10s %8s  %s\n" "query" "type" "serial_ms" "par_ms" "speedup" "status"
echo "-------------------------------------------------------------"

IDX=0
for entry in "${QUERIES[@]}"; do
    IFS='|' read -r QTYPE q <<< "${entry}"
    name=$(basename "$q")
    IDX=$((IDX+1))

    log_s="${OUTPUT_DIR}/${name}.${QTYPE}.serial.log"
    log_p="${OUTPUT_DIR}/${name}.${QTYPE}.parallel.log"

    # Skip if already completed
    if [ -f "${log_s}" ] && grep -q "Incremental Matching:" "${log_s}" 2>/dev/null &&
       [ -f "${log_p}" ] && grep -q "Incremental Matching:" "${log_p}" 2>/dev/null; then
        ms_s=$(parse_time "${log_s}")
        ms_p=$(parse_time "${log_p}")
        pos_s=$(parse_pos "${log_s}")
        pos_p=$(parse_pos "${log_p}")
        neg_s=$(parse_neg "${log_s}")
        neg_p=$(parse_neg "${log_p}")
        status="CACHED"
    else
        # Run serial
        ec_s=0; run_one "graphflow" 1 "$q" "$log_s" || ec_s=$?
        ms_s=$(parse_time "${log_s}")
        pos_s=$(parse_pos "${log_s}")
        neg_s=$(parse_neg "${log_s}")

        # Run parallel
        ec_p=0; run_one "parallel_graphflow" "${THREADS}" "$q" "$log_p" || ec_p=$?
        ms_p=$(parse_time "${log_p}")
        pos_p=$(parse_pos "${log_p}")
        neg_p=$(parse_neg "${log_p}")

        if [[ $ec_s -ne 0 || $ec_p -ne 0 ]]; then
            status="TIMEOUT(s=${ec_s},p=${ec_p})"
        elif [[ "${pos_s}" == "${pos_p}" && "${neg_s}" == "${neg_p}" ]]; then
            status="PASS"
        else
            status="MISMATCH"
        fi
    fi

    speedup="N/A"
    if [[ -n "${ms_s}" && -n "${ms_p}" && "${ms_p}" != "0" ]]; then
        speedup=$(awk "BEGIN{printf \"%.2f\", ${ms_s}/${ms_p}}")
    fi

    printf "[%d/%d] %-12s %-6s %10s %10s %8s  %s\n" "${IDX}" "${TOTAL}" "${name}" "${QTYPE}" "${ms_s:-?}" "${ms_p:-?}" "${speedup}" "${status}"
    echo "${name},${QTYPE},${ms_s:-?},${ms_p:-?},${speedup},${pos_s:-?},${pos_p:-?},${neg_s:-?},${neg_p:-?},${status}" >> "${SUMMARY}"
done

echo ""
echo "================================================================"
echo "Results with speedup > 1.5x:"
awk -F, 'NR>1 && $5+0 > 1.5 {printf "  %-12s %-6s  serial=%-12s  parallel=%-12s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}' "${SUMMARY}"
echo ""
echo "Results with speedup > 2.0x:"
awk -F, 'NR>1 && $5+0 > 2.0 {printf "  %-12s %-6s  serial=%-12s  parallel=%-12s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}' "${SUMMARY}"
echo ""
echo "All results sorted by speedup (descending):"
head -1 "${SUMMARY}"
tail -n +2 "${SUMMARY}" | sort -t, -k5 -rn | head -20
echo ""
echo "Summary CSV: ${SUMMARY}"
echo "Done."
