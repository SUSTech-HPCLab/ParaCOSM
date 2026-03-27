#!/bin/bash
# =============================================================================
# Fast screen: pick a sample of queries from each type, run serial vs parallel
# =============================================================================
set -euo pipefail

CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"
QUERY_BASE="${HOME}/haibin/paracosm_data/lj/livejournal/random_walk/7_self"

THREADS=8
TIME_LIMIT=180
RUN_TIMEOUT=480

OUTPUT_DIR="${1:-logs_txt/screen_lj_7self_fast}"
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

# Sample queries from each type
SAMPLE_SPARSE="Q_1 Q_3 Q_10 Q_22 Q_42 Q_60 Q_80 Q_100"
SAMPLE_TREE="Q_2 Q_5 Q_16 Q_18 Q_28 Q_30 Q_34 Q_36"
SAMPLE_DENSE="Q_9 Q_24 Q_25 Q_40 Q_54 Q_70 Q_85 Q_92"

echo ""
echo "================================================================"
echo " Serial graphflow (1T) vs Parallel_graphflow (${THREADS}T)"
echo " Dataset: LiveJournal 7_self  |  mode: single (no classify)"
echo "================================================================"
printf "%-12s %-6s %12s %12s %8s  %s\n" "query" "type" "serial_ms" "par_ms" "speedup" "status"
echo "--------------------------------------------------------------"

do_query() {
    local QTYPE="$1" name="$2"
    local q="${QUERY_BASE}/${QTYPE}/${name}"
    [ -e "${q}" ] || { printf "%-12s %-6s %12s %12s %8s  %s\n" "${name}" "${QTYPE}" "-" "-" "-" "MISSING"; return; }

    local log_s="${OUTPUT_DIR}/${name}.${QTYPE}.serial.log"
    local log_p="${OUTPUT_DIR}/${name}.${QTYPE}.parallel.log"

    # Serial
    local ec_s=0; run_one "graphflow" 1 "$q" "$log_s" || ec_s=$?
    local ms_s=$(parse_time "${log_s}")
    local pos_s=$(parse_pos "${log_s}")
    local neg_s=$(parse_neg "${log_s}")

    # Parallel
    local ec_p=0; run_one "parallel_graphflow" "${THREADS}" "$q" "$log_p" || ec_p=$?
    local ms_p=$(parse_time "${log_p}")
    local pos_p=$(parse_pos "${log_p}")
    local neg_p=$(parse_neg "${log_p}")

    local status="PASS"
    if [[ $ec_s -ne 0 || $ec_p -ne 0 ]]; then
        status="TIMEOUT(s=${ec_s},p=${ec_p})"
    elif [[ "${pos_s:-x}" != "${pos_p:-y}" || "${neg_s:-x}" != "${neg_p:-y}" ]]; then
        status="MISMATCH"
    fi

    local speedup="N/A"
    if [[ -n "${ms_s}" && -n "${ms_p}" && "${ms_p}" != "0" ]]; then
        speedup=$(awk "BEGIN{printf \"%.2f\", ${ms_s}/${ms_p}}")
    fi

    printf "%-12s %-6s %12s %12s %8s  %s\n" "${name}" "${QTYPE}" "${ms_s:-?}" "${ms_p:-?}" "${speedup}" "${status}"
    echo "${name},${QTYPE},${ms_s:-?},${ms_p:-?},${speedup},${pos_s:-?},${pos_p:-?},${neg_s:-?},${neg_p:-?},${status}" >> "${SUMMARY}"
}

for name in ${SAMPLE_SPARSE}; do do_query sparse "${name}"; done
for name in ${SAMPLE_TREE};   do do_query tree   "${name}"; done
for name in ${SAMPLE_DENSE};  do do_query dense  "${name}"; done

echo ""
echo "================================================================"
echo "Top speedup queries:"
tail -n +2 "${SUMMARY}" | sort -t, -k5 -rn | head -10 | \
    awk -F, '{printf "  %-12s %-6s  serial=%-12s  parallel=%-12s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}'
echo ""
echo "Summary CSV: ${SUMMARY}"
echo "Done."
