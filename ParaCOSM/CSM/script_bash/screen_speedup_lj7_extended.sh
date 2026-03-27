#!/bin/bash
# =============================================================================
# Extended screen: longer timeout to get serial times for heavy queries,
# plus dense queries that were missed.
# =============================================================================
set -euo pipefail

CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"
QUERY_BASE="${HOME}/haibin/paracosm_data/lj/livejournal/random_walk/7_self"

THREADS=8
TIME_LIMIT=600
RUN_TIMEOUT=900

OUTPUT_DIR="${1:-logs_txt/screen_lj_7self_extended}"
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

# Heavy tree queries + more sparse with known high speedup + all dense
declare -A QUERIES
QUERIES[tree]="Q_5 Q_16 Q_18 Q_28 Q_30 Q_34 Q_36"
QUERIES[dense]="Q_9 Q_24 Q_25 Q_40 Q_54 Q_70 Q_85 Q_92"
QUERIES[sparse]="Q_46 Q_55 Q_65 Q_75 Q_90 Q_95 Q_99"

echo ""
echo "================================================================"
echo " Extended screen: serial graphflow (1T) vs parallel (${THREADS}T)"
echo " TIME_LIMIT=${TIME_LIMIT}s  RUN_TIMEOUT=${RUN_TIMEOUT}s"
echo "================================================================"
printf "%-12s %-6s %12s %12s %8s  %s\n" "query" "type" "serial_ms" "par_ms" "speedup" "status"
echo "--------------------------------------------------------------"

for QTYPE in tree sparse dense; do
    for name in ${QUERIES[${QTYPE}]}; do
        q="${QUERY_BASE}/${QTYPE}/${name}"
        [ -e "${q}" ] || { printf "%-12s %-6s %12s %12s %8s  %s\n" "${name}" "${QTYPE}" "-" "-" "-" "MISSING"; continue; }

        log_s="${OUTPUT_DIR}/${name}.${QTYPE}.serial.log"
        log_p="${OUTPUT_DIR}/${name}.${QTYPE}.parallel.log"

        ec_s=0; run_one "graphflow" 1 "$q" "$log_s" || ec_s=$?
        ms_s=$(parse_time "${log_s}")
        pos_s=$(parse_pos "${log_s}")
        neg_s=$(parse_neg "${log_s}")

        ec_p=0; run_one "parallel_graphflow" "${THREADS}" "$q" "$log_p" || ec_p=$?
        ms_p=$(parse_time "${log_p}")
        pos_p=$(parse_pos "${log_p}")
        neg_p=$(parse_neg "${log_p}")

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

        printf "%-12s %-6s %12s %12s %8s  %s\n" "${name}" "${QTYPE}" "${ms_s:-?}" "${ms_p:-?}" "${speedup}" "${status}"
        echo "${name},${QTYPE},${ms_s:-?},${ms_p:-?},${speedup},${pos_s:-?},${pos_p:-?},${neg_s:-?},${neg_p:-?},${status}" >> "${SUMMARY}"
    done
done

echo ""
echo "================================================================"
echo "Top speedups:"
tail -n +2 "${SUMMARY}" | sort -t, -k5 -rn | head -10 | \
    awk -F, '{printf "  %-12s %-6s  serial=%-12s  parallel=%-12s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}'
echo ""
echo "Summary: ${SUMMARY}"
