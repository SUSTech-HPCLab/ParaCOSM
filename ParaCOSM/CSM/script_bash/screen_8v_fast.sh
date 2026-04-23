#!/bin/bash
# =============================================================================
# Fast screen: run parallel first (short timeout), then serial only if parallel
# finished. This skips truly intractable queries quickly.
# =============================================================================
set -euo pipefail

CSM_BIN="./build/bin/csm"
DATA_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/data.graph"
INSERT_GRAPH="${HOME}/haibin/paracosm_data/lj/livejournal/data_graph/insertion.graph"

QUERY_DIR="${1:?Usage: $0 <query_dir> [output_dir]}"
OUTPUT_DIR="${2:-logs_txt/screen_8v_fast}"
mkdir -p "${OUTPUT_DIR}"

THREADS=8
# Short limits: parallel must finish in 120s or we skip
PAR_TIME_LIMIT=${PAR_TIME_LIMIT:-60}
PAR_TIMEOUT=${PAR_TIMEOUT:-180}
# Serial gets more time since it's slower
SER_TIME_LIMIT=${SER_TIME_LIMIT:-180}
SER_TIMEOUT=${SER_TIMEOUT:-480}

SUMMARY="${OUTPUT_DIR}/summary.csv"
if [ ! -s "${SUMMARY}" ]; then
    echo "query,edges,serial_ms,parallel_ms,speedup,pos_serial,pos_parallel,neg_serial,neg_parallel,status" > "${SUMMARY}"
fi

parse_time() { awk '/^Incremental Matching:/ {gsub(/ms/,"",$3); print $3; exit}' "$1"; }
parse_pos()  { awk '$2=="positive" && $3=="matches." {print $1; exit}' "$1"; }
parse_neg()  { awk '$2=="negative" && $3=="matches." {print $1; exit}' "$1"; }

QUERIES=( $(ls "${QUERY_DIR}"/Q_gen_* 2>/dev/null | sort) )
TOTAL=${#QUERIES[@]}

echo ""
echo "================================================================"
echo " Fast screen: parallel first, then serial if feasible"
echo " Queries: ${TOTAL}  PAR_LIMIT=${PAR_TIME_LIMIT}s  SER_LIMIT=${SER_TIME_LIMIT}s"
echo "================================================================"
printf "%-14s %3s %10s %10s %7s  %s\n" "query" "E" "serial" "parallel" "speedup" "status"
echo "--------------------------------------------------------------"

DONE=0
SKIPPED=0
GOOD=0

for idx in "${!QUERIES[@]}"; do
    q="${QUERIES[$idx]}"
    name=$(basename "$q")
    edges=$(grep -c '^e ' "$q")

    # Skip if already in summary
    if grep -q "^${name}," "${SUMMARY}" 2>/dev/null; then
        DONE=$((DONE+1))
        continue
    fi

    log_p="${OUTPUT_DIR}/${name}.parallel.log"
    log_s="${OUTPUT_DIR}/${name}.serial.log"

    # Step 1: Run parallel with short timeout
    ec_p=0
    timeout --kill-after=5 "${PAR_TIMEOUT}" \
        "${CSM_BIN}" -a parallel_graphflow \
        -d "${DATA_GRAPH}" -u "${INSERT_GRAPH}" -q "$q" \
        --time-limit "${PAR_TIME_LIMIT}" --report-initial off \
        -t "${THREADS}" --auto-tuning 0 -m single \
        > "$log_p" 2>&1 || ec_p=$?

    ms_p=$(parse_time "$log_p")
    pos_p=$(parse_pos "$log_p")
    neg_p=$(parse_neg "$log_p")

    # If parallel timed out, skip serial entirely
    if [[ $ec_p -ne 0 || -z "$ms_p" ]]; then
        SKIPPED=$((SKIPPED+1))
        printf "[%d/%d] %-14s %3d %10s %10s %7s  %s\n" "$((idx+1))" "$TOTAL" "$name" "$edges" "skip" "${ms_p:-timeout}" "N/A" "PAR_TIMEOUT"
        echo "${name},${edges},skip,${ms_p:-timeout},N/A,skip,${pos_p:-?},skip,${neg_p:-?},PAR_TIMEOUT" >> "${SUMMARY}"
        continue
    fi

    # Step 2: Run serial with longer timeout
    ec_s=0
    timeout --kill-after=5 "${SER_TIMEOUT}" \
        "${CSM_BIN}" -a graphflow \
        -d "${DATA_GRAPH}" -u "${INSERT_GRAPH}" -q "$q" \
        --time-limit "${SER_TIME_LIMIT}" --report-initial off \
        -t 1 --auto-tuning 0 -m single \
        > "$log_s" 2>&1 || ec_s=$?

    ms_s=$(parse_time "$log_s")
    pos_s=$(parse_pos "$log_s")
    neg_s=$(parse_neg "$log_s")

    status="PASS"
    if [[ $ec_s -ne 0 ]]; then
        status="SER_TIMEOUT"
    elif [[ "${pos_s:-x}" != "${pos_p:-y}" || "${neg_s:-x}" != "${neg_p:-y}" ]]; then
        status="MISMATCH"
    fi

    speedup="N/A"
    if [[ -n "${ms_s}" && -n "${ms_p}" && "${ms_p}" != "0" ]]; then
        speedup=$(awk "BEGIN{printf \"%.2f\", ${ms_s}/${ms_p}}")
        if (( $(awk "BEGIN{print (${ms_s}/${ms_p} > 3.0)}") )); then
            GOOD=$((GOOD+1))
        fi
    fi

    DONE=$((DONE+1))
    printf "[%d/%d] %-14s %3d %10s %10s %7s  %s\n" "$((idx+1))" "$TOTAL" "$name" "$edges" "${ms_s:-timeout}" "${ms_p}" "${speedup}" "${status}"
    echo "${name},${edges},${ms_s:-timeout},${ms_p},${speedup},${pos_s:-?},${pos_p},${neg_s:-?},${neg_p},${status}" >> "${SUMMARY}"
done

echo ""
echo "================================================================"
echo "Completed: done=${DONE} skipped=${SKIPPED} good(>3x)=${GOOD}"
echo ""
echo "Queries with speedup > 3.0x:"
awk -F, 'NR>1 && $5+0 > 3.0 {printf "  %-14s edges=%-2s  serial=%-10s  parallel=%-10s  speedup=%sx\n", $1, $2, $3, $4, $5}' "${SUMMARY}"
echo ""
echo "Top 20 by speedup:"
tail -n +2 "${SUMMARY}" | sort -t, -k5 -rn | head -20 | \
    awk -F, '{printf "  %-14s edges=%-2s  serial=%-10s  parallel=%-10s  speedup=%sx  %s\n", $1, $2, $3, $4, $5, $10}'
echo ""
echo "Summary: ${SUMMARY}"
