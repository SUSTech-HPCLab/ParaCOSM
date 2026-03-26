#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXECUTABLE="${EXECUTABLE:-${REPO_ROOT}/matching_executor/Parallel_NewSP/build/csm}"
DATA_ROOT="${DATA_ROOT:-/home/v-haibinlai/haibin/paracosm_data}"
DATA_GRAPH="${DATA_GRAPH:-${DATA_ROOT}/6/data_graph/data.graph}"
UPDATE_GRAPH="${UPDATE_GRAPH:-${DATA_ROOT}/6/data_graph/insertion.graph}"
QUERY_DIR="${QUERY_DIR:-${DATA_ROOT}/8_self/sparse}"

BASELINE_THREADS="${BASELINE_THREADS:-1}"
TEST_THREADS="${TEST_THREADS:-8}"
TIME_LIMIT="${TIME_LIMIT:-60}"
INITIAL_TIME_LIMIT="${INITIAL_TIME_LIMIT:-60}"
RUN_TIMEOUT="${RUN_TIMEOUT:-180}"
REPORT_INITIAL="${REPORT_INITIAL:-1}"
KEEP_LOGS="${KEEP_LOGS:-1}"
LOG_DIR="${LOG_DIR:-${REPO_ROOT}/logs_txt/parallel_newsp_correctness}"

QUERY_IDS_INPUT="${QUERY_IDS:-3 4 5}"
read -r -a QUERY_IDS <<< "${QUERY_IDS_INPUT}"

mkdir -p "${LOG_DIR}"

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "[ERROR] executable not found or not executable: ${EXECUTABLE}" >&2
    exit 1
fi

if [[ ! -f "${DATA_GRAPH}" ]]; then
    echo "[ERROR] data graph not found: ${DATA_GRAPH}" >&2
    exit 1
fi

if [[ ! -f "${UPDATE_GRAPH}" ]]; then
    echo "[ERROR] update graph not found: ${UPDATE_GRAPH}" >&2
    exit 1
fi

extract_initial_matches() {
    local log_file="$1"
    awk 'BEGIN {found = 0} /^[0-9]+ matches[.][[:space:]]*$/ {print $1; found = 1; exit} END {if (found == 0) print "NA"}' "${log_file}"
}

extract_metric() {
    local log_file="$1"
    local label="$2"
    awk -v label="$label" 'BEGIN {found = 0} $2 == label && $3 == "matches." {print $1; found = 1; exit} END {if (found == 0) print "NA"}' "${log_file}"
}

extract_time_ms() {
    local log_file="$1"
    awk '/^#Time:/ {print $2; exit}' "${log_file}"
}

run_case() {
    local threads="$1"
    local query_path="$2"
    local log_file="$3"

    (
        cd "$(dirname "${EXECUTABLE}")"
        export OMP_NUM_THREADS="${threads}"
        timeout --kill-after=5 "${RUN_TIMEOUT}" \
            "${EXECUTABLE}" \
            -d "${DATA_GRAPH}" \
            -u "${UPDATE_GRAPH}" \
            -q "${query_path}" \
            --report-initial "${REPORT_INITIAL}" \
            --initial-time-limit "${INITIAL_TIME_LIMIT}" \
            --time-limit "${TIME_LIMIT}"
    ) >"${log_file}" 2>&1
}

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

echo "Parallel_NewSP correctness check"
echo "  executable      : ${EXECUTABLE}"
echo "  data graph      : ${DATA_GRAPH}"
echo "  update graph    : ${UPDATE_GRAPH}"
echo "  query dir       : ${QUERY_DIR}"
echo "  threads         : ${BASELINE_THREADS} vs ${TEST_THREADS}"
echo "  query ids       : ${QUERY_IDS[*]}"
echo "  log dir         : ${LOG_DIR}"
echo

for query_id in "${QUERY_IDS[@]}"; do
    query_path="${QUERY_DIR}/Q_${query_id}"
    if [[ ! -f "${query_path}" ]]; then
        echo "[SKIP] missing query: ${query_path}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    TOTAL=$((TOTAL + 1))

    base_log="${LOG_DIR}/Q_${query_id}.t${BASELINE_THREADS}.log"
    test_log="${LOG_DIR}/Q_${query_id}.t${TEST_THREADS}.log"

    echo "[RUN ] Q_${query_id} baseline (${BASELINE_THREADS} thread)"
    if ! run_case "${BASELINE_THREADS}" "${query_path}" "${base_log}"; then
        echo "[FAIL] Q_${query_id} baseline run failed"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "[RUN ] Q_${query_id} parallel (${TEST_THREADS} threads)"
    if ! run_case "${TEST_THREADS}" "${query_path}" "${test_log}"; then
        echo "[FAIL] Q_${query_id} parallel run failed"
        FAILED=$((FAILED + 1))
        continue
    fi

    base_initial="$(extract_initial_matches "${base_log}")"
    test_initial="$(extract_initial_matches "${test_log}")"
    base_positive="$(extract_metric "${base_log}" "positive")"
    test_positive="$(extract_metric "${test_log}" "positive")"
    base_negative="$(extract_metric "${base_log}" "negative")"
    test_negative="$(extract_metric "${test_log}" "negative")"
    base_time="$(extract_time_ms "${base_log}")"
    test_time="$(extract_time_ms "${test_log}")"

    if [[ "${base_initial}" == "${test_initial}" && "${base_positive}" == "${test_positive}" && "${base_negative}" == "${test_negative}" ]]; then
        echo "[PASS] Q_${query_id} initial=${base_initial} positive=${base_positive} negative=${base_negative} time_ms=${base_time}/${test_time}"
        PASSED=$((PASSED + 1))
        if [[ "${KEEP_LOGS}" == "0" ]]; then
            rm -f "${base_log}" "${test_log}"
        fi
    else
        echo "[FAIL] Q_${query_id} mismatch"
        echo "       baseline: initial=${base_initial} positive=${base_positive} negative=${base_negative} time_ms=${base_time}"
        echo "       parallel: initial=${test_initial} positive=${test_positive} negative=${test_negative} time_ms=${test_time}"
        FAILED=$((FAILED + 1))
    fi
done

echo
echo "Summary: total=${TOTAL} passed=${PASSED} failed=${FAILED} skipped=${SKIPPED}"

if [[ "${FAILED}" -ne 0 ]]; then
    exit 1
fi