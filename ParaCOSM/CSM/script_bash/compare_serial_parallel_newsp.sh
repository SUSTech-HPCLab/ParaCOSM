#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SUBPROJECT_ROOT="${REPO_ROOT}/matching_executor/Parallel_NewSP"

DATA_ROOT="${DATA_ROOT:-/home/v-haibinlai/haibin/paracosm_data}"
DATA_GRAPH="${DATA_GRAPH:-${DATA_ROOT}/6/data_graph/data.graph}"
UPDATE_GRAPH="${UPDATE_GRAPH:-${DATA_ROOT}/6/data_graph/insertion.graph}"
QUERY_DIR="${QUERY_DIR:-${DATA_ROOT}/8_self/sparse}"

PARALLEL_EXE="${PARALLEL_EXE:-${SUBPROJECT_ROOT}/build/csm}"
SERIAL_EXE="${SERIAL_EXE:-${SUBPROJECT_ROOT}/build/csm_serial}"

SERIAL_THREADS="${SERIAL_THREADS:-1}"
PARALLEL_THREADS="${PARALLEL_THREADS:-8}"
TIME_LIMIT="${TIME_LIMIT:-60}"
INITIAL_TIME_LIMIT="${INITIAL_TIME_LIMIT:-60}"
RUN_TIMEOUT="${RUN_TIMEOUT:-180}"
REPORT_INITIAL="${REPORT_INITIAL:-1}"

QUERY_IDS_INPUT="${QUERY_IDS:-3 4 5}"
read -r -a QUERY_IDS <<< "${QUERY_IDS_INPUT}"

LOG_DIR="${LOG_DIR:-${REPO_ROOT}/logs_txt/serial_vs_parallel_newsp}"
mkdir -p "${LOG_DIR}"
mkdir -p "${SUBPROJECT_ROOT}/build/obj"

build_serial_reference() {
    echo "[BUILD] serial_newsp reference executable"
    (
        cd "${SUBPROJECT_ROOT}"
        g++ -std=c++17 -O2 -fopenmp -DSERIAL_NEWSP_ROOT_DISPATCH -I. -I./tbb/include \
            -c matching/CSMPP.cpp -o build/obj/CSMPP.serial.o
        g++ -std=c++17 -O2 -fopenmp -I. -I./tbb/include \
            graph/graph.cpp \
            matching/matching.cpp \
            matching/main.cpp \
            build/obj/CSMPP.serial.o \
            utils/globals.cpp \
            DecisionMakingSystem/DecisionMakingSystem.cpp \
            -o build/csm_serial
    )
}

build_parallel_if_missing() {
    if [[ -x "${PARALLEL_EXE}" ]]; then
        return
    fi

    echo "[BUILD] parallel_newsp executable"
    (
        cd "${SUBPROJECT_ROOT}"
        g++ -std=c++17 -O2 -fopenmp -I. -I./tbb/include \
            -c matching/CSMPP.cpp -o build/obj/CSMPP.o
        g++ -std=c++17 -O2 -fopenmp -I. -I./tbb/include \
            graph/graph.cpp \
            matching/matching.cpp \
            matching/main.cpp \
            build/obj/CSMPP.o \
            utils/globals.cpp \
            DecisionMakingSystem/DecisionMakingSystem.cpp \
            -o build/csm
    )
}

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
    local exe="$1"
    local threads="$2"
    local query_path="$3"
    local log_file="$4"

    (
        cd "$(dirname "${exe}")"
        export OMP_NUM_THREADS="${threads}"
        timeout --kill-after=5 "${RUN_TIMEOUT}" \
            "${exe}" \
            -d "${DATA_GRAPH}" \
            -u "${UPDATE_GRAPH}" \
            -q "${query_path}" \
            --report-initial "${REPORT_INITIAL}" \
            --initial-time-limit "${INITIAL_TIME_LIMIT}" \
            --time-limit "${TIME_LIMIT}"
    ) >"${log_file}" 2>&1
}

if [[ ! -f "${DATA_GRAPH}" || ! -f "${UPDATE_GRAPH}" ]]; then
    echo "[ERROR] data/update graph path invalid" >&2
    exit 1
fi

build_parallel_if_missing
build_serial_reference

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

echo "serial_newsp vs parallel_newsp"
echo "  serial exe      : ${SERIAL_EXE}"
echo "  parallel exe    : ${PARALLEL_EXE}"
echo "  data graph      : ${DATA_GRAPH}"
echo "  update graph    : ${UPDATE_GRAPH}"
echo "  query dir       : ${QUERY_DIR}"
echo "  threads         : ${SERIAL_THREADS} vs ${PARALLEL_THREADS}"
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
    serial_log="${LOG_DIR}/Q_${query_id}.serial.log"
    parallel_log="${LOG_DIR}/Q_${query_id}.parallel.log"

    echo "[RUN ] Q_${query_id} serial_newsp"
    if ! run_case "${SERIAL_EXE}" "${SERIAL_THREADS}" "${query_path}" "${serial_log}"; then
        echo "[FAIL] Q_${query_id} serial_newsp run failed"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "[RUN ] Q_${query_id} parallel_newsp"
    if ! run_case "${PARALLEL_EXE}" "${PARALLEL_THREADS}" "${query_path}" "${parallel_log}"; then
        echo "[FAIL] Q_${query_id} parallel_newsp run failed"
        FAILED=$((FAILED + 1))
        continue
    fi

    serial_initial="$(extract_initial_matches "${serial_log}")"
    parallel_initial="$(extract_initial_matches "${parallel_log}")"
    serial_positive="$(extract_metric "${serial_log}" "positive")"
    parallel_positive="$(extract_metric "${parallel_log}" "positive")"
    serial_negative="$(extract_metric "${serial_log}" "negative")"
    parallel_negative="$(extract_metric "${parallel_log}" "negative")"
    serial_time="$(extract_time_ms "${serial_log}")"
    parallel_time="$(extract_time_ms "${parallel_log}")"

    if [[ "${serial_initial}" == "${parallel_initial}" && "${serial_positive}" == "${parallel_positive}" && "${serial_negative}" == "${parallel_negative}" ]]; then
        echo "[PASS] Q_${query_id} initial=${serial_initial} positive=${serial_positive} negative=${serial_negative} time_ms=${serial_time}/${parallel_time}"
        PASSED=$((PASSED + 1))
    else
        echo "[FAIL] Q_${query_id} mismatch"
        echo "       serial   : initial=${serial_initial} positive=${serial_positive} negative=${serial_negative} time_ms=${serial_time}"
        echo "       parallel : initial=${parallel_initial} positive=${parallel_positive} negative=${parallel_negative} time_ms=${parallel_time}"
        FAILED=$((FAILED + 1))
    fi
done

echo
echo "Summary: total=${TOTAL} passed=${PASSED} failed=${FAILED} skipped=${SKIPPED}"

if [[ "${FAILED}" -ne 0 ]]; then
    exit 1
fi