#!/usr/bin/env bash
# ------------------------------------------------------------
# ParaCOSM Test Runner (tidied & reusable)
# ------------------------------------------------------------
# Defaults are set to match the original snippet:
#   - algorithm: parallel_symbi
#   - dataset : livejournal
#   - suffixes: (6)
#   - threads : 32
#   - query i : from 3 to 3
#   - timeouts: main 3600s with --kill-after=5
#
# Usage examples:
#   ./run_paracosm.sh
#   ./run_paracosm.sh -A parallel_symbi -D livejournal -s 6 -t 32 -q 3-3 \
#       -o logs_txt/amazon/Paracosm
# ------------------------------------------------------------

set -Eeuo pipefail
IFS=$'\n\t'

# ---------- Defaults (can be overridden via CLI) ----------
export LD_LIBRARY_PATH="/home/cc/intel/oneapi/tbb/latest/lib:${LD_LIBRARY_PATH:-}"

BIN="./build/bin/csm"
ALGORITHM="parallel_symbi"
DATASET_NAME="livejournal"
SUFFIXES=("6")
THREADS=32
QUERY_RANGE="3-3"          # inclusive range: START-END
TIME_LIMIT=1800            # csm's --time-limit
TIMEOUT_SEC=3600           # outer timeout
KILL_AFTER=5               # timeout --kill-after
OUTPUT_DIR="logs_txt/${DATASET_NAME}/Graphflow"  # default log dir
LOGGING=true               # toggle with --no-log

# Data roots
DIR="/home/cc/haibin2"
DATA_GRAPH="/home/cc/haibin2/paracosm/${DATASET_NAME}/data_graph/data.graph"
INSERT_GRAPH="/home/cc/haibin2/paracosm/${DATASET_NAME}/data_graph/insertion.graph"

# Query configs (type path). You can add 'dense' and 'tree' similarly.
# 注意：若你想强制覆盖查询目录，可在下面“OVERRIDE_QUERY_DIR”里设置。
CONFIGS=(
  "sparse ${DIR}/livejournal/6/query_graph/sparse_6"
)

# 如果需要与原脚本一致地“硬覆盖”查询图目录，请在此赋值；留空则使用 CONFIGS 中的目录
OVERRIDE_QUERY_DIR="/home/cc/haibin2/paracosm/livejournal/random_walk/6_self/sparse"

# ---------- Helpers ----------
die() { echo "Error: $*" >&2; exit 1; }

usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -A, --algorithm NAME       Algorithm (default: ${ALGORITHM})
  -D, --dataset  NAME       Dataset name (default: ${DATASET_NAME})
  -s, --suffixes  CSV       Comma-separated suffix list (default: ${SUFFIXES[*]})
  -t, --threads   N         Threads (default: ${THREADS})
  -q, --queries   A-B       Query index inclusive range (default: ${QUERY_RANGE})
  -o, --outdir    DIR       Output directory for logs (default: ${OUTPUT_DIR})
      --bin       PATH      Path to csm binary (default: ${BIN})
      --time-limit SEC      In-app time limit (default: ${TIME_LIMIT})
      --timeout    SEC      Outer timeout seconds (default: ${TIMEOUT_SEC})
      --kill-after SEC      timeout --kill-after (default: ${KILL_AFTER})
      --no-log               Do not write logs to files (stdout only)
  -h, --help

Examples:
  $0 -A parallel_symbi -D livejournal -s 6 -t 32 -q 3-3 -o logs_txt/livejournal/Graphflow
EOF
}

parse_range() {
  local rg="$1"
  [[ "$rg" =~ ^([0-9]+)-([0-9]+)$ ]] || die "Bad range: $rg (expect A-B)"
  echo "${BASH_REMATCH[1]} ${BASH_REMATCH[2]}"
}

# ---------- Parse CLI ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -A|--algorithm) ALGORITHM="$2"; shift 2;;
    -D|--dataset)   DATASET_NAME="$2"; shift 2;;
    -s|--suffixes)  IFS=',' read -r -a SUFFIXES <<< "$2"; shift 2;;
    -t|--threads)   THREADS="${2}"; shift 2;;
    -q|--queries)   QUERY_RANGE="$2"; shift 2;;
    -o|--outdir)    OUTPUT_DIR="$2"; shift 2;;
    --bin)          BIN="$2"; shift 2;;
    --time-limit)   TIME_LIMIT="$2"; shift 2;;
    --timeout)      TIMEOUT_SEC="$2"; shift 2;;
    --kill-after)   KILL_AFTER="$2"; shift 2;;
    --no-log)       LOGGING=false; shift 1;;
    -h|--help)      usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# Update paths that depend on DATASET_NAME
DATA_GRAPH="/home/cc/haibin2/paracosm/${DATASET_NAME}/data_graph/data.graph"
INSERT_GRAPH="/home/cc/haibin2/paracosm/${DATASET_NAME}/data_graph/insertion.graph"

# ---------- Validation ----------
[[ -x "$BIN" ]] || die "Binary not found or not executable: $BIN"
[[ -f "$DATA_GRAPH" ]] || die "Data graph not found: $DATA_GRAPH"
[[ -f "$INSERT_GRAPH" ]] || die "Insertion graph not found: $INSERT_GRAPH"

# ensure output dir (if logging)
if $LOGGING; then
  mkdir -p "$OUTPUT_DIR"
fi

tmp="$(parse_range "$QUERY_RANGE")"
IFS=' ' read -r Q_START Q_END <<< "$tmp"
(( Q_START <= Q_END )) || die "Query range invalid (start > end): $QUERY_RANGE"

# ---------- Runner ----------
run_one() {
  local type="$1"
  local qdir="$2"
  local suffix="$3"

  # override query dir if requested (to mimic your original override behavior)
  if [[ -n "${OVERRIDE_QUERY_DIR:-}" ]]; then
    qdir="$OVERRIDE_QUERY_DIR"
  fi

  [[ -d "$qdir" ]] || die "Query dir not found: $qdir"

  local data_set="${DATASET_NAME}/${suffix}"
  local out_file="${ALGORITHM}_${DATASET_NAME^}_${suffix}_P2.txt"  # e.g., parallel_symbi_Livejournal_6_P2.txt
  local target_dir="${OUTPUT_DIR}/${type}_${suffix}"
  $LOGGING && mkdir -p "$target_dir"

  echo "==> Start: type=${type}, suffix=${suffix}, dataset=${data_set}, qdir=${qdir}"

  for (( i=Q_START; i<=Q_END; ++i )); do
    local qpath="${qdir}/Q_${i}"
    [[ -d "$qpath" || -f "$qpath" ]] || die "Query path not found: ${qpath}"

    echo "--> Running query: ${qpath}"
    # Compose command
    cmd=( timeout --kill-after="${KILL_AFTER}" "${TIMEOUT_SEC}" "${BIN}"
        -a "${ALGORITHM}"
        -d "${DATA_GRAPH}"
        -u "${INSERT_GRAPH}"
        -q "${qpath}"
        --time-limit "${TIME_LIMIT}"
        --report-initial off
        -t "${THREADS}"
        --auto-tuning 0
    )

    if $LOGGING; then
      # Log to file and stdout
      "${cmd[@]}" | tee -a "${target_dir}/${out_file}"
    else
      "${cmd[@]}"
    fi

    echo "<-- Finished query: ${qpath}"
  done

  echo "==> Done: type=${type}, suffix=${suffix}"
}

# ---------- Main loop ----------
for suffix in "${SUFFIXES[@]}"; do
  for cfg in "${CONFIGS[@]}"; do
    # split "type path"
    read -r type qdir <<<"$cfg"
    run_one "$type" "$qdir" "$suffix"
  done
done

echo "All tests completed for algorithm=${ALGORITHM}, suffixes=(${SUFFIXES[*]})."
