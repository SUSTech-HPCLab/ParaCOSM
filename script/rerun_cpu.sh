#!/usr/bin/env bash
# 重跑 CPU 失败 (FAIL_129) 任务: 直接读 summary.csv 提取 cpu/FAIL_129 行
# 用 setsid 完全脱离控制终端, 防止 SIGHUP
set -u

OUT=${OUT:-/home/superbench/haibin/paracosm/bench_lj_8v_20260602_174246}
BIN=${CSM_BIN:-/home/superbench/haibin/paracosm/ParaCOSM-gpu/ParaCOSM/CSM/build/bin/csm}
ROOT=/home/superbench/haibin/paracosm
DATA=$ROOT/LJ/data_graph/data.graph
UPD=$ROOT/LJ/data_graph/insertion.graph
QDIR=$ROOT/LJ/8v
TIMEOUT_S=${TIMEOUT_S:-1800}
PARALLEL=${PARALLEL:-50}

JOBS=$OUT/cpu_rerun.tsv
SUMMARY=$OUT/summary_rerun.csv

[[ -f $JOBS ]] || { echo "no $JOBS"; exit 1; }
echo "mode,algo,query,status,wall_ms,positive,negative,update_ms" > "$SUMMARY"

run_one() {
    local algo=$1 q=$2
    local qpath="$QDIR/$q"
    local log="$OUT/cpu/$algo/$q.log"
    local t0 t1 rc wall
    t0=$(date +%s%3N)
    timeout --foreground "${TIMEOUT_S}s" \
      env CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 \
      "$BIN" -q "$qpath" -d "$DATA" -u "$UPD" \
      -a "$algo" -m single --time-limit "$TIMEOUT_S" \
      -t 1 --auto-tuning 0 --report-initial 0 \
      >"$log" 2>&1
    rc=$?
    t1=$(date +%s%3N); wall=$((t1 - t0))
    local status
    case $rc in
        0)   status=OK ;;
        124) status=TIMEOUT ;;
        137) status=OOM_OR_KILL ;;
        139) status=SEGV ;;
        *)   status=FAIL_$rc ;;
    esac
    local pos neg upd_ms
    pos=$(grep -aE 'positive_num_results|positive results' "$log" | tail -1 | grep -oE '[0-9]+' | tail -1)
    neg=$(grep -aE 'negative_num_results|negative results' "$log" | tail -1 | grep -oE '[0-9]+' | tail -1)
    upd_ms=$(grep -aE 'Update Time|Incremental .* Time|Total .* time' "$log" | tail -1 | grep -oE '[0-9]+(\.[0-9]+)?' | tail -1)
    printf 'cpu,%s,%s,%s,%d,%s,%s,%s\n' \
        "$algo" "$q" "$status" "$wall" "${pos:-}" "${neg:-}" "${upd_ms:-}" >> "$SUMMARY"
    printf '[cpu-rerun] %-19s %-7s -> %-10s %6dms\n' "$algo" "$q" "$status" "$wall"
}
export -f run_one
export BIN DATA UPD QDIR OUT SUMMARY TIMEOUT_S

echo "=== rerun $(wc -l <"$JOBS") cpu jobs, ${PARALLEL}-way parallel ==="
xargs -a "$JOBS" -P "$PARALLEL" -n 1 -I{} \
    bash -c 'IFS=$'"'"'\t'"'"' read -r a q <<<"{}"; run_one "$a" "$q"'

echo "=== done. summary: $SUMMARY ==="
awk -F, 'NR>1{c[$4]++} END{for(k in c) print k": "c[k]}' "$SUMMARY" | sort
