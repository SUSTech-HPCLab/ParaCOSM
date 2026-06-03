#!/usr/bin/env bash
# 多 GPU 并行重跑 GPU 阶段剩余 + FAIL_129 任务
# 用法: GPUS="0 2 3" bash rerun_gpu_multi.sh
set -u

OUT=${OUT:-/home/superbench/haibin/paracosm/bench_lj_8v_20260602_174246}
BIN=${CSM_BIN:-/home/superbench/haibin/paracosm/ParaCOSM-gpu/ParaCOSM/CSM/build/bin/csm}
ROOT=/home/superbench/haibin/paracosm
DATA=$ROOT/LJ/data_graph/data.graph
UPD=$ROOT/LJ/data_graph/insertion.graph
QDIR=$ROOT/LJ/8v
TIMEOUT_S=${TIMEOUT_S:-1800}
GPUS=${GPUS:-"0 2 3"}
GPU_MODE=${GPU_MODE:-gpu_bfs_versioned}

PENDING=$OUT/gpu_pending.txt
SUMMARY=$OUT/summary_gpu_rerun.csv
[[ -f $PENDING ]] || { echo "missing $PENDING"; exit 1; }

echo "mode,algo,query,status,wall_ms,positive,negative,update_ms" > "$SUMMARY"

# 把 pending 按 GPU 数轮转切片
declare -a GPU_ARR=($GPUS)
N=${#GPU_ARR[@]}
for g in "${GPU_ARR[@]}"; do : > "$OUT/gpu_pending_g${g}.tsv"; done

i=0
while IFS=/ read -r algo q; do
    g=${GPU_ARR[$((i % N))]}
    printf '%s\t%s\n' "$algo" "$q" >> "$OUT/gpu_pending_g${g}.tsv"
    i=$((i+1))
done < "$PENDING"

echo "=== plan ==="
for g in "${GPU_ARR[@]}"; do
    echo "  GPU $g : $(wc -l <"$OUT/gpu_pending_g${g}.tsv") jobs"
done

run_one() {
    local algo=$1 q=$2 gpu=$3
    local qpath="$QDIR/$q"
    local log="$OUT/gpu/$algo/$q.log"
    mkdir -p "$(dirname "$log")"
    local t0 t1 rc wall
    t0=$(date +%s%3N)
    timeout --foreground "${TIMEOUT_S}s" \
      env CUDA_VISIBLE_DEVICES="$gpu" \
      "$BIN" -q "$qpath" -d "$DATA" -u "$UPD" \
      -a "$algo" -m "$GPU_MODE" --time-limit "$TIMEOUT_S" \
      --report-initial 0 \
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
    {
      flock -x 9
      printf 'gpu,%s,%s,%s,%d,%s,%s,%s\n' \
          "$algo" "$q" "$status" "$wall" "${pos:-}" "${neg:-}" "${upd_ms:-}" >> "$SUMMARY"
    } 9>"$SUMMARY.lock"
    printf '[gpu%s] %-19s %-7s -> %-10s %6dms\n' "$gpu" "$algo" "$q" "$status" "$wall"
}
export -f run_one
export BIN DATA UPD QDIR OUT SUMMARY TIMEOUT_S GPU_MODE

# 每张卡一个串行 worker
PIDS=()
for g in "${GPU_ARR[@]}"; do
    (
      while IFS=$'\t' read -r a q; do
          [[ -z "$a" ]] && continue
          run_one "$a" "$q" "$g"
      done < "$OUT/gpu_pending_g${g}.tsv"
    ) &
    PIDS+=($!)
    echo "  started worker for GPU $g, pid=${PIDS[-1]}"
done

echo "=== waiting on workers: ${PIDS[*]} ==="
wait "${PIDS[@]}"
echo "=== all gpu workers done ==="
awk -F, 'NR>1{c[$4]++} END{for(k in c) print k": "c[k]}' "$SUMMARY" | sort
