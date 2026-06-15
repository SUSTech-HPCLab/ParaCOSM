#!/usr/bin/env bash
# ============================================================================
# bench_lj_8v.sh
#
# 在 LiveJournal 数据集上对 8v 查询进行 CPU vs GPU 基准测试
#
# - 数据集:  /home/superbench/haibin/paracosm/LJ
# - 查询集:  LJ/8v/Q_*  (取前 NUM_QUERIES 个)
# - 算法:    parallel_graphflow / parallel_symbi / parallel_turboflux
#            parallel_calig    / parallel_newsp                (5 个)
# - CPU:     -m single  -t 1   (单线程)，并行 20 个进程
# - GPU:     -m gpu_bfs_versioned    (单 GPU 进程串行)
# - timeout: 每个 query 1800s
#
# 输出:
#   $OUT_DIR/cpu/<algo>/<query>.log     (单条日志)
#   $OUT_DIR/gpu/<algo>/<query>.log
#   $OUT_DIR/summary.csv                (汇总: algo,query,mode,status,time_ms,...)
#
# 用法:
#   bash bench_lj_8v.sh [DRY]
#     DRY=1 只打印命令，不执行
# ============================================================================
set -u
shopt -s nullglob

# ------------------------------ 配置 ----------------------------------------
ROOT=/home/superbench/haibin/paracosm
CSM=$ROOT/ParaCOSM/ParaCOSM/CSM
BIN=${CSM_BIN:-$CSM/build/csm}

DATA=$ROOT/LJ/data_graph/data.graph
UPD=$ROOT/LJ/data_graph/insertion.graph
QDIR=$ROOT/LJ/8v

NUM_QUERIES=${NUM_QUERIES:-50}
TIMEOUT_S=${TIMEOUT_S:-1800}
CPU_PARALLEL=${CPU_PARALLEL:-50}
GPU_DEVICE=${GPU_DEVICE:-0}

ALGOS=(parallel_graphflow parallel_symbi parallel_turboflux parallel_calig parallel_newsp)
CPU_MODE=single                # CPU 单线程
GPU_MODE=gpu_bfs_versioned     # GPU BFS versioned (inner-update 语义)

STAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR=${OUT_DIR:-$ROOT/bench_lj_8v_$STAMP}
SUMMARY=$OUT_DIR/summary.csv

DRY=${DRY:-0}

# ------------------------------ 校验 ----------------------------------------
err() { echo "[ERR] $*" >&2; exit 1; }
[[ -x "$BIN" ]]      || err "binary not found / not executable: $BIN  (先 build 出 csm)"
[[ -f "$DATA" ]]     || err "data graph not found: $DATA"
[[ -f "$UPD"  ]]     || err "update stream not found: $UPD"
[[ -d "$QDIR" ]]     || err "query dir not found: $QDIR"
command -v timeout >/dev/null || err "'timeout' (coreutils) required"

# 选 query: 字典序前 NUM_QUERIES 个 (Q_1, Q_10, Q_100, Q_11, ...).
# 若想按数字顺序 Q_1..Q_20，可改成: ls Q_* | sort -t_ -k2,2n
mapfile -t QUERIES < <(ls "$QDIR" | sort -t_ -k2,2n | head -n "$NUM_QUERIES")
[[ ${#QUERIES[@]} -gt 0 ]] || err "no queries selected from $QDIR"

mkdir -p "$OUT_DIR"
for a in "${ALGOS[@]}"; do
    mkdir -p "$OUT_DIR/cpu/$a" "$OUT_DIR/gpu/$a"
done

cat <<EOF | tee "$OUT_DIR/run_info.txt"
=== bench_lj_8v.sh ===
binary       : $BIN
data         : $DATA
update       : $UPD
qdir         : $QDIR
queries (${#QUERIES[@]}): ${QUERIES[*]}
algos        : ${ALGOS[*]}
cpu mode     : $CPU_MODE  (1 thread, ${CPU_PARALLEL} concurrent procs)
gpu mode     : $GPU_MODE  (1 process at a time, GPU $GPU_DEVICE)
timeout      : ${TIMEOUT_S}s
out          : $OUT_DIR
EOF

# ------------------------------ Warmup -------------------------------------
# 把数据图 + 插入流先读一遍喂进 OS page cache, 后续每次 LoadFromFile 走内存
# 若装了 vmtouch 还可用 'vmtouch -t' 强制 touch+锁定
echo
echo "================== Warmup (preload data into page cache) =================="
warm_start=$(date +%s)
if command -v vmtouch >/dev/null 2>&1; then
    vmtouch -t "$DATA" "$UPD" 2>&1 | tail -5
    vmtouch    "$DATA" "$UPD" 2>&1 | tail -5
else
    # 并行读两个文件到 /dev/null
    { time { cat "$DATA" >/dev/null & cat "$UPD" >/dev/null & wait; }; } 2>&1 | sed 's/^/  /'
fi
warm_end=$(date +%s)
free -h | sed 's/^/  /'
echo "  warmup took $((warm_end - warm_start))s"

# CSV 表头
echo "mode,algo,query,status,wall_ms,positive,negative,update_ms" > "$SUMMARY"

# ------------------------------ 单次运行 ------------------------------------
# run_one MODE_TAG ALGO QUERY UPDATE_MODE [GPU_VISIBLE]
run_one() {
    local mode_tag=$1 algo=$2 q=$3 upd_mode=$4 gpu_vis=${5:-}
    local qpath="$QDIR/$q"
    local log="$OUT_DIR/$mode_tag/$algo/$q.log"
    local extra_env=()

    if [[ "$mode_tag" == "cpu" ]]; then
        extra_env=(env CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1)
        local extra_args=(-t 1 --auto-tuning 0 --report-initial 0)
    else
        extra_env=(env CUDA_VISIBLE_DEVICES="$gpu_vis")
        local extra_args=(--report-initial 0)
    fi

    local cmd=(timeout --foreground "${TIMEOUT_S}s" "${extra_env[@]}" \
        "$BIN" -q "$qpath" -d "$DATA" -u "$UPD" \
        -a "$algo" -m "$upd_mode" \
        --time-limit "$TIMEOUT_S" \
        "${extra_args[@]}")

    if [[ "$DRY" == "1" ]]; then
        echo "[DRY] $mode_tag $algo $q :: ${cmd[*]} > $log"
        return 0
    fi

    local t0 t1 rc
    t0=$(date +%s%3N)
    "${cmd[@]}" >"$log" 2>&1
    rc=$?
    t1=$(date +%s%3N)
    local wall=$((t1 - t0))

    local status
    if   [[ $rc -eq 0   ]]; then status=OK
    elif [[ $rc -eq 124 ]]; then status=TIMEOUT
    elif [[ $rc -eq 137 ]]; then status=OOM_OR_KILL
    elif [[ $rc -eq 139 ]]; then status=SEGV
    else                          status=FAIL_$rc
    fi

    local pos neg upd_ms
    pos=$(grep -aE 'positive_num_results|positive results' "$log" | tail -1 \
            | grep -oE '[0-9]+' | tail -1)
    neg=$(grep -aE 'negative_num_results|negative results' "$log" | tail -1 \
            | grep -oE '[0-9]+' | tail -1)
    upd_ms=$(grep -aE 'Update Time|Incremental .* Time|Total .* time' "$log" \
            | tail -1 | grep -oE '[0-9]+(\.[0-9]+)?' | tail -1)

    printf '%s,%s,%s,%s,%d,%s,%s,%s\n' \
        "$mode_tag" "$algo" "$q" "$status" "$wall" \
        "${pos:-}" "${neg:-}" "${upd_ms:-}" >> "$SUMMARY"

    printf '[%s] %-19s %-7s %-22s -> %-10s %6dms\n' \
        "$mode_tag" "$algo" "$q" "$status" "$wall"
}

export -f run_one
export BIN DATA UPD QDIR OUT_DIR SUMMARY TIMEOUT_S DRY

# ------------------------------ CPU + GPU 阶段并发 --------------------------
echo
echo "================== CPU PHASE (${CPU_PARALLEL}-way) + GPU PHASE (serial dev $GPU_DEVICE) concurrently =================="
cpu_jobs=$(mktemp)
trap 'rm -f "$cpu_jobs"' EXIT
for a in "${ALGOS[@]}"; do
    for q in "${QUERIES[@]}"; do
        printf '%s\t%s\n' "$a" "$q"
    done
done > "$cpu_jobs"

# CPU 后台跑 (xargs 调用 run_one cpu <algo> <query> single)
( xargs -a "$cpu_jobs" -P "$CPU_PARALLEL" -n 1 -I{} \
      bash -c 'IFS=$'"'"'\t'"'"' read -r a q <<<"{}"; run_one cpu "$a" "$q" '"$CPU_MODE" ) &
CPU_PID=$!
echo "  CPU phase started in background, pid=$CPU_PID"

# GPU 串行前台跑
for a in "${ALGOS[@]}"; do
    for q in "${QUERIES[@]}"; do
        run_one gpu "$a" "$q" "$GPU_MODE" "$GPU_DEVICE"
    done
done
echo "  GPU phase done, waiting for CPU phase (pid=$CPU_PID) ..."
wait "$CPU_PID" || true
echo "  CPU phase done"

# ------------------------------ 汇总 ----------------------------------------
echo
echo "================== Summary =================="
echo "summary csv: $SUMMARY"
awk -F, 'NR>1 {n[$1"|"$2]++; t[$1"|"$2]+=$5; ok[$1"|"$2]+=($4=="OK")}
         END { printf "%-4s %-22s %5s %5s %12s\n","mode","algo","n","ok","avg_ms";
               for (k in n) { split(k,a,"|");
                 printf "%-4s %-22s %5d %5d %12.0f\n",a[1],a[2],n[k],ok[k],t[k]/n[k]; } }' \
    "$SUMMARY" | sort
