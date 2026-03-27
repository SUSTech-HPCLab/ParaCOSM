#!/usr/bin/env python3
"""Screen queries to find ones with large 2-level expansion (vv_size) and measure speedup."""

import os
import re
import subprocess
import sys
from pathlib import Path

VV_SIZE_RE = re.compile(r"\[PFM2\] order=\d+ vv_size=(\d+)")
INITIAL_TIME_RE = re.compile(r"^Initial Matching:\s+([0-9.]+)ms\s*$", re.MULTILINE)
INCR_TIME_RE = re.compile(r"^Incremental Matching:\s+([0-9.]+)ms\s*$", re.MULTILINE)
INCR_TIME_RE2 = re.compile(r"^#Time:\s+([0-9.]+)\s+ms\s*$", re.MULTILINE)

def extract(pattern, text):
    m = pattern.search(text)
    return float(m.group(1)) if m else None

def run_query(exe, data, update, query, algo, threads, time_limit, run_timeout):
    cmd = [str(exe), "-a", algo, "-d", str(data), "-u", str(update),
           "-q", str(query), "--report-initial", "on",
           "--initial-time-limit", str(time_limit), "--time-limit", str(time_limit),
           "-t", str(threads)]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=run_timeout, env=env, check=False)
        output = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    except subprocess.TimeoutExpired:
        output = ""
    
    init_ms = extract(INITIAL_TIME_RE, output)
    incr_ms = extract(INCR_TIME_RE, output)
    if incr_ms is None:
        incr_ms = extract(INCR_TIME_RE2, output)
    
    # Extract all vv_size values
    vv_sizes = [int(m) for m in VV_SIZE_RE.findall(output)]
    total_vv = sum(vv_sizes)
    max_vv = max(vv_sizes) if vv_sizes else 0
    
    return init_ms, incr_ms, total_vv, max_vv, len(vv_sizes)


def main():
    exe = Path("/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/build/bin/csm")
    data = Path("/home/v-haibinlai/haibin/paracosm_data/6/data_graph/data.graph")
    update = Path("/home/v-haibinlai/haibin/paracosm_data/6/data_graph/insertion.graph")
    
    # Collect all query sources
    query_dirs = [
        ("8v_seed", Path("/home/v-haibinlai/haibin/paracosm_data/8_self/sparse")),
        ("9v_seed", Path("/home/v-haibinlai/haibin/paracosm_data/9_self/sparse")),
        ("8v_gen", Path("/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/generated_slow_queries/8v/accepted_queries")),
        ("9v_gen", Path("/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/generated_slow_queries/9v/accepted_queries")),
    ]
    
    time_limit = 30  # short limit for screening
    run_timeout = 90
    
    results = []
    
    for label, qdir in query_dirs:
        queries = sorted(qdir.glob("Q_*"))
        for qpath in queries:
            tag = f"{label}/{qpath.name}"
            
            # Run parallel first to get vv_size
            print(f"[SCREEN] {tag} parallel_graphflow(6T)...", end="", flush=True)
            p_init, p_incr, total_vv, max_vv, n_calls = run_query(
                exe, data, update, qpath, "parallel_graphflow", 6, time_limit, run_timeout)
            
            p_total = None
            if p_init is not None and p_incr is not None:
                p_total = p_init + p_incr
            
            print(f" vv={total_vv} max={max_vv} calls={n_calls}", end="", flush=True)
            
            # Only run single-thread if there's meaningful vv_size and parallel didn't fully timeout
            if total_vv > 50 and (p_init is not None or p_incr is not None):
                print(f" → graphflow(1T)...", end="", flush=True)
                s_init, s_incr, _, _, _ = run_query(
                    exe, data, update, qpath, "graphflow", 1, time_limit, run_timeout)
                s_total = None
                if s_init is not None and s_incr is not None:
                    s_total = s_init + s_incr
            else:
                s_init = s_incr = s_total = None
            
            speedup = None
            if s_total and p_total and p_total > 0:
                speedup = s_total / p_total
            
            results.append((tag, total_vv, max_vv, n_calls, s_init, s_incr, s_total, p_init, p_incr, p_total, speedup))
            
            sp_str = f" speedup={speedup:.2f}x" if speedup else ""
            s_str = f" s={s_total/1000:.1f}s" if s_total else ""
            p_str = f" p={p_total/1000:.1f}s" if p_total else ""
            print(f"{s_str}{p_str}{sp_str}")
    
    # Sort by total_vv descending
    results.sort(key=lambda r: r[1], reverse=True)
    
    print("\n" + "=" * 130)
    print(f"{'Query':<25} {'TotalVV':>9} {'MaxVV':>8} {'Calls':>6} {'Single(s)':>10} {'Para(s)':>10} {'Speedup':>8}")
    print("-" * 130)
    for tag, tvv, mvv, nc, si, sn, st, pi, pn, pt, sp in results:
        st_s = f"{st/1000:.2f}" if st else "TO/NA"
        pt_s = f"{pt/1000:.2f}" if pt else "TO/NA"
        sp_s = f"{sp:.2f}x" if sp else "-"
        print(f"{tag:<25} {tvv:>9} {mvv:>8} {nc:>6} {st_s:>10} {pt_s:>10} {sp_s:>8}")
    print("=" * 130)


if __name__ == "__main__":
    main()
