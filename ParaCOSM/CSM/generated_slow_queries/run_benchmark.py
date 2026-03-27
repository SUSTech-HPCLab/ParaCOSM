#!/usr/bin/env python3
"""Benchmark generated slow queries with graphflow (1 thread) and parallel_graphflow (6 threads)
on both amazon and lj datasets."""

import csv
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

INITIAL_TIME_RE = re.compile(r"^Initial Matching:\s+([0-9.]+)ms\s*$", re.MULTILINE)
# Parallel_NewSP uses "#Time: XXX ms", main binary uses "Incremental Matching: XXXms"
INCREMENTAL_TIME_RE = re.compile(r"^#Time:\s+([0-9.]+)\s+ms\s*$", re.MULTILINE)
INCREMENTAL_TIME_RE2 = re.compile(r"^Incremental Matching:\s+([0-9.]+)ms\s*$", re.MULTILINE)
INITIAL_MATCHES_RE = re.compile(r"^([0-9]+) matches\.\s*$", re.MULTILINE)
POSITIVE_MATCHES_RE = re.compile(r"^([0-9]+) positive matches\.\s*$", re.MULTILINE)
NEGATIVE_MATCHES_RE = re.compile(r"^([0-9]+) negative matches\.\s*$", re.MULTILINE)
EDGE_UPDATES_RE = re.compile(r"^([0-9]+) edge updates\.\s*$", re.MULTILINE)
TIMEOUT_RE = re.compile(r"Timeout\s+([0-9]+)s")


@dataclass
class BenchResult:
    query_name: str
    vertex_count: int
    edge_count: int
    dataset: str
    algorithm: str
    threads: int
    initial_ms: float | None
    incremental_ms: float | None
    total_ms: float | None
    initial_matches: int | None
    positive_matches: int | None
    negative_matches: int | None
    edge_updates: int | None
    timed_out: bool
    return_code: int


def extract_first(pattern, text):
    m = pattern.search(text)
    if not m:
        return None
    v = m.group(1)
    return float(v) if "." in v else int(v)


def count_shape(path: Path) -> tuple[int, int]:
    vc = ec = 0
    for line in path.read_text().splitlines():
        if line.startswith("v "):
            vc += 1
        elif line.startswith("e "):
            ec += 1
    return vc, ec


def run_bench(
    exe: Path,
    data_graph: Path,
    update_graph: Path,
    query_path: Path,
    algorithm: str,
    threads: int,
    time_limit: int,
    run_timeout: int,
    dataset_name: str,
) -> BenchResult:
    vc, ec = count_shape(query_path)
    cmd = [
        str(exe),
        "-a", algorithm,
        "-d", str(data_graph),
        "-u", str(update_graph),
        "-q", str(query_path),
        "--report-initial", "on",
        "--initial-time-limit", str(time_limit),
        "--time-limit", str(time_limit),
        "-t", str(threads),
    ]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=run_timeout, env=env, check=False
        )
        output = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        output = ""
        rc = -1

    init_ms = extract_first(INITIAL_TIME_RE, output)
    incr_ms = extract_first(INCREMENTAL_TIME_RE, output)
    if incr_ms is None:
        incr_ms = extract_first(INCREMENTAL_TIME_RE2, output)
    total = None
    if init_ms is not None and incr_ms is not None:
        total = init_ms + incr_ms
    elif init_ms is not None:
        total = init_ms
    elif incr_ms is not None:
        total = incr_ms

    return BenchResult(
        query_name=query_path.name,
        vertex_count=vc,
        edge_count=ec,
        dataset=dataset_name,
        algorithm=algorithm,
        threads=threads,
        initial_ms=init_ms,
        incremental_ms=incr_ms,
        total_ms=total,
        initial_matches=extract_first(INITIAL_MATCHES_RE, output),
        positive_matches=extract_first(POSITIVE_MATCHES_RE, output),
        negative_matches=extract_first(NEGATIVE_MATCHES_RE, output),
        edge_updates=extract_first(EDGE_UPDATES_RE, output),
        timed_out=bool(TIMEOUT_RE.search(output)),
        return_code=rc,
    )


def main():
    exe = Path("/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/build/bin/csm")
    data_root = Path("/home/v-haibinlai/haibin/paracosm_data")
    gen_root = Path("/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/generated_slow_queries")
    output_dir = gen_root / "benchmark_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = {
        "amazon": {
            "data": data_root / "6/data_graph/data.graph",
            "update": data_root / "6/data_graph/insertion.graph",
        },
        "lj": {
            "data": data_root / "lj/livejournal/data_graph/data.graph",
            "update": data_root / "lj/livejournal/data_graph/insertion.graph",
        },
    }

    configs = [
        ("graphflow", 1),
        ("parallel_graphflow", 6),
    ]

    # Collect all query paths
    queries_8v = sorted((gen_root / "8v/accepted_queries").glob("Q_gen_*"))
    queries_9v = sorted((gen_root / "9v/accepted_queries").glob("Q_gen_*"))
    all_queries = [("8v", q) for q in queries_8v] + [("9v", q) for q in queries_9v]

    time_limit = 120  # seconds
    run_timeout = 360  # outer timeout

    results: list[BenchResult] = []
    total_runs = len(all_queries) * len(datasets) * len(configs)
    run_idx = 0

    for ds_name, ds_paths in datasets.items():
        for algo, threads in configs:
            for size_label, qpath in all_queries:
                run_idx += 1
                tag = f"[{run_idx}/{total_runs}]"
                print(f"{tag} {ds_name} | {algo}(t={threads}) | {size_label}/{qpath.name}", flush=True)
                t0 = time.time()
                r = run_bench(
                    exe, ds_paths["data"], ds_paths["update"], qpath,
                    algo, threads, time_limit, run_timeout, ds_name,
                )
                elapsed = time.time() - t0
                init_s = f"{r.initial_ms/1000:.2f}s" if r.initial_ms else "NA"
                incr_s = f"{r.incremental_ms/1000:.2f}s" if r.incremental_ms else "NA"
                to = " TIMEOUT" if r.timed_out else ""
                print(f"       => initial={init_s}  incremental={incr_s}{to}  ({elapsed:.1f}s wall)", flush=True)
                results.append(r)

    # Write CSV
    csv_path = output_dir / "benchmark_results.csv"
    fields = [
        "query_name", "vertex_count", "edge_count", "dataset", "algorithm", "threads",
        "initial_ms", "incremental_ms", "total_ms",
        "initial_matches", "positive_matches", "negative_matches", "edge_updates",
        "timed_out", "return_code",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow(r.__dict__)

    print(f"\nResults written to {csv_path}")

    # Print summary table
    print("\n" + "=" * 120)
    print(f"{'Query':<14} {'V':>2} {'E':>2} {'Dataset':<8} {'Algorithm':<20} {'T':>2} {'Initial(s)':>12} {'Incr(s)':>12} {'Total(s)':>12} {'Timeout':>8}")
    print("-" * 120)
    for r in results:
        init = f"{r.initial_ms/1000:.2f}" if r.initial_ms else "NA"
        incr = f"{r.incremental_ms/1000:.2f}" if r.incremental_ms else "NA"
        total = f"{r.total_ms/1000:.2f}" if r.total_ms else "NA"
        to = "YES" if r.timed_out else "no"
        print(f"{r.query_name:<14} {r.vertex_count:>2} {r.edge_count:>2} {r.dataset:<8} {r.algorithm:<20} {r.threads:>2} {init:>12} {incr:>12} {total:>12} {to:>8}")
    print("=" * 120)


if __name__ == "__main__":
    main()
