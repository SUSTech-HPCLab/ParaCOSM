#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


INITIAL_MATCHES_RE = re.compile(r"^([0-9]+) matches\.\s*$", re.MULTILINE)
POSITIVE_MATCHES_RE = re.compile(r"^([0-9]+) positive matches\.\s*$", re.MULTILINE)
NEGATIVE_MATCHES_RE = re.compile(r"^([0-9]+) negative matches\.\s*$", re.MULTILINE)
TIME_RE = re.compile(r"^(?:#Time:\s+([0-9.]+)\s+ms|Incremental Matching:\s+([0-9.]+)ms)\s*$", re.MULTILINE)
TIMEOUT_RE = re.compile(r"Timeout\s+([0-9]+)s")
EDGE_UPDATES_RE = re.compile(r"^([0-9]+) edge updates\.\s*$", re.MULTILINE)


@dataclass
class RunStats:
    return_code: int
    initial_matches: int | None
    positive_matches: int | None
    negative_matches: int | None
    incremental_ms: float | None
    edge_updates: int | None
    timed_out: bool
    log_path: str


@dataclass
class CompareResult:
    query_name: str
    query_path: str
    serial_algorithm: str
    parallel_algorithm: str
    serial_initial: int | None
    parallel_initial: int | None
    serial_positive: int | None
    parallel_positive: int | None
    serial_negative: int | None
    parallel_negative: int | None
    serial_time_ms: float | None
    parallel_time_ms: float | None
    serial_edge_updates: int | None
    parallel_edge_updates: int | None
    serial_timed_out: bool
    parallel_timed_out: bool
    serial_return_code: int
    parallel_return_code: int
    match_status: str
    serial_log: str
    parallel_log: str


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    data_root = Path("/home/v-haibinlai/haibin/paracosm_data")

    parser = argparse.ArgumentParser(description="Compare two top-level CSM algorithms on all query files in a directory.")
    parser.add_argument("--query-dir", type=Path, required=True)
    parser.add_argument("--data-graph", type=Path, default=data_root / "6/data_graph/data.graph")
    parser.add_argument("--update-graph", type=Path, default=data_root / "6/data_graph/insertion.graph")
    parser.add_argument("--exe", type=Path, default=repo_root / "build/bin/csm")
    parser.add_argument("--serial-algorithm", required=True)
    parser.add_argument("--parallel-algorithm", required=True)
    parser.add_argument("--serial-threads", type=int, default=1)
    parser.add_argument("--parallel-threads", type=int, default=8)
    parser.add_argument("--serial-auto-tuning", type=int, default=0)
    parser.add_argument("--parallel-auto-tuning", type=int, default=0)
    parser.add_argument("--time-limit", type=int, default=180)
    parser.add_argument("--initial-time-limit", type=int, default=180)
    parser.add_argument("--run-timeout", type=int, default=420)
    parser.add_argument("--report-initial", choices=("on", "off", "1", "0"), default="on")
    parser.add_argument("--max-queries", type=int)
    parser.add_argument("--output-dir", type=Path, default=repo_root / "logs_txt/compare_csm_algorithms_batch")
    return parser.parse_args()


def extract_first_int(pattern: re.Pattern[str], text: str) -> int | None:
    match = pattern.search(text)
    return int(match.group(1)) if match else None


def extract_time_ms(text: str) -> float | None:
    match = TIME_RE.search(text)
    if not match:
        return None
    for group in match.groups():
        if group is not None:
            return float(group)
    return None


def run_one(
    exe: Path,
    algorithm: str,
    threads: int,
    auto_tuning: int,
    query_path: Path,
    args: argparse.Namespace,
    log_path: Path,
) -> RunStats:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(exe.resolve()),
        "-a",
        algorithm,
        "-d",
        str(args.data_graph.resolve()),
        "-u",
        str(args.update_graph.resolve()),
        "-q",
        str(query_path.resolve()),
        "--report-initial",
        args.report_initial,
        "--initial-time-limit",
        str(args.initial_time_limit),
        "--time-limit",
        str(args.time_limit),
        "-t",
        str(threads),
        "--auto-tuning",
        str(auto_tuning),
    ]

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=args.run_timeout,
        env=env,
        cwd=str(exe.resolve().parent),
        check=False,
    )
    combined_output = completed.stdout if not completed.stderr else f"{completed.stdout}\n[stderr]\n{completed.stderr}"
    log_path.write_text(combined_output, encoding="utf-8")

    return RunStats(
        return_code=completed.returncode,
        initial_matches=extract_first_int(INITIAL_MATCHES_RE, combined_output),
        positive_matches=extract_first_int(POSITIVE_MATCHES_RE, combined_output),
        negative_matches=extract_first_int(NEGATIVE_MATCHES_RE, combined_output),
        incremental_ms=extract_time_ms(combined_output),
        edge_updates=extract_first_int(EDGE_UPDATES_RE, combined_output),
        timed_out=bool(TIMEOUT_RE.search(combined_output)),
        log_path=str(log_path),
    )


def compare_stats(args: argparse.Namespace, query_name: str, query_path: Path, serial: RunStats, parallel: RunStats) -> CompareResult:
    same_counts = (
        serial.initial_matches == parallel.initial_matches
        and serial.positive_matches == parallel.positive_matches
        and serial.negative_matches == parallel.negative_matches
    )
    same_updates = serial.edge_updates == parallel.edge_updates

    if same_counts and same_updates:
        status = "pass"
    elif serial.timed_out or parallel.timed_out:
        status = "timeout-mismatch-progress" if not same_updates else "timeout-mismatch"
    elif serial.positive_matches == parallel.positive_matches and serial.negative_matches == parallel.negative_matches:
        status = "initial-mismatch"
    else:
        status = "fail"

    return CompareResult(
        query_name=query_name,
        query_path=str(query_path),
        serial_algorithm=args.serial_algorithm,
        parallel_algorithm=args.parallel_algorithm,
        serial_initial=serial.initial_matches,
        parallel_initial=parallel.initial_matches,
        serial_positive=serial.positive_matches,
        parallel_positive=parallel.positive_matches,
        serial_negative=serial.negative_matches,
        parallel_negative=parallel.negative_matches,
        serial_time_ms=serial.incremental_ms,
        parallel_time_ms=parallel.incremental_ms,
        serial_edge_updates=serial.edge_updates,
        parallel_edge_updates=parallel.edge_updates,
        serial_timed_out=serial.timed_out,
        parallel_timed_out=parallel.timed_out,
        serial_return_code=serial.return_code,
        parallel_return_code=parallel.return_code,
        match_status=status,
        serial_log=serial.log_path,
        parallel_log=parallel.log_path,
    )


def main() -> int:
    args = parse_args()
    args.output_dir = args.output_dir.resolve()
    if not args.query_dir.is_dir():
        print(f"Query directory not found: {args.query_dir}", file=sys.stderr)
        return 2
    if not args.exe.is_file():
        print(f"Executable not found: {args.exe}", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = args.output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    queries = sorted(path for path in args.query_dir.iterdir() if path.is_file())
    if args.max_queries is not None:
        queries = queries[: args.max_queries]

    results: list[CompareResult] = []
    for query_path in queries:
        query_name = query_path.name
        serial_log = logs_dir / f"{query_name}.serial.log"
        parallel_log = logs_dir / f"{query_name}.parallel.log"
        print(f"[RUN ] {query_name} {args.serial_algorithm}")
        serial_stats = run_one(args.exe, args.serial_algorithm, args.serial_threads, args.serial_auto_tuning, query_path, args, serial_log)
        print(f"[RUN ] {query_name} {args.parallel_algorithm}")
        parallel_stats = run_one(args.exe, args.parallel_algorithm, args.parallel_threads, args.parallel_auto_tuning, query_path, args, parallel_log)
        result = compare_stats(args, query_name, query_path, serial_stats, parallel_stats)
        results.append(result)
        print(
            f"[{result.match_status.upper()}] {query_name} "
            f"serial={result.serial_time_ms} parallel={result.parallel_time_ms} "
            f"updates={result.serial_edge_updates}/{result.parallel_edge_updates}"
        )

    summary_path = args.output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(results[0]).keys()) if results else list(CompareResult.__annotations__.keys()))
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))

    counts: dict[str, int] = {}
    for result in results:
        counts[result.match_status] = counts.get(result.match_status, 0) + 1

    print()
    print(f"Summary: total={len(results)} statuses={counts}")
    print(f"Summary CSV: {summary_path}")
    return 0 if counts.get("fail", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())