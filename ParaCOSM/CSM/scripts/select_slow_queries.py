#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


INITIAL_TIME_RE = re.compile(r"^Initial Matching:\s+([0-9.]+)ms\s*$", re.MULTILINE)
INCREMENTAL_TIME_RE = re.compile(r"^#Time:\s+([0-9.]+)\s+ms\s*$", re.MULTILINE)
INITIAL_MATCHES_RE = re.compile(r"^([0-9]+) matches\.\s*$", re.MULTILINE)
POSITIVE_MATCHES_RE = re.compile(r"^([0-9]+) positive matches\.\s*$", re.MULTILINE)
NEGATIVE_MATCHES_RE = re.compile(r"^([0-9]+) negative matches\.\s*$", re.MULTILINE)
EDGE_UPDATES_RE = re.compile(r"^([0-9]+) edge updates\.\s*$", re.MULTILINE)
TIMEOUT_RE = re.compile(r"Timeout\s+([0-9]+)s")


@dataclass
class QueryRunResult:
    query_name: str
    query_path: str
    vertex_count: int
    edge_count: int
    return_code: int
    selected_metric: float | None
    initial_ms: float | None
    incremental_ms: float | None
    initial_matches: int | None
    positive_matches: int | None
    negative_matches: int | None
    edge_updates: int | None
    timed_out: bool
    qualifies: bool
    command: list[str]
    log_path: str | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CSM on query graphs, keep only the cases whose runtime exceeds a threshold, "
            "and optionally copy those query files to an output directory."
        )
    )
    repo_root = Path(__file__).resolve().parents[1]
    data_root = Path("/home/v-haibinlai/haibin/paracosm_data")

    parser.add_argument(
        "--exe",
        type=Path,
        default=repo_root / "matching_executor/Parallel_NewSP/build/csm",
        help="Path to the executable that evaluates a query graph.",
    )
    parser.add_argument(
        "--data-graph",
        type=Path,
        default=data_root / "6/data_graph/data.graph",
        help="Path to the initial data graph.",
    )
    parser.add_argument(
        "--update-graph",
        type=Path,
        default=data_root / "6/data_graph/insertion.graph",
        help="Path to the update stream.",
    )
    parser.add_argument(
        "--query-dir",
        type=Path,
        action="append",
        required=True,
        help="Directory that contains query files such as Q_21. Repeat this option for multiple directories.",
    )
    parser.add_argument(
        "--query-ids",
        type=int,
        nargs="*",
        help="Optional query ids. If omitted, the script scans all Q_* files in the directory.",
    )
    parser.add_argument(
        "--vertex-count",
        type=int,
        default=8,
        help="Only keep query graphs with exactly this number of vertices.",
    )
    parser.add_argument(
        "--min-seconds",
        type=float,
        default=30.0,
        help="Threshold in seconds. Queries whose selected metric is at least this value will be selected.",
    )
    parser.add_argument(
        "--max-seconds",
        type=float,
        help="Optional upper bound in seconds for the selected metric.",
    )
    parser.add_argument(
        "--metric",
        choices=("initial", "incremental", "total", "max"),
        default="incremental",
        help="Runtime metric used for threshold filtering.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="OMP_NUM_THREADS used when running the executable.",
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=60,
        help="Value passed to --time-limit of the executable, in seconds.",
    )
    parser.add_argument(
        "--initial-time-limit",
        type=int,
        default=60,
        help="Value passed to --initial-time-limit of the executable, in seconds.",
    )
    parser.add_argument(
        "--run-timeout",
        type=int,
        default=180,
        help="Outer subprocess timeout in seconds.",
    )
    parser.add_argument(
        "--report-initial",
        choices=("on", "off"),
        default="on",
        help="Value passed to --report-initial.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "logs_txt/slow_queries_8v",
        help="Directory used for logs, summary files, and copied query graphs.",
    )
    parser.add_argument(
        "--copy-queries",
        action="store_true",
        help="Copy selected query files into output-dir/queries.",
    )
    parser.add_argument(
        "--keep-all-logs",
        action="store_true",
        help="Keep logs for all runs. By default only selected queries keep their logs.",
    )
    parser.add_argument(
        "--require-no-timeout",
        action="store_true",
        help="Only select query graphs whose run does not report a timeout.",
    )
    return parser.parse_args()


def count_query_shape(query_path: Path) -> tuple[int, int]:
    vertex_count = 0
    edge_count = 0
    with query_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("v "):
                vertex_count += 1
            elif line.startswith("e "):
                edge_count += 1
    return vertex_count, edge_count


def iter_queries(query_dirs: Iterable[Path], query_ids: list[int] | None) -> list[Path]:
    query_paths: list[Path] = []
    for query_dir in query_dirs:
        if query_ids:
            for query_id in query_ids:
                query_path = query_dir / f"Q_{query_id}"
                if query_path.is_file():
                    query_paths.append(query_path)
        else:
            query_paths.extend(sorted(path for path in query_dir.glob("Q_*") if path.is_file()))
    return sorted(query_paths, key=lambda path: path.name)


def extract_first_number(pattern: re.Pattern[str], text: str) -> float | int | None:
    match = pattern.search(text)
    if not match:
        return None
    value = match.group(1)
    if "." in value:
        return float(value)
    return int(value)


def select_metric(metric: str, initial_ms: float | None, incremental_ms: float | None) -> float | None:
    if metric == "initial":
        return initial_ms
    if metric == "incremental":
        return incremental_ms
    if metric == "total":
        if initial_ms is None and incremental_ms is None:
            return None
        return (initial_ms or 0.0) + (incremental_ms or 0.0)
    values = [value for value in (initial_ms, incremental_ms) if value is not None]
    return max(values) if values else None


def run_query(args: argparse.Namespace, query_path: Path, logs_dir: Path) -> QueryRunResult:
    exe_path = args.exe.resolve()
    data_graph_path = args.data_graph.resolve()
    update_graph_path = args.update_graph.resolve()
    query_graph_path = query_path.resolve()

    command = [
        str(exe_path),
        "-d",
        str(data_graph_path),
        "-u",
        str(update_graph_path),
        "-q",
        str(query_graph_path),
        "--report-initial",
        args.report_initial,
        "--initial-time-limit",
        str(args.initial_time_limit),
        "--time-limit",
        str(args.time_limit),
    ]

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.threads)

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=args.run_timeout,
        env=env,
        cwd=str(exe_path.parent),
        check=False,
    )

    stdout = completed.stdout
    stderr = completed.stderr
    combined_output = stdout if not stderr else f"{stdout}\n[stderr]\n{stderr}"

    initial_ms = extract_first_number(INITIAL_TIME_RE, combined_output)
    incremental_ms = extract_first_number(INCREMENTAL_TIME_RE, combined_output)
    selected_ms = select_metric(args.metric, initial_ms, incremental_ms)
    timed_out = bool(TIMEOUT_RE.search(combined_output))

    vertex_count, edge_count = count_query_shape(query_path)
    qualifies = selected_ms is not None and selected_ms >= args.min_seconds * 1000.0
    if qualifies and args.max_seconds is not None:
        qualifies = selected_ms <= args.max_seconds * 1000.0
    if qualifies and args.require_no_timeout:
        qualifies = not timed_out

    log_path = logs_dir / f"{query_path.name}.log"
    if args.keep_all_logs or qualifies:
        log_path.write_text(combined_output, encoding="utf-8")
        log_path_str: str | None = str(log_path)
    else:
        log_path_str = None

    return QueryRunResult(
        query_name=query_path.name,
        query_path=str(query_path),
        vertex_count=vertex_count,
        edge_count=edge_count,
        return_code=completed.returncode,
        selected_metric=selected_ms,
        initial_ms=initial_ms,
        incremental_ms=incremental_ms,
        initial_matches=extract_first_number(INITIAL_MATCHES_RE, combined_output),
        positive_matches=extract_first_number(POSITIVE_MATCHES_RE, combined_output),
        negative_matches=extract_first_number(NEGATIVE_MATCHES_RE, combined_output),
        edge_updates=extract_first_number(EDGE_UPDATES_RE, combined_output),
        timed_out=timed_out,
        qualifies=qualifies,
        command=command,
        log_path=log_path_str,
    )


def write_csv(results: list[QueryRunResult], path: Path) -> None:
    fieldnames = list(asdict(results[0]).keys()) if results else [
        "query_name",
        "query_path",
        "vertex_count",
        "edge_count",
        "return_code",
        "selected_metric",
        "initial_ms",
        "incremental_ms",
        "initial_matches",
        "positive_matches",
        "negative_matches",
        "edge_updates",
        "timed_out",
        "qualifies",
        "command",
        "log_path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = asdict(result)
            row["command"] = " ".join(result.command)
            writer.writerow(row)


def write_json(results: list[QueryRunResult], path: Path) -> None:
    payload = []
    for result in results:
        item = asdict(result)
        item["command"] = result.command
        payload.append(item)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def print_summary(results: list[QueryRunResult], min_seconds: float, metric: str) -> None:
    selected = [result for result in results if result.qualifies]
    print(f"Scanned {len(results)} query graphs.")
    print(f"Selected {len(selected)} queries with {metric} runtime >= {min_seconds:.1f}s.")
    if not selected:
        return

    print()
    print("Selected queries:")
    for result in sorted(selected, key=lambda item: item.selected_metric or 0.0, reverse=True):
        metric_seconds = (result.selected_metric or 0.0) / 1000.0
        initial_seconds = (result.initial_ms or 0.0) / 1000.0 if result.initial_ms is not None else None
        incremental_seconds = (result.incremental_ms or 0.0) / 1000.0 if result.incremental_ms is not None else None
        print(
            f"- {result.query_name}: selected={metric_seconds:.3f}s "
            f"initial={initial_seconds if initial_seconds is not None else 'NA'} "
            f"incremental={incremental_seconds if incremental_seconds is not None else 'NA'} "
            f"timeout={'yes' if result.timed_out else 'no'}"
        )


def copy_selected_queries(results: list[QueryRunResult], output_dir: Path) -> None:
    queries_dir = output_dir / "queries"
    queries_dir.mkdir(parents=True, exist_ok=True)
    for result in results:
        if result.qualifies:
            shutil.copy2(result.query_path, queries_dir / result.query_name)


def main() -> int:
    args = parse_args()

    if not args.exe.is_file():
        print(f"Executable not found: {args.exe}", file=sys.stderr)
        return 2
    if not args.data_graph.is_file():
        print(f"Data graph not found: {args.data_graph}", file=sys.stderr)
        return 2
    if not args.update_graph.is_file():
        print(f"Update graph not found: {args.update_graph}", file=sys.stderr)
        return 2

    for query_dir in args.query_dir:
        if not query_dir.is_dir():
            print(f"Query directory not found: {query_dir}", file=sys.stderr)
            return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = args.output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    results: list[QueryRunResult] = []
    query_paths = iter_queries(args.query_dir, args.query_ids)
    for query_path in query_paths:
        vertex_count, _ = count_query_shape(query_path)
        if vertex_count != args.vertex_count:
            continue
        print(f"[RUN ] {query_path}")
        try:
            result = run_query(args, query_path, logs_dir)
        except subprocess.TimeoutExpired:
            print(f"[FAIL] {query_path} exceeded outer timeout of {args.run_timeout}s", file=sys.stderr)
            continue
        results.append(result)
        selected_ms = "NA" if result.selected_metric is None else f"{result.selected_metric / 1000.0:.3f}s"
        status = "SELECT" if result.qualifies else "skip"
        print(f"[{status}] {query_path.name} metric={selected_ms} rc={result.return_code}")

    if args.copy_queries:
        copy_selected_queries(results, args.output_dir)

    write_csv(results, args.output_dir / "summary.csv")
    write_json(results, args.output_dir / "summary.json")
    print_summary(results, args.min_seconds, args.metric)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())