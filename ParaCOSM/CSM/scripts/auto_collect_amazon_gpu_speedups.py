#!/usr/bin/env python3
"""Continuously collect distinct, full-stream Amazon GPU-speedup queries."""

from __future__ import annotations

import argparse
import csv
import fcntl
import math
import os
import random
import signal
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

import networkx as nx


FIELDS = [
    "query", "vertices", "edges", "positive_matches", "cpu_16core_ms",
    "gpu_v100_ms", "speedup", "counts_equal", "source",
]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    workspace = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vertex-counts", type=int, nargs="+",
                        default=[8, 9, 10, 11, 12])
    parser.add_argument("--per-size-count", type=int, default=10)
    parser.add_argument("--candidate-count", type=int, default=10)
    parser.add_argument("--per-round-target", type=int, default=3)
    parser.add_argument("--target-speedup", type=float, default=20.0)
    parser.add_argument(
        "--sampling-strategy", choices=("random", "high-volume-pattern"),
        default="random",
    )
    parser.add_argument("--cpu-threads", type=int, default=16)
    parser.add_argument("--cpu-numa-node", type=int, default=0)
    parser.add_argument("--gpu-id", type=int, default=0)
    parser.add_argument("--gpu-timeout", type=float, default=180.0)
    parser.add_argument("--cpu-timeout", type=float, default=600.0)
    parser.add_argument("--base-seed", type=int, default=20260912)
    parser.add_argument("--max-rounds", type=int, default=0,
                        help="Zero means continue until target-count is reached.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Prepare state and print the next round without executing it.")
    parser.add_argument(
        "--bootstrap-dir", type=Path,
        default=root / "generated_gpu_speedup_queries/amazon_verified_unique",
    )
    parser.add_argument("--seed-root", type=Path,
                        default=workspace / "datasets/amazon/AZ")
    parser.add_argument(
        "--output-dir", type=Path,
        default=root / "generated_gpu_speedup_queries/amazon_auto_50",
    )
    parser.add_argument("--log-file", type=Path,
                        help="Redirect this process and all child output here.")
    return parser.parse_args()


def redirect_output(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    stream = path.open("a", buffering=1, encoding="utf-8")
    os.dup2(stream.fileno(), 1)
    os.dup2(stream.fileno(), 2)
    sys.stdout = os.fdopen(1, "w", buffering=1, closefd=False)
    sys.stderr = os.fdopen(2, "w", buffering=1, closefd=False)


def load_graph(path: Path) -> nx.Graph:
    graph = nx.Graph()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v":
                graph.add_node(int(parts[1]), label=int(parts[2]))
            elif parts[0] == "e":
                graph.add_edge(int(parts[1]), int(parts[2]), label=int(parts[3]))
    return graph


NODE_MATCH = nx.algorithms.isomorphism.categorical_node_match("label", None)
EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("label", None)


def is_duplicate(graph: nx.Graph, known: list[nx.Graph]) -> bool:
    return any(
        graph.number_of_nodes() == prior.number_of_nodes()
        and graph.number_of_edges() == prior.number_of_edges()
        and nx.is_isomorphic(graph, prior, node_match=NODE_MATCH,
                             edge_match=EDGE_MATCH)
        for prior in known
    )


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in FIELDS} for row in rows)
    temporary.replace(path)


def add_query(source: Path, result: dict[str, str], rows: list[dict[str, str]],
              known: list[nx.Graph], query_dir: Path, origin: str) -> bool:
    graph = load_graph(source)
    if is_duplicate(graph, known):
        print(f"[AUTO] isomorphic duplicate skipped: {source}", flush=True)
        return False
    index = len(rows) + 1
    name = f"A{index:03d}_{graph.number_of_nodes()}v_{graph.number_of_edges()}e.graph"
    destination = query_dir / name
    shutil.copy2(source, destination)
    rows.append({
        "query": name,
        "vertices": str(graph.number_of_nodes()),
        "edges": str(graph.number_of_edges()),
        "positive_matches": result.get("gpu_positive_matches", result.get("positive_matches", "")),
        "cpu_16core_ms": result.get("cpu_incremental_ms", result.get("cpu_16core_ms", "")),
        "gpu_v100_ms": result.get("gpu_incremental_ms", result.get("gpu_v100_ms", "")),
        "speedup": result.get("speedup", ""),
        "counts_equal": "True",
        "source": origin,
    })
    known.append(graph)
    print(f"[AUTO] accepted unique {name}: {rows[-1]['speedup']}x", flush=True)
    return True


def bootstrap(args: argparse.Namespace, rows: list[dict[str, str]],
              known: list[nx.Graph], query_dir: Path) -> None:
    if rows:
        for row in rows:
            known.append(load_graph(query_dir / row["query"]))
        return
    source_manifest = args.bootstrap_dir / "manifest.csv"
    for row in read_manifest(source_manifest):
        source = args.bootstrap_dir / "queries" / row["query"]
        add_query(source, row, rows, known, query_dir, "bootstrap")


def mutation_steps(round_index: int) -> tuple[int, int]:
    return ((1, 1), (1, 2), (2, 3), (1, 3))[(round_index - 1) % 4]


def size_counts(rows: list[dict[str, str]], sizes: list[int]) -> dict[int, int]:
    return {size: sum(int(row["vertices"]) == size for row in rows)
            for size in sizes}


def high_volume_seed_score(graph: nx.Graph, row: dict[str, str]) -> float:
    """Rank verified seeds using the structural signals from the 100-query study."""
    degrees = [degree for _, degree in graph.degree()]
    leaves = sum(degree == 1 for degree in degrees)
    bridges = sum(1 for _ in nx.bridges(graph))
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + 1
    triangles = sum(nx.triangles(graph).values()) // 3
    diameter = nx.diameter(graph)
    matches = max(1, int(row.get("positive_matches") or 1))
    speedup = float(row.get("speedup") or 0.0)
    leaf_score = 3.0 if 2 <= leaves <= 3 else -2.0 * abs(leaves - 2)
    return (
        3.0 * math.log10(matches)
        + min(speedup, 80.0) / 20.0
        + leaf_score + bridges + 0.5 * diameter
        - 2.0 * cycle_rank - triangles
    )


def make_12v_seed(args: argparse.Namespace, output: Path,
                  rows: list[dict[str, str]], query_dir: Path,
                  round_index: int) -> Path:
    seed_dir = output / "seeds"
    seed_dir.mkdir(parents=True, exist_ok=True)
    destination = seed_dir / f"Q_seed_12v_round_{round_index:04d}"
    if destination.is_file():
        return destination

    verified_11v = [query_dir / row["query"] for row in rows
                    if int(row["vertices"]) == 11]
    if verified_11v:
        source = verified_11v[(round_index - 1) % len(verified_11v)]
    else:
        original_11v = sorted((args.seed_root / "11v").glob("Q_*"))
        if not original_11v:
            raise RuntimeError("cannot synthesize 12v seed without an 11v query")
        source = original_11v[(round_index - 1) % len(original_11v)]
    graph = load_graph(source)
    new_vertex = graph.number_of_nodes()
    rng = random.Random(args.base_seed + round_index * 7919)
    # Amazon labels are nearly uniform. Add a two- or three-edge constrained
    # vertex to a verified high-speed 11v graph; a leaf made the intermediate
    # frontier explode on every tested 12v candidate.
    graph.add_node(new_vertex, label=round_index % 6)
    attachment_count = 2 + (round_index % 2)
    for target in rng.sample(list(graph.nodes)[:-1], attachment_count):
        graph.add_edge(target, new_vertex, label=0)
    with destination.open("w", encoding="utf-8") as handle:
        for vertex in sorted(graph.nodes):
            handle.write(f"v {vertex} {graph.nodes[vertex]['label']}\n")
        for u, v, data in sorted(graph.edges(data=True)):
            handle.write(f"e {min(u, v)} {max(u, v)} {data['label']}\n")
    return destination


def choose_focus(args: argparse.Namespace, rows: list[dict[str, str]],
                 query_dir: Path, output: Path, size: int,
                 round_index: int) -> Path:
    size_rows = [row for row in rows if int(row["vertices"]) == size]
    if args.sampling_strategy == "high-volume-pattern" and size_rows:
        size_rows.sort(
            key=lambda row: high_volume_seed_score(
                load_graph(query_dir / row["query"]), row
            ),
            reverse=True,
        )
        # Rotate among several strong families rather than cloning one topology.
        size_rows = size_rows[:min(8, len(size_rows))]
    verified = [query_dir / row["query"] for row in size_rows]
    if verified:
        return verified[(round_index - 1) % len(verified)]
    seed_dir = args.seed_root / f"{size}v"
    seeds = sorted(seed_dir.glob("Q_*")) if seed_dir.is_dir() else []
    if seeds:
        return seeds[(round_index - 1) % len(seeds)]
    if size == 12:
        return make_12v_seed(args, output, rows, query_dir, round_index)
    raise RuntimeError(f"no seed query available for {size} vertices")


def main() -> int:
    args = parse_args()
    if args.log_file:
        redirect_output(args.log_file.resolve())

    root = Path(__file__).resolve().parents[1]
    collector = root / "scripts/collect_gpu_speedup_queries.py"
    output = args.output_dir.resolve()
    query_dir = output / "queries"
    rounds_dir = output / "rounds"
    query_dir.mkdir(parents=True, exist_ok=True)
    rounds_dir.mkdir(parents=True, exist_ok=True)

    lock_handle = (output / "collector.lock").open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"collector already running: {output}", file=sys.stderr)
        return 2
    pid_path = output / "collector.pid"
    pid_path.write_text(f"{os.getpid()}\n", encoding="utf-8")

    rows = read_manifest(output / "manifest.csv")
    known: list[nx.Graph] = []
    bootstrap(args, rows, known, query_dir)
    write_manifest(output / "manifest.csv", rows)
    target_total = len(args.vertex_counts) * args.per_size_count
    counts = size_counts(rows, args.vertex_counts)
    print(f"[AUTO] start host={socket.gethostname()} pid={os.getpid()} "
          f"verified={len(rows)}/{target_total} by_size={counts}", flush=True)

    existing = sorted(rounds_dir.glob("round_*"))
    round_index = len(existing) + 1
    completed_this_run = 0
    active_child: subprocess.Popen[bytes] | None = None

    def stop_child(signum: int, _frame: object) -> None:
        nonlocal active_child
        print(f"[AUTO] received signal {signum}; stopping active round", flush=True)
        if active_child is not None and active_child.poll() is None:
            active_child.terminate()
            try:
                active_child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                active_child.kill()
                active_child.wait()
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, stop_child)
    signal.signal(signal.SIGINT, stop_child)
    try:
        while any(value < args.per_size_count for value in counts.values()):
            if args.max_rounds and completed_this_run >= args.max_rounds:
                break
            pending_sizes = [size for size in args.vertex_counts
                             if counts[size] < args.per_size_count]
            size = pending_sizes[(round_index - 1) % len(pending_sizes)]
            focus = choose_focus(args, rows, query_dir, output, size, round_index)
            low, high = mutation_steps(round_index)
            round_output = rounds_dir / f"round_{round_index:04d}"
            command = [
                sys.executable, str(collector),
                "--gpu-ids", str(args.gpu_id),
                "--cpu-threads", str(args.cpu_threads),
                "--cpu-numa-node", str(args.cpu_numa_node),
                "--focus-query", str(focus),
                "--candidate-count", str(args.candidate_count),
                "--target-count", str(args.per_round_target),
                "--target-speedup", str(args.target_speedup),
                "--sampling-strategy", args.sampling_strategy,
                "--mutation-steps", str(low), str(high),
                "--random-seed", str(args.base_seed + round_index),
                "--update-limit", "0",
                "--gpu-timeout", str(args.gpu_timeout),
                "--cpu-timeout", str(args.cpu_timeout),
                "--output-dir", str(round_output),
            ]
            print(f"[AUTO] round={round_index} focus={focus.name} "
                  f"size={size}v mutation={low}..{high} "
                  f"verified={len(rows)}/{target_total} by_size={counts}",
                  flush=True)
            if args.dry_run:
                print("[AUTO] command: " + " ".join(command), flush=True)
                break
            active_child = subprocess.Popen(command)
            returncode = active_child.wait()
            active_child = None
            results_path = round_output / "screening_results.csv"
            if results_path.is_file():
                for result in read_manifest(results_path):
                    if result.get("accepted") != "True":
                        continue
                    if counts[size] >= args.per_size_count:
                        break
                    add_query(Path(result["query_path"]), result, rows, known,
                              query_dir, str(results_path))
                    counts = size_counts(rows, args.vertex_counts)
                write_manifest(output / "manifest.csv", rows)
            print(f"[AUTO] round={round_index} exit={returncode} "
                  f"verified={len(rows)}/{target_total} by_size={counts}", flush=True)
            round_index += 1
            completed_this_run += 1
            time.sleep(1)
    finally:
        if pid_path.exists() and pid_path.read_text().strip() == str(os.getpid()):
            pid_path.unlink()
    finished = all(value >= args.per_size_count for value in counts.values())
    print(f"[AUTO] finished verified={len(rows)}/{target_total} by_size={counts}",
          flush=True)
    return 0 if finished else 1


if __name__ == "__main__":
    raise SystemExit(main())
