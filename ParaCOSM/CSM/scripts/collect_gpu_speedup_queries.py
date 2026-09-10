#!/usr/bin/env python3
"""Sample and collect Amazon queries with large single-GPU speedups.

The workflow is intentionally correctness-first:

1. Mutate connected sparse queries from the supplied Amazon query set.
2. Screen candidates concurrently on the selected GPUs.
3. Verify GPU-heavy candidates with CPU versioned execution.
4. Accept only queries whose CPU/GPU counts match and CPU_ms/GPU_ms reaches
   the target. CPU timeouts are kept as promising, never as verified speedups.

Example:
  python3 scripts/collect_gpu_speedup_queries.py \
    --gpu-ids 0 --cpu-threads 16 --cpu-numa-node 0 \
    --target-speedup 20 --target-count 10
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import queue
import random
import re
import shutil
import socket
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import networkx as nx


INCREMENTAL_RE = re.compile(r"Incremental Matching:\s*([0-9.]+)ms")
POSITIVE_RE = re.compile(r"^([0-9]+) positive matches\.", re.MULTILINE)
NEGATIVE_RE = re.compile(r"^([0-9]+) negative matches\.", re.MULTILINE)
GPU_SEARCH_RE = re.compile(r"BFS search:\s*([0-9.]+) ms")
CUDA_ERROR_RE = re.compile(
    r"CUDA error|invalid device function|out of memory|Segmentation fault|std::bad_alloc",
    re.IGNORECASE,
)
OVERFLOW_RE = re.compile(r"overflow|retry chunk", re.IGNORECASE)


@dataclass(frozen=True)
class QueryGraph:
    labels: tuple[int, ...]
    edges: tuple[tuple[int, int, int], ...]


@dataclass
class RunResult:
    status: str
    incremental_ms: float | None
    positive_matches: int | None
    negative_matches: int | None
    gpu_search_ms: float | None
    wall_ms: float
    overflow_seen: bool
    log_path: str


@dataclass
class CandidateResult:
    query_name: str
    query_path: str
    vertex_count: int
    edge_count: int
    gpu_id: int | None = None
    gpu_status: str = "not_run"
    gpu_incremental_ms: float | None = None
    gpu_search_ms: float | None = None
    gpu_positive_matches: int | None = None
    gpu_overflow_seen: bool = False
    cpu_status: str = "not_run"
    cpu_incremental_ms: float | None = None
    cpu_positive_matches: int | None = None
    counts_equal: bool = False
    speedup: float | None = None
    accepted: bool = False
    classification: str = "not_screened"
    gpu_log: str = ""
    cpu_log: str = ""


def parse_args() -> argparse.Namespace:
    csm_root = Path(__file__).resolve().parents[1]
    workspace = Path(__file__).resolve().parents[4]
    amazon = workspace / "datasets/amazon/AZ"
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--exe", type=Path,
                        default=csm_root / "build-v100-gcc10/bin/csm")
    parser.add_argument("--data-graph", type=Path,
                        default=amazon / "data_graph/data.graph")
    parser.add_argument("--update-graph", type=Path,
                        default=amazon / "data_graph/insertion.graph")
    parser.add_argument("--seed-root", type=Path, default=amazon)
    parser.add_argument(
        "--focus-query", type=Path,
        help="Mutate only this query instead of sampling the seed directories.",
    )
    parser.add_argument("--vertex-counts", type=int, nargs="+", default=[8, 9, 10])
    parser.add_argument("--candidate-count", type=int, default=200)
    parser.add_argument("--target-count", type=int, default=10)
    parser.add_argument("--target-speedup", type=float, default=20.0)
    parser.add_argument("--random-seed", type=int, default=20260910)
    parser.add_argument("--min-extra-edges", type=int, default=0)
    parser.add_argument("--max-extra-edges", type=int, default=3)
    parser.add_argument("--mutation-steps", type=int, nargs=2, default=[1, 4],
                        metavar=("MIN", "MAX"))
    parser.add_argument("--include-seeds", action="store_true")
    parser.add_argument("--update-limit", type=int, default=50_000,
                        help="Number of update lines used for screening; 0 means full.")
    parser.add_argument("--gpu-ids", default="0")
    parser.add_argument("--gpu-timeout", type=float, default=120.0)
    parser.add_argument("--cpu-timeout", type=float, default=300.0)
    parser.add_argument("--cpu-threads", type=int, default=16)
    parser.add_argument("--cpu-numa-node", type=int, default=0)
    parser.add_argument("--max-cpu-validations", type=int, default=80)
    parser.add_argument("--min-gpu-matches", type=int, default=1_000_000)
    parser.add_argument("--gpu-buffer-matches", type=int, default=0,
                        help="Set GPU_BFS_MAX_MATCHES; 0 uses auto sizing.")
    parser.add_argument("--numa-bind", action=argparse.BooleanOptionalAction,
                        default=True)
    parser.add_argument("--output-dir", type=Path,
                        default=csm_root / "generated_gpu_speedup_queries/amazon")
    parser.add_argument("--generate-only", action="store_true")
    return parser.parse_args()


def load_query(path: Path) -> QueryGraph:
    labels: dict[int, int] = {}
    edges: list[tuple[int, int, int]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            parts = raw.split()
            if not parts:
                continue
            if parts[0] == "v":
                labels[int(parts[1])] = int(parts[2])
            elif parts[0] == "e":
                u, v, label = map(int, parts[1:4])
                if u > v:
                    u, v = v, u
                edges.append((u, v, label))
    if not labels or sorted(labels) != list(range(len(labels))):
        raise ValueError(f"query vertex IDs must be contiguous from zero: {path}")
    return QueryGraph(
        tuple(labels[index] for index in range(len(labels))),
        tuple(sorted(edges)),
    )


def write_query(graph: QueryGraph, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for vertex, label in enumerate(graph.labels):
            handle.write(f"v {vertex} {label}\n")
        for u, v, label in graph.edges:
            handle.write(f"e {u} {v} {label}\n")


def to_networkx(graph: QueryGraph) -> nx.Graph:
    result = nx.Graph()
    for vertex, label in enumerate(graph.labels):
        result.add_node(vertex, label=label)
    for u, v, label in graph.edges:
        result.add_edge(u, v, label=label)
    return result


def is_labeled_isomorphic(left: QueryGraph, right: QueryGraph) -> bool:
    if len(left.labels) != len(right.labels) or len(left.edges) != len(right.edges):
        return False
    return nx.is_isomorphic(
        to_networkx(left), to_networkx(right),
        node_match=nx.algorithms.isomorphism.categorical_node_match("label", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("label", None),
    )


def is_connected(vertex_count: int, edges: list[tuple[int, int, int]]) -> bool:
    adjacency = [[] for _ in range(vertex_count)]
    for u, v, _ in edges:
        adjacency[u].append(v)
        adjacency[v].append(u)
    seen = {0}
    pending = [0]
    while pending:
        u = pending.pop()
        for v in adjacency[u]:
            if v not in seen:
                seen.add(v)
                pending.append(v)
    return len(seen) == vertex_count


def mutate_query(seed: QueryGraph, rng: random.Random, label_pool: list[int],
                 min_edges: int, max_edges: int,
                 steps: tuple[int, int]) -> QueryGraph:
    labels = list(seed.labels)
    edges = list(seed.edges)
    vertex_count = len(labels)

    for _ in range(rng.randint(*steps)):
        operation = rng.choice(("rewire", "swap", "label", "remove", "add"))
        pairs = {(u, v) for u, v, _ in edges}
        if operation == "rewire" and edges:
            index = rng.randrange(len(edges))
            edge_label = edges[index][2]
            remaining = edges[:index] + edges[index + 1:]
            remaining_pairs = {(u, v) for u, v, _ in remaining}
            missing = [(u, v) for u in range(vertex_count)
                       for v in range(u + 1, vertex_count)
                       if (u, v) not in remaining_pairs]
            rng.shuffle(missing)
            for u, v in missing:
                trial = remaining + [(u, v, edge_label)]
                if is_connected(vertex_count, trial):
                    edges = trial
                    break
        elif operation == "swap" and vertex_count >= 2:
            u, v = rng.sample(range(vertex_count), 2)
            labels[u], labels[v] = labels[v], labels[u]
        elif operation == "label":
            labels[rng.randrange(vertex_count)] = rng.choice(label_pool)
        elif operation == "remove" and len(edges) > min_edges:
            indices = list(range(len(edges)))
            rng.shuffle(indices)
            for index in indices:
                trial = edges[:index] + edges[index + 1:]
                if is_connected(vertex_count, trial):
                    edges = trial
                    break
        elif operation == "add" and len(edges) < max_edges:
            missing = [(u, v) for u in range(vertex_count)
                       for v in range(u + 1, vertex_count) if (u, v) not in pairs]
            if missing:
                u, v = rng.choice(missing)
                edges.append((u, v, rng.choice(edges)[2] if edges else 0))

    while len(edges) > max_edges:
        removable = [index for index in range(len(edges))
                     if is_connected(vertex_count, edges[:index] + edges[index + 1:])]
        if not removable:
            break
        edges.pop(rng.choice(removable))
    while len(edges) < min_edges:
        pairs = {(u, v) for u, v, _ in edges}
        missing = [(u, v) for u in range(vertex_count)
                   for v in range(u + 1, vertex_count) if (u, v) not in pairs]
        if not missing:
            break
        u, v = rng.choice(missing)
        edges.append((u, v, rng.choice(edges)[2] if edges else 0))
    return QueryGraph(tuple(labels), tuple(sorted(edges)))


def create_update_prefix(source: Path, destination: Path, limit: int) -> int:
    count = 0
    with source.open(encoding="utf-8") as src, destination.open(
        "w", encoding="utf-8"
    ) as dst:
        for line in src:
            if limit and count >= limit:
                break
            dst.write(line)
            count += 1
    return count


def run_matcher(args: argparse.Namespace, query_path: Path, update_path: Path,
                mode: str, threads: int, timeout: float, log_path: Path,
                gpu_id: int | None = None) -> RunResult:
    command = [
        str(args.exe), "-q", str(query_path), "-d", str(args.data_graph),
        "-u", str(update_path), "-a", "parallel_graphflow", "-m", mode,
        "-t", str(threads), "--auto-tuning", "0", "--report-initial", "0",
    ]
    env = dict(os.environ)
    env.update({
        "OMP_NUM_THREADS": str(threads),
        "OMP_DYNAMIC": "FALSE",
        "OMP_PROC_BIND": "spread",
        "OMP_PLACES": "cores",
    })
    if args.gpu_buffer_matches:
        env["GPU_BFS_MAX_MATCHES"] = str(args.gpu_buffer_matches)
    if gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        if args.numa_bind and shutil.which("numactl"):
            command = [
                "numactl", f"--cpunodebind={gpu_id}", f"--membind={gpu_id}"
            ] + command
    elif args.numa_bind and shutil.which("numactl"):
        command = [
            "numactl", f"--cpunodebind={args.cpu_numa_node}",
            f"--membind={args.cpu_numa_node}",
        ] + command

    started = time.monotonic()
    try:
        process = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=timeout, env=env, check=False,
        )
        output = process.stdout
        status = "ok" if process.returncode == 0 else "failed"
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        output = partial.decode(errors="replace") if isinstance(partial, bytes) else partial
        output += f"\n[collector] outer timeout after {timeout:.1f}s\n"
        status = "timeout"
    wall_ms = (time.monotonic() - started) * 1000.0
    log_path.write_text(
        "command: " + " ".join(command) + "\n\n" + output, encoding="utf-8"
    )

    incremental = INCREMENTAL_RE.search(output)
    positive = POSITIVE_RE.search(output)
    negative = NEGATIVE_RE.search(output)
    gpu_search = GPU_SEARCH_RE.search(output)
    if CUDA_ERROR_RE.search(output):
        status = "cuda_error"
    elif status == "ok" and (incremental is None or positive is None):
        status = "parse_error"
    return RunResult(
        status=status,
        incremental_ms=float(incremental.group(1)) if incremental else None,
        positive_matches=int(positive.group(1)) if positive else None,
        negative_matches=int(negative.group(1)) if negative else None,
        gpu_search_ms=float(gpu_search.group(1)) if gpu_search else None,
        wall_ms=wall_ms,
        overflow_seen=bool(OVERFLOW_RE.search(output)),
        log_path=str(log_path),
    )


def collect_seeds(args: argparse.Namespace) -> dict[int, list[tuple[Path, QueryGraph]]]:
    if args.focus_query is not None:
        graph = load_query(args.focus_query)
        return {len(graph.labels): [(args.focus_query, graph)]}

    result: dict[int, list[tuple[Path, QueryGraph]]] = {}
    for vertex_count in args.vertex_counts:
        seed_dir = args.seed_root / f"{vertex_count}v"
        seeds = [(path, load_query(path)) for path in sorted(seed_dir.glob("Q_*"))
                 if path.is_file()]
        seeds = [(path, graph) for path, graph in seeds
                 if len(graph.labels) == vertex_count]
        if not seeds:
            raise RuntimeError(f"no {vertex_count}-vertex seeds found in {seed_dir}")
        result[vertex_count] = seeds
    return result


def generate_candidates(args: argparse.Namespace, generated_dir: Path) -> list[Path]:
    rng = random.Random(args.random_seed)
    seeds_by_size = collect_seeds(args)
    label_pool = [label for seeds in seeds_by_size.values()
                  for _, graph in seeds for label in graph.labels]
    known: set[QueryGraph] = set()
    known_by_size: dict[tuple[int, int], list[QueryGraph]] = {}
    candidates: list[Path] = []

    def remember_if_unique(graph: QueryGraph) -> bool:
        if graph in known:
            return False
        shape = (len(graph.labels), len(graph.edges))
        if any(is_labeled_isomorphic(graph, prior)
               for prior in known_by_size.get(shape, [])):
            return False
        known.add(graph)
        known_by_size.setdefault(shape, []).append(graph)
        return True

    if args.include_seeds:
        for vertex_count in seeds_by_size:
            for seed_path, graph in seeds_by_size[vertex_count]:
                if len(candidates) >= args.candidate_count:
                    break
                if not remember_if_unique(graph):
                    continue
                destination = generated_dir / f"seed_{vertex_count}v_{seed_path.name}"
                shutil.copy2(seed_path, destination)
                candidates.append(destination)
    else:
        # The seed is a search origin, not a new result. Register its whole
        # isomorphism class so renumbered/no-op mutations are not collected.
        for seeds in seeds_by_size.values():
            for _, graph in seeds:
                remember_if_unique(graph)

    attempts = 0
    while len(candidates) < args.candidate_count and attempts < args.candidate_count * 30:
        attempts += 1
        vertex_count = rng.choice(list(seeds_by_size))
        _, seed = rng.choice(seeds_by_size[vertex_count])
        minimum = vertex_count - 1 + args.min_extra_edges
        maximum = vertex_count - 1 + args.max_extra_edges
        graph = mutate_query(
            seed, rng, label_pool, minimum, maximum, tuple(args.mutation_steps)
        )
        if (not is_connected(vertex_count, list(graph.edges))
                or not remember_if_unique(graph)):
            continue
        path = generated_dir / f"Q_gpu_{len(candidates) + 1:04d}_{vertex_count}v"
        write_query(graph, path)
        candidates.append(path)
    if len(candidates) < args.candidate_count:
        raise RuntimeError(f"generated only {len(candidates)} unique candidates")
    return candidates


def write_results(path: Path, results: list[CandidateResult]) -> None:
    if not results:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(asdict(results[0])), lineterminator="\n"
        )
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))


def validate_args(args: argparse.Namespace) -> list[int]:
    for path, label in (
        (args.exe, "executable"), (args.data_graph, "data graph"),
        (args.update_graph, "update graph"), (args.seed_root, "seed root"),
    ):
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
    if not os.access(args.exe, os.X_OK):
        raise PermissionError(f"executable is not runnable: {args.exe}")
    if args.focus_query is not None and not args.focus_query.is_file():
        raise FileNotFoundError(f"focus query not found: {args.focus_query}")
    gpu_ids = [int(item) for item in args.gpu_ids.split(",") if item.strip()]
    if not gpu_ids:
        raise ValueError("--gpu-ids must contain at least one device")
    if args.max_extra_edges < args.min_extra_edges:
        raise ValueError("max-extra-edges must be >= min-extra-edges")
    return gpu_ids


def main() -> int:
    args = parse_args()
    gpu_ids = validate_args(args)
    output = args.output_dir.resolve()
    generated_dir = output / "generated"
    verified_dir = output / "accepted_verified"
    promising_dir = output / "promising"
    logs_dir = output / "logs"
    for directory in (generated_dir, verified_dir, promising_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    update_path = output / (
        f"amazon_updates_{args.update_limit}.graph"
        if args.update_limit else "amazon_updates_full.graph"
    )
    update_count = create_update_prefix(args.update_graph, update_path, args.update_limit)
    candidates = generate_candidates(args, generated_dir)
    config = vars(args).copy()
    config.update({
        "exe": str(args.exe.resolve()),
        "data_graph": str(args.data_graph.resolve()),
        "update_graph": str(args.update_graph.resolve()),
        "seed_root": str(args.seed_root.resolve()),
        "output_dir": str(output),
        "hostname": socket.gethostname(),
        "actual_update_count": update_count,
        "gpu_ids_parsed": gpu_ids,
    })
    (output / "run_config.json").write_text(
        json.dumps(config, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print(f"Generated {len(candidates)} candidates; updates={update_count}; output={output}")
    if args.generate_only:
        return 0

    available_gpus: queue.Queue[int] = queue.Queue()
    for gpu_id in gpu_ids:
        available_gpus.put(gpu_id)

    def gpu_screen(path: Path) -> tuple[Path, int, RunResult]:
        gpu_id = available_gpus.get()
        try:
            log = logs_dir / f"{path.name}.gpu{gpu_id}.log"
            result = run_matcher(
                args, path, update_path, "gpu_bfs_versioned", 1,
                args.gpu_timeout, log, gpu_id,
            )
            return path, gpu_id, result
        finally:
            available_gpus.put(gpu_id)

    records: dict[Path, CandidateResult] = {}
    with ThreadPoolExecutor(max_workers=len(gpu_ids)) as executor:
        futures = [executor.submit(gpu_screen, path) for path in candidates]
        for done, future in enumerate(as_completed(futures), 1):
            path, gpu_id, gpu = future.result()
            graph = load_query(path)
            record = CandidateResult(
                query_name=path.name,
                query_path=str(path),
                vertex_count=len(graph.labels),
                edge_count=len(graph.edges),
                gpu_id=gpu_id,
                gpu_status=gpu.status,
                gpu_incremental_ms=gpu.incremental_ms,
                gpu_search_ms=gpu.gpu_search_ms,
                gpu_positive_matches=gpu.positive_matches,
                gpu_overflow_seen=gpu.overflow_seen,
                gpu_log=gpu.log_path,
            )
            if gpu.status != "ok":
                record.classification = "gpu_failed"
            elif (gpu.positive_matches or 0) < args.min_gpu_matches:
                record.classification = "gpu_too_little_work"
            else:
                record.classification = "awaiting_cpu"
            records[path] = record
            print(
                f"[GPU {done:4d}/{len(candidates)}] {path.name} dev={gpu_id} "
                f"status={gpu.status} matches={gpu.positive_matches} "
                f"incr_ms={gpu.incremental_ms} overflow={gpu.overflow_seen}",
                flush=True,
            )

    cpu_queue = [
        records[path] for path in candidates
        if records[path].classification == "awaiting_cpu"
    ]
    for record in cpu_queue:
        if ((record.gpu_incremental_ms or float("inf")) * args.target_speedup
                > args.cpu_timeout * 1000.0):
            record.classification = "gpu_too_slow_for_cpu_timeout"
    cpu_queue = [record for record in cpu_queue
                 if record.classification == "awaiting_cpu"]
    cpu_queue.sort(
        key=lambda item: (
            item.gpu_incremental_ms or float("inf"),
            -((item.gpu_positive_matches or 0)
              / (item.gpu_incremental_ms or float("inf"))),
        )
    )
    cpu_queue = cpu_queue[:args.max_cpu_validations]

    accepted = 0
    for index, record in enumerate(cpu_queue, 1):
        path = Path(record.query_path)
        cpu_log = logs_dir / f"{path.name}.cpu{args.cpu_threads}.log"
        cpu = run_matcher(
            args, path, update_path, "versioned", args.cpu_threads,
            args.cpu_timeout, cpu_log,
        )
        record.cpu_status = cpu.status
        record.cpu_incremental_ms = cpu.incremental_ms
        record.cpu_positive_matches = cpu.positive_matches
        record.cpu_log = cpu.log_path
        if cpu.status == "timeout":
            record.classification = "promising_cpu_timeout"
            shutil.copy2(path, promising_dir / path.name)
        elif cpu.status != "ok":
            record.classification = "cpu_failed"
        else:
            record.counts_equal = (
                cpu.positive_matches == record.gpu_positive_matches
                and cpu.negative_matches == 0
            )
            if not record.counts_equal:
                record.classification = "count_mismatch"
            elif record.gpu_incremental_ms and cpu.incremental_ms:
                record.speedup = cpu.incremental_ms / record.gpu_incremental_ms
                if record.speedup >= args.target_speedup:
                    record.accepted = True
                    record.classification = "verified"
                    shutil.copy2(path, verified_dir / path.name)
                    accepted += 1
                else:
                    record.classification = "verified_below_target"
        print(
            f"[CPU {index:3d}/{len(cpu_queue)}] {path.name} status={cpu.status} "
            f"counts_equal={record.counts_equal} speedup={record.speedup} "
            f"class={record.classification}",
            flush=True,
        )
        write_results(output / "screening_results.csv", list(records.values()))
        if accepted >= args.target_count:
            break

    ordered = sorted(
        records.values(),
        key=lambda item: item.speedup if item.speedup is not None else -1,
        reverse=True,
    )
    write_results(output / "screening_results.csv", ordered)
    print(f"Verified >= {args.target_speedup:.1f}x: {accepted}")
    print(f"Results: {output / 'screening_results.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
