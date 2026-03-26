#!/usr/bin/env python3

import argparse
import csv
import random
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from select_slow_queries import QueryRunResult, count_query_shape, run_query


@dataclass(frozen=True)
class QueryGraph:
    labels: tuple[int, ...]
    edges: tuple[tuple[int, int, int], ...]


@dataclass
class AcceptedQuery:
    query_name: str
    seed_name: str
    mutation: str
    query_path: str
    vertex_count: int
    edge_count: int
    selected_metric: float | None
    initial_ms: float | None
    incremental_ms: float | None
    initial_matches: int | None
    positive_matches: int | None
    negative_matches: int | None
    edge_updates: int | None
    timed_out: bool
    log_path: str | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate new 8-vertex query graphs by mutating seed queries, then validate them "
            "with the matcher and keep only the slow ones."
        )
    )
    repo_root = Path(__file__).resolve().parents[1]
    data_root = Path("/home/v-haibinlai/haibin/paracosm_data")

    parser.add_argument(
        "--seed-dir",
        type=Path,
        action="append",
        required=True,
        help="Directory that contains seed query graphs. Repeat for multiple directories.",
    )
    parser.add_argument(
        "--seed-query-ids",
        type=int,
        nargs="*",
        help="Optional subset of seed query ids to use.",
    )
    parser.add_argument(
        "--exe",
        type=Path,
        default=repo_root / "matching_executor/Parallel_NewSP/build/csm",
        help="Matcher executable used for validation.",
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
        "--target-count",
        type=int,
        default=5,
        help="Stop after accepting this many new slow queries.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=100,
        help="Maximum number of generated candidates to validate.",
    )
    parser.add_argument(
        "--vertex-count",
        type=int,
        default=8,
        help="Seed queries must have exactly this many vertices.",
    )
    parser.add_argument(
        "--min-seconds",
        type=float,
        default=30.0,
        help="Keep only queries whose selected runtime is at least this threshold in seconds.",
    )
    parser.add_argument(
        "--metric",
        choices=("initial", "incremental", "total", "max"),
        default="max",
        help="Metric used to test whether a generated query is slow enough.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="OMP_NUM_THREADS used when validating a generated query.",
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=60,
        help="Value passed to --time-limit during validation.",
    )
    parser.add_argument(
        "--initial-time-limit",
        type=int,
        default=60,
        help="Value passed to --initial-time-limit during validation.",
    )
    parser.add_argument(
        "--run-timeout",
        type=int,
        default=180,
        help="Outer timeout for one validation run.",
    )
    parser.add_argument(
        "--report-initial",
        choices=("on", "off"),
        default="on",
        help="Forwarded to the matcher.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=20260326,
        help="Random seed for reproducible generation.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "logs_txt/generated_slow_queries_8v",
        help="Output directory for generated queries, logs, and summary files.",
    )
    parser.add_argument(
        "--keep-all-logs",
        action="store_true",
        help="Keep logs for all generated candidates, not only the accepted ones.",
    )
    return parser.parse_args()


def load_query(query_path: Path) -> QueryGraph:
    labels: dict[int, int] = {}
    edges: list[tuple[int, int, int]] = []
    with query_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if parts[0] == "v":
                labels[int(parts[1])] = int(parts[2])
            elif parts[0] == "e":
                u = int(parts[1])
                v = int(parts[2])
                edge_label = int(parts[3])
                if u > v:
                    u, v = v, u
                edges.append((u, v, edge_label))
    ordered_labels = tuple(labels[index] for index in sorted(labels))
    edges.sort()
    return QueryGraph(labels=ordered_labels, edges=tuple(edges))


def write_query(query: QueryGraph, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index, label in enumerate(query.labels):
            handle.write(f"v {index} {label}\n")
        for u, v, edge_label in query.edges:
            handle.write(f"e {u} {v} {edge_label}\n")


def is_connected(query: QueryGraph) -> bool:
    size = len(query.labels)
    adjacency = {index: set() for index in range(size)}
    for u, v, _ in query.edges:
        adjacency[u].add(v)
        adjacency[v].add(u)

    stack = [0]
    seen = {0}
    while stack:
        node = stack.pop()
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == size


def edge_set(query: QueryGraph) -> set[tuple[int, int, int]]:
    return set(query.edges)


def graph_signature(query: QueryGraph) -> tuple[tuple[int, ...], tuple[tuple[int, int, int], ...]]:
    return query.labels, tuple(sorted(query.edges))


def build_label_pool(seed_graphs: list[QueryGraph]) -> list[int]:
    pool: list[int] = []
    for graph in seed_graphs:
        pool.extend(graph.labels)
    return pool


def mutate_query(seed_query: QueryGraph, label_pool: list[int], rng: random.Random) -> tuple[QueryGraph, str]:
    mutation_type = rng.choice(("rewire", "swap_labels", "change_label", "rewire_and_label"))
    labels = list(seed_query.labels)
    edges = list(seed_query.edges)
    size = len(labels)

    def rewire_once() -> bool:
        existing = {(u, v) for u, v, _ in edges}
        removable = list(edges)
        rng.shuffle(removable)
        for removed_u, removed_v, removed_label in removable:
            candidate_edges = [(u, v, label) for u, v, label in edges if (u, v, label) != (removed_u, removed_v, removed_label)]
            candidate_existing = {(u, v) for u, v, _ in candidate_edges}
            all_missing = [(u, v) for u in range(size) for v in range(u + 1, size) if (u, v) not in candidate_existing]
            rng.shuffle(all_missing)
            for add_u, add_v in all_missing:
                candidate_query = QueryGraph(labels=tuple(labels), edges=tuple(sorted(candidate_edges + [(add_u, add_v, removed_label)])))
                if is_connected(candidate_query):
                    edges.clear()
                    edges.extend(candidate_query.edges)
                    return True
        return False

    if mutation_type in {"rewire", "rewire_and_label"}:
        if not rewire_once():
            mutation_type = "swap_labels"

    if mutation_type in {"swap_labels", "rewire_and_label"}:
        left, right = rng.sample(range(size), 2)
        labels[left], labels[right] = labels[right], labels[left]

    if mutation_type == "change_label":
        position = rng.randrange(size)
        current = labels[position]
        candidates = [label for label in label_pool if label != current]
        if candidates:
            labels[position] = rng.choice(candidates)
        else:
            mutation_type = "swap_labels"
            left, right = rng.sample(range(size), 2)
            labels[left], labels[right] = labels[right], labels[left]

    if mutation_type == "rewire_and_label":
        position = rng.randrange(size)
        candidates = [label for label in label_pool if label != labels[position]]
        if candidates:
            labels[position] = rng.choice(candidates)

    mutated = QueryGraph(labels=tuple(labels), edges=tuple(sorted(edges)))
    return mutated, mutation_type


def iter_seed_paths(seed_dirs: list[Path], seed_query_ids: list[int] | None) -> list[Path]:
    paths: list[Path] = []
    for seed_dir in seed_dirs:
        if seed_query_ids:
            for query_id in seed_query_ids:
                candidate = seed_dir / f"Q_{query_id}"
                if candidate.is_file():
                    paths.append(candidate)
        else:
            paths.extend(sorted(path for path in seed_dir.glob("Q_*") if path.is_file()))
    return sorted(paths, key=lambda path: path.name)


def write_summary(path: Path, accepted: list[AcceptedQuery]) -> None:
    fieldnames = [
        "query_name",
        "seed_name",
        "mutation",
        "query_path",
        "vertex_count",
        "edge_count",
        "selected_metric",
        "initial_ms",
        "incremental_ms",
        "initial_matches",
        "positive_matches",
        "negative_matches",
        "edge_updates",
        "timed_out",
        "log_path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in accepted:
            writer.writerow(item.__dict__)


def make_validator_args(args: argparse.Namespace) -> SimpleNamespace:
    return SimpleNamespace(
        exe=args.exe,
        data_graph=args.data_graph,
        update_graph=args.update_graph,
        metric=args.metric,
        min_seconds=args.min_seconds,
        threads=args.threads,
        time_limit=args.time_limit,
        initial_time_limit=args.initial_time_limit,
        run_timeout=args.run_timeout,
        report_initial=args.report_initial,
        keep_all_logs=args.keep_all_logs,
    )


def main() -> int:
    args = parse_args()
    rng = random.Random(args.random_seed)

    for seed_dir in args.seed_dir:
        if not seed_dir.is_dir():
            raise SystemExit(f"Seed directory not found: {seed_dir}")
    if not args.exe.is_file():
        raise SystemExit(f"Executable not found: {args.exe}")

    seed_paths = iter_seed_paths(args.seed_dir, args.seed_query_ids)
    seed_entries: list[tuple[Path, QueryGraph]] = []
    for seed_path in seed_paths:
        vertex_count, _ = count_query_shape(seed_path)
        if vertex_count == args.vertex_count:
            seed_entries.append((seed_path, load_query(seed_path)))

    if not seed_entries:
        raise SystemExit("No matching seed queries found.")

    label_pool = build_label_pool([graph for _, graph in seed_entries])
    known_signatures = {graph_signature(graph) for _, graph in seed_entries}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    generated_dir = args.output_dir / "generated_queries"
    generated_dir.mkdir(parents=True, exist_ok=True)
    accepted_dir = args.output_dir / "accepted_queries"
    accepted_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = args.output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    validator_args = make_validator_args(args)
    accepted: list[AcceptedQuery] = []
    attempts = 0
    generated_index = 1

    while attempts < args.max_attempts and len(accepted) < args.target_count:
        seed_path, seed_query = rng.choice(seed_entries)
        candidate_query, mutation = mutate_query(seed_query, label_pool, rng)
        signature = graph_signature(candidate_query)
        attempts += 1
        if signature in known_signatures:
            continue
        known_signatures.add(signature)

        candidate_name = f"Q_gen_{generated_index:03d}"
        generated_index += 1
        candidate_path = generated_dir / candidate_name
        write_query(candidate_query, candidate_path)

        print(f"[TRY ] {candidate_name} from {seed_path.name} using {mutation}")
        result: QueryRunResult = run_query(validator_args, candidate_path, logs_dir)
        metric_text = "NA" if result.selected_metric is None else f"{result.selected_metric / 1000.0:.3f}s"
        if result.qualifies:
            final_path = accepted_dir / candidate_name
            candidate_path.replace(final_path)
            accepted.append(
                AcceptedQuery(
                    query_name=candidate_name,
                    seed_name=seed_path.name,
                    mutation=mutation,
                    query_path=str(final_path),
                    vertex_count=result.vertex_count,
                    edge_count=result.edge_count,
                    selected_metric=result.selected_metric,
                    initial_ms=result.initial_ms,
                    incremental_ms=result.incremental_ms,
                    initial_matches=result.initial_matches,
                    positive_matches=result.positive_matches,
                    negative_matches=result.negative_matches,
                    edge_updates=result.edge_updates,
                    timed_out=result.timed_out,
                    log_path=result.log_path,
                )
            )
            print(f"[KEEP] {candidate_name} metric={metric_text}")
        else:
            print(f"[DROP] {candidate_name} metric={metric_text}")

    write_summary(args.output_dir / "accepted_summary.csv", accepted)
    print()
    print(f"Attempts: {attempts}")
    print(f"Accepted: {len(accepted)} / {args.target_count}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())