#!/usr/bin/env python3
"""
find_high_speedup_queries.py — Automated workflow to discover 9-vertex query graphs
that achieve ≥4x parallel speedup with parallel_graphflow.

Workflow:
  1. Generate random 9v query graphs by mutating seeds (varied edge counts)
  2. Screen each: run single-thread and 6-thread, measure incremental speedup
  3. Keep only queries with speedup ≥ target (default 4.0x)
  4. Output accepted queries + summary CSV

Usage:
  python3 find_high_speedup_queries.py [--target-speedup 4.0] [--target-count 5]
"""

import argparse
import csv
import os
import random
import re
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path


# ── Regex for parsing CSM output ──
INITIAL_TIME_RE = re.compile(r"Initial Matching:\s+([0-9.]+)ms", re.MULTILINE)
INCR_TIME_RE = re.compile(r"Incremental Matching:\s+([0-9.]+)ms", re.MULTILINE)
INCR_TIME_RE2 = re.compile(r"#Time:\s+([0-9.]+)\s+ms", re.MULTILINE)
POS_MATCHES_RE = re.compile(r"^([0-9]+) positive matches\.", re.MULTILINE)
TIMEOUT_RE = re.compile(r"Timeout\s+([0-9]+)s")


@dataclass
class QueryGraph:
    labels: tuple
    edges: tuple

@dataclass
class ScreenResult:
    query_name: str
    query_path: str
    vertex_count: int
    edge_count: int
    single_init_ms: float
    single_incr_ms: float
    single_total_ms: float
    parallel_init_ms: float
    parallel_incr_ms: float
    parallel_total_ms: float
    speedup_incr: float
    speedup_total: float
    positive_matches: int
    timed_out: bool
    accepted: bool


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    repo = Path(__file__).resolve().parents[1]
    data = Path("/home/v-haibinlai/haibin/paracosm_data")

    p.add_argument("--exe", type=Path,
                   default=repo / "build/bin/csm")
    p.add_argument("--data-graph", type=Path,
                   default=data / "6/data_graph/data.graph")
    p.add_argument("--update-graph", type=Path,
                   default=data / "6/data_graph/insertion.graph")
    p.add_argument("--seed-dir", type=Path, action="append",
                   default=[data / "9_self/sparse"])
    p.add_argument("--vertex-count", type=int, default=9)
    p.add_argument("--target-speedup", type=float, default=4.0,
                   help="Minimum incremental speedup to accept a query.")
    p.add_argument("--target-count", type=int, default=5,
                   help="Stop after finding this many accepted queries.")
    p.add_argument("--max-attempts", type=int, default=300,
                   help="Maximum number of candidates to generate and screen.")
    p.add_argument("--min-single-ms", type=float, default=5000,
                   help="Minimum single-thread incremental time (ms) to consider.")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--time-limit", type=int, default=60,
                   help="CSM --time-limit for screening runs.")
    p.add_argument("--run-timeout", type=int, default=180,
                   help="Outer subprocess timeout.")
    p.add_argument("--random-seed", type=int, default=20260327)
    p.add_argument("--output-dir", type=Path,
                   default=repo / "generated_slow_queries/high_speedup_9v")
    p.add_argument("--edge-range", type=int, nargs=2, default=[8, 11],
                   help="Target edge count range (inclusive). Sparser = better speedup.")
    return p.parse_args()


# ── Query I/O ──

def load_query(path: Path) -> QueryGraph:
    labels, edges = {}, []
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts: continue
        if parts[0] == "v":
            labels[int(parts[1])] = int(parts[2])
        elif parts[0] == "e":
            u, v, el = int(parts[1]), int(parts[2]), int(parts[3])
            if u > v: u, v = v, u
            edges.append((u, v, el))
    ordered = tuple(labels[i] for i in sorted(labels))
    edges.sort()
    return QueryGraph(labels=ordered, edges=tuple(edges))


def write_query(q: QueryGraph, path: Path):
    with path.open("w") as f:
        for i, l in enumerate(q.labels):
            f.write(f"v {i} {l}\n")
        for u, v, el in q.edges:
            f.write(f"e {u} {v} {el}\n")


def count_shape(path: Path):
    vc = ec = 0
    for line in path.read_text().splitlines():
        if line.startswith("v "): vc += 1
        elif line.startswith("e "): ec += 1
    return vc, ec


def is_connected(q: QueryGraph) -> bool:
    n = len(q.labels)
    adj = {i: set() for i in range(n)}
    for u, v, _ in q.edges:
        adj[u].add(v); adj[v].add(u)
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for nb in adj[node]:
            if nb not in seen:
                seen.add(nb); stack.append(nb)
    return len(seen) == n


def graph_sig(q: QueryGraph):
    return (q.labels, tuple(sorted(q.edges)))


# ── Mutation ──

def mutate(seed: QueryGraph, label_pool: list, rng: random.Random,
           target_edges: range) -> QueryGraph:
    """Mutate seed, optionally adjusting edge count toward target range."""
    labels = list(seed.labels)
    edges = list(seed.edges)
    n = len(labels)
    cur_edges = len(edges)

    # Decide mutation type based on edge count vs target
    if cur_edges > target_edges.stop:
        # Too many edges — try to remove one
        mutation = "remove_edge"
    elif cur_edges < target_edges.start:
        # Too few edges — try to add one
        mutation = "add_edge"
    else:
        mutation = rng.choice(["rewire", "swap_labels", "change_label",
                               "remove_edge", "add_edge", "rewire_and_label"])

    existing = {(u, v) for u, v, _ in edges}

    if mutation == "remove_edge" and len(edges) > n - 1:
        rng.shuffle(edges)
        for idx, (eu, ev, el) in enumerate(edges):
            candidate = edges[:idx] + edges[idx+1:]
            cq = QueryGraph(labels=tuple(labels), edges=tuple(sorted(candidate)))
            if is_connected(cq):
                return cq
        mutation = "swap_labels"  # fallback

    if mutation == "add_edge":
        missing = [(u, v) for u in range(n) for v in range(u+1, n) if (u,v) not in existing]
        if missing:
            u, v = rng.choice(missing)
            edges.append((u, v, 0))
            return QueryGraph(labels=tuple(labels), edges=tuple(sorted(edges)))
        mutation = "swap_labels"

    if mutation in ("rewire", "rewire_and_label"):
        removable = list(edges)
        rng.shuffle(removable)
        for ru, rv, rl in removable:
            cand_edges = [(u,v,l) for u,v,l in edges if (u,v,l) != (ru,rv,rl)]
            cand_set = {(u,v) for u,v,_ in cand_edges}
            missing = [(u,v) for u in range(n) for v in range(u+1,n) if (u,v) not in cand_set]
            rng.shuffle(missing)
            for au, av in missing:
                cq = QueryGraph(labels=tuple(labels), edges=tuple(sorted(cand_edges + [(au,av,rl)])))
                if is_connected(cq):
                    if mutation == "rewire_and_label":
                        ll = list(cq.labels)
                        pos = rng.randrange(n)
                        cands = [l for l in label_pool if l != ll[pos]]
                        if cands: ll[pos] = rng.choice(cands)
                        return QueryGraph(labels=tuple(ll), edges=cq.edges)
                    return cq
        mutation = "swap_labels"

    if mutation == "swap_labels":
        a, b = rng.sample(range(n), 2)
        labels[a], labels[b] = labels[b], labels[a]

    if mutation == "change_label":
        pos = rng.randrange(n)
        cands = [l for l in label_pool if l != labels[pos]]
        if cands:
            labels[pos] = rng.choice(cands)
        else:
            a, b = rng.sample(range(n), 2)
            labels[a], labels[b] = labels[b], labels[a]

    return QueryGraph(labels=tuple(labels), edges=tuple(sorted(edges)))


# ── CSM Runner ──

def extract_first(pattern, text):
    m = pattern.search(text)
    return float(m.group(1)) if m else None


def run_csm(exe, data, update, query, algo, threads, time_limit, run_timeout):
    cmd = [str(exe), "-a", algo, "-d", str(data), "-u", str(update),
           "-q", str(query), "--report-initial", "on",
           "--initial-time-limit", str(time_limit),
           "--time-limit", str(time_limit), "-t", str(threads)]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=run_timeout, env=env, check=False)
        out = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    except subprocess.TimeoutExpired:
        out = ""

    init = extract_first(INITIAL_TIME_RE, out)
    incr = extract_first(INCR_TIME_RE, out)
    if incr is None:
        incr = extract_first(INCR_TIME_RE2, out)
    pos = extract_first(POS_MATCHES_RE, out)
    to = bool(TIMEOUT_RE.search(out))
    return init, incr, pos, to


def screen_query(args, query_path: Path) -> ScreenResult:
    """Run single-thread and parallel, compute speedup."""
    vc, ec = count_shape(query_path)
    name = query_path.name

    # Single-thread
    s_init, s_incr, s_pos, s_to = run_csm(
        args.exe, args.data_graph, args.update_graph, query_path,
        "graphflow", 1, args.time_limit, args.run_timeout)

    # If single-thread is too fast, not interesting
    if s_incr is not None and s_incr < args.min_single_ms:
        return ScreenResult(
            query_name=name, query_path=str(query_path),
            vertex_count=vc, edge_count=ec,
            single_init_ms=s_init or 0, single_incr_ms=s_incr or 0,
            single_total_ms=(s_init or 0) + (s_incr or 0),
            parallel_init_ms=0, parallel_incr_ms=0, parallel_total_ms=0,
            speedup_incr=0, speedup_total=0,
            positive_matches=int(s_pos) if s_pos else 0,
            timed_out=s_to, accepted=False)

    # Parallel
    p_init, p_incr, p_pos, p_to = run_csm(
        args.exe, args.data_graph, args.update_graph, query_path,
        "parallel_graphflow", args.threads, args.time_limit, args.run_timeout)

    s_total = ((s_init or 0) + (s_incr or 0)) if s_incr else 0
    p_total = ((p_init or 0) + (p_incr or 0)) if p_incr else 0

    sp_incr = (s_incr / p_incr) if (s_incr and p_incr and p_incr > 0) else 0
    sp_total = (s_total / p_total) if (s_total and p_total and p_total > 0) else 0

    accepted = sp_incr >= args.target_speedup and not p_to

    return ScreenResult(
        query_name=name, query_path=str(query_path),
        vertex_count=vc, edge_count=ec,
        single_init_ms=s_init or 0, single_incr_ms=s_incr or 0,
        single_total_ms=s_total,
        parallel_init_ms=p_init or 0, parallel_incr_ms=p_incr or 0,
        parallel_total_ms=p_total,
        speedup_incr=sp_incr, speedup_total=sp_total,
        positive_matches=int(s_pos or p_pos or 0),
        timed_out=p_to, accepted=accepted)


# ── Main ──

def main():
    args = parse_args()
    rng = random.Random(args.random_seed)

    # Validate paths
    for p, name in [(args.exe, "exe"), (args.data_graph, "data-graph"),
                    (args.update_graph, "update-graph")]:
        if not p.is_file():
            sys.exit(f"{name} not found: {p}")

    # Load seeds
    seeds = []
    for sd in args.seed_dir:
        for sp in sorted(sd.glob("Q_*")):
            vc, ec = count_shape(sp)
            if vc == args.vertex_count:
                seeds.append((sp, load_query(sp)))
    if not seeds:
        sys.exit("No matching seed queries found.")
    print(f"Loaded {len(seeds)} seed queries")

    label_pool = []
    for _, g in seeds:
        label_pool.extend(g.labels)
    label_pool = list(set(label_pool))

    known = {graph_sig(g) for _, g in seeds}

    # Setup output
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gen_dir = args.output_dir / "generated"
    acc_dir = args.output_dir / "accepted"
    gen_dir.mkdir(exist_ok=True)
    acc_dir.mkdir(exist_ok=True)

    target_edge_range = range(args.edge_range[0], args.edge_range[1] + 1)

    accepted: list[ScreenResult] = []
    all_results: list[ScreenResult] = []
    gen_idx = 1
    attempts = 0

    print(f"Target: {args.target_count} queries with ≥{args.target_speedup}x speedup")
    print(f"Edge range: {args.edge_range[0]}-{args.edge_range[1]}, "
          f"min single-thread: {args.min_single_ms}ms")
    print()

    while attempts < args.max_attempts and len(accepted) < args.target_count:
        _, seed_q = rng.choice(seeds)
        candidate = mutate(seed_q, label_pool, rng, target_edge_range)
        sig = graph_sig(candidate)
        attempts += 1

        if sig in known:
            continue
        known.add(sig)

        ec = len(candidate.edges)
        if ec not in target_edge_range:
            continue

        cand_name = f"Q_hsp_{gen_idx:03d}"
        gen_idx += 1
        cand_path = gen_dir / cand_name
        write_query(candidate, cand_path)

        print(f"[{attempts:>3}/{args.max_attempts}] {cand_name} (e={ec}) ", end="", flush=True)

        t0 = time.time()
        result = screen_query(args, cand_path)
        elapsed = time.time() - t0
        all_results.append(result)

        si = f"{result.single_incr_ms/1000:.1f}s" if result.single_incr_ms else "NA"
        pi = f"{result.parallel_incr_ms/1000:.1f}s" if result.parallel_incr_ms else "NA"

        if result.accepted:
            final = acc_dir / cand_name
            cand_path.rename(final)
            result.query_path = str(final)
            accepted.append(result)
            print(f"1T={si} 6T={pi} speedup={result.speedup_incr:.2f}x ★ ACCEPTED ({elapsed:.0f}s)")
        elif result.speedup_incr > 0:
            print(f"1T={si} 6T={pi} speedup={result.speedup_incr:.2f}x ({elapsed:.0f}s)")
        else:
            reason = "too fast" if result.single_incr_ms and result.single_incr_ms < args.min_single_ms else "TO/err"
            print(f"1T={si} 6T={pi} [{reason}] ({elapsed:.0f}s)")

    # Also screen existing seeds
    if len(accepted) < args.target_count:
        print(f"\nScreening {len(seeds)} existing seed queries...")
        for sp, _ in seeds:
            if len(accepted) >= args.target_count:
                break
            ec = count_shape(sp)[1]
            if ec not in target_edge_range:
                continue
            print(f"[SEED] {sp.name} (e={ec}) ", end="", flush=True)
            t0 = time.time()
            result = screen_query(args, sp)
            elapsed = time.time() - t0
            all_results.append(result)

            si = f"{result.single_incr_ms/1000:.1f}s" if result.single_incr_ms else "NA"
            pi = f"{result.parallel_incr_ms/1000:.1f}s" if result.parallel_incr_ms else "NA"

            if result.accepted:
                import shutil
                shutil.copy2(sp, acc_dir / sp.name)
                result.query_path = str(acc_dir / sp.name)
                accepted.append(result)
                print(f"1T={si} 6T={pi} speedup={result.speedup_incr:.2f}x ★ ACCEPTED ({elapsed:.0f}s)")
            elif result.speedup_incr > 0:
                print(f"1T={si} 6T={pi} speedup={result.speedup_incr:.2f}x ({elapsed:.0f}s)")
            else:
                print(f"1T={si} 6T={pi} [skip] ({elapsed:.0f}s)")

    # Write CSV
    csv_path = args.output_dir / "screening_results.csv"
    fields = list(asdict(all_results[0]).keys()) if all_results else []
    if fields:
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in sorted(all_results, key=lambda x: x.speedup_incr, reverse=True):
                w.writerow(asdict(r))

    # Summary
    print(f"\n{'='*80}")
    print(f"Attempts: {attempts}, Generated: {gen_idx-1}, Accepted: {len(accepted)}/{args.target_count}")
    print(f"Output: {args.output_dir}")
    if accepted:
        print(f"\nAccepted queries (speedup ≥ {args.target_speedup}x):")
        print(f"  {'Name':<18} {'E':>3} {'1T_incr':>10} {'6T_incr':>10} {'Speedup':>8}")
        print(f"  {'-'*55}")
        for r in sorted(accepted, key=lambda x: x.speedup_incr, reverse=True):
            print(f"  {r.query_name:<18} {r.edge_count:>3}"
                  f"  {r.single_incr_ms/1000:>9.2f}s"
                  f"  {r.parallel_incr_ms/1000:>9.2f}s"
                  f"  {r.speedup_incr:>7.2f}x")

    # Also show top-10 by speedup regardless of threshold
    screened = [r for r in all_results if r.speedup_incr > 0]
    if screened:
        print(f"\nTop-10 by speedup (all screened):")
        print(f"  {'Name':<18} {'E':>3} {'1T_incr':>10} {'6T_incr':>10} {'Speedup':>8}")
        print(f"  {'-'*55}")
        for r in sorted(screened, key=lambda x: x.speedup_incr, reverse=True)[:10]:
            mark = "★" if r.accepted else " "
            print(f"  {r.query_name:<18} {r.edge_count:>3}"
                  f"  {r.single_incr_ms/1000:>9.2f}s"
                  f"  {r.parallel_incr_ms/1000:>9.2f}s"
                  f"  {r.speedup_incr:>7.2f}x {mark}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
