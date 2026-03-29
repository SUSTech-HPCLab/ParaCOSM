#!/usr/bin/env python3
"""
Generate random connected 8-vertex query graphs with sparse edges (7..23).
All edge labels are 0. Vertex labels drawn from LJ's label set (0..29).
Outputs one file per query in the specified directory.
"""
import argparse
import os
import random


def gen_random_tree(n):
    """Generate a random labeled tree on n vertices (Prüfer sequence)."""
    if n == 1:
        return []
    if n == 2:
        return [(0, 1)]
    seq = [random.randint(0, n - 1) for _ in range(n - 2)]
    degree = [1] * n
    for v in seq:
        degree[v] += 1
    edges = []
    for v in seq:
        for u in range(n):
            if degree[u] == 1:
                edges.append((min(u, v), max(u, v)))
                degree[u] -= 1
                degree[v] -= 1
                break
    # last edge
    last = [u for u in range(n) if degree[u] == 1]
    edges.append((min(last[0], last[1]), max(last[0], last[1])))
    return edges


def gen_connected_graph(n, num_edges, max_labels=30):
    """Generate a connected graph with n vertices and num_edges edges."""
    assert num_edges >= n - 1, "Need at least n-1 edges for connectivity"
    max_possible = n * (n - 1) // 2
    assert num_edges <= max_possible

    # Start with a random spanning tree
    edges = set()
    tree_edges = gen_random_tree(n)
    for e in tree_edges:
        edges.add(e)

    # Add random edges until we reach num_edges
    all_possible = [(i, j) for i in range(n) for j in range(i + 1, n)]
    random.shuffle(all_possible)
    for e in all_possible:
        if len(edges) >= num_edges:
            break
        edges.add(e)

    # Assign random vertex labels
    vlabels = [random.randint(0, max_labels - 1) for _ in range(n)]

    return vlabels, sorted(edges)


def write_query(path, vlabels, edges, elabel=0):
    with open(path, 'w') as f:
        for i, lbl in enumerate(vlabels):
            f.write(f"v {i} {lbl}\n")
        for u, v in edges:
            f.write(f"e {u} {v} {elabel}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-vertices", type=int, default=8)
    parser.add_argument("--min-edges", type=int, default=7,
                        help="Minimum edges (>= n-1 for connectivity)")
    parser.add_argument("--max-edges", type=int, default=23,
                        help="Maximum edges (< n*(n-1)/2 = 28 for 8 vertices)")
    parser.add_argument("--max-labels", type=int, default=30,
                        help="Vertex labels in [0, max_labels)")
    parser.add_argument("--count", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260327)
    parser.add_argument("--output-dir", type=str, required=True)
    args = parser.parse_args()

    random.seed(args.seed)
    os.makedirs(args.output_dir, exist_ok=True)

    n = args.num_vertices
    assert args.min_edges >= n - 1

    for i in range(args.count):
        num_edges = random.randint(args.min_edges, args.max_edges)
        vlabels, edges = gen_connected_graph(n, num_edges, args.max_labels)
        name = f"Q_gen_{i+1:03d}"
        write_query(os.path.join(args.output_dir, name), vlabels, edges)

    print(f"Generated {args.count} queries in {args.output_dir}")
    # Print edge distribution
    from collections import Counter
    edge_counts = []
    for i in range(args.count):
        name = f"Q_gen_{i+1:03d}"
        with open(os.path.join(args.output_dir, name)) as f:
            ec = sum(1 for line in f if line.startswith('e '))
        edge_counts.append(ec)
    ctr = Counter(edge_counts)
    print("Edge count distribution:")
    for k in sorted(ctr):
        print(f"  {k} edges: {ctr[k]} queries")


if __name__ == "__main__":
    main()
