#!/usr/bin/env python3
"""
Generate GPU profiling charts for ParaCOSM PPT presentation.
Compares GPU DFS vs GPU BFS on NCU metrics + performance data.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import csv
import os

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    'font.size': 13,
    'figure.dpi': 200,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
    'axes.titlesize': 15,
    'axes.labelsize': 13,
})

COLORS = {
    'dfs': '#e74c3c',
    'bfs': '#2ecc71',
    'neutral': '#3498db',
    'gray': '#bdc3c7',
}

# ================================================================
# Chart 1: NCU Key Metrics Comparison (bar chart)
# ================================================================
def chart_ncu_comparison():
    metrics = [
        'Achieved\nOccupancy (%)',
        'Active Threads\nPer Warp',
        'Registers\nPer Thread',
        'Active Warps\nPer Scheduler',
    ]
    dfs_vals = [5.65, 1.02, 72, 1.97]
    bfs_vals = [16.50, 2.07, 40, 4.59]

    x = np.arange(len(metrics))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 5))
    bars1 = ax.bar(x - width/2, dfs_vals, width, label='GPU DFS', color=COLORS['dfs'], edgecolor='white', linewidth=0.5)
    bars2 = ax.bar(x + width/2, bfs_vals, width, label='GPU BFS', color=COLORS['bfs'], edgecolor='white', linewidth=0.5)

    # Add value labels
    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{bar.get_height():.1f}' if bar.get_height() < 10 else f'{bar.get_height():.0f}',
                ha='center', va='bottom', fontsize=11, fontweight='bold', color=COLORS['dfs'])
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{bar.get_height():.1f}' if bar.get_height() < 10 else f'{bar.get_height():.0f}',
                ha='center', va='bottom', fontsize=11, fontweight='bold', color=COLORS['bfs'])

    ax.set_ylabel('Value')
    ax.set_title('NCU Profiling: GPU DFS vs GPU BFS (8v/Q_001, Amazon, A100)')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend(fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, max(max(dfs_vals), max(bfs_vals)) * 1.2)

    fig.savefig(os.path.join(OUT_DIR, 'ncu_metrics_comparison.png'))
    plt.close(fig)
    print('  -> ncu_metrics_comparison.png')


# ================================================================
# Chart 2: Local Memory Traffic (dramatic comparison)
# ================================================================
def chart_local_memory():
    fig, ax = plt.subplots(figsize=(6, 5))

    categories = ['GPU DFS', 'GPU BFS']
    # local memory load + store in GB
    load_vals = [58.75, 0]
    store_vals = [39.65, 0]

    x = np.arange(len(categories))
    width = 0.5

    bars1 = ax.bar(x, load_vals, width, label='Local Load', color=COLORS['dfs'], alpha=0.85)
    bars2 = ax.bar(x, store_vals, width, bottom=load_vals, label='Local Store', color='#c0392b', alpha=0.85)

    ax.text(0, 98.4 + 2, '98.4 GB', ha='center', va='bottom', fontsize=14, fontweight='bold', color=COLORS['dfs'])
    ax.text(1, 1.5, '0 bytes', ha='center', va='bottom', fontsize=14, fontweight='bold', color=COLORS['bfs'])

    ax.set_ylabel('Local Memory Traffic (GB)')
    ax.set_title('Recursive Stack Spill: Local Memory Traffic')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=13)
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, 120)

    # Add annotation
    ax.annotate('DFS recursion spills\n72 registers/thread\nto slow local memory',
                xy=(0, 60), fontsize=10, ha='center', color='white', fontweight='bold')

    fig.savefig(os.path.join(OUT_DIR, 'local_memory_traffic.png'))
    plt.close(fig)
    print('  -> local_memory_traffic.png')


# ================================================================
# Chart 3: Kernel Time Breakdown (BFS pipeline)
# ================================================================
def chart_bfs_pipeline():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: DFS single kernel
    ax = axes[0]
    ax.barh(['gpu_batch_search\n_kernel'], [4572], color=COLORS['dfs'], height=0.4)
    ax.set_xlabel('Duration (ms)')
    ax.set_title('GPU DFS: Single Monolithic Kernel')
    ax.text(4572 + 50, 0, '4,572 ms', va='center', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 5500)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Right: BFS multi-stage pipeline
    ax = axes[1]
    stages = ['Init', 'Expand (d2->3)', 'Expand (d3->4)', 'Expand (d4->5)',
              'Expand (d5->6)', 'Expand (d6->7)', 'Count']
    times = [0.089, 2.76, 2.97, 4.33, 13.06, 26.64, 21.00]
    colors = ['#27ae60'] + ['#2ecc71'] * 5 + ['#1abc9c']

    bars = ax.barh(stages, times, color=colors, height=0.6)
    ax.set_xlabel('Duration (ms)')
    ax.set_title('GPU BFS: Multi-Stage Pipeline')
    for bar, t in zip(bars, times):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                f'{t:.1f} ms', va='center', fontsize=10)
    ax.set_xlim(0, 35)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add total annotation
    total_bfs = sum(times)
    fig.text(0.75, 0.02, f'BFS total: {total_bfs:.1f} ms  |  DFS: 4,572 ms  |  Speedup: {4572/total_bfs:.0f}x',
             ha='center', fontsize=12, fontweight='bold', color=COLORS['neutral'])

    fig.suptitle('Kernel Duration Comparison (8v/Q_001)', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'kernel_time_breakdown.png'))
    plt.close(fig)
    print('  -> kernel_time_breakdown.png')


# ================================================================
# Chart 4: Warp Utilization (thread activity per warp)
# ================================================================
def chart_warp_utilization():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for ax, label, active, pred, color in [
        (axes[0], 'GPU DFS', 1.02, 0.94, COLORS['dfs']),
        (axes[1], 'GPU BFS', 2.07, 1.87, COLORS['bfs']),
    ]:
        # Draw 32 thread slots in a 4x8 grid
        active_int = int(round(active))
        pred_int = int(round(pred))
        for i in range(32):
            row, col = i // 8, i % 8
            if i < pred_int:
                c = color
                alpha = 1.0
            elif i < active_int:
                c = color
                alpha = 0.4
            else:
                c = COLORS['gray']
                alpha = 0.3
            rect = plt.Rectangle((col, 3 - row), 0.9, 0.9, facecolor=c, alpha=alpha,
                                  edgecolor='white', linewidth=1.5)
            ax.add_patch(rect)

        ax.set_xlim(-0.2, 8.2)
        ax.set_ylim(-0.5, 4.5)
        ax.set_aspect('equal')
        ax.set_title(f'{label}\nAvg Active: {active}/32 ({active/32*100:.1f}%)', fontsize=12)
        ax.axis('off')

    fig.suptitle('Average Active Threads Per Warp', fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'warp_utilization.png'))
    plt.close(fig)
    print('  -> warp_utilization.png')


# ================================================================
# Chart 5: Speedup across algorithms (from results.csv)
# ================================================================
def chart_speedup_by_algo():
    csv_path = os.path.join(os.path.dirname(OUT_DIR), '..', 'logs', 'az_full_0517', 'results.csv')
    csv_path = os.path.normpath(csv_path)

    if not os.path.exists(csv_path):
        print(f'  [SKIP] results.csv not found at {csv_path}')
        return

    # Read CSV
    data = {}  # algo -> list of speedups
    with open(csv_path) as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 10:
                continue
            dataset, size, query, algo, single_ms, gpu_ms, single_pos, gpu_pos, correct, speedup = row[:10]
            if size != '8v':
                continue
            if correct not in ('PASS', 'PASS_vs_TF'):
                continue
            try:
                sp = float(speedup)
                if sp > 0:
                    data.setdefault(algo, []).append(sp)
            except ValueError:
                continue

    if not data:
        print('  [SKIP] No 8v data found')
        return

    algo_order = ['turboflux', 'graphflow', 'symbi', 'calig', 'newsp']
    algo_labels = ['TurboFlux', 'GraphFlow', 'SymBi', 'CaLiG', 'NewSP']
    algo_colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']

    fig, ax = plt.subplots(figsize=(10, 5.5))

    positions = []
    for i, (algo, label, color) in enumerate(zip(algo_order, algo_labels, algo_colors)):
        if algo not in data:
            continue
        vals = sorted(data[algo])
        pos = i
        positions.append(pos)

        bp = ax.boxplot([vals], positions=[pos], widths=0.5, patch_artist=True,
                        boxprops=dict(facecolor=color, alpha=0.6),
                        medianprops=dict(color='black', linewidth=2),
                        whiskerprops=dict(linewidth=1.2),
                        capprops=dict(linewidth=1.2),
                        flierprops=dict(marker='o', markersize=3, alpha=0.4))

        # Add median and max labels
        median = np.median(vals)
        maxval = max(vals)
        ax.text(pos, median + maxval * 0.02, f'med={median:.1f}x', ha='center', fontsize=9, fontweight='bold')

    ax.set_xticks(range(len(algo_order)))
    ax.set_xticklabels(algo_labels, fontsize=12)
    ax.set_ylabel('Speedup (single / GPU BFS versioned)')
    ax.set_title('GPU BFS Speedup Distribution (Amazon, 8-vertex queries, 100 queries per algorithm)')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{y:.0f}x' if y >= 1 else f'{y:.1f}x'))
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(4.6, 1.1, 'breakeven', color='gray', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    fig.savefig(os.path.join(OUT_DIR, 'speedup_8v_by_algo.png'))
    plt.close(fig)
    print('  -> speedup_8v_by_algo.png')


# ================================================================
# Chart 6: End-to-end time breakdown (BFS pipeline stages)
# ================================================================
def chart_e2e_breakdown():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # DFS
    dfs_stages = ['Classify', 'Add Edges', 'Build CSR', 'GPU Search']
    dfs_times = [10.6, 61.4, 367.7, 3509.1]
    dfs_colors = ['#f1c40f', '#e67e22', '#e74c3c', '#c0392b']

    ax = axes[0]
    wedges, texts, autotexts = ax.pie(dfs_times, labels=dfs_stages, colors=dfs_colors,
                                       autopct='%1.0f%%', startangle=90,
                                       textprops={'fontsize': 10})
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title(f'GPU DFS\nTotal: {sum(dfs_times):.0f} ms', fontsize=12, fontweight='bold')

    # BFS
    bfs_stages = ['Classify', 'Add Edges', 'Build CSR', 'BFS Search']
    bfs_times = [38.8, 106.1, 365.1, 57.5]
    bfs_colors = ['#f1c40f', '#e67e22', '#2ecc71', '#27ae60']

    ax = axes[1]
    wedges, texts, autotexts = ax.pie(bfs_times, labels=bfs_stages, colors=bfs_colors,
                                       autopct='%1.0f%%', startangle=90,
                                       textprops={'fontsize': 10})
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title(f'GPU BFS Versioned\nTotal: {sum(bfs_times):.0f} ms', fontsize=12, fontweight='bold')

    fig.suptitle('End-to-End Time Breakdown (8v/Q_001)', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'e2e_time_breakdown.png'))
    plt.close(fig)
    print('  -> e2e_time_breakdown.png')


# ================================================================
# Chart 7: Occupancy comparison (visual gauge)
# ================================================================
def chart_occupancy():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    for ax, label, theoretical, achieved, color in [
        (axes[0], 'GPU DFS', 43.75, 5.65, COLORS['dfs']),
        (axes[1], 'GPU BFS', 75.00, 16.50, COLORS['bfs']),
    ]:
        # Draw stacked horizontal bar
        ax.barh([0], [achieved], height=0.5, color=color, alpha=0.9, label=f'Achieved: {achieved}%')
        ax.barh([0], [theoretical - achieved], height=0.5, left=[achieved],
                color=color, alpha=0.25, label=f'Theoretical: {theoretical}%')
        ax.barh([0], [100 - theoretical], height=0.5, left=[theoretical],
                color=COLORS['gray'], alpha=0.15)

        ax.set_xlim(0, 105)
        ax.set_title(f'{label}', fontsize=13, fontweight='bold')
        ax.set_yticks([])
        ax.set_xlabel('Occupancy (%)')
        ax.legend(loc='upper right', fontsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    fig.suptitle('SM Occupancy: Theoretical vs Achieved', fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'occupancy_comparison.png'))
    plt.close(fig)
    print('  -> occupancy_comparison.png')


# ================================================================
# Main
# ================================================================
if __name__ == '__main__':
    print('Generating charts...')
    chart_ncu_comparison()
    chart_local_memory()
    chart_bfs_pipeline()
    chart_warp_utilization()
    chart_speedup_by_algo()
    chart_e2e_breakdown()
    chart_occupancy()
    print('Done!')
