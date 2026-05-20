#!/usr/bin/env python3
"""Plot time growth for FlipDist vs Java from parity CSV."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import math


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot time growth from parity CSV.")
    p.add_argument("--input", required=True, help="Parity CSV path")
    p.add_argument("--output", required=True, help="Output PNG path")
    p.add_argument("--bruteforce-input", default=None, help="Optional brute-force CSV with time_ms_bruteforce column")
    p.add_argument("--agg", choices=["max", "mean", "median"], default="max")
    p.add_argument("--per-point", action="store_true", help="Plot each (n,seed) row instead of aggregating per n.")
    p.add_argument("--label-values", action="store_true", help="Label points with time_ms values (per-point mode).")
    p.add_argument("--small-threshold-ms", type=float, default=10.0, help="Treat times below this as small for emphasis.")
    p.add_argument("--short-labels", action="store_true", help="Use compact x labels like n12s3.")
    p.add_argument("--sort-by-time", action="store_true", help="Sort per-point rows by time (min->max).")
    p.add_argument("--zoom-quantile", type=float, default=1.0, help="Zoom y-axis to this quantile of times (linear only).")
    p.add_argument("--logy", action="store_true", help="Use log scale for y")
    p.add_argument("--title", default="FlipDist vs Java Time Growth")
    return p.parse_args()


def aggregate(values: list[float], mode: str) -> float:
    if not values:
        return float("nan")
    if mode == "max":
        return max(values)
    if mode == "mean":
        return sum(values) / len(values)
    if mode == "median":
        vals = sorted(values)
        mid = len(vals) // 2
        return vals[mid] if len(vals) % 2 == 1 else (vals[mid - 1] + vals[mid]) / 2.0
    raise ValueError(mode)


def main() -> int:
    cfg = parse_args()

    per_n_flip: dict[int, list[float]] = defaultdict(list)
    per_n_java: dict[int, list[float]] = defaultdict(list)
    per_n_brute: dict[int, list[float]] = defaultdict(list)
    per_point: list[tuple[int, int, float, float, float]] = []

    brute_map: dict[tuple[int, int, str], float] = {}
    if cfg.bruteforce_input:
        with open(cfg.bruteforce_input, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    n = int(row["n"])
                    seed = int(row["seed"])
                    direction = row.get("direction", "")
                    tb = row.get("time_ms_bruteforce")
                    if tb in (None, ""):
                        continue
                    brute_map[(n, seed, direction)] = float(tb)
                    per_n_brute[n].append(float(tb))
                except Exception:
                    continue

    with open(cfg.input, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                n = int(row["n"])
                seed = int(row["seed"])
            except Exception:
                continue
            # time_ms_flipdist/time_ms_java can be empty on timeouts
            tf = row.get("time_ms_flipdist")
            tj = row.get("time_ms_java")
            direction = row.get("direction", "")
            if tf not in (None, ""):
                per_n_flip[n].append(float(tf))
            if tj not in (None, ""):
                per_n_java[n].append(float(tj))
            tb = brute_map.get((n, seed, direction), float("nan"))
            if tb == tb:
                per_n_brute[n].append(tb)
            if tf not in (None, "") or tj not in (None, "") or tb == tb:
                per_point.append(
                    (n, seed,
                     float(tf) if tf not in (None, "") else float("nan"),
                     float(tj) if tj not in (None, "") else float("nan"),
                     float(tb))
                )

    ns = sorted(set(per_n_flip) | set(per_n_java))
    if not ns and not per_point:
        raise SystemExit("No rows to plot.")

    flip_y = [aggregate(per_n_flip.get(n, []), cfg.agg) for n in ns]
    java_y = [aggregate(per_n_java.get(n, []), cfg.agg) for n in ns]
    brute_y = [aggregate(per_n_brute.get(n, []), cfg.agg) for n in ns]

    # Build a short description of seeds per n (if present).
    seed_sets: dict[int, set[int]] = defaultdict(set)
    with open(cfg.input, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                n = int(row["n"])
                seed = int(row["seed"])
            except Exception:
                continue
            seed_sets[n].add(seed)

    def seed_desc(n: int) -> str:
        seeds = sorted(seed_sets.get(n, []))
        if not seeds:
            return "seeds=?"
        if len(seeds) == 1:
            return f"seed={seeds[0]}"
        return f"seeds={seeds[0]}..{seeds[-1]} ({len(seeds)})"

    subtitle = " | ".join(f"n={n} {seed_desc(n)}" for n in ns)

    if cfg.sort_by_time and per_point:
        def sort_key(row: tuple[int, int, float, float, float]) -> float:
            _, _, tf, tj, tb = row
            vals = [v for v in (tf, tj, tb) if v == v]
            return max(vals) if vals else float("inf")
        per_point = sorted(per_point, key=sort_key)

    def quantile(vals: list[float], q: float) -> float | None:
        clean = sorted(v for v in vals if v == v and v > 0)
        if not clean:
            return None
        if q <= 0:
            return clean[0]
        if q >= 1:
            return clean[-1]
        idx = int((len(clean) - 1) * q)
        return clean[idx]

    out_path = Path(cfg.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Lazy import so script can be inspected without matplotlib.
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 5))
        if cfg.per_point:
            if cfg.short_labels:
                xlabels = [f"n{n}s{seed}" for n, seed, _, _, _ in per_point]
            else:
                xlabels = [f"n={n} s={seed}" for n, seed, _, _, _ in per_point]
            x = list(range(len(per_point)))
            flip_vals = [tf for _, _, tf, _, _ in per_point]
            java_vals = [tj for _, _, _, tj, _ in per_point]
            brute_vals = [tb for _, _, _, _, tb in per_point]
            flip_sizes = [40 if (v == v and v < cfg.small_threshold_ms) else 26 for v in flip_vals]
            java_sizes = [28 if (v == v and v < cfg.small_threshold_ms) else 18 for v in java_vals]
            brute_sizes = [28 if (v == v and v < cfg.small_threshold_ms) else 18 for v in brute_vals]
            x_flip = [xi + 0.1 for xi in x]
            x_java = [xi - 0.1 for xi in x]
            x_brute = x
            # Trajectory lines + points
            plt.plot(x_flip, flip_vals, color="#1f77b4", linewidth=1.6, alpha=0.9, zorder=2)
            plt.scatter(x_flip, flip_vals, label="FlipDist", s=flip_sizes, clip_on=True, zorder=3)
            plt.plot(x_java, java_vals, color="#ff7f0e", linewidth=1.4, alpha=0.8, zorder=1)
            plt.scatter(x_java, java_vals, label="Java BFS", s=java_sizes, clip_on=True, zorder=2)
            if any(v == v for v in brute_vals):
                plt.plot(x_brute, brute_vals, color="#2ca02c", linewidth=1.4, alpha=0.8, zorder=1)
                plt.scatter(x_brute, brute_vals, label="Brute Force", s=brute_sizes, clip_on=True, zorder=2)
            # draw a faint line between the solvers for direct comparison
            for xi, tf, tj, tb in zip(x, flip_vals, java_vals, brute_vals):
                pts = [v for v in (tf, tj, tb) if v == v]
                if len(pts) >= 2:
                    plt.plot([xi, xi], [min(pts), max(pts)], color="#999", linewidth=0.6, alpha=0.5, zorder=0)
            plt.xticks(x, xlabels, rotation=60, ha="right", fontsize=8)
            if cfg.label_values:
                vals = [v for v in flip_vals + java_vals + brute_vals if v == v]
                span = max(vals) - min(vals) if vals else 1.0
                if span <= 0:
                    span = max(vals) if vals else 1.0
                for idx, ((n, seed, tf, tj, tb), xi) in enumerate(zip(per_point, x)):
                    y_jitter = ((idx % 3) - 1) * 0.02 * span
                    if tf == tf and tf > 0:
                        fsize = 9 if tf < cfg.small_threshold_ms else 7
                        label = f"n{n}s{seed} {tf:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tf:.0f} ms"
                        plt.text(xi + 0.16, tf + y_jitter, label, fontsize=fsize, ha="left", va="bottom", clip_on=True)
                    if tj == tj and tj > 0:
                        fsize = 9 if tj < cfg.small_threshold_ms else 7
                        label = f"n{n}s{seed} {tj:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tj:.0f} ms"
                        plt.text(xi - 0.16, tj - y_jitter, label, fontsize=fsize, ha="right", va="top", clip_on=True)
                    if tb == tb and tb > 0:
                        fsize = 9 if tb < cfg.small_threshold_ms else 7
                        label = f"n{n}s{seed} {tb:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tb:.0f} ms"
                        plt.text(xi, tb + y_jitter * 0.5, label, fontsize=fsize, ha="center", va="bottom", clip_on=True)
        else:
            plt.plot(ns, flip_y, marker="o", label="FlipDist")
            plt.plot(ns, java_y, marker="o", label="Java BFS")
            if any(v == v for v in brute_y):
                plt.plot(ns, brute_y, marker="o", label="Brute Force")
        plt.title(cfg.title)
        if not cfg.per_point:
            plt.suptitle(subtitle, fontsize=9, y=0.99)
        plt.xlabel("n, seed" if cfg.per_point else "n")
        plt.ylabel(f"time_ms ({cfg.agg})")
        if cfg.logy:
            plt.yscale("log")
        plt.grid(True, linestyle="--", alpha=0.4)
        if not cfg.logy and cfg.zoom_quantile and 0 < cfg.zoom_quantile < 1:
            all_times = []
            if cfg.per_point:
                all_times.extend([v for v in flip_vals if v == v])
                all_times.extend([v for v in java_vals if v == v])
                all_times.extend([v for v in brute_vals if v == v])
            else:
                all_times.extend([v for v in flip_y if v == v])
                all_times.extend([v for v in java_y if v == v])
                all_times.extend([v for v in brute_y if v == v])
            qv = quantile(all_times, cfg.zoom_quantile)
            if qv is not None and qv > 0:
                plt.ylim(bottom=0, top=qv * 1.15)
        plt.legend(loc="upper left", bbox_to_anchor=(0.01, 0.99), frameon=True, framealpha=0.6, fontsize=9)
        if cfg.per_point:
            plt.subplots_adjust(bottom=0.33, left=0.12, right=0.96, top=0.9)
        else:
            plt.subplots_adjust(left=0.12, right=0.96, top=0.9)
        plt.savefig(out_path)
        return 0
    except ModuleNotFoundError:
        # Simple SVG fallback (no external deps).
        width, height = 800, 500
        margin = 60
        plot_w = width - 2 * margin
        plot_h = height - 2 * margin

        def safe_log(v: float) -> float:
            return math.log10(v) if v > 0 else float("-inf")

        ys = flip_y + java_y + brute_y
        if cfg.logy:
            ys = [safe_log(v) for v in ys if v > 0]
        else:
            ys = [v for v in ys if v == v]
        if not cfg.logy and cfg.zoom_quantile and 0 < cfg.zoom_quantile < 1:
            qv = quantile(ys, cfg.zoom_quantile)
            if qv is not None and qv > 0:
                ys = [v for v in ys if v <= qv * 1.15]
        y_min = min(ys) if ys else 0.0
        y_max = max(ys) if ys else 1.0
        if y_min == y_max:
            y_max = y_min + 1.0

        x_min = min(ns)
        x_max = max(ns)
        if x_min == x_max:
            x_max = x_min + 1

        def x_map(x: int) -> float:
            return margin + (x - x_min) / (x_max - x_min) * plot_w

        def y_map(y: float) -> float:
            yv = safe_log(y) if cfg.logy else y
            return margin + (1.0 - (yv - y_min) / (y_max - y_min)) * plot_h

        def points_for(values: list[float]) -> str:
            pts = []
            for x, y in zip(ns, values):
                if y != y or y <= 0:
                    continue
                pts.append(f"{x_map(x):.2f},{y_map(y):.2f}")
            return " ".join(pts)

        def points_for_per_point(values: list[float], x_offset: float = 0.0) -> str:
            pts = []
            denom = max(1, len(values) - 1)
            for idx, y in enumerate(values):
                if y != y or y <= 0:
                    continue
                x = margin + (idx / denom) * plot_w + x_offset
                pts.append(f"{x:.2f},{y_map(y):.2f}")
            return " ".join(pts)

        if cfg.per_point:
            flip_pts = points_for_per_point([tf for _, _, tf, _, _ in per_point], x_offset=8)
            java_pts = points_for_per_point([tj for _, _, _, tj, _ in per_point], x_offset=-8)
            brute_pts = points_for_per_point([tb for _, _, _, _, tb in per_point], x_offset=0)
        else:
            flip_pts = points_for(flip_y)
            java_pts = points_for(java_y)
            brute_pts = points_for(brute_y)

        clip_id = "plot-clip"
        svg = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
            f'<rect x="0" y="0" width="{width}" height="{height}" fill="white"/>',
            f'<text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="16">{cfg.title}</text>',
            f'<text x="{width/2:.1f}" y="42" text-anchor="middle" font-size="10">{subtitle}</text>',
            f'<text x="{width/2:.1f}" y="{height-10}" text-anchor="middle" font-size="12">{"n, seed" if cfg.per_point else "n"}</text>',
            f'<text x="15" y="{height/2:.1f}" text-anchor="middle" font-size="12" transform="rotate(-90 15 {height/2:.1f})">time_ms ({cfg.agg}){" log10" if cfg.logy else ""}</text>',
            f'<defs><clipPath id="{clip_id}"><rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}"/></clipPath></defs>',
            f'<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#444"/>',
        ]

        if flip_pts:
            svg.append(f'<polyline points="{flip_pts}" fill="none" stroke="#1f77b4" stroke-width="3" clip-path="url(#{clip_id})"/>')
        if java_pts:
            svg.append(f'<polyline points="{java_pts}" fill="none" stroke="#ff7f0e" stroke-width="2" clip-path="url(#{clip_id})"/>')
        if brute_pts:
            svg.append(f'<polyline points="{brute_pts}" fill="none" stroke="#2ca02c" stroke-width="2" clip-path="url(#{clip_id})"/>')

        if cfg.per_point and cfg.label_values:
            denom = max(1, len(per_point) - 1)
            span = y_max - y_min if y_max > y_min else 1.0
            for idx, (n, seed, tf, tj, tb) in enumerate(per_point):
                x = margin + (idx / denom) * plot_w
                y_jitter = ((idx % 3) - 1) * 0.02 * span
                if tf == tf and tf > 0:
                    fsize = 10 if tf < cfg.small_threshold_ms else 8
                    label = f"n{n}s{seed} {tf:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tf:.0f} ms"
                    svg.append(f'<text x="{x+10:.2f}" y="{y_map(tf + y_jitter)-2:.2f}" font-size="{fsize}" text-anchor="start" clip-path="url(#{clip_id})">{label}</text>')
                if tj == tj and tj > 0:
                    fsize = 10 if tj < cfg.small_threshold_ms else 8
                    label = f"n{n}s{seed} {tj:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tj:.0f} ms"
                    svg.append(f'<text x="{x-10:.2f}" y="{y_map(tj - y_jitter)+10:.2f}" font-size="{fsize}" text-anchor="end" clip-path="url(#{clip_id})">{label}</text>')
                if tb == tb and tb > 0:
                    fsize = 10 if tb < cfg.small_threshold_ms else 8
                    label = f"n{n}s{seed} {tb:.0f}ms" if cfg.short_labels else f"n {n}  seed {seed}  {tb:.0f} ms"
                    svg.append(f'<text x="{x:.2f}" y="{y_map(tb + y_jitter * 0.5)-2:.2f}" font-size="{fsize}" text-anchor="middle" clip-path="url(#{clip_id})">{label}</text>')

        # Legend (place outside the plot, left side)
        legend_x = margin + 10
        legend_y = margin
        svg.append(f'<rect x="{legend_x}" y="{legend_y}" width="150" height="64" fill="white" stroke="#ccc"/>')
        svg.append(f'<line x1="{legend_x+10}" y1="{legend_y+16}" x2="{legend_x+40}" y2="{legend_y+16}" stroke="#1f77b4" stroke-width="2"/>')
        svg.append(f'<text x="{legend_x+50}" y="{legend_y+20}" font-size="12">FlipDist</text>')
        svg.append(f'<line x1="{legend_x+10}" y1="{legend_y+34}" x2="{legend_x+40}" y2="{legend_y+34}" stroke="#ff7f0e" stroke-width="2"/>')
        svg.append(f'<text x="{legend_x+50}" y="{legend_y+38}" font-size="12">Java BFS</text>')
        svg.append(f'<line x1="{legend_x+10}" y1="{legend_y+52}" x2="{legend_x+40}" y2="{legend_y+52}" stroke="#2ca02c" stroke-width="2"/>')
        svg.append(f'<text x="{legend_x+50}" y="{legend_y+56}" font-size="12">Brute Force</text>')

        svg.append("</svg>")
        out_path.write_text("\\n".join(svg))
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
