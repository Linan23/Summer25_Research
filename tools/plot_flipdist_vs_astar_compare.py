#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot aggregated FlipDist vs A* runtime growth from compare_flipdist_vs_astar CSV"
    )
    parser.add_argument("--input", required=True, help="Input comparison CSV")
    parser.add_argument("--output", required=True, help="Output SVG path")
    parser.add_argument("--summary-output", default=None, help="Optional aggregated summary CSV path")
    parser.add_argument("--logy", action="store_true", help="Use log scaling on the y-axis")
    parser.add_argument(
        "--agg",
        choices=["median", "mean", "max"],
        default="median",
        help="Aggregation across seeds for each n",
    )
    parser.add_argument(
        "--title",
        default="FlipDist vs A* Runtime Growth",
        help="Chart title",
    )
    parser.add_argument(
        "--subtitle",
        default="",
        help="Optional subtitle shown below the title",
    )
    return parser.parse_args()


def aggregate(values: list[float], mode: str) -> float:
    if not values:
        return float("nan")
    if mode == "median":
        return float(statistics.median(values))
    if mode == "mean":
        return float(statistics.fmean(values))
    if mode == "max":
        return float(max(values))
    raise ValueError(mode)


def finite_positive(values: list[float]) -> list[float]:
    return [v for v in values if math.isfinite(v) and v > 0.0]


def svg_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def main() -> int:
    cfg = parse_args()
    input_path = Path(cfg.input)
    output_path = Path(cfg.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    by_n: dict[int, dict[str, list[float]]] = {}
    with input_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                n = int(row["n"])
            except Exception:
                continue
            bucket = by_n.setdefault(
                n,
                {
                    "flipdist": [],
                    "astar_simple": [],
                    "astar_combined": [],
                },
            )

            if row.get("flipdist_ok_both", "").lower() in {"true", "1"}:
                try:
                    bucket["flipdist"].append(float(row["flipdist_max_ms"]))
                except Exception:
                    pass

            for src_col, dst_key in (
                ("astar_simple_mean_ms", "astar_simple"),
                ("astar_combined_mean_ms", "astar_combined"),
            ):
                raw = row.get(src_col)
                if raw in (None, ""):
                    continue
                try:
                    value = float(raw)
                except ValueError:
                    continue
                if value > 0.0:
                    bucket[dst_key].append(value)

    ns = sorted(by_n.keys())
    if not ns:
        raise SystemExit("No rows found in input CSV")

    summary_rows: list[dict[str, object]] = []
    flip_vals: list[float] = []
    astar_simple_vals: list[float] = []
    astar_combined_vals: list[float] = []

    for n in ns:
        flip = aggregate(by_n[n]["flipdist"], cfg.agg)
        simple = aggregate(by_n[n]["astar_simple"], cfg.agg)
        combined = aggregate(by_n[n]["astar_combined"], cfg.agg)
        flip_vals.append(flip)
        astar_simple_vals.append(simple)
        astar_combined_vals.append(combined)
        summary_rows.append(
            {
                "n": n,
                "flipdist_ms": f"{flip:.3f}" if math.isfinite(flip) else "",
                "astar_simple_ms": f"{simple:.3f}" if math.isfinite(simple) else "",
                "astar_combined_ms": f"{combined:.3f}" if math.isfinite(combined) else "",
                "flipdist_samples": len(by_n[n]["flipdist"]),
                "astar_simple_samples": len(by_n[n]["astar_simple"]),
                "astar_combined_samples": len(by_n[n]["astar_combined"]),
            }
        )

    if cfg.summary_output:
        summary_path = Path(cfg.summary_output)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "n",
                    "flipdist_ms",
                    "astar_simple_ms",
                    "astar_combined_ms",
                    "flipdist_samples",
                    "astar_simple_samples",
                    "astar_combined_samples",
                ],
            )
            writer.writeheader()
            writer.writerows(summary_rows)

    width = 920
    height = 540
    left = 90
    right = 40
    top = 70
    bottom = 90
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_y = finite_positive(flip_vals) + finite_positive(astar_simple_vals) + finite_positive(astar_combined_vals)
    if not all_y:
        raise SystemExit("No numeric series to plot")

    if cfg.logy:
        transformed = [math.log10(v) for v in all_y]
        y_min = min(transformed)
        y_max = max(transformed)
        if y_min == y_max:
            y_max = y_min + 1.0

        def y_to_px(value: float) -> float:
            if not math.isfinite(value) or value <= 0.0:
                return top + plot_h
            t = math.log10(value)
            frac = (t - y_min) / (y_max - y_min)
            return top + plot_h - frac * plot_h

        tick_vals = []
        low = int(math.floor(y_min))
        high = int(math.ceil(y_max))
        for exp in range(low, high + 1):
            tick_vals.append(10 ** exp)
    else:
        y_min = 0.0
        y_max = max(all_y)
        if y_max <= 0.0:
            y_max = 1.0
        y_max *= 1.1

        def y_to_px(value: float) -> float:
            if not math.isfinite(value):
                return top + plot_h
            frac = (value - y_min) / (y_max - y_min)
            return top + plot_h - frac * plot_h

        tick_vals = [y_max * i / 5.0 for i in range(6)]

    x_min = min(ns)
    x_max = max(ns)
    x_span = max(1, x_max - x_min)

    def x_to_px(n: int) -> float:
        frac = (n - x_min) / x_span
        return left + frac * plot_w

    def line_points(series: list[float]) -> str:
        pts = []
        for n, value in zip(ns, series):
            if math.isfinite(value) and value > 0.0:
                pts.append(f"{x_to_px(n):.2f},{y_to_px(value):.2f}")
        return " ".join(pts)

    def circles(series: list[float], color: str) -> list[str]:
        out = []
        for n, value in zip(ns, series):
            if math.isfinite(value) and value > 0.0:
                x = x_to_px(n)
                y = y_to_px(value)
                out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" fill="{color}" />')
                out.append(
                    f'<text x="{x:.2f}" y="{y - 10:.2f}" font-size="11" text-anchor="middle" fill="{color}">{value:.1f}</text>'
                )
        return out

    grid_lines = []
    y_labels = []
    for tick in tick_vals:
        y = y_to_px(tick)
        label = f"{tick:.0f}" if tick >= 10 else f"{tick:.2f}".rstrip("0").rstrip(".")
        grid_lines.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#d8d8d8" stroke-dasharray="4 4" />'
        )
        y_labels.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" font-size="11" text-anchor="end" fill="#444">{svg_escape(label)}</text>'
        )

    x_labels = []
    for n in ns:
        x = x_to_px(n)
        x_labels.append(f'<line x1="{x:.2f}" y1="{top + plot_h}" x2="{x:.2f}" y2="{top + plot_h + 6}" stroke="#444" />')
        x_labels.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 24}" font-size="12" text-anchor="middle" fill="#222">n={n}</text>'
        )

    legend_x = left + 10
    legend_y = 20
    legend = [
        f'<rect x="{legend_x}" y="{legend_y}" width="360" height="56" rx="6" fill="#ffffff" stroke="#cccccc" />',
        f'<line x1="{legend_x + 12}" y1="{legend_y + 14}" x2="{legend_x + 40}" y2="{legend_y + 14}" stroke="#1f77b4" stroke-width="2.4" />',
        f'<text x="{legend_x + 48}" y="{legend_y + 18}" font-size="12" fill="#1f77b4">FlipDist ({cfg.agg} over seeds, max direction time)</text>',
        f'<line x1="{legend_x + 12}" y1="{legend_y + 30}" x2="{legend_x + 40}" y2="{legend_y + 30}" stroke="#d95f02" stroke-width="2.0" />',
        f'<text x="{legend_x + 48}" y="{legend_y + 34}" font-size="12" fill="#d95f02">A* simple ({cfg.agg} of aggregated runtime values)</text>',
        f'<line x1="{legend_x + 12}" y1="{legend_y + 46}" x2="{legend_x + 40}" y2="{legend_y + 46}" stroke="#2ca02c" stroke-width="2.0" />',
        f'<text x="{legend_x + 48}" y="{legend_y + 50}" font-size="12" fill="#2ca02c">A* combined ({cfg.agg} of aggregated runtime values)</text>',
    ]

    subtitle = ""
    if cfg.subtitle:
        subtitle = (
            f'<text x="{width / 2:.1f}" y="36" font-size="12" text-anchor="middle" fill="#555">'
            f'{svg_escape(cfg.subtitle)}</text>'
        )

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="{width}" height="{height}" fill="#ffffff" />
<text x="{width / 2:.1f}" y="18" font-size="18" text-anchor="middle" fill="#111">{svg_escape(cfg.title)}</text>
{subtitle}
{''.join(legend)}
<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#fafafa" stroke="#444" />
{''.join(grid_lines)}
{''.join(y_labels)}
{''.join(x_labels)}
<polyline fill="none" stroke="#1f77b4" stroke-width="2.4" points="{line_points(flip_vals)}" />
<polyline fill="none" stroke="#d95f02" stroke-width="2.0" points="{line_points(astar_simple_vals)}" />
<polyline fill="none" stroke="#2ca02c" stroke-width="2.0" points="{line_points(astar_combined_vals)}" />
{''.join(circles(flip_vals, "#1f77b4"))}
{''.join(circles(astar_simple_vals, "#d95f02"))}
{''.join(circles(astar_combined_vals, "#2ca02c"))}
<text x="{left + plot_w / 2:.2f}" y="{height - 22}" font-size="13" text-anchor="middle" fill="#222">Problem size n</text>
<text x="22" y="{top + plot_h / 2:.2f}" font-size="13" text-anchor="middle" fill="#222" transform="rotate(-90 22,{top + plot_h / 2:.2f})">Runtime (ms)</text>
</svg>
'''

    output_path.write_text(svg, encoding="utf-8")
    print(f"Wrote {output_path}")
    if cfg.summary_output:
        print(f"Wrote {cfg.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
