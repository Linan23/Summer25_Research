#!/usr/bin/env python3
"""Sweep parity checks between FlipDist and Java BFS across n/seeds.

This wraps run_flipdist_java_parity.py and aggregates outputs into a single CSV,
while enforcing the Java safety cap (n<=15 by default).
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Sweep FlipDist vs Java BFS parity across n/seeds.")
    p.add_argument("--case", dest="case_type", choices=["comb", "random"], default="random")
    p.add_argument("--n-min", type=int, required=True)
    p.add_argument("--n-max", type=int, required=True)
    p.add_argument("--seed-min", type=int, default=0)
    p.add_argument("--seed-max", type=int, default=0)
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--cpp-binary", default="./build/flipdist")
    p.add_argument("--max-k", type=int, default=None, help="Max k for FlipDistMinK (defaults to 3n+10).")
    p.add_argument("--bfs-cap", type=int, default=1)
    p.add_argument("--java-binary", default="java")
    p.add_argument("--java-out", default="triangulation/out")
    p.add_argument("--java-lib", default="triangulation/lib/acm.jar")
    p.add_argument("--java-time-limit", type=float, default=120.0)
    p.add_argument("--java-visited-cap", type=int, default=15_000_000)
    p.add_argument("--java-queue-cap", type=int, default=15_000_000)
    p.add_argument("--max-java-n", type=int, default=15)
    p.add_argument("--allow-java-n-over", action="store_true")
    p.add_argument(
        "--output",
        default=None,
        help="CSV output path. If omitted, the script asks whether to save a CSV.",
    )
    p.add_argument("--print", action="store_true", help="Print each parity run summary to stdout.")
    return p.parse_args()


def resolve_output_path(cfg) -> Path | None:
    if cfg.output:
        return Path(cfg.output)
    if not sys.stdin.isatty():
        return None

    default_path = (
        f"results/parity_flipdist_vs_java_{cfg.case_type}_n{cfg.n_min}_{cfg.n_max}_"
        f"seeds{cfg.seed_min}_{cfg.seed_max}.csv"
    )
    try:
        resp = input("Save CSV results? (y/n): ").strip().lower()
    except EOFError:
        return None
    if resp not in {"y", "yes"}:
        return None
    try:
        path = input(f"CSV path [default: {default_path}]: ").strip()
    except EOFError:
        return Path(default_path)
    return Path(path or default_path)


def run_one(cfg, n: int, seed: int, out_csv: Path) -> int:
    max_k = cfg.max_k if cfg.max_k is not None else max(1, 3 * n + 10)
    cmd = [
        sys.executable,
        "scripts/run_flipdist_java_parity.py",
        "--case",
        cfg.case_type,
        "--n",
        str(n),
        "--count",
        str(cfg.count),
        "--seed",
        str(seed),
        "--cpp-binary",
        cfg.cpp_binary,
        "--max-k",
        str(max_k),
        "--bfs-cap",
        str(cfg.bfs_cap),
        "--java-binary",
        cfg.java_binary,
        "--java-out",
        cfg.java_out,
        "--java-lib",
        cfg.java_lib,
        "--java-time-limit",
        str(cfg.java_time_limit),
        "--java-visited-cap",
        str(cfg.java_visited_cap),
        "--java-queue-cap",
        str(cfg.java_queue_cap),
        "--max-java-n",
        str(cfg.max_java_n),
    ]
    if cfg.allow_java_n_over:
        cmd.append("--allow-java-n-over")
    cmd.extend(["--output", str(out_csv)])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if cfg.print:
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
    return proc.returncode


def main() -> int:
    cfg = parse_args()
    if cfg.n_max < cfg.n_min:
        raise SystemExit("--n-max must be >= --n-min")
    if cfg.seed_max < cfg.seed_min:
        raise SystemExit("--seed-max must be >= --seed-min")

    out_path = resolve_output_path(cfg)
    out_dir = out_path.parent if out_path else Path("results")
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    failures = 0
    for n in range(cfg.n_min, cfg.n_max + 1):
        seeds = range(cfg.seed_min, cfg.seed_max + 1) if cfg.case_type == "random" else [-1]
        for seed in seeds:
            sys.stderr.write(f"Parity sweep: n={n} seed={seed} case={cfg.case_type}\n")
            tmp_csv = out_dir / f".parity_tmp_n{n}_s{seed}.csv"
            rc = run_one(cfg, n, seed, tmp_csv)
            if rc != 0:
                failures += 1
                continue
            if not tmp_csv.exists():
                failures += 1
                continue
            with tmp_csv.open() as f:
                reader = csv.DictReader(f)
                all_rows.extend(list(reader))
            tmp_csv.unlink(missing_ok=True)

    if not all_rows:
        sys.stderr.write("No parity rows produced.\n")
        return 1

    if out_path:
        with out_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_rows[0].keys())
            writer.writeheader()
            writer.writerows(all_rows)

    mismatches = [
        row
        for row in all_rows
        if row.get("distance_java") is not None and row.get("distance_flipdist") != row.get("distance_java")
    ]
    if mismatches:
        sys.stderr.write(f"WARNING: {len(mismatches)} distance mismatches detected\n")
        return 2

    if failures:
        sys.stderr.write(f"Completed with {failures} failed runs (timeouts or errors).\n")
        return 3

    if out_path:
        sys.stderr.write(f"All distances match. Wrote {out_path} ({len(all_rows)} rows).\n")
    else:
        sys.stderr.write(f"All distances match. CSV not saved ({len(all_rows)} rows).\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
