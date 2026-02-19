#!/usr/bin/env python3
"""Parity harness: C++ FlipDist vs Java BFS (triangulation oracle).

This runs the C++ `flipdist_fast` CLI to generate random/comb tree pairs and
computes the exact flip distance using the Java BFS implementation
(`TriangulationMetricsCli`). It then reports any distance mismatches.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str], *, input_text: str | None = None) -> str:
    proc = subprocess.run(
        cmd,
        input=input_text,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc.stdout


def run_cpp(cfg) -> list[dict]:
    cmd = [
        cfg.cpp_binary,
        "--case",
        cfg.case_type,
        "--n",
        str(cfg.n),
        "--count",
        str(cfg.count),
        "--bfs-cap",
        str(cfg.bfs_cap),
        "--max-k",
        str(cfg.max_k),
    ]
    if cfg.case_type == "random":
        cmd.extend(["--seed", str(cfg.seed)])
    out = run_cmd(cmd)
    rows = [json.loads(line) for line in out.splitlines() if line.strip()]
    return rows


def run_java(cfg, forward_rows: list[dict]) -> list[dict]:
    if not forward_rows:
        return []

    java_classpath = f"{cfg.java_out}:{cfg.java_lib}"
    cmd = [
        cfg.java_binary,
        "-cp",
        java_classpath,
        "TriangulationMetricsCli",
    ]

    payload_lines: list[str] = []
    for row in forward_rows:
        payload = {
            "case_type": row.get("case_type", "unknown"),
            "seed": row.get("seed", -1),
            "tree_a": row["tree_a"],
            "tree_b": row["tree_b"],
            "time_limit": cfg.java_time_limit,
            "visited_cap": cfg.java_visited_cap,
            "queue_cap": cfg.java_queue_cap,
        }
        payload_lines.append(json.dumps(payload))

    input_data = "\n".join(payload_lines) + "\n"
    out = run_cmd(cmd, input_text=input_data)
    rows = [json.loads(line) for line in out.splitlines() if line.strip()]
    return rows


def merge(cpp_rows: list[dict], java_rows: list[dict]) -> list[dict]:
    java_map: dict[tuple, dict] = {}
    for row in java_rows:
        key = (
            row.get("case_type"),
            row.get("seed"),
            row.get("tree_a"),
            row.get("tree_b"),
            row.get("direction"),
        )
        java_map[key] = row

    merged: list[dict] = []
    for row in cpp_rows:
        key = (
            row.get("case_type"),
            row.get("seed"),
            row.get("tree_a"),
            row.get("tree_b"),
            row.get("direction"),
        )
        merged.append(
            {
                "case_type": row.get("case_type"),
                "n": row.get("n"),
                "seed": row.get("seed"),
                "direction": row.get("direction"),
                "distance_flipdist": row.get("distance_flipdist", row.get("distance")),
                "status_flipdist": row.get("status_flipdist", row.get("status")),
                "time_ms_flipdist": row.get("time_ms_flipdist", row.get("time_ms")),
                "distance_java": java_map.get(key, {}).get("distance"),
                "status_java": java_map.get(key, {}).get("status"),
                "time_ms_java": java_map.get(key, {}).get("time_ms"),
                "tree_a": row.get("tree_a"),
                "tree_b": row.get("tree_b"),
            }
        )
    return merged


def write_csv(path: str, rows: list[dict]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_type",
        "n",
        "seed",
        "direction",
        "distance_flipdist",
        "distance_java",
        "status_flipdist",
        "status_java",
        "time_ms_flipdist",
        "time_ms_java",
        "tree_a",
        "tree_b",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def maybe_prompt_output(cfg) -> None:
    if cfg.output is not None:
        return
    if not sys.stdin.isatty():
        return
    default_path = f"results/parity_flipdist_vs_java_{cfg.case_type}_n{cfg.n}_seed{cfg.seed}.csv"
    try:
        resp = input("Save CSV? (y/n): ").strip().lower()
    except EOFError:
        return
    if resp not in {"y", "yes"}:
        return
    try:
        path = input(f"CSV path [default: {default_path}]: ").strip()
    except EOFError:
        return
    cfg.output = path or default_path


def parse_args():
    p = argparse.ArgumentParser(description="Verify C++ flipdist against Java BFS oracle.")
    p.add_argument("--case", dest="case_type", choices=["comb", "random"], default="random")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--cpp-binary", default="./build/flipdist")
    p.add_argument("--max-k", type=int, default=None, help="Max k for FlipDistMinK")
    p.add_argument("--bfs-cap", type=int, default=1, help="Keep small (we use Java as oracle)")
    p.add_argument("--java-binary", default="java")
    p.add_argument("--java-out", default="triangulation/out")
    p.add_argument("--java-lib", default="triangulation/lib/acm.jar")
    p.add_argument("--java-time-limit", type=float, default=120.0)
    p.add_argument("--java-visited-cap", type=int, default=15_000_000)
    p.add_argument("--java-queue-cap", type=int, default=15_000_000)
    p.add_argument(
        "--max-java-n",
        type=int,
        default=15,
        help="Safety cap for Java BFS (default: 15).",
    )
    p.add_argument(
        "--allow-java-n-over",
        action="store_true",
        help="Allow n > --max-java-n (may be slow or halt).",
    )
    p.add_argument("--output", default=None, help="CSV output path; omit to skip writing a file.")
    p.add_argument("--print", action="store_true", help="Print merged comparison rows to stdout.")
    return p.parse_args()


def main() -> int:
    cfg = parse_args()
    if cfg.n > cfg.max_java_n and not cfg.allow_java_n_over:
        raise SystemExit(
            f"Refusing to run Java BFS at n={cfg.n} (max {cfg.max_java_n}). "
            "Pass --allow-java-n-over to override."
        )
    if cfg.max_k is None:
        cfg.max_k = max(1, 3 * cfg.n + 10)

    sys.stderr.write(f"Parity run: n={cfg.n} seed={cfg.seed} case={cfg.case_type}\n")

    cpp_rows = run_cpp(cfg)
    forward_rows = [row for row in cpp_rows if row.get("direction") == "a->b"]
    java_rows = run_java(cfg, forward_rows)
    merged = merge(cpp_rows, java_rows)
    maybe_prompt_output(cfg)
    if cfg.output:
        write_csv(cfg.output, merged)

    if cfg.print or not cfg.output:
        # Pretty-print merged rows to stdout for terminal use.
        cols = ["n", "seed", "direction", "distance_flipdist", "distance_java",
                "status_flipdist", "status_java", "time_ms_flipdist", "time_ms_java"]
        header = " | ".join(f"{c:>10}" for c in cols)
        print(header)
        print("-" * len(header))
        for row in merged:
            print(" | ".join(f"{str(row.get(c,''))[:10]:>10}" for c in cols))

    mismatches = [
        row
        for row in merged
        if row["distance_java"] is not None and row["distance_flipdist"] != row["distance_java"]
    ]
    if merged:
        times_flip = [float(r["time_ms_flipdist"]) for r in merged if r.get("time_ms_flipdist") not in ("", None)]
        times_java = [float(r["time_ms_java"]) for r in merged if r.get("time_ms_java") not in ("", None)]
        max_flip = max(times_flip) if times_flip else float("nan")
        max_java = max(times_java) if times_java else float("nan")
        sys.stderr.write(
            f"Max times: flipdist={max_flip:.3f}ms java={max_java:.3f}ms "
            f"(rows={len(merged)})\n"
        )
    if mismatches:
        sys.stderr.write(f"WARNING: {len(mismatches)} distance mismatches detected\n")
        return 2
    sys.stderr.write("All distances match between FlipDist and Java BFS.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
