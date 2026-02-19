#!/usr/bin/env python3
"""Parity harness: brute-force bf_bst vs Java BFS (triangulation oracle)."""
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
        "--seed",
        str(cfg.seed),
        "--count",
        str(cfg.count),
    ]
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
                "distance_bruteforce": row.get("distance"),
                "status_bruteforce": row.get("status"),
                "time_ms_bruteforce": row.get("time_ms"),
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
        "distance_bruteforce",
        "distance_java",
        "status_bruteforce",
        "status_java",
        "time_ms_bruteforce",
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
    default_path = f"results/bruteforce_vs_java_{cfg.case_type}_n{cfg.n}_seed{cfg.seed}.csv"
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


def print_table(rows: list[dict]) -> None:
    cols = [
        "n",
        "seed",
        "direction",
        "distance_bruteforce",
        "distance_java",
        "time_ms_bruteforce",
        "time_ms_java",
        "status_bruteforce",
        "status_java",
    ]
    header = " | ".join(f"{c:>12}" for c in cols)
    print(header)
    print("-" * len(header))
    for row in rows:
        print(" | ".join(f"{str(row.get(c, ''))[:12]:>12}" for c in cols))


def parse_args():
    p = argparse.ArgumentParser(description="Compare bf_bst brute-force vs Java BFS.")
    p.add_argument("--case", dest="case_type", choices=["random"], default="random")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--count", type=int, default=1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--cpp-binary", default="./build/bf_bst")
    p.add_argument("--java-binary", default="java")
    p.add_argument("--java-out", default="triangulation/out")
    p.add_argument("--java-lib", default="triangulation/lib/acm.jar")
    p.add_argument("--java-time-limit", type=float, default=120.0)
    p.add_argument("--java-visited-cap", type=int, default=15_000_000)
    p.add_argument("--java-queue-cap", type=int, default=15_000_000)
    p.add_argument("--max-java-n", type=int, default=15)
    p.add_argument("--allow-java-n-over", action="store_true")
    p.add_argument("--output", default=None, help="CSV output path; omit to skip writing a file.")
    return p.parse_args()


def main() -> int:
    cfg = parse_args()
    if cfg.n > cfg.max_java_n and not cfg.allow_java_n_over:
        raise SystemExit(
            f"Refusing to run Java BFS at n={cfg.n} (max {cfg.max_java_n}). "
            "Pass --allow-java-n-over to override."
        )

    cpp_rows = run_cpp(cfg)
    forward_rows = [row for row in cpp_rows if row.get("direction") == "a->b"]
    java_rows = run_java(cfg, forward_rows)
    merged = merge(cpp_rows, java_rows)

    print_table(merged)
    maybe_prompt_output(cfg)
    if cfg.output:
        write_csv(cfg.output, merged)

    mismatches = [
        row
        for row in merged
        if row["distance_java"] is not None and row["distance_bruteforce"] != row["distance_java"]
    ]
    if mismatches:
        sys.stderr.write(f"WARNING: {len(mismatches)} distance mismatches detected\n")
        return 2

    sys.stderr.write("All distances match between brute-force and Java BFS.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
