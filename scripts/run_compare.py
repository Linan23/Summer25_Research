#!/usr/bin/env python3
"""Comparison harness for C++ vs Java rotation solvers.

The script orchestrates parity runs between `./test_asan` and the Java
Triangulation program. It accepts the same knobs used by the C++ CLI (case
type, size, limits) and optionally replays timeouts with the hashed
bidirectional solver so both sides report a distance. Passing `--no-auto-bidir`
captures the raw BFS behaviour instead. CSV output includes solver metadata,
timings, queue statistics, and canonical tree strings so any mismatch is fully
reproducible.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace


def run_cpp_cli(cfg):
    time_limit = cfg.time_limit
    if time_limit is None and getattr(cfg, "case_type", None) == "random" \
            and not getattr(cfg, "use_bidir", False) \
            and not getattr(cfg, "prefer_bidir", False):
        time_limit = 0.5
    cmd = [
        cfg.cpp_binary,
        "--case", cfg.case_type,
        "--n", str(cfg.n),
        "--count", str(cfg.count),
        "--program", cfg.cpp_program,
    ]
    if cfg.case_type == "random":
        cmd.extend(["--seed", str(cfg.seed)])
    if time_limit is not None:
        cmd.extend(["--time-limit", str(time_limit)])
    if cfg.visited_cap is not None:
        cmd.extend(["--visited-cap", str(cfg.visited_cap)])
    if cfg.queue_cap is not None:
        cmd.extend(["--queue-cap", str(cfg.queue_cap)])
    if cfg.use_bidir:
        cmd.append("--fallback-bidir")
    if cfg.bidir_cap is not None:
        cmd.extend(["--bidir-cap", str(cfg.bidir_cap)])
    if getattr(cfg, "prefer_bidir", False):
        cmd.append("--prefer-bidir")

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"C++ compare_cli failed with code {result.returncode}")

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    rows = [json.loads(line) for line in lines]
    return rows


def run_java_cli(cfg, forward_rows):
    if not forward_rows:
        return []

    java_classpath = f"{cfg.java_out}:{cfg.java_lib}"
    cmd = [
        cfg.java_binary,
        "-cp", java_classpath,
        "TriangulationMetricsCli",
    ]

    payload_lines = []
    for row in forward_rows:
        payload = {
            "case_type": row.get("case_type", "unknown"),
            "tree_a": row["tree_a"],
            "tree_b": row["tree_b"],
            "seed": row.get("seed", -1),
        }
        if cfg.time_limit is not None:
            payload["time_limit"] = cfg.time_limit
        if cfg.visited_cap is not None:
            payload["visited_cap"] = cfg.visited_cap
        if cfg.queue_cap is not None:
            payload["queue_cap"] = cfg.queue_cap
        payload_lines.append(json.dumps(payload))

    input_data = "\n".join(payload_lines) + "\n"
    result = subprocess.run(cmd, input=input_data, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"Java TriangulationMetricsCli failed with code {result.returncode}")

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    rows = [json.loads(line) for line in lines]
    return rows

def merge_results(cpp_rows, java_rows):
    java_map = {}
    for row in java_rows:
        key = (
            row.get("case_type", "unknown"),
            row.get("seed", -1),
            row.get("tree_a"),
            row.get("tree_b"),
            row.get("direction"),
        )
        java_map[key] = row

    merged = []
    for cpp_row in cpp_rows:
        key = (
            cpp_row.get("case_type", "unknown"),
            cpp_row.get("seed", -1),
            cpp_row.get("tree_a"),
            cpp_row.get("tree_b"),
            cpp_row.get("direction"),
        )
        java_row = java_map.get(key)
        merged.append((cpp_row, java_row))
    return merged


def write_csv(cfg, merged_rows):
    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    with open(cfg.output, "w", newline="") as fh:
        fieldnames = [
            "case_type",
            "n",
            "seed",
            "direction",
            "distance_cpp",
            "distance_cpp_bfs",
            "distance_java",
            "time_ms_cpp",
            "time_ms_cpp_flipdist",
            "time_ms_cpp_bfs",
            "time_ms_java",
            "expanded_cpp",
            "expanded_java",
            "enqueued_cpp",
            "enqueued_java",
            "visited_cpp",
            "visited_java",
            "max_queue_cpp",
            "max_queue_java",
            "duplicates_cpp",
            "duplicates_java",
            "status_cpp",
            "status_cpp_flipdist",
            "status_cpp_bfs",
            "status_java",
            "max_k_cpp",
            "solver_cpp",
            "solver_java",
            "tree_a",
            "tree_b",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for cpp_row, java_row in merged_rows:
            out = {
                "case_type": cpp_row.get("case_type"),
                "n": cpp_row.get("n"),
                "seed": cpp_row.get("seed"),
                "direction": cpp_row.get("direction"),
                "distance_cpp": cpp_row.get("distance"),
                "time_ms_cpp": cpp_row.get("time_ms"),
                "expanded_cpp": cpp_row.get("expanded"),
                "enqueued_cpp": cpp_row.get("enqueued"),
                "visited_cpp": cpp_row.get("visited"),
                "max_queue_cpp": cpp_row.get("max_queue"),
                "duplicates_cpp": cpp_row.get("duplicates"),
                "status_cpp": cpp_row.get("status"),
                "solver_cpp": cpp_row.get("solver"),
                "tree_a": cpp_row.get("tree_a"),
                "tree_b": cpp_row.get("tree_b"),
            }
            if "distance_bfs" in cpp_row:
                out["distance_cpp_bfs"] = cpp_row.get("distance_bfs")
            if "time_ms_flipdist" in cpp_row:
                out["time_ms_cpp_flipdist"] = cpp_row.get("time_ms_flipdist")
            if "time_ms_bfs" in cpp_row:
                out["time_ms_cpp_bfs"] = cpp_row.get("time_ms_bfs")
            if "status_flipdist" in cpp_row:
                out["status_cpp_flipdist"] = cpp_row.get("status_flipdist")
            if "status_bfs" in cpp_row:
                out["status_cpp_bfs"] = cpp_row.get("status_bfs")
            if "max_k" in cpp_row:
                out["max_k_cpp"] = cpp_row.get("max_k")
            if java_row:
                out.update({
                    "distance_java": java_row.get("distance"),
                    "time_ms_java": java_row.get("time_ms"),
                    "expanded_java": java_row.get("expanded"),
                    "enqueued_java": java_row.get("enqueued"),
                    "visited_java": java_row.get("visited"),
                    "max_queue_java": java_row.get("max_queue"),
                    "duplicates_java": java_row.get("duplicates"),
                    "status_java": java_row.get("status"),
                    "solver_java": java_row.get("solver"),
                })
            writer.writerow(out)


def auto_fill_fallback(cfg, cpp_rows):
    if getattr(cfg, "use_bidir", False):
        return cpp_rows
    if getattr(cfg, "disable_auto_bidir", False):
        return cpp_rows
    if getattr(cfg, "cpp_program", "") == "flipdist_asan":
        return cpp_rows
    if any(row.get("solver") == "flipdist" for row in cpp_rows):
        return cpp_rows

    pending = {
        (row.get("case_type"), row.get("n"), row.get("seed"))
        for row in cpp_rows
        if row.get("status") not in {"ok", "fallback_ok"}
    }
    if not pending:
        return cpp_rows

    for case_type, n, seed in sorted(pending):
        fb_cfg = SimpleNamespace(
            cpp_binary=cfg.cpp_binary,
            case_type=case_type,
            n=n,
            count=1,
            seed=seed,
            cpp_program=cfg.cpp_program,
            time_limit=None,
            visited_cap=cfg.visited_cap,
            queue_cap=cfg.queue_cap,
            use_bidir=True,
            bidir_cap=cfg.bidir_cap,
            disable_auto_bidir=True,
            prefer_bidir=True,
        )
        fb_rows = run_cpp_cli(fb_cfg)
        rows_by_dir = {
            (case_type, n, row.get("seed"), row.get("direction")): row
            for row in fb_rows
        }
        for idx, row in enumerate(cpp_rows):
            key = (
                row.get("case_type"),
                row.get("n"),
                row.get("seed"),
                row.get("direction"),
            )
            if key in rows_by_dir:
                cpp_rows[idx] = rows_by_dir[key]
    return cpp_rows


def parse_args():
    parser = argparse.ArgumentParser(description="Run C++ and Java rotation distance comparisons.")
    parser.add_argument("--case", dest="case_type", choices=["comb", "random"], default="comb")
    parser.add_argument("--n", type=int, required=True, help="Number of internal nodes")
    parser.add_argument("--count", type=int, default=1, help="Number of instances to generate")
    parser.add_argument("--seed", type=int, default=12345, help="Base seed for random cases")
    parser.add_argument("--output", required=True, help="Path to write CSV output")
    parser.add_argument("--cpp-binary", default="cmake-build/compare_cli", help="Path to compare_cli executable")
    parser.add_argument("--cpp-program", default="test_asan", help="Program label for C++ rows")
    parser.add_argument("--java-binary", default="java", help="Java binary")
    parser.add_argument("--java-out", default="triangulation/out", help="Java compiled classes directory")
    parser.add_argument("--java-lib", default="triangulation/lib/acm.jar", help="Java ACM library jar")
    parser.add_argument("--time-limit", type=float, default=None, help="Time limit (seconds)")
    parser.add_argument("--visited-cap", type=int, default=None, help="Visited cap")
    parser.add_argument("--queue-cap", type=int, default=None, help="Queue cap")
    parser.add_argument("--use-bidir", action="store_true", help="Enable fallback bidirectional BFS in C++ runs")
    parser.add_argument("--bidir-cap", type=int, default=None, help="State cap for bidirectional fallback")
    parser.add_argument("--no-auto-bidir", dest="disable_auto_bidir", action="store_true",
                        help="Do not automatically rerun timed-out instances with bidirectional search")
    return parser.parse_args()


def main():
    cfg = parse_args()
    cpp_rows = run_cpp_cli(cfg)
    forward_rows = [
        row for row in cpp_rows
        if row.get("direction") == "a->b"
    ]
    java_rows = run_java_cli(cfg, forward_rows)
    cpp_rows = auto_fill_fallback(cfg, cpp_rows)
    merged = merge_results(cpp_rows, java_rows)
    write_csv(cfg, merged)

    mismatches = [
        (cpp_row, java_row)
        for cpp_row, java_row in merged
        if java_row and cpp_row.get("distance") != java_row.get("distance")
    ]
    if mismatches:
        sys.stderr.write(f"WARNING: {len(mismatches)} distance mismatches detected\n")
    else:
        sys.stderr.write("All distances match between C++ and Java runs.\n")


if __name__ == "__main__":
    main()
