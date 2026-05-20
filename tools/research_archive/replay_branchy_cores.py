#!/usr/bin/env python3
"""Replay exported branchy-core cases through FlipDist and optional Java exact."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay exported branchy-core cases through FlipDist and Java exact",
    )
    parser.add_argument("--manifest", required=True, help="Manifest CSV from export_branchy_replay_corpus.py")
    parser.add_argument("--cpp-binary", default="./build/flipdist")
    parser.add_argument("--bfs-cap", type=int, default=1)
    parser.add_argument("--max-k", type=int, default=0, help="Override max-k for FlipDist")
    parser.add_argument(
        "--flipdist-env",
        action="append",
        default=[],
        help="Extra environment variable in KEY=VALUE form for FlipDist",
    )
    parser.add_argument("--java", action="store_true", help="Also run Java exact on the same replay cases")
    parser.add_argument("--java-binary", default="java")
    parser.add_argument("--java-out", default="oracle/java/out")
    parser.add_argument("--java-lib", default="oracle/java/lib/acm.jar")
    parser.add_argument("--java-time-limit", type=float, default=120.0)
    parser.add_argument("--java-visited-cap", type=int, default=15_000_000)
    parser.add_argument("--java-queue-cap", type=int, default=15_000_000)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--print", action="store_true", dest="do_print")
    return parser.parse_args()


def run_cmd(cmd: list[str], *, env: dict[str, str] | None = None, input_text: str | None = None) -> str:
    proc = subprocess.run(
        cmd,
        input=input_text,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc.stdout


def parse_n_from_canonical(canonical: str) -> int:
    if not canonical.startswith("P:"):
        raise ValueError("Canonical traversal must start with P:")
    split = canonical.find(";I:")
    if split < 0:
        raise ValueError("Canonical traversal missing ;I:")
    inorder = canonical[split + 3 :].strip()
    if not inorder:
        return 0
    return len([part for part in inorder.split(",") if part.strip()])


def compute_max_k(row: dict, override: int) -> int:
    if override > 0:
        return override
    n = parse_n_from_canonical(row["tree_a"])
    return max(1, 3 * n + 10)


def flipdist_env(extra_vars: list[str]) -> dict[str, str]:
    env = os.environ.copy()
    for item in extra_vars:
        if "=" not in item:
            raise ValueError(f"Invalid --flipdist-env value: {item}")
        key, value = item.split("=", 1)
        env[key] = value
    return env


def run_flipdist(cfg: argparse.Namespace, row: dict, env: dict[str, str]) -> list[dict]:
    cmd = [
        cfg.cpp_binary,
        "--tree-a-file",
        row["tree_a_file"],
        "--tree-b-file",
        row["tree_b_file"],
        "--max-k",
        str(compute_max_k(row, cfg.max_k)),
        "--bfs-cap",
        str(cfg.bfs_cap),
    ]
    out = run_cmd(cmd, env=env)
    return [json.loads(line) for line in out.splitlines() if line.strip().startswith("{")]


def run_java(cfg: argparse.Namespace, row: dict) -> list[dict]:
    payload = {
        "case_type": row["case_id"],
        "seed": 0,
        "tree_a": row["tree_a"],
        "tree_b": row["tree_b"],
        "time_limit": cfg.java_time_limit,
        "visited_cap": cfg.java_visited_cap,
        "queue_cap": cfg.java_queue_cap,
    }
    java_classpath = f"{cfg.java_out}:{cfg.java_lib}"
    cmd = [
        cfg.java_binary,
        "-cp",
        java_classpath,
        "TriangulationMetricsCli",
    ]
    out = run_cmd(cmd, input_text=json.dumps(payload) + "\n")
    return [json.loads(line) for line in out.splitlines() if line.strip().startswith("{")]


def merge_case(manifest_row: dict, flipdist_rows: list[dict], java_rows: list[dict]) -> list[dict]:
    java_map = {
        (row.get("tree_a"), row.get("tree_b"), row.get("direction")): row
        for row in java_rows
    }
    merged: list[dict] = []
    for row in flipdist_rows:
        key = (row.get("tree_a"), row.get("tree_b"), row.get("direction"))
        jrow = java_map.get(key, {})
        merged.append(
            {
                "case_id": manifest_row["case_id"],
                "count": manifest_row.get("count", ""),
                "motif": manifest_row.get("motif", ""),
                "reduced_edges": manifest_row.get("reduced_edges", ""),
                "start_branching_nodes": manifest_row.get("start_branching_nodes", ""),
                "target_branching_nodes": manifest_row.get("target_branching_nodes", ""),
                "legal_children": manifest_row.get("legal_children", ""),
                "plateau_buckets": manifest_row.get("plateau_buckets", ""),
                "direction": row.get("direction", ""),
                "distance_flipdist": row.get("distance"),
                "status_flipdist": row.get("status"),
                "time_ms_flipdist": row.get("time_ms"),
                "distance_java": jrow.get("distance", ""),
                "status_java": jrow.get("status", ""),
                "time_ms_java": jrow.get("time_ms", ""),
                "tree_a": row.get("tree_a"),
                "tree_b": row.get("tree_b"),
            }
        )
    return merged


def load_manifest(path: Path, limit: int) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if limit > 0:
        rows = rows[:limit]
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "case_id",
        "count",
        "motif",
        "reduced_edges",
        "start_branching_nodes",
        "target_branching_nodes",
        "legal_children",
        "plateau_buckets",
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
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    cfg = parse_args()
    manifest_rows = load_manifest(Path(cfg.manifest), cfg.limit)
    env = flipdist_env(cfg.flipdist_env)

    merged_rows: list[dict] = []
    mismatches = 0
    failures = 0

    for idx, row in enumerate(manifest_rows, start=1):
        if cfg.do_print:
            print(
                f"Replay case {idx}/{len(manifest_rows)}: {row['case_id']} "
                f"edges={row.get('reduced_edges','')} motif={row.get('motif','')}"
            )
        try:
            flipdist_rows = run_flipdist(cfg, row, env)
            java_rows = run_java(cfg, row) if cfg.java else []
            case_rows = merge_case(row, flipdist_rows, java_rows)
            for merged in case_rows:
                if cfg.java and merged["distance_java"] != "" and merged["distance_flipdist"] != merged["distance_java"]:
                    mismatches += 1
            merged_rows.extend(case_rows)
        except Exception as exc:
            failures += 1
            merged_rows.append(
                {
                    "case_id": row["case_id"],
                    "count": row.get("count", ""),
                    "motif": row.get("motif", ""),
                    "reduced_edges": row.get("reduced_edges", ""),
                    "start_branching_nodes": row.get("start_branching_nodes", ""),
                    "target_branching_nodes": row.get("target_branching_nodes", ""),
                    "legal_children": row.get("legal_children", ""),
                    "plateau_buckets": row.get("plateau_buckets", ""),
                    "direction": "",
                    "distance_flipdist": "",
                    "distance_java": "",
                    "status_flipdist": f"error:{exc}",
                    "status_java": "",
                    "time_ms_flipdist": "",
                    "time_ms_java": "",
                    "tree_a": row.get("tree_a", ""),
                    "tree_b": row.get("tree_b", ""),
                }
            )

    output_path = Path(cfg.output)
    write_csv(output_path, merged_rows)

    if cfg.do_print:
        print(f"Wrote {len(merged_rows)} rows to {output_path}")
        if cfg.java:
            if mismatches:
                print(f"WARNING: {mismatches} distance mismatches detected")
            else:
                print("All replayed distances match Java exact")
        if failures:
            print(f"Completed with {failures} failed replay cases")

    return 0 if mismatches == 0 and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
