#!/usr/bin/env python3
"""Replay captured plateau decision states with their original budget k."""

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
        description="Replay plateau decision states through FlipDist with captured k",
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--cpp-binary", default="./build/flipdist")
    parser.add_argument("--bfs-cap", type=int, default=1)
    parser.add_argument("--k-offset", type=int, default=0, help="Additive offset applied to captured k")
    parser.add_argument("--flipdist-env", action="append", default=[])
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--print", action="store_true", dest="do_print")
    return parser.parse_args()


def flipdist_env(extra_vars: list[str]) -> dict[str, str]:
    env = os.environ.copy()
    for item in extra_vars:
        if "=" not in item:
            raise ValueError(f"Invalid --flipdist-env value: {item}")
        key, value = item.split("=", 1)
        env[key] = value
    return env


def load_manifest(path: Path, limit: int) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if limit > 0:
        rows = rows[:limit]
    return rows


def run_cmd(cmd: list[str], env: dict[str, str]) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc.stdout


def run_flipdist(cfg: argparse.Namespace, row: dict, env: dict[str, str]) -> list[dict]:
    k = int(row["k"]) + cfg.k_offset
    cmd = [
        cfg.cpp_binary,
        "--tree-a-file",
        row["tree_a_file"],
        "--tree-b-file",
        row["tree_b_file"],
        "--max-k",
        str(max(1, k)),
        "--bfs-cap",
        str(cfg.bfs_cap),
    ]
    out = run_cmd(cmd, env)
    return [json.loads(line) for line in out.splitlines() if line.strip().startswith("{")]


def merge_rows(manifest_row: dict, flipdist_rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for row in flipdist_rows:
        out.append(
            {
                "case_id": manifest_row["case_id"],
                "count": manifest_row.get("count", ""),
                "motif": manifest_row.get("motif", ""),
                "reduced_edges": manifest_row.get("reduced_edges", ""),
                "conflicts": manifest_row.get("conflicts", ""),
                "k": manifest_row.get("k", ""),
                "start_branching_nodes": manifest_row.get("start_branching_nodes", ""),
                "target_branching_nodes": manifest_row.get("target_branching_nodes", ""),
                "legal_children": manifest_row.get("legal_children", ""),
                "plateau_buckets": manifest_row.get("plateau_buckets", ""),
                "direction": row.get("direction", ""),
                "distance_flipdist": row.get("distance"),
                "status_flipdist": row.get("status"),
                "time_ms_flipdist": row.get("time_ms"),
                "tree_a": row.get("tree_a"),
                "tree_b": row.get("tree_b"),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "case_id",
        "count",
        "motif",
        "reduced_edges",
        "conflicts",
        "k",
        "start_branching_nodes",
        "target_branching_nodes",
        "legal_children",
        "plateau_buckets",
        "direction",
        "distance_flipdist",
        "status_flipdist",
        "time_ms_flipdist",
        "tree_a",
        "tree_b",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    cfg = parse_args()
    rows = load_manifest(Path(cfg.manifest), cfg.limit)
    env = flipdist_env(cfg.flipdist_env)
    merged: list[dict] = []
    failures = 0
    for idx, row in enumerate(rows, start=1):
        if cfg.do_print:
            print(
                f"Replay decision {idx}/{len(rows)}: {row['case_id']} "
                f"edges={row.get('reduced_edges','')} k={row.get('k','')} motif={row.get('motif','')}"
            )
        try:
            merged.extend(merge_rows(row, run_flipdist(cfg, row, env)))
        except Exception as exc:
            failures += 1
            merged.append(
                {
                    "case_id": row["case_id"],
                    "count": row.get("count", ""),
                    "motif": row.get("motif", ""),
                    "reduced_edges": row.get("reduced_edges", ""),
                    "conflicts": row.get("conflicts", ""),
                    "k": row.get("k", ""),
                    "start_branching_nodes": row.get("start_branching_nodes", ""),
                    "target_branching_nodes": row.get("target_branching_nodes", ""),
                    "legal_children": row.get("legal_children", ""),
                    "plateau_buckets": row.get("plateau_buckets", ""),
                    "direction": "",
                    "distance_flipdist": "",
                    "status_flipdist": f"error:{exc}",
                    "time_ms_flipdist": "",
                    "tree_a": row.get("tree_a", ""),
                    "tree_b": row.get("tree_b", ""),
                }
            )
    output_path = Path(cfg.output)
    write_csv(output_path, merged)
    if cfg.do_print:
        print(f"Wrote {len(merged)} rows to {output_path}")
        if failures:
            print(f"Completed with {failures} failed replay cases")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
