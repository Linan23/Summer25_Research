#!/usr/bin/env python3
"""Summarize TreeDistI post-I prune JSONL records."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze TreeDistI post-I prune profile JSONL")
    parser.add_argument("--inputs", nargs="+", required=True, help="JSONL files or directories")
    parser.add_argument("--output", required=True, help="Output CSV")
    parser.add_argument("--top", type=int, default=50)
    return parser.parse_args()


def iter_paths(inputs: list[str]) -> list[Path]:
    out: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_file():
            out.append(path)
        elif path.is_dir():
            out.extend(sorted(p for p in path.rglob("*.jsonl") if p.is_file()))
    return out


def main() -> int:
    cfg = parse_args()
    groups: dict[tuple[str, str, int], dict] = {}
    for path in iter_paths(cfg.inputs):
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if row.get("record_type") != "tdi_post_i_state":
                    continue
                key = (row["start_tree"], row["target_tree"], int(row.get("s_size", 0)))
                agg = groups.setdefault(
                    key,
                    {
                        "count": 0,
                        "recurrence_count": 0,
                        "conflicts": int(row.get("conflicts", 0)),
                        "s_size": int(row.get("s_size", 0)),
                        "min_remaining_k": int(row.get("min_remaining_k", -1)),
                        "max_remaining_k": int(row.get("max_remaining_k", -1)),
                        "min_lower_bound": int(row.get("min_lower_bound", -1)),
                        "max_lower_bound": int(row.get("max_lower_bound", -1)),
                        "start_branching_nodes": int(row.get("start_branching_nodes", 0)),
                        "target_branching_nodes": int(row.get("target_branching_nodes", 0)),
                        "start_degree_histogram": row.get("start_degree_histogram", ""),
                        "target_degree_histogram": row.get("target_degree_histogram", ""),
                        "start_tree": row.get("start_tree", ""),
                        "target_tree": row.get("target_tree", ""),
                    },
                )
                agg["count"] += 1
                agg["recurrence_count"] += int(row.get("recurrence_count", 0))
                agg["min_remaining_k"] = min(agg["min_remaining_k"], int(row.get("min_remaining_k", -1)))
                agg["max_remaining_k"] = max(agg["max_remaining_k"], int(row.get("max_remaining_k", -1)))
                agg["min_lower_bound"] = min(agg["min_lower_bound"], int(row.get("min_lower_bound", -1)))
                agg["max_lower_bound"] = max(agg["max_lower_bound"], int(row.get("max_lower_bound", -1)))

    rows = sorted(
        groups.values(),
        key=lambda r: (
            -r["recurrence_count"],
            -r["count"],
            -r["max_lower_bound"],
            -r["conflicts"],
            r["start_tree"],
            r["target_tree"],
        ),
    )
    if cfg.top > 0:
        rows = rows[: cfg.top]

    for row in rows:
        row["min_gap"] = row["min_lower_bound"] - row["max_remaining_k"]
        row["max_gap"] = row["max_lower_bound"] - row["min_remaining_k"]

    fieldnames = [
        "count",
        "recurrence_count",
        "conflicts",
        "s_size",
        "min_remaining_k",
        "max_remaining_k",
        "min_lower_bound",
        "max_lower_bound",
        "min_gap",
        "max_gap",
        "start_branching_nodes",
        "target_branching_nodes",
        "start_degree_histogram",
        "target_degree_histogram",
        "start_tree",
        "target_tree",
    ]
    out_path = Path(cfg.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
