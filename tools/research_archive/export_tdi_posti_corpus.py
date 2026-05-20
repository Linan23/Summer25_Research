#!/usr/bin/env python3
"""Export clustered TreeDistI post-I states into a replay corpus."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export TreeDistI post-I summary rows into a replayable corpus",
    )
    parser.add_argument("--input", required=True, help="CSV from analyze_tdi_posti_profile.py")
    parser.add_argument("--output-dir", required=True, help="Directory for manifest and tree files")
    parser.add_argument("--manifest", default="manifest.csv")
    parser.add_argument("--top", type=int, default=20, help="Top rows to export")
    parser.add_argument("--min-recurrence", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    cfg = parse_args()
    input_path = Path(cfg.input)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with input_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    rows = [r for r in rows if int(r.get("recurrence_count", 0)) >= cfg.min_recurrence]
    rows = rows[: cfg.top] if cfg.top > 0 else rows

    manifest_path = output_dir / cfg.manifest
    fieldnames = [
        "case_id",
        "count",
        "motif",
        "reduced_edges",
        "start_branching_nodes",
        "target_branching_nodes",
        "legal_children",
        "plateau_buckets",
        "conflicts",
        "s_size",
        "min_remaining_k",
        "max_remaining_k",
        "min_lower_bound",
        "max_lower_bound",
        "start_degree_histogram",
        "target_degree_histogram",
        "tree_a_file",
        "tree_b_file",
        "tree_a",
        "tree_b",
    ]

    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for idx, row in enumerate(rows, start=1):
            case_id = f"tdi_posti_{idx:04d}"
            case_dir = output_dir / case_id
            case_dir.mkdir(parents=True, exist_ok=True)
            tree_a_path = case_dir / "tree_a.txt"
            tree_b_path = case_dir / "tree_b.txt"
            tree_a = row["start_tree"]
            tree_b = row["target_tree"]
            tree_a_path.write_text(tree_a + "\n", encoding="utf-8")
            tree_b_path.write_text(tree_b + "\n", encoding="utf-8")
            writer.writerow(
                {
                    "case_id": case_id,
                    "count": row.get("recurrence_count", ""),
                    "motif": "tdi_post_i",
                    "reduced_edges": row.get("conflicts", ""),
                    "start_branching_nodes": row.get("start_branching_nodes", ""),
                    "target_branching_nodes": row.get("target_branching_nodes", ""),
                    "legal_children": "",
                    "plateau_buckets": "",
                    "conflicts": row.get("conflicts", ""),
                    "s_size": row.get("s_size", ""),
                    "min_remaining_k": row.get("min_remaining_k", ""),
                    "max_remaining_k": row.get("max_remaining_k", ""),
                    "min_lower_bound": row.get("min_lower_bound", ""),
                    "max_lower_bound": row.get("max_lower_bound", ""),
                    "start_degree_histogram": row.get("start_degree_histogram", ""),
                    "target_degree_histogram": row.get("target_degree_histogram", ""),
                    "tree_a_file": str(tree_a_path),
                    "tree_b_file": str(tree_b_path),
                    "tree_a": tree_a,
                    "tree_b": tree_b,
                }
            )

    print(f"Exported {len(rows)} TDI post-I cases to {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
