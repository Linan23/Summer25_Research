#!/usr/bin/env python3
"""Rank recurring family patterns inside a branchy replay manifest."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze branchy replay manifest families and optionally emit a filtered manifest"
    )
    parser.add_argument("--manifest", required=True, help="Replay manifest CSV")
    parser.add_argument("--output-summary", required=True, help="Family summary CSV")
    parser.add_argument(
        "--output-top-manifest",
        help="Optional filtered manifest containing only the top-ranked family",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of families to write in the summary",
    )
    return parser.parse_args()


def family_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row.get("motif", ""),
        row.get("reduced_edges", ""),
        row.get("start_branching_nodes", ""),
        row.get("target_branching_nodes", ""),
        row.get("start_degree_histogram", ""),
        row.get("target_degree_histogram", ""),
        row.get("legal_children", ""),
        row.get("plateau_buckets", ""),
    )


def main() -> int:
    cfg = parse_args()
    manifest_path = Path(cfg.manifest)
    rows = list(csv.DictReader(manifest_path.open(newline="")))

    families: dict[tuple[str, ...], dict[str, object]] = defaultdict(
        lambda: {
            "count": 0,
            "source_file_count_sum": 0,
            "case_ids": [],
            "edge_counter": Counter(),
            "branch_pair": "",
        }
    )

    for row in rows:
        key = family_key(row)
        bucket = families[key]
        bucket["count"] = int(bucket["count"]) + 1
        bucket["source_file_count_sum"] = int(bucket["source_file_count_sum"]) + int(
            row.get("source_file_count", "0") or "0"
        )
        bucket["case_ids"].append(row["case_id"])
        bucket["edge_counter"][row.get("reduced_edges", "")] += 1
        bucket["branch_pair"] = f"{row.get('start_branching_nodes','')},{row.get('target_branching_nodes','')}"

    ranked = sorted(
        families.items(),
        key=lambda kv: (
            -int(kv[1]["count"]),
            -int(kv[1]["source_file_count_sum"]),
            kv[0],
        ),
    )

    summary_path = Path(cfg.output_summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "rank",
                "family_key",
                "count",
                "source_file_count_sum",
                "motif",
                "reduced_edges",
                "branch_pair",
                "start_degree_histogram",
                "target_degree_histogram",
                "legal_children",
                "plateau_buckets",
                "sample_case_ids",
            ],
        )
        writer.writeheader()
        for idx, (key, data) in enumerate(ranked[: cfg.top], start=1):
            motif, reduced_edges, start_bn, target_bn, start_deg, target_deg, legal_children, plateau_buckets = key
            writer.writerow(
                {
                    "rank": idx,
                    "family_key": "|".join(key),
                    "count": data["count"],
                    "source_file_count_sum": data["source_file_count_sum"],
                    "motif": motif,
                    "reduced_edges": reduced_edges,
                    "branch_pair": f"{start_bn},{target_bn}",
                    "start_degree_histogram": start_deg,
                    "target_degree_histogram": target_deg,
                    "legal_children": legal_children,
                    "plateau_buckets": plateau_buckets,
                    "sample_case_ids": ";".join(data["case_ids"][:10]),
                }
            )

    if cfg.output_top_manifest and ranked:
        top_key = ranked[0][0]
        filtered_rows = [row for row in rows if family_key(row) == top_key]
        out_path = Path(cfg.output_top_manifest)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(filtered_rows)

    print(f"Loaded {len(rows)} replay cases")
    print(f"Wrote {summary_path}")
    if cfg.output_top_manifest and ranked:
        print(f"Wrote top-family manifest with {len(filtered_rows)} rows to {cfg.output_top_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
