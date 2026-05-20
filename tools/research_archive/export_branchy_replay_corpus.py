#!/usr/bin/env python3
"""Export a deduped replay corpus from plateau JSONL dumps."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Aggregate:
    count: int = 0
    min_k: int = 10**9
    max_k: int = -1
    min_conflicts: int = 10**9
    max_conflicts: int = -1
    source_files: set[str] | None = None
    sample: dict | None = None

    def __post_init__(self) -> None:
        if self.source_files is None:
            self.source_files = set()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export reduced branchy plateau cores into a replay corpus",
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="JSONL files or directories that contain plateau dump records",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory that will receive tree files and the manifest CSV",
    )
    parser.add_argument(
        "--manifest",
        default="manifest.csv",
        help="Manifest filename relative to --output-dir",
    )
    parser.add_argument(
        "--motifs",
        nargs="+",
        default=["branchy", "branchy_to_broom"],
        help="Motifs to keep",
    )
    parser.add_argument("--min-reduced-edges", type=int, default=11)
    parser.add_argument("--max-reduced-edges", type=int, default=14)
    parser.add_argument(
        "--branch-pairs",
        nargs="*",
        default=["3,3", "3,4", "4,4"],
        help="Allowed branching-node pairs as min,max",
    )
    parser.add_argument(
        "--dead-end-only",
        action="store_true",
        help="Keep only states with legal_children=0 and plateau_buckets=0",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=0,
        help="Keep only the top N aggregated cases by occurrence count",
    )
    return parser.parse_args()


def normalize_pair(text: str) -> tuple[int, int]:
    left, right = text.split(",", 1)
    a = int(left.strip())
    b = int(right.strip())
    return (min(a, b), max(a, b))


def iter_jsonl_paths(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_file():
            paths.append(path)
            continue
        if path.is_dir():
            paths.extend(sorted(p for p in path.rglob("*.jsonl") if p.is_file()))
    return paths


def keep_record(row: dict, motifs: set[str], branch_pairs: set[tuple[int, int]], cfg: argparse.Namespace) -> bool:
    if "start_tree" not in row or "target_tree" not in row:
        return False
    if row.get("motif") not in motifs:
        return False
    reduced_edges = int(row.get("reduced_edges", -1))
    if reduced_edges < cfg.min_reduced_edges or reduced_edges > cfg.max_reduced_edges:
        return False
    pair = (
        min(int(row.get("start_branching_nodes", 0)), int(row.get("target_branching_nodes", 0))),
        max(int(row.get("start_branching_nodes", 0)), int(row.get("target_branching_nodes", 0))),
    )
    if pair not in branch_pairs:
        return False
    if cfg.dead_end_only:
        if int(row.get("legal_children", -1)) != 0:
            return False
        if int(row.get("plateau_buckets", -1)) != 0:
            return False
    return True


def aggregate_records(paths: list[Path], cfg: argparse.Namespace) -> list[tuple[tuple[str, str], Aggregate]]:
    motifs = set(cfg.motifs)
    branch_pairs = {normalize_pair(text) for text in cfg.branch_pairs}
    groups: dict[tuple[str, str], Aggregate] = {}

    for path in paths:
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not keep_record(row, motifs, branch_pairs, cfg):
                    continue
                key = (row["start_tree"], row["target_tree"])
                agg = groups.setdefault(key, Aggregate())
                agg.count += 1
                agg.min_k = min(agg.min_k, int(row.get("k", -1)))
                agg.max_k = max(agg.max_k, int(row.get("k", -1)))
                agg.min_conflicts = min(agg.min_conflicts, int(row.get("conflicts", -1)))
                agg.max_conflicts = max(agg.max_conflicts, int(row.get("conflicts", -1)))
                agg.source_files.add(str(path))
                if agg.sample is None:
                    agg.sample = row

    items = sorted(
        groups.items(),
        key=lambda kv: (
            -kv[1].count,
            -int(kv[1].sample.get("reduced_edges", -1)),
            kv[1].sample.get("motif", ""),
            kv[0][0],
            kv[0][1],
        ),
    )
    if cfg.top and cfg.top > 0:
        items = items[: cfg.top]
    return items


def export_corpus(items: list[tuple[tuple[str, str], Aggregate]], output_dir: Path, manifest_name: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / manifest_name
    fieldnames = [
        "case_id",
        "count",
        "motif",
        "reduced_edges",
        "min_conflicts",
        "max_conflicts",
        "min_k",
        "max_k",
        "start_branching_nodes",
        "target_branching_nodes",
        "start_degree_histogram",
        "target_degree_histogram",
        "legal_children",
        "plateau_buckets",
        "source_file_count",
        "tree_a_file",
        "tree_b_file",
        "tree_a",
        "tree_b",
    ]

    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for idx, ((tree_a, tree_b), agg) in enumerate(items, start=1):
            case_id = f"core_{idx:04d}"
            case_dir = output_dir / case_id
            case_dir.mkdir(parents=True, exist_ok=True)
            tree_a_path = case_dir / "tree_a.txt"
            tree_b_path = case_dir / "tree_b.txt"
            tree_a_path.write_text(tree_a + "\n", encoding="utf-8")
            tree_b_path.write_text(tree_b + "\n", encoding="utf-8")

            row = {
                "case_id": case_id,
                "count": agg.count,
                "motif": agg.sample.get("motif", ""),
                "reduced_edges": agg.sample.get("reduced_edges", ""),
                "min_conflicts": agg.min_conflicts,
                "max_conflicts": agg.max_conflicts,
                "min_k": agg.min_k if agg.min_k < 10**9 else "",
                "max_k": agg.max_k if agg.max_k >= 0 else "",
                "start_branching_nodes": agg.sample.get("start_branching_nodes", ""),
                "target_branching_nodes": agg.sample.get("target_branching_nodes", ""),
                "start_degree_histogram": agg.sample.get("start_degree_histogram", ""),
                "target_degree_histogram": agg.sample.get("target_degree_histogram", ""),
                "legal_children": agg.sample.get("legal_children", ""),
                "plateau_buckets": agg.sample.get("plateau_buckets", ""),
                "source_file_count": len(agg.source_files),
                "tree_a_file": str(tree_a_path),
                "tree_b_file": str(tree_b_path),
                "tree_a": tree_a,
                "tree_b": tree_b,
            }
            writer.writerow(row)

    return manifest_path


def main() -> int:
    cfg = parse_args()
    paths = iter_jsonl_paths(cfg.inputs)
    if not paths:
        raise SystemExit("No JSONL inputs found")
    items = aggregate_records(paths, cfg)
    manifest_path = export_corpus(items, Path(cfg.output_dir), cfg.manifest)
    print(f"Loaded {len(paths)} JSONL files")
    print(f"Exported {len(items)} replay cases to {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
