#!/usr/bin/env python3
"""Export captured plateau decision states into a replay corpus."""

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
    source_files: set[str] | None = None
    sample: dict | None = None

    def __post_init__(self) -> None:
        if self.source_files is None:
            self.source_files = set()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export plateau decision states into a replayable corpus",
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="JSONL files or directories containing plateau dump records",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", default="manifest.csv")
    parser.add_argument("--motifs", nargs="+", default=["branchy", "branchy_to_broom"])
    parser.add_argument("--min-reduced-edges", type=int, default=0)
    parser.add_argument("--max-reduced-edges", type=int, default=10**9)
    parser.add_argument(
        "--branch-pairs",
        nargs="*",
        default=[],
        help="Allowed branching-node pairs as min,max",
    )
    parser.add_argument("--dead-end-only", action="store_true")
    parser.add_argument("--require-k-equals-conflicts", action="store_true")
    parser.add_argument("--target-tree", default="")
    parser.add_argument("--start-tree", default="")
    parser.add_argument("--start-degree-histogram", default="")
    parser.add_argument("--target-degree-histogram", default="")
    parser.add_argument("--top", type=int, default=0)
    return parser.parse_args()


def normalize_pair(text: str) -> tuple[int, int]:
    left, right = text.split(",", 1)
    a = int(left.strip())
    b = int(right.strip())
    return (min(a, b), max(a, b))


def iter_jsonl_paths(inputs: list[str]) -> list[Path]:
    out: list[Path] = []
    for item in inputs:
        path = Path(item)
        if path.is_file():
            out.append(path)
            continue
        if path.is_dir():
            out.extend(sorted(p for p in path.rglob("*.jsonl") if p.is_file()))
    return out


def keep_record(
    row: dict,
    *,
    motifs: set[str],
    branch_pairs: set[tuple[int, int]],
    cfg: argparse.Namespace,
) -> bool:
    if "start_tree" not in row or "target_tree" not in row:
        return False
    if row.get("motif") not in motifs:
        return False
    reduced_edges = int(row.get("reduced_edges", -1))
    if reduced_edges < cfg.min_reduced_edges or reduced_edges > cfg.max_reduced_edges:
        return False
    if cfg.branch_pairs:
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
    if cfg.require_k_equals_conflicts and int(row.get("k", -1)) != int(row.get("conflicts", -2)):
        return False
    if cfg.target_tree and row.get("target_tree") != cfg.target_tree:
        return False
    if cfg.start_tree and row.get("start_tree") != cfg.start_tree:
        return False
    if cfg.start_degree_histogram and row.get("start_degree_histogram") != cfg.start_degree_histogram:
        return False
    if cfg.target_degree_histogram and row.get("target_degree_histogram") != cfg.target_degree_histogram:
        return False
    return True


def aggregate_records(paths: list[Path], cfg: argparse.Namespace) -> list[tuple[tuple[str, str, int], Aggregate]]:
    motifs = set(cfg.motifs)
    branch_pairs = {normalize_pair(text) for text in cfg.branch_pairs}
    groups: dict[tuple[str, str, int], Aggregate] = {}
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
                if not keep_record(row, motifs=motifs, branch_pairs=branch_pairs, cfg=cfg):
                    continue
                key = (row["start_tree"], row["target_tree"], int(row.get("k", -1)))
                agg = groups.setdefault(key, Aggregate())
                agg.count += 1
                agg.source_files.add(str(path))
                if agg.sample is None:
                    agg.sample = row
    items = sorted(
        groups.items(),
        key=lambda kv: (
            -kv[1].count,
            -int(kv[1].sample.get("reduced_edges", -1)),
            kv[1].sample.get("motif", ""),
            kv[0][2],
            kv[0][0],
            kv[0][1],
        ),
    )
    if cfg.top and cfg.top > 0:
        items = items[: cfg.top]
    return items


def export_corpus(items: list[tuple[tuple[str, str, int], Aggregate]], output_dir: Path, manifest_name: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / manifest_name
    fieldnames = [
        "case_id",
        "count",
        "motif",
        "reduced_edges",
        "conflicts",
        "k",
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
        for idx, ((tree_a, tree_b, k), agg) in enumerate(items, start=1):
            case_id = f"decision_{idx:04d}"
            case_dir = output_dir / case_id
            case_dir.mkdir(parents=True, exist_ok=True)
            tree_a_path = case_dir / "tree_a.txt"
            tree_b_path = case_dir / "tree_b.txt"
            tree_a_path.write_text(tree_a + "\n", encoding="utf-8")
            tree_b_path.write_text(tree_b + "\n", encoding="utf-8")
            sample = agg.sample or {}
            writer.writerow(
                {
                    "case_id": case_id,
                    "count": agg.count,
                    "motif": sample.get("motif", ""),
                    "reduced_edges": sample.get("reduced_edges", ""),
                    "conflicts": sample.get("conflicts", ""),
                    "k": k,
                    "start_branching_nodes": sample.get("start_branching_nodes", ""),
                    "target_branching_nodes": sample.get("target_branching_nodes", ""),
                    "start_degree_histogram": sample.get("start_degree_histogram", ""),
                    "target_degree_histogram": sample.get("target_degree_histogram", ""),
                    "legal_children": sample.get("legal_children", ""),
                    "plateau_buckets": sample.get("plateau_buckets", ""),
                    "source_file_count": len(agg.source_files),
                    "tree_a_file": str(tree_a_path),
                    "tree_b_file": str(tree_b_path),
                    "tree_a": tree_a,
                    "tree_b": tree_b,
                }
            )
    return manifest_path


def main() -> int:
    cfg = parse_args()
    paths = iter_jsonl_paths(cfg.inputs)
    items = aggregate_records(paths, cfg)
    manifest = export_corpus(items, Path(cfg.output_dir), cfg.manifest)
    print(f"Exported {len(items)} plateau decision states to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
