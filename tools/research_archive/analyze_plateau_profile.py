#!/usr/bin/env python3
"""Summarize plateau-state JSONL records emitted by FlipDist profile mode."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate plateau_state JSONL records from FlipDist profile output."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Profile output files to read")
    parser.add_argument("--output", default=None, help="Optional CSV summary path")
    parser.add_argument("--top", type=int, default=25, help="Rows to print to stdout")
    return parser.parse_args()


def load_records(paths: list[Path]) -> list[dict]:
    best_records: dict[tuple[str, str, str, str], tuple[tuple[int, int, int, int], dict]] = {}
    direction_has_final: set[tuple[str, str]] = set()
    for path in paths:
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or not line.startswith("{"):
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if row.get("record_type") not in {"plateau_state", "plateau_state_snapshot"}:
                    continue
                row["_source_file"] = str(path)
                if row.get("record_type") == "plateau_state":
                    direction_has_final.add((str(path), str(row.get("direction", ""))))
                key = (
                    str(path),
                    str(row.get("direction", "")),
                    str(row.get("state_hi", "")),
                    str(row.get("state_lo", "")),
                )
                rank = (
                    1 if row.get("record_type") == "plateau_state" else 0,
                    int(row.get("snapshot_seq", 0) or 0),
                    int(row.get("recurrence_count", 0) or 0),
                    int(row.get("success_count", 0) or 0)
                    + int(row.get("fail_count", 0) or 0)
                    + int(row.get("timeout_count", 0) or 0),
                )
                prev = best_records.get(key)
                if prev is None or rank > prev[0]:
                    best_records[key] = (rank, row)
    records: list[dict] = []
    for (path_str, direction, _, _), (_, row) in best_records.items():
        if row.get("record_type") == "plateau_state_snapshot" and (path_str, direction) not in direction_has_final:
            row["timeout_count"] = max(1, int(row.get("timeout_count", 0) or 0))
        records.append(row)
    return records


def aggregate(records: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], dict] = {}
    seen_runs: dict[tuple[str, str], set[tuple[str, int, int, str]]] = defaultdict(set)

    for row in records:
        key = (
            str(row.get("core_signature") or ""),
            f"{row.get('state_hi', '')}:{row.get('state_lo', '')}",
        )
        rec = grouped.setdefault(
            key,
            {
                "aggregate_key": key[0] or key[1],
                "motif": row.get("motif", "unknown"),
                "core_signature": row.get("core_signature", ""),
                "internal_edges_max": 0,
                "reduced_core_internal_edges_min": None,
                "conflicts_max": 0,
                "min_k": None,
                "max_k": None,
                "legal_children_max": 0,
                "plateau_buckets_max": 0,
                "start_branching_nodes_max": 0,
                "target_branching_nodes_max": 0,
                "start_max_branch_subtree_edges_max": 0,
                "target_max_branch_subtree_edges_max": 0,
                "chain_arm_count_max": 0,
                "broom_arm_count_max": 0,
                "branchy_reduction_hits": 0,
                "branchy_forced_cost_max": 0,
                "branchy_residual_edges_min": None,
                "branchy_exact_bound_hits": 0,
                "branchy_exact_decision_hits": 0,
                "branchy_residual_motif_examples": set(),
                "recurrence_count": 0,
                "success_count": 0,
                "fail_count": 0,
                "timeout_count": 0,
                "free_edge_hits": 0,
                "free_edge_misses": 0,
                "start_degree_histogram_examples": set(),
                "target_degree_histogram_examples": set(),
                "run_count": 0,
                "sources": set(),
            },
        )

        rec["internal_edges_max"] = max(rec["internal_edges_max"], int(row.get("internal_edges", 0) or 0))
        reduced_core_internal_edges = int(row.get("reduced_core_internal_edges", 0) or 0)
        if rec["reduced_core_internal_edges_min"] is None:
            rec["reduced_core_internal_edges_min"] = reduced_core_internal_edges
        else:
            rec["reduced_core_internal_edges_min"] = min(rec["reduced_core_internal_edges_min"],
                                                         reduced_core_internal_edges)
        rec["conflicts_max"] = max(rec["conflicts_max"], int(row.get("conflicts", 0) or 0))
        rec["legal_children_max"] = max(rec["legal_children_max"], int(row.get("legal_children", 0) or 0))
        rec["plateau_buckets_max"] = max(rec["plateau_buckets_max"], int(row.get("plateau_buckets", 0) or 0))
        rec["start_branching_nodes_max"] = max(rec["start_branching_nodes_max"],
                                                int(row.get("start_branching_nodes", 0) or 0))
        rec["target_branching_nodes_max"] = max(rec["target_branching_nodes_max"],
                                                int(row.get("target_branching_nodes", 0) or 0))
        rec["start_max_branch_subtree_edges_max"] = max(
            rec["start_max_branch_subtree_edges_max"],
            int(row.get("start_max_branch_subtree_edges", 0) or 0),
        )
        rec["target_max_branch_subtree_edges_max"] = max(
            rec["target_max_branch_subtree_edges_max"],
            int(row.get("target_max_branch_subtree_edges", 0) or 0),
        )
        rec["chain_arm_count_max"] = max(rec["chain_arm_count_max"], int(row.get("chain_arm_count", 0) or 0))
        rec["broom_arm_count_max"] = max(rec["broom_arm_count_max"], int(row.get("broom_arm_count", 0) or 0))
        rec["branchy_reduction_hits"] += int(row.get("branchy_reduction_hits", 0) or 0)
        rec["branchy_forced_cost_max"] = max(
            rec["branchy_forced_cost_max"],
            int(row.get("branchy_forced_cost_max", 0) or 0),
        )
        branchy_residual_edges = int(row.get("branchy_residual_edges_min", -1) or -1)
        if branchy_residual_edges >= 0:
            if rec["branchy_residual_edges_min"] is None:
                rec["branchy_residual_edges_min"] = branchy_residual_edges
            else:
                rec["branchy_residual_edges_min"] = min(
                    rec["branchy_residual_edges_min"], branchy_residual_edges
                )
        rec["branchy_exact_bound_hits"] += int(row.get("branchy_exact_bound_hits", 0) or 0)
        rec["branchy_exact_decision_hits"] += int(row.get("branchy_exact_decision_hits", 0) or 0)
        rec["recurrence_count"] += int(row.get("recurrence_count", 0) or 0)
        rec["success_count"] += int(row.get("success_count", 0) or 0)
        rec["fail_count"] += int(row.get("fail_count", 0) or 0)
        rec["timeout_count"] += int(row.get("timeout_count", 0) or 0)
        rec["free_edge_hits"] += int(row.get("free_edge_hits", 0) or 0)
        rec["free_edge_misses"] += int(row.get("free_edge_misses", 0) or 0)
        rec["sources"].add(row["_source_file"])
        if row.get("start_degree_histogram"):
            rec["start_degree_histogram_examples"].add(str(row["start_degree_histogram"]))
        if row.get("target_degree_histogram"):
            rec["target_degree_histogram_examples"].add(str(row["target_degree_histogram"]))
        if row.get("branchy_residual_motif"):
            rec["branchy_residual_motif_examples"].add(str(row["branchy_residual_motif"]))

        min_k = row.get("min_k")
        max_k = row.get("max_k")
        if min_k is not None:
            min_k = int(min_k)
            rec["min_k"] = min_k if rec["min_k"] is None else min(rec["min_k"], min_k)
        if max_k is not None:
            max_k = int(max_k)
            rec["max_k"] = max_k if rec["max_k"] is None else max(rec["max_k"], max_k)

        run_id = (
            str(row.get("case_type", "")),
            int(row.get("n", -1) or -1),
            int(row.get("seed", -1) or -1),
            str(row.get("direction", "")),
        )
        seen_runs[key].add(run_id)

    out: list[dict] = []
    for key, rec in grouped.items():
        rec["run_count"] = len(seen_runs[key])
        rec["source_count"] = len(rec["sources"])
        rec["sources"] = ";".join(sorted(rec["sources"]))
        rec["start_degree_histogram_examples"] = ";".join(sorted(rec["start_degree_histogram_examples"]))
        rec["target_degree_histogram_examples"] = ";".join(sorted(rec["target_degree_histogram_examples"]))
        rec["branchy_residual_motif_examples"] = ";".join(sorted(rec["branchy_residual_motif_examples"]))
        out.append(rec)

    out.sort(
        key=lambda r: (
            -int(r["fail_count"]),
            -int(r["timeout_count"]),
            -int(r["recurrence_count"]),
            -int(r["run_count"]),
            -int(r["internal_edges_max"]),
        )
    )
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "aggregate_key",
        "motif",
        "core_signature",
        "internal_edges_max",
        "reduced_core_internal_edges_min",
        "conflicts_max",
        "min_k",
        "max_k",
        "legal_children_max",
        "plateau_buckets_max",
        "start_branching_nodes_max",
        "target_branching_nodes_max",
        "start_max_branch_subtree_edges_max",
        "target_max_branch_subtree_edges_max",
        "chain_arm_count_max",
        "broom_arm_count_max",
        "branchy_reduction_hits",
        "branchy_forced_cost_max",
        "branchy_residual_edges_min",
        "branchy_exact_bound_hits",
        "branchy_exact_decision_hits",
        "branchy_residual_motif_examples",
        "recurrence_count",
        "success_count",
        "fail_count",
        "timeout_count",
        "free_edge_hits",
        "free_edge_misses",
        "start_degree_histogram_examples",
        "target_degree_histogram_examples",
        "run_count",
        "source_count",
        "sources",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    cfg = parse_args()
    rows = aggregate(load_records([Path(p) for p in cfg.inputs]))
    if cfg.output:
        write_csv(Path(cfg.output), rows)

    for row in rows[: cfg.top]:
        print(
            f"fail={row['fail_count']} timeout={row['timeout_count']} recur={row['recurrence_count']} "
            f"runs={row['run_count']} edges={row['internal_edges_max']} conflicts={row['conflicts_max']} "
            f"motif={row['motif']} key={row['aggregate_key']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
