#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rerun top ranked branchy seeds and capture plateau tree dumps"
    )
    parser.add_argument(
        "--seed-ranking",
        default="results/branchy_deadend_seed_ranking.csv",
        help="CSV from branchy dead-end seed ranking",
    )
    parser.add_argument(
        "--cpp-binary",
        default="./build/flipdist",
        help="Path to flipdist binary",
    )
    parser.add_argument(
        "--output-dir",
        default="results/branchy_selected_dumps",
        help="Directory for plateau tree dump JSONL files",
    )
    parser.add_argument(
        "--summary-csv",
        default="results/branchy_selected_dump_runs.csv",
        help="CSV summary for capture runs",
    )
    parser.add_argument(
        "--top-seeds",
        type=int,
        default=8,
        help="Number of top multi-n seeds to replay",
    )
    parser.add_argument(
        "--timeout-sec",
        type=float,
        default=12.0,
        help="Per-run timeout in seconds",
    )
    parser.add_argument(
        "--max-k-mult",
        type=int,
        default=2,
        help="Set max-k to n * this multiplier",
    )
    parser.add_argument(
        "--dump-min-edges",
        type=int,
        default=11,
        help="Minimum reduced edge count for tree dumps",
    )
    parser.add_argument(
        "--dump-max-edges",
        type=int,
        default=14,
        help="Maximum reduced edge count for tree dumps",
    )
    parser.add_argument(
        "--dump-limit",
        type=int,
        default=80,
        help="Maximum tree dumps per run",
    )
    parser.add_argument(
        "--branch-pairs",
        default="3,3;3,4;4,4",
        help="Semicolon-separated branch-pair filter for tree dumps",
    )
    parser.add_argument(
        "--seeds",
        nargs="*",
        type=int,
        help="Optional explicit seed list override",
    )
    return parser.parse_args()


def load_selected_seeds(path: Path, top_seeds: int, explicit: list[int] | None) -> list[tuple[int, list[int], int]]:
    rows = list(csv.DictReader(path.open()))
    selected: list[tuple[int, list[int], int]] = []
    explicit_set = set(explicit or [])
    for row in rows:
        seed = int(row["seed"])
        if explicit_set and seed not in explicit_set:
            continue
        n_values = [int(x) for x in row["n_values"].split(";") if x]
        if len(n_values) < 2:
            continue
        selected.append((seed, n_values, int(row["mass"])))
        if not explicit_set and len(selected) >= top_seeds:
            break
    return selected


def run_capture(args: argparse.Namespace) -> int:
    root = Path.cwd()
    seed_ranking = root / args.seed_ranking
    cpp_binary = root / args.cpp_binary
    output_dir = root / args.output_dir
    summary_csv = root / args.summary_csv
    output_dir.mkdir(parents=True, exist_ok=True)

    selected = load_selected_seeds(seed_ranking, args.top_seeds, args.seeds)
    results: list[dict[str, str | int]] = []

    for seed, n_values, mass in selected:
        for n in n_values:
            dump_file = output_dir / f"n{n}_s{seed}.jsonl"
            if dump_file.exists():
                dump_file.unlink()
            env = os.environ.copy()
            env.update(
                {
                    "FLIPDIST_ASTAR_ORDER": "1",
                    "FLIPDIST_PROGRESSIVE_S": "0",
                    "FLIPDIST_INCUMBENT_PRUNE": "0",
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_TREES": "1",
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_MIN_EDGES": str(args.dump_min_edges),
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_MAX_EDGES": str(args.dump_max_edges),
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_LIMIT": str(args.dump_limit),
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_BRANCH_PAIRS": args.branch_pairs,
                    "FLIPDIST_PROFILE_PLATEAU_DUMP_FILE": str(dump_file),
                }
            )
            cmd = [
                str(cpp_binary),
                "--case",
                "random",
                "--n",
                str(n),
                "--seed",
                str(seed),
                "--count",
                "1",
                "--max-k",
                str(n * args.max_k_mult),
                "--bfs-cap",
                "1",
            ]
            start = time.time()
            status = "ok"
            try:
                proc = subprocess.run(
                    cmd,
                    cwd=root,
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout_sec,
                )
                if proc.returncode != 0:
                    status = f"rc={proc.returncode}"
            except subprocess.TimeoutExpired:
                status = "timeout"
            elapsed = time.time() - start
            dump_records = 0
            if dump_file.exists():
                with dump_file.open() as fh:
                    dump_records = sum(1 for line in fh if line.strip())
            results.append(
                {
                    "seed": seed,
                    "n": n,
                    "mass": mass,
                    "status": status,
                    "elapsed_sec": f"{elapsed:.3f}",
                    "dump_file": str(dump_file.relative_to(root)),
                    "dump_records": dump_records,
                }
            )
            print(
                f"n={n} seed={seed} mass={mass} status={status} "
                f"elapsed={elapsed:.2f}s dumps={dump_records}"
            )

    with summary_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "seed",
                "n",
                "mass",
                "status",
                "elapsed_sec",
                "dump_file",
                "dump_records",
            ],
        )
        writer.writeheader()
        writer.writerows(results)
    print(f"wrote {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_capture(parse_args()))
