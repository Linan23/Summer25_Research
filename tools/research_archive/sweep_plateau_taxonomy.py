#!/usr/bin/env python3
"""Run plateau-profile sweeps and summarize recurring empty-S states."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep FlipDist with plateau profiling enabled and rank recurring plateau states."
    )
    parser.add_argument("--n-min", type=int, required=True)
    parser.add_argument("--n-max", type=int, required=True)
    parser.add_argument("--seed-min", type=int, default=0)
    parser.add_argument("--seed-max", type=int, default=0)
    parser.add_argument("--cpp-binary", default="./build/flipdist")
    parser.add_argument(
        "--case",
        dest="case_type",
        choices=["random", "simple", "comb"],
        default="random",
        help="'simple' is the clearer alias for the older 'comb'.",
    )
    parser.add_argument("--timeout-sec", type=float, default=10.0)
    parser.add_argument("--max-k-mult", type=int, default=2)
    parser.add_argument("--raw-dir", required=True, help="Directory for raw profile outputs")
    parser.add_argument("--summary-csv", required=True, help="CSV output from analyze_plateau_profile.py")
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--enable-core-reduce", action="store_true")
    parser.add_argument("--core-db-m", type=int, default=8)
    parser.add_argument("--enable-motif-reduce", action="store_true")
    parser.add_argument("--motif-m", type=int, default=7)
    parser.add_argument("--enable-branchy-core-lb", action="store_true")
    parser.add_argument("--branchy-core-m", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    cfg = parse_args()
    if cfg.case_type == "comb":
        cfg.case_type = "simple"
    raw_dir = Path(cfg.raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    raw_files: list[str] = []
    for n in range(cfg.n_min, cfg.n_max + 1):
        seeds = range(cfg.seed_min, cfg.seed_max + 1) if cfg.case_type == "random" else [-1]
        for seed in seeds:
            out_path = raw_dir / f"plateau_n{n}_s{seed}.jsonl"
            cmd = [
                cfg.cpp_binary,
                "--case",
                cfg.case_type,
                "--n",
                str(n),
                "--count",
                "1",
                "--max-k",
                str(max(1, cfg.max_k_mult * n)),
                "--bfs-cap",
                "1",
            ]
            if cfg.case_type == "random":
                cmd.extend(["--seed", str(seed)])

            env = os.environ.copy()
            env.update(
                {
                    "FLIPDIST_PROFILE": "1",
                    "FLIPDIST_PROFILE_PLATEAU": "1",
                    "FLIPDIST_PROFILE_PLATEAU_FILE": str(out_path),
                    "FLIPDIST_ASTAR_ORDER": "1",
                    "FLIPDIST_PROGRESSIVE_S": "0",
                    "FLIPDIST_INCUMBENT_PRUNE": "0",
                }
            )
            if cfg.enable_core_reduce:
                env["FLIPDIST_EMPTY_S_CORE_REDUCE"] = "1"
                env["FLIPDIST_EMPTY_S_CORE_DB_M"] = str(cfg.core_db_m)
            if cfg.enable_motif_reduce:
                env["FLIPDIST_EMPTY_S_MOTIF_REDUCE"] = "1"
                env["FLIPDIST_EMPTY_S_EXACT_MOTIF_M"] = str(cfg.motif_m)
            if cfg.enable_branchy_core_lb:
                env["FLIPDIST_BRANCHY_CORE_LB"] = "1"
                env["FLIPDIST_BRANCHY_CORE_EXACT_M"] = str(cfg.branchy_core_m)

            sys.stderr.write(f"Plateau sweep: n={n} seed={seed}\n")
            try:
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=cfg.timeout_sec,
                    env=env,
                    check=False,
                )
            except subprocess.TimeoutExpired as ex:
                partial = ex.stdout or ""
                if out_path.exists() and out_path.stat().st_size > 0:
                    raw_files.append(str(out_path))
                elif partial:
                    out_path.write_text(partial)
                    raw_files.append(str(out_path))
                continue
            if out_path.exists() and out_path.stat().st_size > 0:
                raw_files.append(str(out_path))
            else:
                out_path.write_text(proc.stdout)
                raw_files.append(str(out_path))

    if not raw_files:
        sys.stderr.write("No raw plateau outputs produced.\n")
        return 1

    analyze_cmd = [
        sys.executable,
        "tools/research_archive/analyze_plateau_profile.py",
        "--inputs",
        *raw_files,
        "--output",
        cfg.summary_csv,
        "--top",
        str(cfg.top),
    ]
    proc = subprocess.run(analyze_cmd, check=False)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
