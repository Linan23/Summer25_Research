#!/usr/bin/env python3
"""Sweep FlipDist runtime limits across n and seeds.

This is a performance harness (not a correctness oracle). It runs the C++
`flipdist_fast` CLI on random or comb cases across a range of n and seeds,
records per-direction status/time, and prints a per-n summary.

Example:
  python3 scripts/sweep_flipdist_limits.py --case random --n-min 15 --n-max 25 \\
    --seed-min 0 --seed-max 9 --timeout-sec 5 \\
    --output results/limit_sweep_random_n15_25_seeds0_9_t5s.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path


def int_list(arg: str) -> list[int]:
    if not arg:
        return []
    try:
        return [int(x) for x in arg.split(",") if x.strip()]
    except ValueError as exc:  # pragma: no cover - defensive
        raise argparse.ArgumentTypeError("must be a comma-separated list of ints") from exc


def float_list(arg: str) -> list[float]:
    if not arg:
        return []
    try:
        return [float(x) for x in arg.split(",") if x.strip()]
    except ValueError as exc:  # pragma: no cover - defensive
        raise argparse.ArgumentTypeError("must be a comma-separated list of floats") from exc


def parse_args():
    p = argparse.ArgumentParser(description="Sweep FlipDist performance across n/seeds.")
    p.add_argument("--case", dest="case_type", choices=["random", "comb"], default="random")
    p.add_argument("--n-min", type=int, required=True)
    p.add_argument("--n-max", type=int, required=True)
    p.add_argument("--seed-min", type=int, default=0, help="Only used for --case random")
    p.add_argument("--seed-max", type=int, default=0, help="Only used for --case random")
    p.add_argument("--timeout-sec", type=float, default=5.0, help="Base timeout for n <= n_threshold (or all n if no threshold).")
    p.add_argument(
        "--high-timeout-sec",
        type=float,
        default=None,
        help="Base timeout for n > n_threshold (defaults to --timeout-sec if unset).",
    )
    p.add_argument("--n-threshold", type=int, default=None, help="Switch to high budgets above this n.")
    p.add_argument("--cpp-binary", default="./build/flipdist")
    p.add_argument("--bfs-cap", type=int, default=1, help="Keep tiny; this is a FlipDist perf sweep.")
    p.add_argument("--max-k-mult", type=int, default=2, help="Base max_k multiplier for n <= n_threshold (or all n if no threshold).")
    p.add_argument(
        "--high-max-k-mult",
        type=int,
        default=None,
        help="Base max_k multiplier for n > n_threshold (defaults to --max-k-mult if unset).",
    )
    p.add_argument(
        "--retry-max-k-mults",
        type=int_list,
        default=[],
        help="Optional comma-separated multipliers for retry attempts (e.g., 3,4).",
    )
    p.add_argument(
        "--retry-timeout-mult",
        type=float,
        default=2.0,
        help="Multiply --timeout-sec by this for retry attempts.",
    )
    p.add_argument(
        "--retry-timeout-mults",
        type=float_list,
        default=None,
        help="Optional comma-separated timeout multipliers aligned with --retry-max-k-mults.",
    )
    p.add_argument(
        "--retry-bfs-cap",
        type=int,
        default=None,
        help="Override BFS cap on retries; defaults to --bfs-cap.",
    )
    p.add_argument(
        "--output",
        default=None,
        help="CSV output path. If omitted, the script asks whether to save a CSV.",
    )
    return p.parse_args()


def resolve_output_path(cfg) -> Path | None:
    if cfg.output:
        return Path(cfg.output)
    if not sys.stdin.isatty():
        return None

    default_path = (
        f"results/flipdist_limits_{cfg.case_type}_n{cfg.n_min}_{cfg.n_max}_"
        f"seeds{cfg.seed_min}_{cfg.seed_max}.csv"
    )
    try:
        resp = input("Save CSV results? (y/n): ").strip().lower()
    except EOFError:
        return None
    if resp not in {"y", "yes"}:
        return None
    try:
        path = input(f"CSV path [default: {default_path}]: ").strip()
    except EOFError:
        return Path(default_path)
    return Path(path or default_path)


def invoke_flipdist(cfg, n: int, seed: int, max_k: int, timeout_sec: float, bfs_cap: int, attempt: str):
    cmd = [
        cfg.cpp_binary,
        "--case",
        cfg.case_type,
        "--n",
        str(n),
        "--count",
        "1",
        "--max-k",
        str(max_k),
        "--bfs-cap",
        str(bfs_cap),
    ]
    if cfg.case_type == "random":
        cmd.extend(["--seed", str(seed)])

    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec, check=False)
    except subprocess.TimeoutExpired:
        rows = [
            {
                "case_type": cfg.case_type,
                "n": n,
                "seed": seed,
                "direction": direction,
                "distance": "",
                "status": "timeout",
                "time_ms": f"{timeout_sec * 1000.0:.3f}",
                "max_k": max_k,
                "attempt": attempt,
            }
            for direction in ("a->b", "b->a")
        ]
        return rows, "timeout"

    elapsed_ms = (time.time() - t0) * 1000.0
    if proc.returncode != 0:
        rows = [
            {
                "case_type": cfg.case_type,
                "n": n,
                "seed": seed,
                "direction": direction,
                "distance": "",
                "status": f"error(rc={proc.returncode})",
                "time_ms": f"{elapsed_ms:.3f}",
                "max_k": max_k,
                "attempt": attempt,
            }
            for direction in ("a->b", "b->a")
        ]
        return rows, "error"

    out_lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    if len(out_lines) != 2:
        rows = [
            {
                "case_type": cfg.case_type,
                "n": n,
                "seed": seed,
                "direction": direction,
                "distance": "",
                "status": f"bad_output(lines={len(out_lines)})",
                "time_ms": f"{elapsed_ms:.3f}",
                "max_k": max_k,
                "attempt": attempt,
            }
            for direction in ("a->b", "b->a")
        ]
        return rows, "bad_output"

    rows = []
    statuses = set()
    for line in out_lines:
        try:
            recj = json.loads(line)
        except json.JSONDecodeError:
            return [
                {
                    "case_type": cfg.case_type,
                    "n": n,
                    "seed": seed,
                    "direction": direction,
                    "distance": "",
                    "status": "json_error",
                    "time_ms": f"{elapsed_ms:.3f}",
                    "max_k": max_k,
                    "attempt": attempt,
                }
                for direction in ("a->b", "b->a")
            ], "json_error"

        status = recj.get("status", "")
        statuses.add(status)
        rows.append(
            {
                "case_type": recj.get("case_type"),
                "n": recj.get("n"),
                "seed": recj.get("seed"),
                "direction": recj.get("direction"),
                "distance": recj.get("distance"),
                "status": status,
                "time_ms": f"{float(recj.get('time_ms', 0.0)):.3f}",
                "max_k": recj.get("max_k"),
                "attempt": attempt,
            }
        )

    overall = "ok" if statuses == {"ok"} else ",".join(sorted(statuses))
    return rows, overall


def main() -> int:
    cfg = parse_args()
    out_path = resolve_output_path(cfg)
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)

    seeds: list[int]
    if cfg.case_type == "random":
        if cfg.seed_max < cfg.seed_min:
            raise SystemExit("--seed-max must be >= --seed-min")
        seeds = list(range(cfg.seed_min, cfg.seed_max + 1))
    else:
        seeds = [-1]

    fieldnames = ["case_type", "n", "seed", "direction", "distance", "status", "time_ms", "max_k", "attempt"]
    rows: list[dict] = []

    retry_max_k_mults = sorted(set(m for m in cfg.retry_max_k_mults if m > 0))
    if cfg.retry_timeout_mults is not None:
        if len(cfg.retry_timeout_mults) != len(retry_max_k_mults):
            raise SystemExit("--retry-timeout-mults must match the length of --retry-max-k-mults")
        retry_timeout_mults = cfg.retry_timeout_mults
    else:
        retry_timeout_mults = [cfg.retry_timeout_mult] * len(retry_max_k_mults)
    retry_bfs_cap = cfg.retry_bfs_cap if cfg.retry_bfs_cap is not None else cfg.bfs_cap

    start = time.time()
    out_file = out_path.open("w", newline="") if out_path else None
    writer = csv.DictWriter(out_file, fieldnames=fieldnames) if out_file else None
    if writer:
        writer.writeheader()

    try:
        for n in range(cfg.n_min, cfg.n_max + 1):
            base_max_k_mult = cfg.max_k_mult
            base_timeout = cfg.timeout_sec
            if cfg.n_threshold is not None and n > cfg.n_threshold:
                base_max_k_mult = cfg.high_max_k_mult if cfg.high_max_k_mult is not None else cfg.max_k_mult
                base_timeout = cfg.high_timeout_sec if cfg.high_timeout_sec is not None else cfg.timeout_sec

            max_k = base_max_k_mult * n
            n_rows_start = len(rows)

            for seed in seeds:
                attempts = [("base", max_k, base_timeout, cfg.bfs_cap)]
                seen = {(max_k, base_timeout, cfg.bfs_cap)}
                for mult, t_mult in zip(retry_max_k_mults, retry_timeout_mults):
                    attempt_max_k = mult * n
                    attempt_timeout = base_timeout * t_mult
                    key = (attempt_max_k, attempt_timeout, retry_bfs_cap)
                    if key in seen:
                        continue
                    seen.add(key)
                    attempts.append(
                        (
                            f"retry-m{mult}",
                            attempt_max_k,
                            attempt_timeout,
                            retry_bfs_cap,
                        )
                    )

                chosen_rows = []
                chosen_status = ""
                for idx, (attempt_label, attempt_max_k, attempt_timeout, attempt_bfs_cap) in enumerate(attempts):
                    attempt_rows, status = invoke_flipdist(
                        cfg, n, seed, attempt_max_k, attempt_timeout, attempt_bfs_cap, attempt_label
                    )
                    if status == "ok":
                        chosen_rows, chosen_status = attempt_rows, status
                        break
                    chosen_rows, chosen_status = attempt_rows, status
                    if idx + 1 < len(attempts):
                        print(
                            f"retrying n={n} seed={seed} after status={status} "
                            f"with max_k={attempts[idx+1][1]} timeout={attempts[idx+1][2]}s "
                            f"bfs_cap={attempts[idx+1][3]}",
                            flush=True,
                        )

                for rec in chosen_rows:
                    rows.append(rec)
                    if writer:
                        writer.writerow(rec)
                if out_file:
                    out_file.flush()

            # Print a simple per-n summary (max over collected rows for that n).
            r_n = rows[n_rows_start:]
            ok = sum(1 for r in r_n if r["status"] == "ok")
            timeout = sum(1 for r in r_n if r["status"] == "timeout")
            not_found = sum(1 for r in r_n if r["status"] == "not_found")
            times = sorted(float(r["time_ms"]) for r in r_n)
            max_ms = times[-1] if times else float("nan")
            p95 = times[int(0.95 * (len(times) - 1))] if times else float("nan")
            print(f"n={n} ok={ok}/{len(r_n)} timeout={timeout} not_found={not_found} max_ms={max_ms:.1f} p95_ms={p95:.1f}", flush=True)
    finally:
        if out_file:
            out_file.close()

    elapsed = time.time() - start
    if out_path:
        print(f"Wrote {out_path} ({len(rows)} rows), wall={elapsed:.1f}s")
    else:
        print(f"Completed sweep without CSV save ({len(rows)} rows), wall={elapsed:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
