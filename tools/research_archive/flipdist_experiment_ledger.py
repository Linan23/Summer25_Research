#!/usr/bin/env python3
"""
Record FlipDist experiment outcomes and guard against regressions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Append experiment metrics to a ledger and fail on regressions."
    )
    parser.add_argument("--change-id", required=True, help="Unique change label, e.g. astar-order-v2")
    parser.add_argument("--seed-band", required=True, help="Seed scope label, e.g. n22_25_s0_100_t10")
    parser.add_argument("--input", required=True, help="Input CSV from a sweep/parity run")
    parser.add_argument("--ledger", required=True, help="Ledger file (.jsonl or .csv)")
    parser.add_argument(
        "--auto-revert-cmd",
        default="",
        help="Optional shell command to run when regression is detected",
    )
    parser.add_argument(
        "--allow-equal-time",
        action="store_true",
        help="Treat equal p95 as non-regression",
    )
    return parser.parse_args()


def maybe_float(value: Any) -> float | None:
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() == "na":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def detect_time_columns(headers: list[str]) -> list[str]:
    preferred = [
        "time_ms",
        "flipdist_max_ms",
        "flipdist_time_ms",
        "flipdist_runtime_ms",
        "flipdist_a2b_ms",
        "flipdist_b2a_ms",
    ]
    found = [c for c in preferred if c in headers]
    if found:
        return found
    return [c for c in headers if c.endswith("_ms") and "astar" not in c and "java" not in c]


def row_timeout(row: dict[str, str]) -> bool:
    status_cols = [
        "status",
        "flipdist_a2b_status",
        "flipdist_b2a_status",
        "status_flipdist",
    ]
    for col in status_cols:
        v = row.get(col, "")
        if "timeout" in v.lower():
            return True
    return False


def row_mismatch(row: dict[str, str]) -> bool:
    mismatch_cols = [
        "match_flipdist_vs_java",
        "match_flipdist_vs_astar_simple",
        "match_flipdist_vs_astar_combined",
        "match",
    ]
    for col in mismatch_cols:
        if col not in row:
            continue
        v = row.get(col, "").strip()
        if v == "":
            continue
        if v in {"0", "False", "false"}:
            return True
    return False


def summarize_csv(path: Path) -> dict[str, Any]:
    rows = list(csv.DictReader(path.open()))
    if not rows:
        raise SystemExit(f"Empty CSV: {path}")

    headers = list(rows[0].keys())
    time_cols = detect_time_columns(headers)
    times: list[float] = []
    for row in rows:
        row_times = [maybe_float(row.get(c)) for c in time_cols]
        row_times = [t for t in row_times if t is not None]
        if row_times:
            times.append(max(row_times))

    timeout_count = sum(1 for r in rows if row_timeout(r))
    mismatch_count = sum(1 for r in rows if row_mismatch(r))
    median_ms = statistics.median(times) if times else math.nan
    p95_ms = statistics.quantiles(times, n=20)[18] if len(times) >= 20 else (max(times) if times else math.nan)

    return {
        "rows": len(rows),
        "timeout_count": timeout_count,
        "mismatch_count": mismatch_count,
        "median_ms": median_ms,
        "p95_ms": p95_ms,
        "time_columns": ",".join(time_cols),
    }


def load_ledger(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    if path.suffix.lower() == ".jsonl":
        out: list[dict[str, Any]] = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            out.append(json.loads(line))
        return out

    rows = list(csv.DictReader(path.open()))
    parsed: list[dict[str, Any]] = []
    for row in rows:
        parsed.append(
            {
                **row,
                "timeout_count": int(row.get("timeout_count", "0") or 0),
                "mismatch_count": int(row.get("mismatch_count", "0") or 0),
                "median_ms": float(row.get("median_ms", "nan")),
                "p95_ms": float(row.get("p95_ms", "nan")),
            }
        )
    return parsed


def append_ledger(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".jsonl":
        with path.open("a") as fh:
            fh.write(json.dumps(record, sort_keys=True) + "\n")
        return

    exists = path.exists()
    fieldnames = [
        "timestamp_utc",
        "change_id",
        "seed_band",
        "input_csv",
        "rows",
        "timeout_count",
        "mismatch_count",
        "median_ms",
        "p95_ms",
        "time_columns",
    ]
    with path.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow({k: record.get(k, "") for k in fieldnames})


def last_baseline(records: list[dict[str, Any]], seed_band: str, change_id: str) -> dict[str, Any] | None:
    filtered = [r for r in records if r.get("seed_band") == seed_band and r.get("change_id") != change_id]
    if not filtered:
        return None
    return filtered[-1]


def is_regression(current: dict[str, Any], baseline: dict[str, Any], allow_equal_time: bool) -> bool:
    if current["mismatch_count"] > baseline["mismatch_count"]:
        return True
    if current["timeout_count"] > baseline["timeout_count"]:
        return True
    cur_p95 = current["p95_ms"]
    base_p95 = baseline["p95_ms"]
    if math.isfinite(cur_p95) and math.isfinite(base_p95):
        if allow_equal_time:
            return cur_p95 > base_p95
        return cur_p95 >= base_p95
    return False


def main() -> int:
    cfg = parse_args()
    input_path = Path(cfg.input)
    ledger_path = Path(cfg.ledger)
    summary = summarize_csv(input_path)

    record = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "change_id": cfg.change_id,
        "seed_band": cfg.seed_band,
        "input_csv": str(input_path),
        **summary,
    }

    existing = load_ledger(ledger_path)
    baseline = last_baseline(existing, cfg.seed_band, cfg.change_id)

    append_ledger(ledger_path, record)
    print(
        f"Ledger append: change={cfg.change_id} seed_band={cfg.seed_band} "
        f"rows={record['rows']} timeout={record['timeout_count']} "
        f"mismatch={record['mismatch_count']} median_ms={record['median_ms']:.3f} "
        f"p95_ms={record['p95_ms']:.3f}"
    )

    if baseline is None:
        print("No baseline found for this seed band; recorded as baseline.")
        return 0

    if not is_regression(record, baseline, cfg.allow_equal_time):
        print("No regression vs latest baseline.")
        return 0

    print("Regression detected vs latest baseline.")
    print(
        f"Baseline: timeout={baseline['timeout_count']} mismatch={baseline['mismatch_count']} "
        f"p95_ms={float(baseline['p95_ms']):.3f}"
    )
    if cfg.auto_revert_cmd:
        print(f"Running auto-revert command: {cfg.auto_revert_cmd}")
        proc = subprocess.run(cfg.auto_revert_cmd, shell=True, text=True)
        return proc.returncode if proc.returncode != 0 else 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
