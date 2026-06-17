#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import generate_shared_convex_instances as gen


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate shared convex instances and compare FlipDist vs A* on the exact same inputs"
    )
    parser.add_argument(
        "--case",
        choices=["random", "simple", "comb"],
        default="random",
        help="'simple' is the clearer alias for the older 'comb'.",
    )
    parser.add_argument("--n-min", type=int, required=True)
    parser.add_argument("--n-max", type=int, required=True)
    parser.add_argument("--seed-min", type=int, default=0)
    parser.add_argument("--seed-max", type=int, default=0)
    parser.add_argument("--shared-root", default="results/shared_convex_bench")
    parser.add_argument("--flipdist-binary", default="./build/flipdist")
    parser.add_argument("--flipdist-max-k", type=int, default=None)
    parser.add_argument(
        "--case-timeout-sec",
        type=float,
        default=30.0,
        help="Maximum wall-clock time allowed per (n, seed) case across all solvers",
    )
    parser.add_argument("--astar-binary", default="third_party/AStarFlipDistance/build-nogurobi/A_star_for_flipdistance")
    parser.add_argument(
        "--astar-algos",
        default="simple,combined",
        help="Comma-separated A* algorithms to run (supported: simple,combined)",
    )
    parser.add_argument("--java-check-max-n", type=int, default=15)
    parser.add_argument("--java-time-limit-sec", type=float, default=60.0)
    parser.add_argument("--java-out", default="oracle/java/out")
    parser.add_argument(
        "--output",
        default=None,
        help="CSV output path. If omitted, the script asks whether to save a CSV.",
    )
    parser.add_argument("--plot-output", default=None, help="Optional SVG output path")
    parser.add_argument("--summary-output", default=None, help="Optional summary CSV path")
    parser.add_argument("--print", action="store_true", help="Print per-run progress")
    return parser.parse_args()


def resolve_output_path(cfg: argparse.Namespace) -> Path | None:
    if cfg.output:
        return Path(cfg.output)
    if not sys.stdin.isatty():
        return None

    default_path = (
        f"results/shared_convex_flipdist_vs_astar_{cfg.case}_"
        f"n{cfg.n_min}_{cfg.n_max}_seeds{cfg.seed_min}_{cfg.seed_max}.csv"
    )
    try:
        resp = input("Save CSV results? (y/n): ").strip().lower()
    except EOFError:
        return None
    if resp not in {"y", "yes"}:
        return None
    try:
        value = input(f"CSV path [default: {default_path}]: ").strip()
    except EOFError:
        return Path(default_path)
    return Path(value or default_path)


def parse_algos(raw: str) -> list[str]:
    algos = []
    for part in raw.split(","):
        value = part.strip().lower()
        if not value:
            continue
        if value not in {"simple", "combined"}:
            raise SystemExit(f"Unsupported A* algo: {value}")
        algos.append(value)
    if not algos:
        raise SystemExit("No A* algorithms selected")
    return sorted(set(algos))


def canonical_case_type(case_type: str) -> str:
    return "simple" if case_type == "comb" else case_type


def generate_shared_case(root: Path, case_type: str, n: int, seed: int) -> dict[str, Path | str | int]:
    case_type = canonical_case_type(case_type)
    case_dir = root / f"n_{n}" / f"set_{seed}"
    tri_dir = case_dir / "triangulation"
    tri_dir.mkdir(parents=True, exist_ok=True)

    (pre_a, in_a), (pre_b, in_b) = gen.generate_case(case_type, n, seed)
    tree_a = gen.tree_to_canonical(pre_a, in_a)
    tree_b = gen.tree_to_canonical(pre_b, in_b)

    triangles_a = gen.triangles_from_tree(pre_a, in_a)
    triangles_b = gen.triangles_from_tree(pre_b, in_b)
    gen.validate_triangulation(n, triangles_a)
    gen.validate_triangulation(n, triangles_b)

    points = gen.convex_points(n + 2)
    points_path = case_dir / "points"
    tri_a_path = tri_dir / "t_0"
    tri_b_path = tri_dir / "t_1"
    tree_a_path = case_dir / "tree_a.txt"
    tree_b_path = case_dir / "tree_b.txt"
    manifest_path = case_dir / "manifest.json"

    gen.write_points(points_path, points)
    gen.write_triangles(tri_a_path, triangles_a)
    gen.write_triangles(tri_b_path, triangles_b)
    tree_a_path.write_text(tree_a + "\n", encoding="utf-8")
    tree_b_path.write_text(tree_b + "\n", encoding="utf-8")

    manifest = {
        "case_type": case_type,
        "n": n,
        "seed": seed,
        "tree_a": tree_a,
        "tree_b": tree_b,
        "points_path": str(points_path),
        "triangulation_1_path": str(tri_a_path),
        "triangulation_2_path": str(tri_b_path),
        "flipdist_tree_a_file": str(tree_a_path),
        "flipdist_tree_b_file": str(tree_b_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    return {
        "case_dir": case_dir,
        "points_path": points_path,
        "tri_a_path": tri_a_path,
        "tri_b_path": tri_b_path,
        "tree_a_path": tree_a_path,
        "tree_b_path": tree_b_path,
        "manifest_path": manifest_path,
    }


def run_flipdist(
    binary: str,
    tree_a_path: Path,
    tree_b_path: Path,
    max_k: int,
    timeout_sec: float | None = None,
) -> dict[str, object]:
    cmd = [
        binary,
        "--tree-a-file",
        str(tree_a_path),
        "--tree-b-file",
        str(tree_b_path),
        "--max-k",
        str(max_k),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"FlipDist failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    if len(lines) != 2:
        raise RuntimeError(f"FlipDist expected 2 JSON lines, saw {len(lines)}")
    recs = [json.loads(line) for line in lines]
    a = next((r for r in recs if r.get("direction") == "a->b"), None)
    b = next((r for r in recs if r.get("direction") == "b->a"), None)
    if a is None or b is None:
        raise RuntimeError("FlipDist output is missing a direction")
    return {
        "flipdist_a2b_ms": float(a.get("time_ms", 0.0)),
        "flipdist_b2a_ms": float(b.get("time_ms", 0.0)),
        "flipdist_max_ms": max(float(a.get("time_ms", 0.0)), float(b.get("time_ms", 0.0))),
        "flipdist_mean_ms": (float(a.get("time_ms", 0.0)) + float(b.get("time_ms", 0.0))) / 2.0,
        "flipdist_a2b_dist": int(a.get("distance", -1)),
        "flipdist_b2a_dist": int(b.get("distance", -1)),
        "flipdist_a2b_status": str(a.get("status", "")),
        "flipdist_b2a_status": str(b.get("status", "")),
        "flipdist_ok_both": (a.get("status") == "ok" and b.get("status") == "ok"),
        "flipdist_distance_agrees": int(a.get("distance", -1) == b.get("distance", -1)),
    }


def run_astar(
    binary: str,
    points: Path,
    tri_a: Path,
    tri_b: Path,
    algo: str,
    timeout_sec: float | None = None,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory() as tmpdir:
        output_name = f"astar_{algo}.json"
        cmd = [
            binary,
            str(points),
            str(tri_a),
            str(tri_b),
            algo,
            tmpdir,
            output_name,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=timeout_sec)
        if proc.returncode != 0:
            raise RuntimeError(f"A* {algo} failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
        payload = json.loads((Path(tmpdir) / output_name).read_text())
        runtime_ms = float(payload.get("runtime", 0.0))
        runtime_plot_ms = runtime_ms if runtime_ms > 0.0 else 0.001
        return {
            f"astar_{algo}_runtime_ms": runtime_ms,
            f"astar_{algo}_mean_ms": runtime_plot_ms,
            f"astar_{algo}_median_ms": runtime_plot_ms,
            f"astar_{algo}_p95_ms": runtime_plot_ms,
            f"astar_{algo}_pairs": 1,
            f"astar_{algo}_timeout_rows": 0,
            f"astar_{algo}_distance": int(payload.get("flip distance", -1)),
            f"astar_{algo}_closed_nodes": int(payload.get("closed_nodes", -1)),
            f"astar_{algo}_opened_nodes": int(payload.get("opened_nodes", -1)),
        }


def run_java(
    java_out: str,
    points: Path,
    tri_a: Path,
    tri_b: Path,
    time_limit_sec: float,
    timeout_sec: float | None = None,
) -> dict[str, object]:
    cmd = [
        "java",
        "-cp",
        java_out,
        "TriangulationFileDistanceCli",
        str(points),
        str(tri_a),
        str(tri_b),
        "--time-limit-sec",
        str(time_limit_sec),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(f"Java BFS failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError("Java BFS produced no output")
    payload = json.loads(lines[-1])
    return {
        "java_distance": int(payload.get("distance", -1)),
        "java_time_ms": float(payload.get("time_ms", 0.0)),
        "java_status": str(payload.get("status", "")),
    }


def maybe_plot(cfg: argparse.Namespace, csv_path: Path) -> None:
    if not cfg.plot_output and not cfg.summary_output:
        return
    temp_plot_path: Path | None = None
    plot_output = cfg.plot_output
    if not plot_output:
        with tempfile.NamedTemporaryFile(suffix=".svg", delete=False) as tmp:
            temp_plot_path = Path(tmp.name)
            plot_output = str(temp_plot_path)
    cmd = [
        sys.executable,
        "tools/plot_flipdist_vs_astar_compare.py",
        "--input",
        str(csv_path),
        "--output",
        str(plot_output),
        "--title",
        "Shared Convex FlipDist vs A* Runtime Growth",
    ]
    if cfg.summary_output:
        cmd.extend(["--summary-output", str(cfg.summary_output)])
    try:
        run_kwargs = {"check": True}
        if temp_plot_path is not None:
            run_kwargs.update({"capture_output": True, "text": True})
        subprocess.run(cmd, **run_kwargs)
    finally:
        if temp_plot_path is not None:
            temp_plot_path.unlink(missing_ok=True)


def remaining_timeout(deadline: float | None) -> float | None:
    if deadline is None:
        return None
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("per-case deadline exceeded")
    return max(0.001, remaining)


def main() -> int:
    cfg = parse_args()
    cfg.case = canonical_case_type(cfg.case)
    if cfg.n_max < cfg.n_min:
        raise SystemExit("--n-max must be >= --n-min")
    if cfg.seed_max < cfg.seed_min:
        raise SystemExit("--seed-max must be >= --seed-min")

    output_path = resolve_output_path(cfg)
    algos = parse_algos(cfg.astar_algos)
    shared_root = Path(cfg.shared_root).resolve()
    shared_root.mkdir(parents=True, exist_ok=True)
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    if cfg.plot_output:
        Path(cfg.plot_output).parent.mkdir(parents=True, exist_ok=True)
    if cfg.summary_output:
        Path(cfg.summary_output).parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    failures = 0
    timeouts = 0

    for n in range(cfg.n_min, cfg.n_max + 1):
        for seed in range(cfg.seed_min, cfg.seed_max + 1):
            if cfg.print:
                print(f"Shared convex sweep: n={n} seed={seed} case={cfg.case}")
            row: dict[str, object] = {
                "n": n,
                "seed": seed,
                "case_type": cfg.case,
            }
            deadline = None if cfg.case_timeout_sec <= 0 else (time.monotonic() + cfg.case_timeout_sec)
            try:
                case_paths = generate_shared_case(shared_root, cfg.case, n, seed)
                max_k = cfg.flipdist_max_k if cfg.flipdist_max_k is not None else max(1, 3 * n + 10)
                row.update(run_flipdist(cfg.flipdist_binary,
                                        case_paths["tree_a_path"],
                                        case_paths["tree_b_path"],
                                        max_k,
                                        remaining_timeout(deadline)))
                row["shared_case_dir"] = str(case_paths["case_dir"])
                row["points_path"] = str(case_paths["points_path"])
                row["triangulation_1_path"] = str(case_paths["tri_a_path"])
                row["triangulation_2_path"] = str(case_paths["tri_b_path"])

                for algo in ("simple", "combined"):
                    if algo in algos:
                        row.update(run_astar(cfg.astar_binary,
                                             case_paths["points_path"],
                                             case_paths["tri_a_path"],
                                             case_paths["tri_b_path"],
                                             algo,
                                             remaining_timeout(deadline)))
                    else:
                        row[f"astar_{algo}_mean_ms"] = ""
                        row[f"astar_{algo}_median_ms"] = ""
                        row[f"astar_{algo}_p95_ms"] = ""
                        row[f"astar_{algo}_runtime_ms"] = ""
                        row[f"astar_{algo}_pairs"] = ""
                        row[f"astar_{algo}_timeout_rows"] = ""
                        row[f"astar_{algo}_distance"] = ""
                        row[f"astar_{algo}_closed_nodes"] = ""
                        row[f"astar_{algo}_opened_nodes"] = ""

                if n <= cfg.java_check_max_n:
                    row.update(run_java(cfg.java_out,
                                        case_paths["points_path"],
                                        case_paths["tri_a_path"],
                                        case_paths["tri_b_path"],
                                        cfg.java_time_limit_sec,
                                        remaining_timeout(deadline)))
                else:
                    row["java_distance"] = ""
                    row["java_time_ms"] = ""
                    row["java_status"] = ""

                flipdist_dist = row.get("flipdist_a2b_dist", -1)
                for algo in ("simple", "combined"):
                    astar_key = f"astar_{algo}_distance"
                    match_key = f"match_flipdist_vs_astar_{algo}"
                    astar_dist = row.get(astar_key, "")
                    if astar_dist == "":
                        row[match_key] = ""
                    else:
                        row[match_key] = int(flipdist_dist == astar_dist)
                java_dist = row.get("java_distance", "")
                row["match_flipdist_vs_java"] = "" if java_dist == "" else int(flipdist_dist == java_dist)
            except subprocess.TimeoutExpired as exc:
                timeouts += 1
                failures += 1
                cmd0 = Path(exc.cmd[0]).name if exc.cmd else "unknown"
                row["status"] = f"timeout:{cmd0}"
            except TimeoutError as exc:
                timeouts += 1
                failures += 1
                row["status"] = f"timeout:{exc}"
            except Exception as exc:
                failures += 1
                row["status"] = f"error:{exc}"
            rows.append(row)

    if not rows:
        raise SystemExit("No rows produced")

    fieldnames = [
        "n",
        "seed",
        "case_type",
        "shared_case_dir",
        "points_path",
        "triangulation_1_path",
        "triangulation_2_path",
        "flipdist_a2b_ms",
        "flipdist_b2a_ms",
        "flipdist_max_ms",
        "flipdist_mean_ms",
        "flipdist_a2b_dist",
        "flipdist_b2a_dist",
        "flipdist_a2b_status",
        "flipdist_b2a_status",
        "flipdist_ok_both",
        "flipdist_distance_agrees",
        "astar_simple_pairs",
        "astar_simple_runtime_ms",
        "astar_simple_mean_ms",
        "astar_simple_median_ms",
        "astar_simple_p95_ms",
        "astar_simple_timeout_rows",
        "astar_simple_distance",
        "astar_simple_closed_nodes",
        "astar_simple_opened_nodes",
        "astar_combined_pairs",
        "astar_combined_runtime_ms",
        "astar_combined_mean_ms",
        "astar_combined_median_ms",
        "astar_combined_p95_ms",
        "astar_combined_timeout_rows",
        "astar_combined_distance",
        "astar_combined_closed_nodes",
        "astar_combined_opened_nodes",
        "java_distance",
        "java_time_ms",
        "java_status",
        "match_flipdist_vs_astar_simple",
        "match_flipdist_vs_astar_combined",
        "match_flipdist_vs_java",
        "status",
    ]

    if output_path:
        with output_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in fieldnames})
        maybe_plot(cfg, output_path)

    mismatch_simple = sum(1 for row in rows if row.get("match_flipdist_vs_astar_simple") == 0)
    mismatch_combined = sum(1 for row in rows if row.get("match_flipdist_vs_astar_combined") == 0)
    mismatch_java = sum(1 for row in rows if row.get("match_flipdist_vs_java") == 0)

    if cfg.print:
        if output_path:
            print(f"Wrote {output_path} ({len(rows)} rows)")
        else:
            print(f"CSV not saved ({len(rows)} rows)")
        if cfg.plot_output:
            print(f"Wrote {cfg.plot_output}")
        print(
            f"Mismatches: simple={mismatch_simple} combined={mismatch_combined} java={mismatch_java} failures={failures} timeouts={timeouts}"
        )

    if mismatch_simple or mismatch_combined or mismatch_java or failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
