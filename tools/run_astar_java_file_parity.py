#!/usr/bin/env python3
import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare A* single-instance distance against Java exact file-based BFS on the same triangulation files"
    )
    parser.add_argument("--n-min", type=int, required=True)
    parser.add_argument("--n-max", type=int, required=True)
    parser.add_argument("--seed-min", type=int, default=0)
    parser.add_argument("--seed-max", type=int, default=0)
    parser.add_argument("--pair-i", type=int, default=0)
    parser.add_argument("--pair-j", type=int, default=1)
    parser.add_argument("--algo", default="simple", choices=["simple", "combined", "decomposition", "eppstein", "bfs"])
    parser.add_argument(
        "--astar-binary",
        default="third_party/AStarFlipDistance/build-nogurobi/A_star_for_flipdistance",
    )
    parser.add_argument(
        "--data-root",
        default="third_party/AStarFlipDistance/data/random_experiments_paper",
    )
    parser.add_argument(
        "--java-src",
        default="oracle/java/src/TriangulationFileDistanceCli.java",
    )
    parser.add_argument(
        "--java-out",
        default="oracle/java/out",
    )
    parser.add_argument("--java-time-limit-sec", type=float, default=60.0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--print", action="store_true", dest="do_print")
    return parser.parse_args()


def compile_java(java_src: Path, java_out: Path) -> None:
    java_out.mkdir(parents=True, exist_ok=True)
    cmd = ["javac", "-d", str(java_out), str(java_src)]
    subprocess.run(cmd, check=True)


def run_astar(astar_binary: Path, points: Path, tri1: Path, tri2: Path, algo: str) -> dict:
    with tempfile.TemporaryDirectory() as tmpdir:
        output_name = "astar_single.json"
        cmd = [
            str(astar_binary),
            str(points),
            str(tri1),
            str(tri2),
            algo,
            tmpdir,
            output_name,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"A* failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
        payload = json.loads((Path(tmpdir) / output_name).read_text())
        return payload


def run_java(java_out: Path, points: Path, tri1: Path, tri2: Path, time_limit_sec: float) -> dict:
    cmd = [
        "java",
        "-cp",
        str(java_out),
        "TriangulationFileDistanceCli",
        str(points),
        str(tri1),
        str(tri2),
        "--time-limit-sec",
        str(time_limit_sec),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Java failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")
    stdout = proc.stdout.strip().splitlines()
    if not stdout:
        raise RuntimeError("Java produced no output")
    return json.loads(stdout[-1])


def main() -> int:
    cfg = parse_args()

    repo_root = Path.cwd()
    astar_binary = (repo_root / cfg.astar_binary).resolve()
    data_root = (repo_root / cfg.data_root).resolve()
    java_src = (repo_root / cfg.java_src).resolve()
    java_out = (repo_root / cfg.java_out).resolve()
    output_path = (repo_root / cfg.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    compile_java(java_src, java_out)

    rows = []
    mismatches = 0
    failures = 0
    incomplete = 0

    for n in range(cfg.n_min, cfg.n_max + 1):
        for seed in range(cfg.seed_min, cfg.seed_max + 1):
            set_dir = data_root / f"n_{n}" / f"set_{seed}"
            points = set_dir / "points"
            tri1 = set_dir / "triangulation" / f"t_{cfg.pair_i}"
            tri2 = set_dir / "triangulation" / f"t_{cfg.pair_j}"

            row = {
                "n": n,
                "seed": seed,
                "pair": f"t_{cfg.pair_i}->t_{cfg.pair_j}",
                "algo": cfg.algo,
                "astar_distance": "",
                "astar_runtime_ms": "",
                "java_distance": "",
                "java_runtime_ms": "",
                "match": "",
                "status": "ok",
            }

            if cfg.do_print:
                print(f"Parity run: n={n} seed={seed} pair=t_{cfg.pair_i}->t_{cfg.pair_j} algo={cfg.algo}")

            try:
                astar = run_astar(astar_binary, points, tri1, tri2, cfg.algo)
                java = run_java(java_out, points, tri1, tri2, cfg.java_time_limit_sec)
                row["astar_distance"] = astar.get("flip distance", "")
                row["astar_runtime_ms"] = astar.get("runtime", "")
                row["java_distance"] = java.get("distance", "")
                row["java_runtime_ms"] = java.get("time_ms", "")
                java_status = java.get("status", "ok")
                if java_status != "ok":
                    row["status"] = f"java_{java_status}"
                    row["match"] = ""
                    incomplete += 1
                else:
                    row["match"] = int(row["astar_distance"] == row["java_distance"])
                    if not row["match"]:
                        mismatches += 1
                        row["status"] = "distance_mismatch"
            except Exception as exc:
                failures += 1
                incomplete += 1
                row["status"] = f"error:{exc}"

            rows.append(row)

    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "n",
                "seed",
                "pair",
                "algo",
                "astar_distance",
                "astar_runtime_ms",
                "java_distance",
                "java_runtime_ms",
                "match",
                "status",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    if cfg.do_print:
        print(f"Wrote {len(rows)} rows to {output_path}")
        if mismatches:
            print(f"WARNING: {mismatches} distance mismatches detected")
        elif incomplete:
            print(f"Completed with {incomplete} incomplete runs (timeouts or errors)")
        else:
            print("All distances match between A* and Java file BFS")
        if failures:
            print(f"Completed with {failures} failed runs")

    return 0 if mismatches == 0 and incomplete == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
