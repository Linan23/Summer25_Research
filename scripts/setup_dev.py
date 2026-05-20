#!/usr/bin/env python3
"""Set up a local developer environment for FlipDist."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Configure, build, and smoke-test the FlipDist workspace.")
    parser.add_argument("--skip-java", action="store_true", help="Skip Java toolchain checks and oracle compilation.")
    parser.add_argument("--skip-python-env", action="store_true", help="Do not create .venv or install requirements.")
    parser.add_argument("--skip-tests", action="store_true", help="Build only; do not run smoke checks.")
    parser.add_argument("--clean-build", action="store_true", help="Remove build/ before configuring CMake.")
    parser.add_argument(
        "--with-astar",
        metavar="PATH",
        default=None,
        help="Optional path to an AStarFlipDistance checkout or A_star_for_flipdistance binary.",
    )
    return parser.parse_args()


def run(cmd: list[str], *, cwd: Path = ROOT) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def require_tool(name: str, hint: str | None = None) -> str:
    path = shutil.which(name)
    if path:
        print(f"found {name}: {path}")
        return path
    suffix = f" ({hint})" if hint else ""
    raise SystemExit(f"missing required tool: {name}{suffix}")


def require_cpp_compiler() -> None:
    cxx = os.environ.get("CXX")
    candidates = [cxx] if cxx else []
    candidates.extend(["c++", "clang++", "g++"])
    for candidate in candidates:
        if candidate and shutil.which(candidate):
            print(f"found C++ compiler: {candidate}")
            return
    raise SystemExit("missing required C++20 compiler; install clang++ or g++, or set CXX.")


def venv_python() -> Path:
    if sys.platform == "win32":
        return ROOT / ".venv" / "Scripts" / "python"
    return ROOT / ".venv" / "bin" / "python"


def setup_python_env(skip: bool) -> None:
    if skip:
        print("skip python environment")
        return
    venv = ROOT / ".venv"
    if not venv.exists():
        run([sys.executable, "-m", "venv", str(venv)])
    run([str(venv_python()), "-m", "pip", "install", "--upgrade", "pip"])
    run([str(venv_python()), "-m", "pip", "install", "-r", "requirements.txt"])


def build_cpp(clean_build: bool) -> None:
    build_dir = ROOT / "build"
    if clean_build and build_dir.exists():
        shutil.rmtree(build_dir)
    run(["cmake", "-S", ".", "-B", "build", "-DCMAKE_BUILD_TYPE=Release"])
    run(["cmake", "--build", "build", "-j"])


def compile_java(skip: bool) -> None:
    if skip:
        print("skip Java oracle")
        return
    source_dir = ROOT / "oracle" / "java" / "src"
    out_dir = ROOT / "oracle" / "java" / "out"
    jar = ROOT / "oracle" / "java" / "lib" / "acm.jar"
    sources = sorted(str(path) for path in source_dir.glob("*.java"))
    if not sources:
        raise SystemExit(f"no Java sources found in {source_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    run(["javac", "-cp", str(jar), "-d", str(out_dir), *sources])


def check_astar(path_arg: str | None) -> None:
    if not path_arg:
        print("skip AStarFlipDistance integration; pass --with-astar PATH to verify a local checkout or binary")
        return
    path = Path(path_arg).expanduser()
    candidates = [path]
    if path.is_dir():
        candidates.append(path / "build-nogurobi" / "A_star_for_flipdistance")
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            print(f"found optional AStar binary: {candidate}")
            return
    raise SystemExit(
        "could not find an executable AStar binary. Pass the binary itself or a checkout "
        "with build-nogurobi/A_star_for_flipdistance."
    )


def run_smoke_tests(skip_tests: bool, skip_java: bool) -> None:
    if skip_tests:
        print("skip smoke tests")
        return
    run(["./build/bf_bst"])
    run(["./build/flipdist", "--case", "random", "--n", "12", "--seed", "0", "--count", "1", "--max-k", "30"])
    if not skip_java:
        run(
            [
                sys.executable,
                "tools/run_flipdist_java_parity.py",
                "--case",
                "random",
                "--n",
                "12",
                "--seed",
                "0",
                "--count",
                "1",
                "--cpp-binary",
                "./build/flipdist",
                "--max-k",
                "36",
                "--java-out",
                "oracle/java/out",
                "--java-lib",
                "oracle/java/lib/acm.jar",
            ]
        )


def main() -> int:
    cfg = parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    os.chdir(ROOT)

    print(f"workspace: {ROOT}")
    require_tool("python3")
    require_tool("cmake")
    require_cpp_compiler()
    if not cfg.skip_java:
        require_tool("java")
        require_tool("javac")

    setup_python_env(cfg.skip_python_env)
    build_cpp(cfg.clean_build)
    compile_java(cfg.skip_java)
    check_astar(cfg.with_astar)
    run_smoke_tests(cfg.skip_tests, cfg.skip_java)
    print("setup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
