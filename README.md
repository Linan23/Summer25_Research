# FlipDist Research Solver

FlipDist is an exact solver for flip distance on rooted binary trees, using the triangulation-equivalent view for validation and comparison. The C++ implementation preserves the Li-Xia search structure and is optimized around the current hard path in `TreeDistS`, especially partition-driven branching in empty-side subproblems.

The repository is organized as a developer-and-research handoff package: buildable C++ sources, a Java exact oracle, maintained benchmark tools, curated current results, and focused documentation.

## Current Solver Status

Latest requested random benchmark: `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

| Metric | Result |
| --- | ---: |
| First-direction exact solves | `274/303 = 90.4%` |
| First-direction solves under 2s | `271/303 = 89.4%` |
| Directed exact solves | `548/606` |
| Directed timeouts | `58` |

Single-run timing has small variance near the 2s boundary, so the under-2s count can move slightly across machines and runs. The retained CSV is in `benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv`.

## AStarFlipDistance Comparison

AStarFlipDistance is optional and is not vendored in this repository. Existing shared-convex benchmark summaries show A* has faster median runtime on the shared-convex cases it solves, while FlipDist has slightly higher completion coverage in the retained summary. Use `--astar-binary` with the benchmark tools to compare against a local external AStar build.

## Quick Start

From a fresh clone on macOS or Linux with no existing build output:

```bash
./setup.sh
```

The setup script checks Python 3, CMake, a C++20 compiler, Java/Javac, creates `.venv/`, installs `requirements.txt`, builds the C++ binaries, compiles the Java oracle, and runs smoke checks.

Manual build and smoke test:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/bf_bst
./build/flipdist --case random --n 12 --seed 0 --count 1 --max-k 30
```

Small parity check against the Java oracle:

```bash
python3 tools/run_flipdist_java_parity.py \
  --case random --n 12 --seed 0 --count 1 \
  --cpp-binary ./build/flipdist --max-k 36 \
  --java-out oracle/java/out --java-lib oracle/java/lib/acm.jar
```

## Directory Map

| Path | Purpose |
| --- | --- |
| `src/flipdist/` | C++ solver, brute-force validator, CLI, memoization, and helper code. |
| `tools/` | Maintained benchmark, parity, plotting, and comparison scripts. |
| `tools/research_archive/` | Historical one-off analysis tools retained for reference only. |
| `oracle/java/` | Java triangulation oracle source and required `acm.jar`. |
| `benchmarks/` | Curated current benchmark summaries used in documentation. |
| `docs/` | Architecture, benchmark, development, and artifact policy documentation. |
| `scripts/setup_dev.py` | All-in-one developer setup script. |

## Documentation

- `docs/architecture.md`: solver components and Li-Xia flow.
- `docs/benchmarks.md`: maintained benchmark and parity commands.
- `docs/development.md`: build/test workflow and contribution expectations.
- `docs/data-artifacts.md`: retained artifacts, ignored outputs, and regeneration policy.

## Requirements

- CMake
- C++20 compiler (`clang++`, `g++`, or `CXX`)
- Python 3
- Java and Javac for oracle checks
- Python package: `matplotlib>=3.8` for plotting tools

Optional AStar comparison setup is documented in `docs/benchmarks.md`; keep external checkouts under `third_party/` or outside the repo.
