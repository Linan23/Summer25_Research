# FlipDist Research Solver

FlipDist is an exact solver for flip distance on rooted binary trees. The same problem can also be viewed through triangulations, and this repository keeps a Java triangulation oracle for independent checks. The C++ solver preserves the Li-Xia search structure and is optimized around the current hard path in `TreeDistS`, especially repeated partition checks in empty-side subproblems.

The repository is organized as a developer-and-research handoff package: buildable C++ sources, a Java exact oracle, maintained benchmark tools, curated current results, and focused documentation.

## Current Solver Status

Latest requested random benchmark: `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

| Metric | Result |
| --- | ---: |
| First-direction exact solves | `288/303 = 95.0%` |
| First-direction solves under 2s | `287/303 = 94.7%` |
| Directed exact solves | `576/606 = 95.0%` |
| Directed timeouts | `30` |

Per-size directed coverage: `n=23: 198/202 = 98.0%`, `n=24: 186/202 = 92.1%`, `n=25: 192/202 = 95.0%`.

Single-run timing has small variance near the 2s boundary, so the under-2s count can move slightly across machines and runs. The retained CSV is in `benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv`.

## AStarFlipDistance Comparison

AStarFlipDistance is optional and is not vendored in this repository. Older retained shared-convex summaries show A* faster on one prior dataset/build. A newer local no-Gurobi AStar comparison on identical shared-convex inputs (`n=22..30`) shows FlipDist faster on paired solved-case median runtime for each n in that sample. Use `--astar-binary` with the benchmark tools to compare against a local external AStar build.

## Hard Limit Snapshot

The current practical limit evidence is summarized in `docs/hard-limit-analysis.md`. On the wider hard-limit sweep, `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`, coverage drops from `84.2%` at `n=26` to `36.6%` at `n=35`.

| n range | Current takeaway |
| --- | --- |
| `n=23..25` | Stable baseline range; combined directed coverage is `95.0%` at `2.5s`. |
| `n=26..27` | Practical boundary; many instances solve quickly, but full `0..100` coverage is below 90% under `2s`. |
| `n=28+` | Current Li-Xia-structured solver is not reliable under the strict `2s` cap. |

The bottleneck is still `TreeDistS/S.empty()`: hard cases repeatedly explore rotation children and partition-side checks. Local caching and ordering improvements help timing-margin cases, but the full seed sweep shows that higher coverage likely requires a structural change.

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
- `docs/hard-limit-analysis.md`: current hard-limit, hard-case profile, and AStar comparison evidence.
- `docs/development.md`: build/test workflow and contribution expectations.
- `docs/data-artifacts.md`: retained artifacts, ignored outputs, and regeneration policy.

## Requirements

- CMake
- C++20 compiler (`clang++`, `g++`, or `CXX`)
- Python 3
- Java and Javac for oracle checks
- Python package: `matplotlib>=3.8` for plotting tools

Optional AStar comparison setup is documented in `docs/benchmarks.md`; keep external checkouts under `third_party/` or outside the repo.
