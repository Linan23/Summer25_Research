# Exact Flip Distance Research Code

FlipDist is a research implementation of an exact solver for flip distance on rooted binary trees, equivalently rotation distance between binary trees or flip distance between triangulations of a convex polygon. The solver is built on the fixed-parameter algorithmic framework of Li and Xia, "An O(3.82^k) Time FPT Algorithm for Convex Flip Distance" (STACS 2023).

The repository is organized for research handoff and reproducibility. It contains buildable C++ solver sources, a Java triangulation oracle for independent checks, benchmark and parity tools, curated current results, and documentation for future extensions.

See `docs/references.md` for the paper citation and problem background. See `docs/architecture.md` for how the paper-level structure maps onto this C++/Python/Java codebase.

## Current Solver Status

Current retained random benchmark: `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

| Metric | Result |
| --- | ---: |
| First-direction exact solves | `288/303 = 95.0%` |
| First-direction solves under 2s | `287/303 = 94.7%` |
| Directed exact solves | `576/606 = 95.0%` |
| Directed timeouts | `30` |

Per-size directed coverage: `n=23: 198/202 = 98.0%`, `n=24: 186/202 = 92.1%`, `n=25: 192/202 = 95.0%`.

Single-run timing has small variance near the 2s boundary, so the under-2s count can move slightly across machines and runs. The retained CSV is in `benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv`.

## AStarFlipDistance Comparison

AStarFlipDistance is optional and is not vendored in this repository. It is used only as an external comparison point on shared-convex inputs. Older retained summaries show A* faster on one prior dataset/build; the current retained local no-Gurobi comparison on identical shared-convex inputs (`n=22..30`) shows FlipDist faster on paired solved-case median runtime for each n in that sample. Use `--astar-binary` with the benchmark tools to compare against a local external AStar build.

## Hard Limit Snapshot

The current practical limit evidence is summarized in `docs/hard-limit-analysis.md`. The n=23..25 baseline remains below the requested 99% per-size target: `n=23: 198/202`, `n=24: 186/202`, `n=25: 192/202`. Exact-safe probes over the current n=23..25 timeout seeds recovered `0/15` under the 2.5s cap, and only `3/15` solved by 10s. After the latest exact budget-probe and direction-order pass, the retained `n=26..27`, seeds `0..100`, timeout `2s`, `max_k=3n` result stands at `n=26: 174/202 = 86.1%` and `n=27: 176/202 = 87.1%`. On the wider retained `n=26..35` sweep, coverage drops to `36.6%` by `n=35`.

| n range | Current takeaway |
| --- | --- |
| `n=23..25` | Stable baseline range; combined directed coverage is `95.0%` at `2.5s`, but the 99% per-size target was not reached. |
| `n=26..27` | Practical boundary; latest full `0..100` coverage is improved but still below 90% under `2s`. |
| `n=28+` | Current Li-Xia-structured solver is not reliable under the strict `2s` cap. |

The bottleneck is still `TreeDistS/S.empty()`: complex cases repeatedly explore rotation children and partition-side checks. Local exact budget-probe tuning and ordering improvements help timing-margin cases, but the retained sweeps suggest that 90%-plus coverage through n=27 under `2s` would require a deeper structural improvement.

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
tests/smoke.sh
ctest --test-dir build --output-on-failure
```

Small parity check against the Java oracle:

```bash
tests/java_parity.sh
```

New readers may want `docs/terminology.md` first. It defines terms such as simple case, complex case, directed row, combined directed coverage, timeout, and hard-limit analysis.

## Directory Map

| Path | Purpose |
| --- | --- |
| `src/flipdist/` | C++ solver, brute-force validator, CLI, memoization, and helper code. |
| `tools/` | Maintained benchmark, parity, plotting, and comparison scripts. |
| `tools/research_archive/` | Historical one-off analysis tools retained for reference only. |
| `tests/` | Smoke, parity, benchmark-slice, and script-help wrappers. |
| `oracle/java/` | Java triangulation oracle source and required `acm.jar`. |
| `benchmarks/` | Curated current benchmark summaries used in documentation. |
| `results/` | Ignored generated sweeps, plots, profiles, and parity outputs. |
| `third_party/` | Ignored optional local checkouts such as AStarFlipDistance. |
| `docs/` | Architecture, benchmark, development, artifact, and research notes. |
| `scripts/setup_dev.py` | All-in-one developer setup script. |

## Documentation

- `docs/architecture.md`: solver components and Li-Xia flow.
- `docs/references.md`: primary paper citation and problem background.
- `docs/terminology.md`: plain-language definitions for benchmark and solver terms.
- `docs/benchmarks.md`: maintained benchmark and parity commands.
- `docs/hard-limit-analysis.md`: current hard-limit, hard-case profile, and AStar comparison evidence.
- `docs/development.md`: build/test workflow and contribution expectations.
- `docs/data-artifacts.md`: retained artifacts, ignored outputs, and regeneration policy.
- `docs/repository-inventory.md`: file classes and tracked/ignored policy.
- `docs/research-notes.md`: short handoff notes on solver contract and current bottleneck.

## Requirements

- CMake
- C++20 compiler (`clang++`, `g++`, or `CXX`)
- Python 3
- Java and Javac for oracle checks
- Python package: `matplotlib>=3.8` for plotting tools

Optional AStar comparison setup is documented in `docs/benchmarks.md`; keep external checkouts under `third_party/` or outside the repo.
