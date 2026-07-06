# Exact Flip Distance Research Code

FlipDist is a research implementation of an exact solver for flip distance on rooted binary trees, equivalently rotation distance between binary trees or flip distance between triangulations of a convex polygon. The solver is built on the fixed-parameter algorithmic framework of Li and Xia, "An O(3.82^k) Time FPT Algorithm for Convex Flip Distance" (STACS 2023).

The repository is organized for research handoff and reproducibility. It contains buildable C++ solver sources, a Java triangulation oracle for independent checks, benchmark and parity tools, curated current results, and documentation for future extensions.

See `docs/references.md` for the paper citation and problem background. See `docs/architecture.md` for how the paper-level structure maps onto this C++/Python/Java codebase.

The organization follows the spirit of the `gchure/reproducible_research` template: each major area has a clear purpose, generated artifacts are separated from tracked summaries, and reproduction commands are documented from the repository root. The adaptation is described in `docs/reproducible-research.md`.

## Current Solver Status

Current exact random benchmark: `n=23..25`, seeds `0..200`, timeout `5s`,
`max_k=3n`, `bfs_cap=1`, both directions.

| n | Pair solves | Pair timeouts | Directed solves | Median solved-pair max wall time | p95 solved-pair max wall time |
| --- | ---: | ---: | ---: | ---: | ---: |
| `23` | `199/201 = 99.00%` | `2/201 = 1.00%` | `398/402 = 99.00%` | `14.224 ms` | `978.041 ms` |
| `24` | `192/201 = 95.52%` | `9/201 = 4.48%` | `384/402 = 95.52%` | `13.657 ms` | `811.542 ms` |
| `25` | `191/201 = 95.02%` | `10/201 = 4.98%` | `382/402 = 95.02%` | `77.078 ms` | `811.223 ms` |

This meets the current `>=95%` pair-solvability target for each of `n=23`,
`n=24`, and `n=25` under the 5s timeout. Final timeout seeds were
`n=23: 24, 177`, `n=24: 15, 33, 73, 97, 150, 154, 156, 161, 173`, and
`n=25: 14, 60, 73, 97, 150, 153, 154, 156, 161, 193`.

Correctness guardrails remained clean: Java parity had zero distance/status
mismatches on feasible sampled oracle ranges through `n=15`, the full benchmark
had no solved-pair distance mismatches between directions, and there were no
`not_found` or error rows. The generated final CSV is
`results/goal_opt_random_n23_25_s0_200_t5_m3_terminal.csv`.

## AStarFlipDistance Comparison

AStarFlipDistance is optional and is not vendored in this repository. It is used only as an external comparison point on shared-convex inputs. Older retained summaries show A* faster on one prior dataset/build; the current retained local no-Gurobi comparison on identical shared-convex inputs (`n=22..30`) shows FlipDist faster on paired solved-case median runtime for each n in that sample. Use `--astar-binary` with the benchmark tools to compare against a local external AStar build.

## Hard Limit Snapshot

The current practical limit evidence is summarized in `docs/hard-limit-analysis.md`. The
latest exact n=23..25 optimization pass reaches the current `>=95%` per-size
target under a 5s cap over seeds `0..200`, but it does not claim a 99% target
or extend that threshold to larger n. Older strict-cap evidence remains useful
for the boundary: the retained `n=26..27`, seeds `0..100`, timeout `2s`,
`max_k=3n` result stands at `n=26: 174/202 = 86.1%` and
`n=27: 176/202 = 87.1%`. On the wider retained `n=26..35` sweep, coverage drops
to `36.6%` by `n=35`.

| n range | Current takeaway |
| --- | --- |
| `n=23..25` | Current exact solver reaches at least `95%` pair solvability per n under `5s` over seeds `0..200`. |
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

Open the browser visualizer and generate a case from the page:

```bash
python3 tools/visualizer/serve.py
```

New readers may want `docs/terminology.md` first. It defines terms such as simple case, complex case, directed row, combined directed coverage, timeout, and hard-limit analysis.

## Directory Map

| Path | Purpose |
| --- | --- |
| `src/flipdist/` | C++ solver, brute-force validator, CLI, memoization, and helper code. |
| `scripts/` | Setup and workflow automation. |
| `tools/` | Maintained benchmark, parity, plotting, and comparison scripts. |
| `tools/visualizer/` | Local browser visualizer for normal `flipdist` JSON-lines output. |
| `tools/research_archive/` | Historical one-off analysis tools retained for reference only. |
| `tests/` | Smoke, parity, benchmark-slice, and script-help wrappers. |
| `oracle/java/` | Java triangulation oracle source and required `acm.jar`. |
| `benchmarks/` | Curated current benchmark summaries used in documentation. |
| `results/` | Ignored generated sweeps, plots, profiles, and parity outputs. |
| `third_party/` | Ignored optional local checkouts such as AStarFlipDistance. |
| `docs/` | Architecture, benchmark, development, artifact, and research notes. |
| `scripts/setup_dev.py` | All-in-one developer setup script. |

## Documentation

- `docs/reproducible-research.md`: how this repository adapts the research-template layout.
- `docs/architecture.md`: solver components and Li-Xia flow.
- `docs/references.md`: primary paper citation and problem background.
- `docs/terminology.md`: plain-language definitions for benchmark and solver terms.
- `docs/benchmarks.md`: maintained benchmark and parity commands.
- `docs/hard-limit-analysis.md`: current hard-limit, hard-case profile, and AStar comparison evidence.
- `docs/development.md`: build/test workflow and contribution expectations.
- `docs/data-artifacts.md`: retained artifacts, ignored outputs, and regeneration policy.
- `docs/repository-inventory.md`: file classes and tracked/ignored policy.
- `docs/research-notes.md`: short handoff notes on solver contract and current bottleneck.

The visualizer workflow is documented in `tools/visualizer/README.md`. It runs `flipdist` through a local server, then reconstructs small rotation paths in the browser for display only; it does not require `--emit-path` and does not change solver output.

## Requirements

- CMake
- C++20 compiler (`clang++`, `g++`, or `CXX`)
- Python 3
- Java and Javac for oracle checks
- Python package: `matplotlib>=3.8` for plotting tools

Optional AStar comparison setup is documented in `docs/benchmarks.md`; keep external checkouts under `third_party/` or outside the repo.
