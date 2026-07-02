# Reproducible Research Organization

This repository follows the spirit of the `gchure/reproducible_research`
template. The goal of doing so is that a future researcher can easily identify
the solver, reproduce the retained results, and understand which files are
source material versus generated output.

## Template Mapping

| Template idea | FlipDist adaptation | Notes |
| --- | --- | --- |
| `software_module/` | `src/flipdist/` | Maintained C++ solver module and binaries built by CMake. |
| `experiments/processing` | `tools/` and `scripts/` | Maintained Python scripts generate cases, run sweeps, compare oracles, and prepare summaries. |
| `experiments/analysis` | `tools/`, `benchmarks/`, and `docs/hard-limit-analysis.md` | Analysis scripts and retained summaries document the current empirical state. |
| `experiments/figures` | Plot scripts in `tools/`; generated plots in `results/` | Generated figures are local artifacts unless promoted deliberately. |
| `data/` | Not used | The project does not require external input data. Benchmark instances are regenerated from documented seeds. |
| `miscellaneous/protocols` | `docs/benchmarks.md`, `docs/development.md` | These files describe reproducible workflows and validation expectations. |
| `miscellaneous/software details` | `README.md`, `setup.sh`, `scripts/setup_dev.py`, `requirements.txt` | Setup and dependency information live at the repository root and in setup scripts. |
| `tests/` | `tests/` | Thin shell/Python wrappers for smoke, parity, benchmark-slice, and script-help checks. |
| `templates/` | Not currently used | There is no repeated manual experiment form. Add this only if future work needs standardized experiment notes. |

## Research Handoff Principles

1. Keep the root `README.md` as the entry point for purpose, quick start,
   current performance, and directory navigation.
2. Keep detailed explanation in `docs/`, with plain terminology for new
   readers and specific benchmark commands for reproduction.
3. Keep solver source in `src/flipdist/` and preserve the public CLI/output
   contract unless a research task explicitly changes it.
4. Keep maintained experiment and comparison scripts in `tools/`; keep setup
   automation in `scripts/`.
5. Keep retained benchmark summaries in `benchmarks/`; keep generated raw
   sweeps, profiles, plots, and scratch outputs in ignored `results/`.
6. Keep optional external dependencies, such as AStarFlipDistance, outside the
   tracked source surface or under ignored `third_party/`.

## Reproduce The Current State

From a fresh clone:

```bash
./setup.sh
```

For a manual build and validation pass:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
tests/smoke.sh
ctest --test-dir build --output-on-failure
tests/java_parity.sh
tests/benchmark_slice.sh
```

For benchmark regeneration and interpretation, use `docs/benchmarks.md` and
`docs/data-artifacts.md`. Promote only compact, documented summaries from
`results/` into `benchmarks/`.

