# Data And Artifacts

The repository tracks source, maintained tools, oracle source, and a small set of curated benchmark summaries. Large generated sweeps, profiles, temporary parity outputs, build outputs, and external dependency builds are intentionally excluded.

## Retained

`benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv` contains the current requested random sweep for `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

`benchmarks/random_hard_limit_n23_35_summary.csv` contains the compact hard-limit summary across the baseline and larger random samples.

`benchmarks/random_n26_35_seeds0_100_t2_m3_summary.csv` contains the latest full hard-limit sweep for `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`.

`benchmarks/random_n26_35_seeds0_100_t2_m3_instances.csv` contains the pair-level status for that same sweep. Use it to see exactly which seeds solved and which seeds timed out.

`benchmarks/random_n26_27_boundary_summary.csv` contains the focused n=26..27 hard-boundary optimization summary, including strict 2s and 2.5s coverage before and after the direction-order lock.

`benchmarks/random_n27_30_90_feasibility_summary.csv` contains the focused n=27..30 90% feasibility checks, including the strict 2s default, the current 2.5s reference, a 10s timeout probe, and Li-Xia-preserving cache/order experiments.

`benchmarks/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv` contains the retained shared-convex AStar comparison summary.

`benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv` contains the fresh local no-Gurobi AStar comparison summary.

## Ignored

`results/` is the default destination for regenerated sweeps, profiles, parity CSVs, plots, and scratch experiment outputs. It is ignored so reviews do not include large or stale generated files.

`build/`, `build-asan/`, `.venv/`, and `oracle/java/out/` are local setup outputs.

`third_party/` is reserved for optional local checkouts such as AStarFlipDistance. Do not vendor AStar source or build outputs into the tracked repository unless project policy changes.

## Regeneration

Use `docs/benchmarks.md` for the maintained commands. After regenerating a result that should become part of the handoff state, copy only the final compact CSV or summary into `benchmarks/` and update `README.md` if the headline numbers change.
