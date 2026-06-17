# Data And Artifacts

The repository tracks the source code, maintained research tools, oracle source, policy notes, and a compact set of curated benchmark summaries. Large generated sweeps, profiles, temporary parity outputs, build outputs, and external dependency builds are intentionally excluded so that the tracked history remains reviewable.

## Directory Policy

The project does not use a tracked external data directory. The retained random experiments are generated from documented seeds by the solver and benchmark harnesses.

`results/` is ignored except for `results/README.md`. Use it for regenerated sweeps, parity CSVs, plots, profiles, shared-convex generated inputs, and temporary experiment outputs.

`benchmarks/` is tracked and should contain only compact curated CSV summaries that document the empirical handoff state.

`third_party/` is ignored except for `third_party/README.md`. Optional dependencies such as AStarFlipDistance should stay as local external checkouts or live outside the repository.

## Retained

`benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv` contains the current requested random sweep for `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

`benchmarks/random_n23_25_profile_instrumented_validation_summary.csv` contains the fresh n=23..25 strict-2.5s validation run after adding optional reuse-counter instrumentation.

`benchmarks/random_n23_25_timeout_exact_switch_probe_summary.csv` contains the exact-safe switch probe over the current n=23..25 timeout seeds. It records that no timeout seed was recovered under the `2.5s` cap.

`benchmarks/random_n23_25_remaining_timeouts_t10_summary.csv` contains the 10s retest summary for the current n=23..25 timeout seeds.

`benchmarks/random_n23_27_goal_90_99_hard_limit_summary.csv` contains the compact status for the requested n=23..25 99% and n=26..27 90% targets.

`benchmarks/random_hard_limit_n23_35_summary.csv` contains the compact hard-limit summary across the baseline and larger random samples.

`benchmarks/random_n26_35_seeds0_100_t2_m3_summary.csv` contains the latest full hard-limit sweep for `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`.

`benchmarks/random_n26_35_seeds0_100_t2_m3_instances.csv` contains the pair-level status for that same sweep. Use it to see exactly which seeds solved and which seeds timed out.

`benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_summary.csv` contains the latest full n=26..27 boundary refresh after exact budget-probe tuning and shape-based direction ordering.

`benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_instances.csv` contains the pair-level status for that same latest n=26..27 refresh.

`benchmarks/random_n26_27_remaining_timeouts_t10_after_budget_direction_summary.csv` contains the compact 10s evidence filtered to seeds that still time out after the latest strict-2s n=26..27 refresh.

`benchmarks/random_n26_27_profile_instrumented_validation_summary.csv` contains the fresh n=26..27 strict-2s validation run after adding optional reuse-counter instrumentation.

`benchmarks/random_n26_27_bottleneck_exact_cache_probe_summary.csv` contains the exact-cache bottleneck probe for the current n=26..27 boundary. It records the partition cache, split cache, combined cache, and incumbent-pruning checks that were not promoted.

`benchmarks/random_n26_27_profile_reuse_counters_summary.csv` contains profile-only unique/repeated state counters for representative persistent n=26..27 hard seeds.

`benchmarks/random_n26_27_seeds0_100_t2_direction_patch_summary.csv` contains the previous full n=26..27 boundary refresh after shape-based direction ordering.

`benchmarks/random_n26_27_seeds0_100_t2_direction_patch_instances.csv` contains the pair-level status for that same previous n=26..27 refresh.

`benchmarks/random_n26_27_remaining_timeouts_t10_summary.csv` contains the compact 10s retest summary for the seeds that still timed out after the previous strict-2s n=26..27 refresh.

`benchmarks/random_n26_27_boundary_summary.csv` contains the focused n=26..27 hard-boundary optimization summary, including strict 2s and 2.5s coverage before and after the direction-order lock.

`benchmarks/random_n27_30_90_feasibility_summary.csv` contains the focused n=27..30 90% feasibility checks, including the strict 2s default, the current 2.5s reference, a 10s timeout probe, and Li-Xia-preserving cache/order experiments.

`benchmarks/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv` contains the retained shared-convex AStar comparison summary.

`benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv` contains the current retained local no-Gurobi AStar comparison summary.

## Ignored

`results/` is the default destination for regenerated sweeps, profiles, parity CSVs, plots, and scratch experiment outputs. It is ignored so reviews do not include large or stale generated files.

`build/`, `build-asan/`, `.venv/`, and `oracle/java/out/` are local setup outputs.

`third_party/` is reserved for optional local checkouts such as AStarFlipDistance. Do not vendor AStar source or build outputs into the tracked repository unless project policy changes.

## Regeneration

Use `docs/benchmarks.md` for maintained commands. After regenerating a result that should become part of the handoff state, copy only the final compact CSV or summary into `benchmarks/` and update `README.md` if the headline numbers change.
