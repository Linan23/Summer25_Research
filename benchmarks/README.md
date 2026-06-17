# Benchmark Summaries

This directory contains curated benchmark artifacts that document the current empirical state of the FlipDist solver. Raw generated sweeps, plots, profile dumps, and temporary parity outputs belong in ignored `results/`.

Keep files here compact and intentional. A benchmark CSV belongs here only when it is cited by `README.md`, `docs/benchmarks.md`, `docs/hard-limit-analysis.md`, or another research handoff note.

## Files

- `random_n23_25_seeds0_100_t2p5_m3.csv`: current retained random sweep for `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.
- `random_n23_25_profile_instrumented_validation_summary.csv`: fresh n=23..25 strict-2.5s validation after adding optional reuse-counter instrumentation.
- `random_n23_25_timeout_exact_switch_probe_summary.csv`: exact-safe switch probe over the current n=23..25 timeout seeds; no timeout seed was recovered under `2.5s`.
- `random_n23_25_remaining_timeouts_t10_summary.csv`: 10s retest summary for the current n=23..25 timeout seeds.
- `random_n23_27_goal_90_99_hard_limit_summary.csv`: compact status for the requested n=23..25 99% and n=26..27 90% targets.
- `random_hard_limit_n23_35_summary.csv`: curated hard-limit summary across the baseline and larger random samples.
- `random_n26_35_seeds0_100_t2_m3_summary.csv`: full hard-limit sweep summary for `n=26..35`, seeds `0..100`, strict `2s`, `max_k=3n`.
- `random_n26_35_seeds0_100_t2_m3_instances.csv`: pair-level status for the same full hard-limit sweep, so solved and timeout seeds can be inspected directly.
- `random_n26_27_seeds0_100_t2_budget_direction_patch_summary.csv`: latest n=26..27 full boundary refresh after exact budget-probe tuning and shape-based direction ordering.
- `random_n26_27_seeds0_100_t2_budget_direction_patch_instances.csv`: pair-level status for the same latest n=26..27 boundary refresh.
- `random_n26_27_remaining_timeouts_t10_after_budget_direction_summary.csv`: 10s evidence filtered to seeds that still time out after the latest strict-2s n=26..27 boundary refresh.
- `random_n26_27_profile_instrumented_validation_summary.csv`: fresh n=26..27 strict-2s validation after adding optional reuse-counter instrumentation.
- `random_n26_27_bottleneck_exact_cache_probe_summary.csv`: current exact-cache bottleneck probe; tested partition cache, split cache, combined cache, and incumbent pruning, with no promoted default because coverage did not improve and one full n=27 run regressed.
- `random_n26_27_profile_reuse_counters_summary.csv`: profile-only unique/repeated state counters for representative persistent n=26..27 hard seeds.
- `random_n26_27_seeds0_100_t2_direction_patch_summary.csv`: previous n=26..27 full boundary refresh after shape-based direction ordering.
- `random_n26_27_seeds0_100_t2_direction_patch_instances.csv`: pair-level status for the same previous n=26..27 boundary refresh.
- `random_n26_27_remaining_timeouts_t10_summary.csv`: 10s retest summary for the seeds that still timed out after the previous strict-2s n=26..27 boundary refresh.
- `random_n26_27_boundary_summary.csv`: focused n=26..27 hard-boundary optimization and timeout evidence.
- `random_n27_30_90_feasibility_summary.csv`: focused evidence for the n=27..30 90% feasibility check.
- `shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv`: retained shared-convex FlipDist vs AStarFlipDistance summary.
- `shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`: current retained local shared-convex comparison against a no-Gurobi AStar binary on identical inputs, with paired-row medians and win counts.

See `docs/benchmarks.md` for regeneration commands and metric definitions.
