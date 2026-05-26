# Benchmark Summaries

This directory contains curated benchmark artifacts that document the current solver state. Raw generated sweeps and plots belong in ignored `results/`.

## Files

- `random_n23_25_seeds0_100_t2p5_m3.csv`: latest requested random sweep for `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.
- `random_hard_limit_n23_35_summary.csv`: curated hard-limit summary across the baseline and larger random samples.
- `random_n26_27_boundary_summary.csv`: focused n=26..27 hard-boundary optimization and timeout evidence.
- `random_n27_30_90_feasibility_summary.csv`: focused evidence for the n=27..30 90% feasibility check.
- `shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv`: retained shared-convex FlipDist vs AStarFlipDistance summary.
- `shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`: fresh local shared-convex comparison against a no-Gurobi AStar binary on identical inputs, with paired-row medians and win counts.

See `docs/benchmarks.md` for regeneration commands and metric definitions.
