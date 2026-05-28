# Hard-Limit Analysis

This page explains where the current FlipDist solver works well, where it starts to time out, and what the benchmark numbers mean. It is written for both researchers and users who want the current project state without reading the solver code first.

FlipDist is still an exact solver. A timeout does not mean the distance is wrong. It means the solver did not finish within the benchmark time cap.

## How To Read The Tables

- `n` is the tree size used by the random generator.
- Each seed creates one tree pair.
- Each pair is tested in two directions: `a->b` and `b->a`.
- `Directed solves` counts both directions, so seeds `0..100` give `202` directed rows.
- `Solved pairs` counts a seed only when both directions finish.
- `max_k=3n` is the search budget cap. The solver remains complete for the configured `max_k`.
- The full `n=26..35` sweep below uses a strict `2s` process cap. If the process exceeds that cap before both JSON rows return, both rows are recorded as timeout by the harness.

## Current Baseline

The maintained baseline remains `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`.

| n | Directed solves | Solved pairs | Median time on solved pairs |
|---:|---:|---:|---:|
| 23 | 198/202 = 98.0% | 99/101 = 98.0% | 37.2 ms |
| 24 | 186/202 = 92.1% | 93/101 = 92.1% | 22.9 ms |
| 25 | 192/202 = 95.0% | 96/101 = 95.0% | 24.3 ms |

Combined directed coverage is `576/606 = 95.0%`. This is the stable baseline range for the current solver.

## Full n=26..35 Sweep

Latest full hard-limit sweep: random `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`.

| n | Solved pairs | Directed solves | Median time on solved pairs | Timeout seeds |
|---:|---:|---:|---:|---|
| 26 | 85/101 = 84.2% | 170/202 = 84.2% | 31.5 ms | 14, 24, 25, 30, 32, 59, 62, 67, 69, 72, ... (16 total) |
| 27 | 80/101 = 79.2% | 160/202 = 79.2% | 43.4 ms | 5, 9, 11, 14, 32, 40, 44, 47, 48, 54, ... (21 total) |
| 28 | 64/101 = 63.4% | 128/202 = 63.4% | 45.7 ms | 1, 5, 7, 9, 11, 14, 17, 20, 25, 26, ... (37 total) |
| 29 | 57/101 = 56.4% | 114/202 = 56.4% | 33.4 ms | 1, 4, 7, 10, 11, 13, 14, 16, 17, 20, ... (44 total) |
| 30 | 54/101 = 53.5% | 108/202 = 53.5% | 32.4 ms | 1, 2, 7, 9, 11, 13, 14, 16, 17, 20, ... (47 total) |
| 31 | 49/101 = 48.5% | 98/202 = 48.5% | 72.3 ms | 1, 2, 4, 5, 8, 9, 11, 12, 14, 17, ... (52 total) |
| 32 | 45/101 = 44.5% | 90/202 = 44.5% | 83.8 ms | 1, 4, 9, 10, 11, 12, 14, 15, 16, 17, ... (56 total) |
| 33 | 47/101 = 46.5% | 94/202 = 46.5% | 128.4 ms | 0, 1, 4, 7, 8, 12, 13, 14, 18, 26, ... (54 total) |
| 34 | 44/101 = 43.6% | 88/202 = 43.6% | 119.9 ms | 0, 1, 4, 6, 9, 10, 11, 12, 13, 14, ... (57 total) |
| 35 | 37/101 = 36.6% | 74/202 = 36.6% | 91.8 ms | 3, 5, 6, 8, 10, 11, 12, 13, 15, 19, ... (64 total) |

The exact solved and timeout seed lists are retained in:

- `benchmarks/random_n26_35_seeds0_100_t2_m3_summary.csv`
- `benchmarks/random_n26_35_seeds0_100_t2_m3_instances.csv`

## Practical Limit

Under the current 2s cap, the practical limit starts at `n=26..27`.

The solver can still solve many larger instances quickly, but the success rate drops below the target range:

- `n=26`: still useful, but no longer reaches 95% on the wider `0..100` seed set.
- `n=27`: below 80% on the wider seed set.
- `n=28+`: not reliable under the current 2s cap.

This is a practical time limit, not a proof that exact solving is impossible at those sizes.

## Why The Solver Times Out

The main bottleneck is `TreeDistS` when `S` is empty.

In plain terms:

1. The solver reaches a state where many rotations are possible.
2. Many rotations lead to similar-looking subproblems.
3. The solver repeatedly splits those subproblems into partition sides.
4. Hard cases create hundreds of thousands of repeated recursive checks before the 2s cap expires.

The current implementation already uses safe caching, child-state deduplication, lower-bound reuse, and direction ordering. These reduce some repeated work, but they do not remove the core growth in the empty-`S` branch.

## Recent Optimization Evidence

The n=26..27 boundary pass improved the smaller seeds `0..20` slice by changing safe direction ordering:

| Scenario | Timeout | n26 | n27 | Combined |
|---|---:|---:|---:|---:|
| Before direction lock | 2s | 38/42 = 90.5% | 34/42 = 81.0% | 72/84 = 85.7% |
| After direction lock | 2s | 40/42 = 95.2% | 36/42 = 85.7% | 76/84 = 90.5% |
| Empty-S pair-bound propagation | 2s | 40/42 = 95.2% | 36/42 = 85.7% | 76/84 = 90.5% |

That improvement came from solving timing-margin cases earlier. It did not fix the persistent hard seeds.

Persistent hard seeds in the smaller slice:

- `n=26`: seed `14`.
- `n=27`: seeds `5`, `9`, and `14`.

Profile samples still show the same bottleneck. For example, an n=27 hard run can reach more than `370,000` `TreeDistS` calls, more than `140,000` empty-`S` calls, and more than `170,000` partition calls before aborting.

## What Would Be Needed Next

Further large gains likely require a structural change, not only more local caching.

Promising directions include:

- stronger lower bounds that reject impossible branches earlier;
- a shared dynamic program for repeated partition-side checks;
- stronger state equivalence so equivalent rotation paths are not explored repeatedly;
- a different empty-`S` search strategy while preserving exactness.

## AStar Comparison

AStarFlipDistance is optional and benchmarked separately on shared-convex inputs.

The current local no-Gurobi comparison on paired solved rows shows FlipDist faster than A* on median runtime for `n=22..30` in that shared-convex sample. Raw medians can be misleading when one solver times out on different rows, so paired-row medians are the preferred comparison.

The retained comparison summary is:

- `benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`

## Validation

Latest validation for this solver state:

- C++ build completed.
- `./build/bf_bst` passed.
- CLI smoke test passed.
- Java parity passed for random `n=12..13`, seeds `0..5`.
- Full random sweep completed for `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`.
