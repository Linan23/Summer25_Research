# Hard-Limit Analysis

This page summarizes the current empirical limits of the FlipDist implementation. It is written for researchers and collaborators who need to understand what the solver can reproduce, where the retained benchmarks start to fail under fixed time caps, and which bottlenecks remain open. Definitions for benchmark terms are in `docs/terminology.md`.

The implementation is based on the Li-Xia exact FPT framework for Convex Flip Distance; see `docs/references.md`. FlipDist remains an exact solver for the configured `max_k`. A timeout does not imply an incorrect distance; it means the process did not finish within the benchmark cap.

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

A fresh validation run on this build recorded `n=23: 198/202`, `n=24: 184/202`, and `n=25: 192/202`. The only row-level difference from the retained baseline was n=24 seed `24`, which solved in `2395.958 ms` in the retained run and now lands around `2.57s` when forced through the lower-cost direction. This is another strict-time timing-margin case. The compact validation summary is retained in `benchmarks/random_n23_25_profile_instrumented_validation_summary.csv`.

## 99% Baseline Attempt

The requested n=23..25 target was at least `200/202 = 99.0%` directed solves for each n under the same `2.5s` cap. The current solver does not reach that target:

| n | Current directed solves | Target | Gap |
|---:|---:|---:|---:|
| 23 | 198/202 = 98.0% | 200/202 | 2 directed solves |
| 24 | 186/202 = 92.1% | 200/202 | 14 directed solves |
| 25 | 192/202 = 95.0% | 200/202 | 8 directed solves |

Exact-safe probes over the current n=23..25 timeout seeds tested forced direction, exact budget-probe variants, partition side ordering, empty-`S` candidate ordering, partition/split caches, incumbent pruning, tight empty-`S` cache reuse, free-hint gating, articulation-arm reduction, common-edge decomposition, and A*-style pair ordering. None recovered any of the `15` timeout seed pairs under the `2.5s` cap.

Retesting those same timeout seeds at `10s` recovered only `3/15`:

| n | 2.5s timeout seeds retested | Solved by 10s | Still timeout at 10s |
|---:|---:|---:|---:|
| 23 | 2 | 0 | 2 |
| 24 | 8 | 2 | 6 |
| 25 | 5 | 1 | 4 |

The retained summaries are:

- `benchmarks/random_n23_25_timeout_exact_switch_probe_summary.csv`
- `benchmarks/random_n23_25_remaining_timeouts_t10_summary.csv`
- `benchmarks/random_n23_27_goal_90_99_hard_limit_summary.csv`

This makes the 99% n=23..25 target a preserved-structure empirical limit for the current implementation, not just a missing direction rule.

## Current Boundary Refresh

Latest focused boundary refresh: random `n=26..27`, seeds `0..100`, timeout `2s`, `max_k=3n`, after exact budget-probe tuning plus shape-based direction ordering. The Li-Xia search itself is unchanged; the solver now spends less time on expensive too-small budget probes at n=26..27 and starts with the direction that is more likely to finish first for a few hard boundary shapes.

| n | Solved pairs | Directed solves | Median time on solved pairs | Recovered seeds |
|---:|---:|---:|---:|---|
| 26 | 87/101 = 86.1% | 174/202 = 86.1% | 27.7 ms | 25, 92 |
| 27 | 88/101 = 87.1% | 176/202 = 87.1% | 52.5 ms | 11, 32, 40, 44, 47, 63, 83, 96 |

The exact seed-level refresh is retained in:

- `benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_summary.csv`
- `benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_instances.csv`

A fresh validation run after adding optional reuse counters recorded `n=26: 174/202` and `n=27: 174/202`. The row-level difference was seed `40` at n=27: the retained run solved `a->b` in `1922.731 ms`, while repeated direct runs on the current build place the same direction around `2.08s`. This is a timing-margin case at the strict 2s cap, not a distance or correctness change. The compact validation summary is retained in `benchmarks/random_n26_27_profile_instrumented_validation_summary.csv`.

This improves the boundary but does not reach the 90% target. The remaining misses are still dominated by the same empty-`S` partition recursion.

## Remaining Misses At 10s

The remaining strict-2s timeout seeds from the latest n=26..27 refresh were checked against the retained `10s` retest evidence. This answers a practical question: are the misses only slightly above 2s, or are many of them much harder?

| n | 2s timeout seeds retested | Solved by 10s | Still timeout at 10s | Solved-by-10s time range |
|---:|---:|---:|---:|---:|
| 26 | 14 | 5 | 9 | 3.84s-8.25s |
| 27 | 13 | 6 | 7 | 4.05s-9.61s |

The exact seed lists are retained in `benchmarks/random_n26_27_remaining_timeouts_t10_after_budget_direction_summary.csv`.

This means the 90% under-2s target is not just blocked by a few bad first directions. Even with five times the timeout, `16` of the `27` remaining boundary seed pairs still do not finish.

## Full n=26..35 Sweep

Retained full hard-limit sweep before the latest n=26..27 boundary refresh: random `n=26..35`, seeds `0..100`, timeout `2s`, `max_k=3n`.

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

Under the current 2s cap, the practical limit still starts at `n=26..27`.

The solver can still solve many larger instances quickly, but the success rate drops below the target range:

- `n=26`: improved to 86.1% on the wider `0..100` seed set, still below the 90%-99% target range.
- `n=27`: retained best run is 87.1% on the wider `0..100` seed set, with fresh validation at 86.1% because one timing-margin seed pair flipped over the strict 2s cap.
- `n=28+`: not reliable under the current 2s cap.

This is a practical time limit, not a proof that exact solving is impossible at those sizes.
For the current Li-Xia-structured solver, the practical 2s hard boundary is n=26..27. Reaching 90% would require four more n=26 seed pairs and three more n=27 seed pairs beyond the latest result; the retained 10s evidence and exact-switch probes show that many of those seeds are still far outside the current search profile.

## Why The Solver Times Out

The main bottleneck is `TreeDistS` when `S` is empty.

Operationally, the bottleneck appears as follows:

1. The solver reaches a state where many rotations are possible.
2. Many rotations lead to similar-looking subproblems.
3. The solver repeatedly splits those subproblems into partition sides.
4. Complex cases create hundreds of thousands of repeated recursive checks before the 2s cap expires.

The current implementation already uses safe caching, child-state deduplication, lower-bound reuse, and direction ordering. These reduce some repeated work, but they do not remove the core growth in the empty-`S` branch.

## Recent Optimization Evidence

The first n=26..27 boundary pass improved the smaller seeds `0..20` slice by changing safe direction ordering:

| Scenario | Timeout | n26 | n27 | Combined |
|---|---:|---:|---:|---:|
| Before direction lock | 2s | 38/42 = 90.5% | 34/42 = 81.0% | 72/84 = 85.7% |
| After direction lock | 2s | 40/42 = 95.2% | 36/42 = 85.7% | 76/84 = 90.5% |
| Empty-S pair-bound propagation | 2s | 40/42 = 95.2% | 36/42 = 85.7% | 76/84 = 90.5% |

That improvement came from solving timing-margin cases earlier. It did not fix the persistent complex seeds.

The latest full `0..100` boundary pass added a narrower shape-based direction rule. It recovered seven previously timed-out seed pairs without changing exact distances or the Li-Xia search:

- `n=26`: seeds `25` and `92`.
- `n=27`: seeds `11`, `32`, `44`, `63`, and `96`.

The follow-up 90% attempt added exact n=26..27 budget-probe tuning and one additional n=27 direction shape. It recovered three more n=27 seed pairs: `40`, `47`, and `83`. The new full result is:

- `n=26`: `174/202 = 86.1%`, unchanged from the previous boundary pass.
- `n=27`: `176/202 = 87.1%`, up from `170/202 = 84.2%`.

This still misses the 90% target of `182/202` directed solves for both sizes.

The combined n=23..27 target status is retained in `benchmarks/random_n23_27_goal_90_99_hard_limit_summary.csv`.

Profile evidence for `n=27 seed=83` shows what the budget-probe tuning fixed. Simulating the previous probe behavior, `a->b` aborted at 2.2s after about `383,000` `TreeDistS` calls, `186,000` empty-`S` calls, `130,000` partition calls, and `2.49s` of accumulated partition budget-loop time. With the latest default, the same direction solves in about `1.58s` with about `277,000` `TreeDistS` calls, `135,000` empty-`S` calls, `94,000` partition calls, and `1.75s` of accumulated partition budget-loop time. The reverse direction then reuses the exact result immediately.

For example, before the patch, `n=27 seed=32` entered the hard `a->b` direction first and aborted after more than `346,000` `TreeDistS` calls, `160,000` empty-`S` calls, and `128,000` partition calls. After the patch, it starts with `b->a`, solves in about `336 ms`, and the reverse direction reuses the exact result immediately.

Persistent complex seeds in the smaller slice:

- `n=26`: seed `14`.
- `n=27`: seeds `5`, `9`, and `14`.

Profile samples still show the same bottleneck. After the direction patch, persistent complex seeds still abort at 5s with very high recursive pressure. Profile times are accumulated across recursive calls, so a counter can be larger than wall-clock time.

| Profile | Direction | TreeDistS calls | Empty-S calls | Partition calls | Accumulated budget-loop time |
|---|---|---:|---:|---:|---:|
| `n=26 seed=14` | `a->b` | 835,399 | 392,060 | 286,518 | 6.66s |
| `n=26 seed=14` | `b->a` | 752,865 | 383,691 | 233,967 | 5.28s |
| `n=27 seed=5` | `a->b` | 729,401 | 252,748 | 340,816 | 2.16s |
| `n=27 seed=5` | `b->a` | 663,600 | 243,520 | 295,404 | 6.58s |
| `n=27 seed=9` | `a->b` | 699,722 | 250,012 | 317,490 | 2.69s |
| `n=27 seed=9` | `b->a` | 740,565 | 383,907 | 240,678 | 5.86s |

The direction-order patch helps when one direction is dramatically easier. It does not remove the repeated empty-`S` partition branching for persistent hard seeds.

The latest opt-in profile counters add a more detailed view. They separate states the solver has seen before from states that are new in the current search. In 2.2s profile slices, the persistent hard seeds are not just revisiting a small set of identical partition structures. They create a large number of new `TreeDistS`, empty-`S`, partition-side, and split-signature states:

| Profile | Direction | TreeDistS unique/repeated | Partition structures unique/repeated | Budget-loop time |
|---|---|---:|---:|---:|
| `n=26 seed=14` | `a->b` | 194,609 / 95,973 | 94,702 / 17,887 | 2.61s |
| `n=26 seed=14` | `b->a` | 194,825 / 84,839 | 81,619 / 11,131 | 1.83s |
| `n=27 seed=5` | `a->b` | 177,968 / 74,843 | 98,276 / 13,745 | 0.94s |
| `n=27 seed=5` | `b->a` | 175,353 / 63,981 | 95,534 / 8,234 | 3.31s |
| `n=27 seed=9` | `a->b` | 185,670 / 69,353 | 92,879 / 12,676 | 1.70s |
| `n=27 seed=9` | `b->a` | 205,760 / 67,777 | 74,899 / 7,690 | 2.31s |

This explains why simply enabling more exact partition caches did not push coverage to 90%. The cache can avoid some repeated work, but many misses are dominated by new branch growth. The compact counter summary is retained in `benchmarks/random_n26_27_profile_reuse_counters_summary.csv`.

Additional exactness-preserving experiments were checked and rejected as defaults for this goal:

- For n=23..25 timeout seeds, forced direction, exact budget-probe variants, partition side ordering, empty-`S` ordering, partition/split caches, incumbent pruning, tight empty-`S` cache reuse, free-hint gating, articulation-arm reduction, common-edge decomposition, and A*-style pair ordering recovered `0/15` under the `2.5s` cap.
- Short direction probes at 7ms and 20ms did not improve the n=26..27 seeds `0..20` slice; 7ms regressed it.
- Forcing partition or split caches at n=26..27 regressed prior small-slice checks.
- A fresh n=26..27 exact-cache probe retested forced partition cache, forced split cache, combined partition/split cache, and incumbent pruning on seeds `0..20`. None improved coverage. Split cache regressed n=26 from `40/42` to `38/42`, and combined partition/split cache regressed n=27 from `36/42` to `34/42`.
- Pair-incumbent pruning looked mildly faster on the small slice, but making it a default for n>=26 regressed the full n=27 seeds `0..100` benchmark from `176/202` to `174/202`. That code change was reverted.
- Native compiler flags did not recover the near-boundary n=27 seeds that solved just over 2s.
- Exact-safe switches for tight empty-`S` child-cache reuse, free-hint gating, articulation-arm reduction, common-edge decomposition, inflight pruning, and alternate empty-`S` orderings did not recover additional current timeout seeds.
- Exact budget-probe step `3` was retested on current `n=27`, seeds `0..100`, timeout `2s`; it regressed the full result to `172/202`, so the existing default was kept.

The compact summary for the latest exact-cache probe is retained in `benchmarks/random_n26_27_bottleneck_exact_cache_probe_summary.csv`.

## What Would Be Needed Next

Further large gains likely require a structural change, not only more local caching.

Promising directions include:

- stronger lower bounds that reject impossible branches earlier;
- a shared dynamic program for partition-side checks that can also merge equivalent new states;
- stronger state equivalence so equivalent rotation paths are not explored repeatedly;
- a different empty-`S` search strategy while preserving exactness.

## AStar Comparison

AStarFlipDistance is benchmarked separately on shared-convex inputs.

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
