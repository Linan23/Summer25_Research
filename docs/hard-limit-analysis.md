# Hard-Limit Analysis

This note summarizes the current evidence for the practical limit of the Li-Xia-structured FlipDist solver. The solver remains exact and search-complete for the configured `max_k`; the measurements below describe wall-clock coverage under fixed benchmark timeouts.

## Current Baseline

The retained baseline is random rooted binary tree instances with `max_k=3n`, timeout `2.5s`, and seeds `0..100` for `n=23..25`.

- Directed exact solves: `576/606 = 95.05%`.
- Per-instance first-direction coverage: `288/303 = 95.05%`.
- Directed coverage by n: `n=23: 198/202`, `n=24: 186/202`, `n=25: 192/202`.

This preserves the target `>=95%` baseline coverage. Timing near the `2.5s` boundary has measurable run-to-run variance.

## Larger Random Instances

The random hard-limit sweep uses the same `max_k=3n` and `2.5s` timeout.

| n | Seeds | Directed solves | Instance solves | Median solved pair max |
|---|---:|---:|---:|---:|
| 26 | 0..20 | 40/42 = 95.24% | 20/21 = 95.24% | 25.300 ms |
| 27 | 0..20 | 36/42 = 85.71% | 18/21 = 85.71% | 43.280 ms |
| 28 | 0..20 | 26/42 = 61.90% | 13/21 = 61.90% | 51.089 ms |
| 29 | 0..20 | 24/42 = 57.14% | 12/21 = 57.14% | 112.539 ms |
| 30 | 0..20 | 24/42 = 57.14% | 12/21 = 57.14% | 82.502 ms |
| 31 | 0..10 | 10/22 = 45.45% | 5/11 = 45.45% | 39.252 ms |
| 32 | 0..10 | 14/22 = 63.64% | 7/11 = 63.64% | 137.275 ms |
| 33 | 0..10 | 10/22 = 45.45% | 5/11 = 45.45% | 25.470 ms |
| 34 | 0..10 | 10/22 = 45.45% | 5/11 = 45.45% | 192.874 ms |
| 35 | 0..10 | 10/22 = 45.45% | 5/11 = 45.45% | 30.529 ms |

The practical limit is not a clean n cutoff. The `n=26`, seeds `0..20` slice still reaches the `>=95%` target, but hard shape families become common enough by `n=27` that 2.5s coverage again drops below the baseline target. A direction-order probe found that the low-millisecond reverse-direction probe used for the `n=23..25` baseline was counterproductive on several `n>=26` hard cases, so the default now keeps that probe only for `n=23..25` while preserving the `FLIPDIST_DIRECTION_PROBE_THRESHOLD_MS` experiment override. The current n>=29 pass adds dynamic partition-cache defaults, conflict-first empty-S ordering, direct-index target child-range lookup, and targeted direction rules for compact reverse-direction cases. This recovered one hard pair each at `n=29` and `n=30`, improved `n=28` median runtime, and preserved the `n=23..26` coverage guardrails without changing solver semantics.

At the current 2.5s timeout, remaining hard seeds are:

- `n=26`: 40/42 directed solves; timeout seed `14`.
- `n=27`: 36/42 directed solves; timeout seeds `5, 9, 14`.
- `n=28`: 26/42 directed solves; timeout seeds `1, 5, 7, 9, 11, 14, 17, 20`.
- `n=29`: 24/42 directed solves; timeout seeds `4, 7, 10, 11, 13, 14, 16, 17, 20`.
- `n=30`: 24/42 directed solves; timeout seeds `1, 2, 9, 11, 13, 14, 16, 17, 20`.

## n=26..27 Boundary Optimization

The current n=26..27 pass added a narrow direction-order lock for cases where a deeper target with an extreme low root would otherwise be overridden by a later reverse-direction shape rule. This is Li-Xia-preserving: it only changes which directed solve is attempted first, and leaves distance semantics, `max_k`, and search completeness unchanged.

| Scenario | Timeout | n26 | n27 | Combined |
|---|---:|---:|---:|---:|
| Before direction lock | 2s | 38/42 = 90.48% | 34/42 = 80.95% | 72/84 = 85.71% |
| After direction lock | 2s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| After direction lock | 2.5s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| Side-budget cache tightening | 2s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| Side-budget cache tightening | 2.5s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| Empty-S pair-bound propagation | 2s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| Empty-S pair-bound propagation | 2.5s | 40/42 = 95.24% | 36/42 = 85.71% | 76/84 = 90.48% |
| Forced partition cache + conflict-first empty-S order | 2s | 40/42 = 95.24% | 34/42 = 80.95% | 74/84 = 88.10% |
| Forced partition split cache | 2s | 40/42 = 95.24% | 34/42 = 80.95% | 74/84 = 88.10% |

The side-budget cache tightening records smaller known feasible side budgets when pair incumbents or existing success bounds prove a lower feasible value than the current queried budget. Empty-S pair-bound propagation additionally records exact empty-S `TreeDistS` outcomes into the pair-distance bounds cache so repeated partition side checks can reuse them through the cheaper pair key. Both changes are exact-safe and preserve coverage, but representative hard-seed profiles did not show a material enough reduction in `S.empty()` call volume, partition call volume, or partition budget-loop time to recover a persistent n=27 timeout seed. The recovered 2s cases remain timing-margin direction-order cases, not new search-space pruning. Persistent timeout seeds remain:

- `n=26`: seed `14`.
- `n=27`: seeds `5, 9, 14`.

Fresh abort profiles on those persistent seeds show no `TreeDistI` work and the same `TreeDistS/S.empty()` partition wall: hundreds of thousands of `S.empty()` calls, hundreds of thousands of partition calls, and partition budget-loop accumulated time in the multi-second range inside each abort. For example, before side-cache tightening, `n=26 seed=14` reached `849,169` `TreeDistS` calls, `397,484` `S.empty()` calls, `291,997` partition calls, and `6,650.647 ms` accumulated partition budget-loop time in the first profiled direction; after tightening, the comparable first direction remained in the same range at `855,923` `TreeDistS` calls, `400,220` `S.empty()` calls, `294,490` partition calls, and `6,637.117 ms` accumulated budget-loop time. In the current n=27 pair-bound pass, seed `5` still aborts at `373,014` `TreeDistS` calls, `140,159` `S.empty()` calls, `170,224` partition calls, and `1,765.462 ms` partition budget-loop time in `a->b`; seed `14` still reaches `440,385` `TreeDistS` calls and `211,360` `S.empty()` calls in `a->b`.

## 90% Feasibility Check

The current goal tested whether `n=27..30`, seeds `0..20`, can reach at least 90% exact solve coverage with Li-Xia-preserving optimization. Under the strict `2s`, `max_k=3n` benchmark, 90% requires at least `38/42` directed solves per n, or `152/168` directed solves across the combined slice.

Current strict-2s default coverage is `110/168 = 65.48%` across `n=27..30`. This matches the 2.5s coverage, so the misses are not concentrated in the 2.0-2.5s band. A higher-timeout `10s` probe reaches `128/168 = 76.19%`, still 24 directed solves short of 90%. Existing Li-Xia-preserving toggles did not recover additional strict-2s coverage:

| Scenario | Timeout | n27 | n28 | n29 | n30 | Combined |
|---|---:|---:|---:|---:|---:|---:|
| Current default | 2s | 36/42 | 26/42 | 24/42 | 24/42 | 110/168 = 65.48% |
| Current default | 2.5s | 36/42 | 26/42 | 24/42 | 24/42 | 110/168 = 65.48% |
| Current default | 10s | 36/42 | 32/42 | 30/42 | 30/42 | 128/168 = 76.19% |
| Forced partition cache + conflict-first empty-S order | 2s | 34/42 | 26/42 | 24/42 | 24/42 | 108/168 = 64.29% |
| Forced partition split cache | 2s | 34/42 | 26/42 | 22/42 | 24/42 | 106/168 = 63.10% |
| Forced partition cache + conflict-first empty-S order | 2.5s | 36/42 | 26/42 | 24/42 | 24/42 | 110/168 = 65.48% |
| Forced partition split cache | 2.5s | 36/42 | 26/42 | 24/42 | 24/42 | 110/168 = 65.48% |

This marks the current Li-Xia-structured solver as reaching its practical hard limit for `n=27..30` under the established 2s benchmark constraints. The gap is too large to plausibly close with local ordering/cache changes alone: `n=27` is close but has persistent timeouts even at 10s, while `n=28..30` remain far below 90% at 2s, 2.5s, and 10s.

## Hard-Case Profiles

Bounded 5s profiles of persistent timeout seeds `n=28 seed=7` and `n=30 seed=1` show the same failure mode:

- `TreeDistI` calls: `0`.
- Search time is dominated by `TreeDistS` with empty `S`.
- Each 5s abort sees hundreds of thousands of `S.empty()` calls and partition calls.
- Partition budget-loop accumulated time remains high because it is nested under many recursive calls, even when iteration counts are reduced by side-budget caches.
- On `n=24 seed=24`, the direct-index target child-range lookup reduced a representative run from roughly `2.58s` to roughly `2.40s` by cutting hot conflict/free-edge lookup overhead.
- Partition side-order probes (`side1`, `side2`, lower-bound, conflict, size, pair-count modes) did not solve the representative persistent hard seeds within 10s.

Representative 5s abort counters:

| Case | Direction | `TreeDistS` calls | `S.empty()` calls | Partition calls | Budget-loop accumulated time | Duplicate child states |
|---|---|---:|---:|---:|---:|---:|
| `n=28 seed=7` | `a->b` | 877,144 | 350,183 | 356,766 | 10,508.634 ms | 227,898 |
| `n=28 seed=7` | `b->a` | 726,591 | 276,175 | 307,205 | 7,328.334 ms | 56,608 |
| `n=30 seed=1` | `a->b` | 858,874 | 264,672 | 408,979 | 11,572.052 ms | 153,989 |
| `n=30 seed=1` | `b->a` | 737,175 | 316,399 | 268,491 | 10,387.712 ms | 176,668 |

Current evidence points to state-space explosion in `TreeDistS/S.empty()` plus repeated partition-driven branching as the practical wall. Further large gains likely require stronger admissible pruning, stronger state equivalence, or a deeper Li-Xia-preserving decomposition of empty-S branch families.

Concretely, the Li-Xia-preserving implementation is limited by the empty-S branch structure: it still enumerates too many rotation children, repeatedly re-enters partition-driven feasibility checks, and cannot derive enough global impossibility information from the current local lower bounds, side-budget caches, or split signatures. Reaching 90%+ coverage on this slice likely requires changing the preserved structure in one of these ways:

- replace local empty-S rotation enumeration with a stronger global search state or canonical state-space traversal;
- add a materially stronger admissible lower bound that reasons across both partition sides before recursive branching;
- change partition handling from repeated side-budget feasibility recursion into a shared decomposition DAG or equivalent global dynamic program;
- introduce stronger state equivalence/dominance than the current pair keys, child signatures, and incumbent bounds can express.

## AStar Comparison

The retained historical shared-convex summary still shows AStarFlipDistance faster on that prior dataset/build. A fresh local comparison was also run using an ignored no-Gurobi build of `A_star_for_flipdistance` on identical shared-convex inputs, case timeout `10s`.

| n | FlipDist median | A* simple median | A* combined median |
|---|---:|---:|---:|
| 22 | 9.764 ms | 168.500 ms | 163.000 ms |
| 23 | 11.517 ms | 125.500 ms | 103.000 ms |
| 24 | 23.542 ms | 185.000 ms | 128.000 ms |
| 25 | 17.944 ms | 40.500 ms | 38.500 ms |
| 26 | 27.463 ms | 327.500 ms | 86.000 ms |
| 27 | 85.647 ms | 350.500 ms | 343.000 ms |
| 28 | 323.901 ms | 475.500 ms | 260.500 ms |
| 29 | 72.993 ms | 759.000 ms | 779.500 ms |
| 30 | 16.869 ms | 149.000 ms | 160.000 ms |

The raw solved-case median is affected by different timeout subsets; for example, at `n=28` the A* combined median is lower because it solved fewer rows. On paired rows solved by both solvers, FlipDist has faster median runtime than both A* modes for every `n=22..30` in this local sample. The paired counts, status counts, and win counts are retained in `benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`.

A current focused `n=28`, seeds `0..10` shared-convex refresh had FlipDist raw solved median `323.901 ms`, A* simple `475.500 ms`, and A* combined `260.500 ms`; A* timed out on several rows. On comparable paired solved rows, FlipDist median remained faster than both modes: `53.855 ms` vs A* simple `475.500 ms`, and `34.169 ms` vs A* combined `260.500 ms`.

## Validation

Latest validation for the tracked solver state:

- `cmake --build build -j`
- `./build/bf_bst`
- `./build/flipdist --case random --n 12 --seed 0 --count 1 --max-k 30 --bfs-cap 1`
- Java parity: random `n=12..13`, seeds `0..5`
- Full random baseline: `n=23..25`, seeds `0..100`, timeout `2.5s`, `max_k=3n`

Curated CSV summaries are in `benchmarks/random_hard_limit_n23_35_summary.csv` and `benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`.
