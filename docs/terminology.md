# Terminology

This page defines project terms in plain language. It is intended for readers who want to run the solver or read benchmark summaries before studying the implementation.

## Problem Terms

| Term | Meaning |
| --- | --- |
| Flip distance | The minimum number of local moves needed to transform one convex polygon triangulation into another. |
| Rotation distance | The same problem viewed on rooted binary trees. A flip in a triangulation corresponds to a rotation in a binary tree. |
| Rooted binary tree | The tree representation used by the C++ solver. Each generated case contains two such trees. |
| Triangulation | A way to split a convex polygon into triangles using non-crossing diagonals. The Java oracle uses this view. |
| Exact solver | A solver that returns the true distance when it finishes within the configured search budget. |
| `max_k` | The largest distance value the exact search is asked to consider. If the true distance is above `max_k`, the run may report `not_found`. |

## Case Terms

| Term | Meaning |
| --- | --- |
| Case | One generated or file-provided pair of trees. |
| Random case | A case generated from two seeded random binary-search-tree insertion orders. Most retained benchmarks use random cases. |
| Simple case | A structured generated case using opposite chain-shaped trees. The CLI also accepts the older name `comb` for compatibility. |
| Complex case | A case that is empirically difficult for the current implementation, usually because it causes many recursive branches before the time cap. |
| Seed | The integer used to regenerate the same random case. |
| Direction | The solver checks each tree pair both ways: `a->b` and `b->a`. These can take different amounts of time even though the exact distance is symmetric. |

## Benchmark Terms

| Term | Meaning |
| --- | --- |
| Directed row | One result for one direction of one case. A random seed normally produces two directed rows: `a->b` and `b->a`. |
| Directed solves | The number of directed rows with `status=ok`. |
| Combined directed coverage | The directed solves across several sizes or benchmarks divided by all directed rows in that group. For example, `576/606 = 95.0%` means 576 of 606 directional runs finished exactly. |
| Solved pair | A seed where both directions finish with `status=ok`. |
| Timeout | The harness stopped the process at the configured wall-clock cap. A timeout is not an incorrect distance. |
| Timing-margin case | A case that finishes very close to a benchmark timeout, so it may cross the cap on another machine or run. |
| Hard-limit analysis | The empirical study of where this implementation stops meeting target coverage under fixed time caps. It is not a mathematical impossibility proof. |
| Curated benchmark summary | A compact CSV retained in `benchmarks/` because it supports a documented project claim. Raw generated outputs stay in ignored `results/`. |

## External Comparison Terms

| Term | Meaning |
| --- | --- |
| AStarFlipDistance | Optional external solver used only for comparison. It is not part of the tracked project source. |
| A* `simple` / `combined` | Algorithm labels used by AStarFlipDistance. They are external names and are kept unchanged in comparison scripts and CSV columns. |
| Paired-row median | A median computed only on rows where both compared solvers produced usable results. This avoids misleading comparisons when timeout subsets differ. |
