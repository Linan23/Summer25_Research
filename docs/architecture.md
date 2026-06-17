# Architecture

FlipDist is a C++ exact solver for rooted binary-tree flip distance. The implementation follows the algorithmic structure of Li and Xia's STACS 2023 FPT algorithm for Convex Flip Distance while using binary trees as the working representation. The Java triangulation oracle is separate and is used only to check C++ answers on small enough cases.

For the primary reference and links, see `docs/references.md`.

## Research Basis

Convex polygon triangulations and rooted full binary trees are equivalent representations for this problem: a flip in a convex triangulation corresponds to a rotation in the binary tree. The Li-Xia algorithm solves the parameterized decision problem "is the flip distance at most `k`?" using structural properties of optimal flip sequences. This repository keeps that exact-search contract and studies practical performance, validation, and bottlenecks in an implementation setting.

## Component Map

| Component | What it does |
|---|---|
| `flipdist` | Main exact solver binary used by benchmarks. |
| `bf_bst` | Small brute-force checker for validation. |
| Java oracle | Independent triangulation-based checker for feasible sizes. |
| Benchmark tools | Python scripts for sweeps, parity checks, plotting, and external comparisons. |
| Curated CSVs | Retained benchmark summaries that document the current empirical state. |

## Source Layout

`src/flipdist/main.cpp` is the `flipdist` CLI entrypoint. It delegates argument parsing, random/comb generation, file input, and JSON-lines output to `testing_cli.*`.

`src/flipdist/bf_bst.*` contains the binary-tree representation, rotation helpers, random/comb construction, canonical serialization, and the `bf_bst` brute-force validator binary.

`src/flipdist/algorithm.*` contains the main exact solver flow. The important parts are:

- `FlipDistTree`: internal tree representation and structural helpers.
- `TreeDistI`: interval subproblem solver.
- `TreeDistS`: side-subproblem solver. The current bottleneck is `S.empty()`, where many recursive choices can appear.
- Partition-driven branching: the solver splits a problem into sides, assigns search budget to each side, and checks whether both sides can fit.

`src/flipdist/memoization.*` owns reusable keying and cache layers used by the recursive solver. `src/flipdist/helpers.*` provides supporting utilities shared by the algorithm and CLI.

## Solver Flow

The public `flipdist` binary emits two directed rows for each case: `a->b` and `b->a`. Each row reports `case_type`, `n`, `seed`, `direction`, `distance`, `time_ms`, `status`, tree encodings, and `max_k`.

The search remains complete for the configured `max_k`. Optimization work should preserve the distance definition, output fields, status meanings, and Li-Xia decomposition. Current bottleneck work is concentrated in `TreeDistS`, where repeated split exploration and partition budget loops dominate hard random cases.

## Validation Layers

`bf_bst` provides brute-force validation for small binary-tree instances.

`oracle/java/` contains a separate triangulation BFS oracle. Python parity tools compare C++ distances against Java for feasible `n` ranges.

Benchmark tools in `tools/` run performance sweeps, Java parity checks, shared-convex AStar comparisons, and plots. Generated outputs are intentionally kept out of git unless they are curated into `benchmarks/`.
