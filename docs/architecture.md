# Architecture

FlipDist is a C++ exact solver for rooted binary-tree flip distance. The implementation keeps the Li-Xia decomposition structure intact and uses the Java triangulation oracle only for independent validation on feasible sizes.

## Source Layout

`src/flipdist/main.cpp` is the `flipdist` CLI entrypoint. It delegates argument parsing, random/comb generation, file input, and JSON-lines output to `testing_cli.*`.

`src/flipdist/bf_bst.*` contains the binary-tree representation, rotation helpers, random/comb construction, canonical serialization, and the `bf_bst` brute-force validator binary.

`src/flipdist/algorithm.*` contains the main exact solver flow. The important phases are:

- `FlipDistTree`: internal tree representation and structural helpers.
- `TreeDistI`: interval subproblem solver.
- `TreeDistS`: side-subproblem solver, including the current performance-sensitive `S.empty()` path.
- Partition-driven branching: split enumeration and recursive budget exploration.

`src/flipdist/memoization.*` owns reusable keying and cache layers used by the recursive solver. `src/flipdist/helpers.*` provides supporting utilities shared by the algorithm and CLI.

## Solver Flow

The public `flipdist` binary emits two directed rows for each case: `a->b` and `b->a`. Each row reports `case_type`, `n`, `seed`, `direction`, `distance`, `time_ms`, `status`, tree encodings, and `max_k`.

The search remains complete for the configured `max_k`. Optimization work should preserve distance semantics, status semantics, and the Li-Xia decomposition. Current bottleneck work is concentrated in `TreeDistS` where repeated split exploration and partition budget loops dominate hard random cases, especially some `b->a` directions.

## Validation Layers

`bf_bst` provides brute-force validation for small binary-tree instances.

`oracle/java/` contains a separate triangulation BFS oracle. Python parity tools compare C++ distances against Java for feasible `n` ranges.

Benchmark tools in `tools/` run performance sweeps, Java parity checks, shared-convex AStar comparisons, and plots. Generated outputs are intentionally kept out of git unless curated into `benchmarks/`.
