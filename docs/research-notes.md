# Research Notes

These notes summarize the current research state without replacing the detailed benchmark records.

## Primary Reference

The implementation builds on Haohong Li and Ge Xia, "An O(3.82^k) Time FPT Algorithm for Convex Flip Distance" (STACS 2023). The paper gives the algorithmic foundation for the exact fixed-parameter search; this repository focuses on an implementation, validation against an independent oracle, and empirical complex-case analysis.

See `docs/references.md` for citation links.

## Solver Contract

FlipDist is an exact solver for rooted binary-tree flip distance. Work in this repository should preserve:

- the Li-Xia decomposition structure;
- exact distance semantics for the configured `max_k`;
- `flipdist` JSON-lines CLI fields and status meanings;
- directed benchmark rows and metric definitions.

## Current Bottleneck

The persistent complex cases concentrate in `TreeDistS` when `S` is empty. Timing-margin improvements from direction ordering and exact budget-probe tuning help selected seeds, but the retained n=23..27 evidence shows many remaining misses still create large amounts of new branch growth.

See `docs/hard-limit-analysis.md` for the current hard-limit evidence and `docs/benchmarks.md` for regeneration commands.

## Research Archive

`tools/research_archive/` contains one-off analysis and replay scripts from earlier optimization passes. They are retained for context, but maintained workflows should live in `tools/` or `tests/`.
