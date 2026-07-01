# FlipDist Solver Source

This directory contains the C++ implementation of the exact FlipDist solver,
the brute-force BST validator, the command-line interface, memoization support,
and generated-case helpers.

Important entry points:

- `main.cpp`: `flipdist` binary entry point.
- `testing_cli.*`: CLI parsing, generated cases, custom tree input, and
  JSON-lines output.
- `algorithm.*`: Li-Xia-structured exact distance implementation.
- `memoization.*`: caches and optional profile counters.
- `bf_bst.*`: binary-tree rotation helpers and brute-force validator binary.

Research contract:

- Preserve exact distance semantics.
- Preserve directed JSON-lines output fields unless an explicit task changes
  the public format.
- Preserve the Li-Xia decomposition unless the research direction deliberately
  changes.
- Validate solver changes with smoke tests and Java parity where feasible.

