# Research Tools

This directory contains maintained Python tools for reproducing and inspecting
the computational experiments around FlipDist.

Main tool classes:

- Benchmark sweeps over generated cases.
- Java oracle parity checks.
- AStarFlipDistance comparison helpers.
- Plotting and summary generation.
- Local browser visualization in `tools/visualizer/`.

Expectations:

- Keep tools non-interactive when `--output` is provided.
- Keep `--help` useful for new researchers.
- Write generated outputs to `results/` by default.
- Promote only compact, documented summaries to `benchmarks/`.
- Put historical one-off scripts in `tools/research_archive/` when they are
  useful for audit but not part of the main workflow.

