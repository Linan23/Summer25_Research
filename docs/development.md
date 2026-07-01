# Development

This project is a research implementation and benchmark package for exact convex flip distance, based on the Li-Xia FPT framework summarized in `docs/references.md`. Keep solver behavior unchanged unless a task explicitly asks for a semantic change. Because FlipDist is an exact solver, performance work must preserve the distance definition, directed output rows, and search-complete behavior for the configured `max_k`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Use `./setup.sh` for full environment setup, including Python dependencies and Java oracle compilation.

## Smoke Tests

```bash
tests/smoke.sh
ctest --test-dir build --output-on-failure
```

Run Java parity after changes that touch solver behavior:

```bash
tests/java_parity.sh
```

Run a small benchmark slice after path or harness changes:

```bash
tests/benchmark_slice.sh
```

Check the maintained script help surface:

```bash
tests/script_help.py
```

## Optimization Guidelines

The current implementation bottleneck is `TreeDistS` in `S.empty()`. Complex cases spend most of their time repeating partition-side checks and budget-loop work. Changes in this area should be paired with:

- Java parity on feasible oracle ranges.
- A random sweep sanity run for `n=23..25`, seeds `0..20`.
- Profile counters or timing summaries showing reduced partition budget-loop calls or time.
- Wall-time comparison against the retained current benchmark.

Prefer narrow changes that fit the existing Li-Xia decomposition, cache structure, and CLI output contract. Avoid algorithm rewrites during optimization passes unless the research direction explicitly changes.

## Before Handing Off A Change

Use this checklist for code changes:

- Build in Release mode.
- Run `tests/smoke.sh`.
- Run Java parity if solver behavior changed.
- Run `tests/benchmark_slice.sh` before a full benchmark.
- Update retained benchmark files only when the run is meant to document project state.
- Update docs when a result changes the current performance or hard-limit story.

## Maintained Script Expectations

Supported scripts live in `tools/`. They should remain non-interactive when `--output` is supplied and should keep a useful `--help` surface. Experimental one-off analysis belongs in `tools/research_archive/` or outside the repository.

Test wrappers live in `tests/`. They should be thin, reproducible wrappers over the public commands documented in `README.md` and `docs/benchmarks.md`.

Generated benchmark outputs belong in `results/`, which is ignored. Curate only stable summaries into `benchmarks/`.

## Browser Visualizer

`tools/visualizer/` is a local browser tool for generating and inspecting small `flipdist` cases. Its Python launcher serves the page and exposes a local endpoint that runs the existing `flipdist` binary. The page displays Tree A and Tree B from the solver's normal `tree_a` and `tree_b` fields, then reconstructs small rotation paths in the browser for visualization only.

Launch it from the repository root:

```bash
python3 tools/visualizer/serve.py
```

Do not treat browser-side path reconstruction as a benchmark or correctness oracle. The C++ solver remains the source of truth for exact distance, status, and runtime.
