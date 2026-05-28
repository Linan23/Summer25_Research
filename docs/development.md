# Development

This project is both a solver implementation and a research benchmark package. Keep solver behavior unchanged unless a task explicitly asks for a semantic change. Because FlipDist is an exact solver, performance work must preserve the distance definition, directed output rows, and search-complete behavior for the configured `max_k`.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Use `./setup.sh` for full environment setup, including Python dependencies and Java oracle compilation.

## Smoke Tests

```bash
./build/bf_bst
./build/flipdist --case random --n 12 --seed 0 --count 1 --max-k 30
```

Run Java parity after changes that touch solver behavior:

```bash
python3 tools/run_flipdist_java_parity.py \
  --case random --n 12 --seed 0 --count 5 \
  --cpp-binary ./build/flipdist --max-k 36 \
  --java-out oracle/java/out --java-lib oracle/java/lib/acm.jar \
  --output results/parity_flipdist_java_n12_seeds0_5.csv
```

## Optimization Guidelines

The current bottleneck is `TreeDistS` in `S.empty()`. Hard cases spend most of their time repeating partition-side checks and budget-loop work. Changes in this area should be paired with:

- Java parity on feasible oracle ranges.
- A random sweep sanity run for `n=23..25`, seeds `0..20`.
- Profile counters or timing summaries showing reduced partition budget-loop calls or time.
- Wall-time comparison against the retained current benchmark.

Prefer narrow changes that fit the existing Li-Xia decomposition, cache structure, and CLI output contract. Avoid algorithm rewrites during optimization passes unless the project explicitly changes direction.

## Before Handing Off A Change

Use this checklist for code changes:

- Build in Release mode.
- Run the two smoke tests.
- Run Java parity if solver behavior changed.
- Run a small random sweep before a full benchmark.
- Update retained benchmark files only when the run is meant to document project state.
- Update docs when a result changes the current performance or hard-limit story.

## Maintained Script Expectations

Supported scripts live in `tools/`. They should remain non-interactive when `--output` is supplied and should keep a useful `--help` surface. Experimental one-off analysis belongs in `tools/research_archive/` or outside the repository.

Generated benchmark outputs belong in `results/`, which is ignored. Curate only stable summaries into `benchmarks/`.
