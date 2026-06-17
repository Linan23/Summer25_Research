# Repository Inventory

This inventory classifies the research handoff surface and the local generated surface. It is meant to help new contributors distinguish reproducible project assets from local experiment output.

## Tracked Project Surface

| Class | Paths | Policy |
| --- | --- | --- |
| C++ solver source | `src/flipdist/`, `CMakeLists.txt` | Tracked. Preserve solver semantics, Li-Xia structure, and CLI/output compatibility unless a task explicitly changes them. |
| Setup scripts | `setup.sh`, `scripts/setup_dev.py` | Tracked. Keep one-command setup working from the repository root. |
| Maintained tools | `tools/*.py` | Tracked. Keep non-interactive `--output` behavior and useful `--help` output. |
| Historical research tools | `tools/research_archive/*.py` | Tracked for context. Do not promote new one-off experiments here unless they are useful for future audit. |
| Tests and wrappers | `tests/` | Tracked. Keep wrappers thin and reproducible. |
| Java oracle | `oracle/java/src/`, `oracle/java/lib/acm.jar` | Tracked. `oracle/java/out/` is generated and ignored. |
| Documentation | `README.md`, `docs/`, `results/README.md`, `third_party/README.md` | Tracked. Keep the README as the GitHub entry point and detailed notes in `docs/`. |
| Curated benchmark summaries | `benchmarks/*.csv`, `benchmarks/README.md` | Tracked only when compact and cited by handoff docs. |

## Ignored Local Surface

| Class | Paths | Policy |
| --- | --- | --- |
| Build outputs | `build/`, `build-*`, `cmake-build-*`, CMake generated files | Ignored. Regenerate with CMake or `./setup.sh`. |
| Raw results | `results/*` except `results/README.md` | Ignored. Copy only compact final summaries to `benchmarks/`. |
| Java build output | `oracle/java/out/`, `*.class` | Ignored. Regenerate with `javac` or `./setup.sh`. |
| Python local state | `.venv/`, caches, pytest/mypy caches | Ignored. Regenerate locally. |
| Optional external dependencies | `third_party/*` except `third_party/README.md` | Ignored. Keep AStarFlipDistance optional and external. |
| Scratch files | `a`, `b`, `tmp/`, `temp/`, `*.tmp`, editor swap files | Ignored. Do not rely on them for reproducibility. |

## Removed During Handoff Cleanup

- `.vscode/tasks.json`: stale editor-specific active-file build task that used a non-CMake C++17 workflow.
- Top-level scratch files `a` and `b`.
