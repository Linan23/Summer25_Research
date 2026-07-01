# Tests

This directory contains thin validation wrappers for the research codebase.

Current wrappers:

- `smoke.sh`: builds confidence that the C++ binaries run and preserve basic
  expected behavior.
- `java_parity.sh`: compares feasible random cases against the Java oracle.
- `benchmark_slice.sh`: runs a small benchmark slice to check script paths and
  output format.
- `script_help.py`: verifies maintained scripts expose a usable help surface.

Run from the repository root after building:

```bash
tests/smoke.sh
ctest --test-dir build --output-on-failure
tests/java_parity.sh
tests/benchmark_slice.sh
tests/script_help.py
```

