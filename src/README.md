# Source Code

This directory contains maintained implementation code. The active solver code
lives in `src/flipdist/`.

Build from the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Do not place generated benchmark output, scratch experiments, or external
dependency checkouts in this directory.

