# Summer25_Research

## Current Progress

- Accuracy parity is strong in tested ranges: n=12..15, seeds 0..4 matched Java BFS distances on all rows (results/parity_flipdist_vs_java_random_n12_15_seeds0_4.csv)
- Performance is good through most n=23..24 random seeds: 40/40 directions solved in the latest n=23..24, seeds 0..9 sweep (results/flipdist_limits_random_n23_24_seeds0_9_t60_90.csv)
- n=25 is still mixed on hard seeds: some runs solve in low seconds, others take tens of seconds or hit timeout caps
- Current problem: TreeDists still blows up on hard cases when partner set S is empty and recursion falls back to expensive generic branching

## What Is In This Repo

- `build/flipdist`: main C++ FlipDist solver
- `build/bf_bst`: brute-force BST solver (exact but very slow for larger n)
- `triangulation/`: Java BFS oracle
- `scripts/`: parity, sweep, and plotting tools

## Setup

Requirements:
- CMake + C++20 compiler
- Python 3
- Java (for oracle/parity scripts)

Optional for plots:
```bash
python3 -m pip install matplotlib
```

Build:
```bash
cmake -S . -B build
cmake --build build --target flipdist bf_bst
```

## Core Commands

FlipDist single run:
```bash
./build/flipdist --case random --n 12 --seed 0 --count 1 --max-k 30 --bfs-cap 1
```

FlipDist single run with rotation-by-rotation tree visualization:
```bash
./build/flipdist --case random --n 8 --seed 0 --count 1 --max-k 20 --bfs-cap 1 --emit-path
```

Brute-force single run:
```bash
./build/bf_bst --case random --n 12 --seed 0 --count 1
```

FlipDist performance sweep:
```bash
python3 scripts/sweep_flipdist_limits.py \
  --case random \
  --n-min 15 --n-max 16 \
  --seed-min 0 --seed-max 1 \
  --timeout-sec 10 \
  --max-k-mult 3 \
  --bfs-cap 1 \
  --output results/flipdist_limits_random_n15_16_seeds0_1_smoke.csv
```
If you omit `--output`, the script asks `Save CSV results? (y/n)`.

FlipDist vs Java parity:
```bash
python3 scripts/run_flipdist_java_parity_sweep.py \
  --case random \
  --n-min 12 --n-max 12 \
  --seed-min 0 --seed-max 1 \
  --cpp-binary ./build/flipdist \
  --bfs-cap 1 \
  --output results/parity_flipdist_vs_java_random_n12_seed0_1_smoke.csv \
  --print
```
If you omit `--output`, the script asks `Save CSV results? (y/n)`.

Brute-force vs Java parity:
```bash
python3 scripts/run_bruteforce_java_parity_sweep.py \
  --case random \
  --n-min 12 --n-max 12 \
  --seed-min 0 --seed-max 0 \
  --cpp-binary ./build/bf_bst \
  --output results/parity_bf_vs_java_random_n12_seed0_smoke.csv \
  --print
```

## Flag Quick Reference

`flipdist` (`./build/flipdist`):
- `--case random|comb`: random seeded trees or comb trees
- `--n <int>`: node count
- `--seed <int>`: seed for random case
- `--count <int>`: number of pairs
- `--max-k <int>`: max distance budget to try
- `--bfs-cap <int>`: local BFS depth cap used as a bounded fallback (usually keep at `1`)
- `--print-trees`: print source/target tree encodings for both `a->b` and `b->a`
- `--emit-path`: print each rotation step with tree visualization
  - includes start tree, each post-rotation tree, final resulting tree, and target tree
  - also includes a per-step table: node parent left right range

Example:
```bash
./build/flipdist --case random --n 12 --seed 0 --count 1 --max-k 30 --bfs-cap 1
```
Example with path/tree visualization:
```bash
./build/flipdist --case random --n 8 --seed 0 --count 1 --max-k 20 --bfs-cap 1 --emit-path
```

`bf_bst` (`./build/bf_bst`):
- `--case random`: only random currently
- `--n <int>`, `--seed <int>`, `--count <int>`
Example:
```bash
./build/bf_bst --case random --n 12 --seed 0 --count 1
```

`sweep_flipdist_limits.py`:
- `--timeout-sec`: base timeout per run
- `--high-timeout-sec`: timeout above `--n-threshold`
- `--n-threshold`: switch point for high timeout/max-k settings
- `--max-k-mult`: base `max_k = n * mult`
- `--high-max-k-mult`: multiplier above threshold
- `--retry-max-k-mults`: optional retry `max_k` multipliers
- `--retry-timeout-mult` or `--retry-timeout-mults`: retry timeout scaling
- `--output`: optional CSV path (if omitted, script prompts yes/no to save)

Example:
```bash
python3 scripts/sweep_flipdist_limits.py \
  --case random \
  --n-min 15 --n-max 16 \
  --seed-min 0 --seed-max 1 \
  --timeout-sec 10 \
  --max-k-mult 3 \
  --bfs-cap 1 \
  --output results/flipdist_limits_random_n15_16_seeds0_1_smoke.csv
```

`run_flipdist_java_parity.py`:
- `--output` is optional
- `--java-time-limit` is important for heavy seeds

Example:
```bash
python3 scripts/run_flipdist_java_parity.py \
  --case random \
  --n 12 \
  --seed 0 \
  --count 1 \
  --cpp-binary ./build/flipdist \
  --bfs-cap 1 \
  --output results/parity_flipdist_vs_java_random_n12_seed0.csv \
  --print
```

`run_bruteforce_java_parity.py`:
- `--output` is optional
- `--java-time-limit` is important for heavy seeds

Example:
```bash
python3 scripts/run_bruteforce_java_parity.py \
  --case random \
  --n 12 \
  --seed 0 \
  --count 1 \
  --cpp-binary ./build/bf_bst \
  --output results/parity_bf_vs_java_random_n12_seed0.csv
```

`run_flipdist_java_parity_sweep.py`:
- `--output` is optional (if omitted, script prompts yes/no to save)
- `--print` shows progress per `(n, seed)`

Example:
```bash
python3 scripts/run_flipdist_java_parity_sweep.py \
  --case random \
  --n-min 12 --n-max 12 \
  --seed-min 0 --seed-max 1 \
  --cpp-binary ./build/flipdist \
  --bfs-cap 1 \
  --output results/parity_flipdist_vs_java_random_n12_seed0_1_smoke.csv \
  --print
```

`run_bruteforce_java_parity_sweep.py`:
- `--output` is optional (if omitted, script prompts yes/no to save)
- `--print` shows each parity run output

Example:
```bash
python3 scripts/run_bruteforce_java_parity_sweep.py \
  --case random \
  --n-min 12 --n-max 12 \
  --seed-min 0 --seed-max 0 \
  --cpp-binary ./build/bf_bst \
  --output results/parity_bf_vs_java_random_n12_seed0_smoke.csv \
  --print
```

## Useful Debug Env Vars

- `FLIPDIST_PROFILE=1`: prints phase-level counters/timing
- `FLIPDIST_PROFILE_ABORT_MS=<ms>`: cuts profile runs early to inspect bottlenecks
- `FLIPDIST_DEBUG=1`: verbose debug logging
