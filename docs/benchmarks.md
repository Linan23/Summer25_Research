# Benchmarks

This page records the benchmark and parity commands used to evaluate the FlipDist research solver. The retained metrics are intended to make the implementation reproducible and to separate exactness checks from performance measurements.

Run all commands from the repository root after building with CMake or `./setup.sh`.

## Fast Validation Wrappers

Use these wrappers for routine handoff checks before changing retained benchmark claims:

```bash
tests/smoke.sh
tests/java_parity.sh
tests/benchmark_slice.sh
tests/script_help.py
```

They call the maintained tools below and write generated CSVs under ignored `results/`.

## Random Solver Sweep

Small sanity sweep:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 23 --n-max 25 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 4 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n23_25_seeds0_20_t4_m3.csv
```

Full retained baseline command:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 23 --n-max 25 \
  --seed-min 0 --seed-max 100 \
  --timeout-sec 2.5 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n23_25_seeds0_100_t2p5_m3.csv
```

The retained headline summary is copied to `benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv`.

The current n=23..25 timeout seed probes are retained as compact summaries:

- `benchmarks/random_n23_25_timeout_exact_switch_probe_summary.csv`
- `benchmarks/random_n23_25_remaining_timeouts_t10_summary.csv`

These support the n=23..25 99% hard-limit analysis in `docs/hard-limit-analysis.md`.

Full hard-limit sweep:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 35 \
  --seed-min 0 --seed-max 100 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n26_35_seeds0_100_t2_m3_raw.csv
```

The retained summaries are:

- `benchmarks/random_hard_limit_n23_35_summary.csv`
- `benchmarks/random_n26_35_seeds0_100_t2_m3_summary.csv`
- `benchmarks/random_n26_35_seeds0_100_t2_m3_instances.csv`

Use the instance file when you need to inspect which seeds solved or timed out.

Smaller hard-limit sanity sweeps:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 30 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n26_30_s0_20_t2_m3.csv

python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 31 --n-max 35 \
  --seed-min 0 --seed-max 10 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n31_35_s0_10_t2_m3.csv
```

Interpretation is in `docs/hard-limit-analysis.md`.

Focused `n=26..27` boundary checks:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n26_27_s0_20_t2_m3.csv

python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2.5 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n26_27_s0_20_t2p5_m3.csv
```

The focused boundary summary is retained in `benchmarks/random_n26_27_boundary_summary.csv`.

Full n=26..27 boundary refresh:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 100 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n26_27_seeds0_100_t2_budget_direction_patch.csv
```

The retained current boundary refresh is:

- `benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_summary.csv`
- `benchmarks/random_n26_27_seeds0_100_t2_budget_direction_patch_instances.csv`

The remaining strict-2s timeout seeds from that refresh were checked against the retained `10s` evidence. The compact retained summary is:

- `benchmarks/random_n26_27_remaining_timeouts_t10_after_budget_direction_summary.csv`

To regenerate that summary, retest the timeout seed lists from the n=26..27 instance file with `timeout=10s` and `max_k=3n`, then keep only the compact solved/still-timeout seed summary in `benchmarks/`.

Current exact-cache bottleneck probe:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/probe_default_n26_27_s0_20_t2.csv

FLIPDIST_DEBUG_ENABLE_PARTITION_CACHE=1 \
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/probe_partition_cache_n26_27_s0_20_t2.csv

FLIPDIST_DEBUG_ENABLE_PARTITION_SPLIT_CACHE=1 \
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 26 --n-max 27 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/probe_split_cache_n26_27_s0_20_t2.csv
```

The compact retained result is `benchmarks/random_n26_27_bottleneck_exact_cache_probe_summary.csv`. It records that the exact-cache switches did not improve coverage enough to promote as defaults.

Profile reuse-counter samples:

```bash
cmake -S . -B build-profile -DCMAKE_BUILD_TYPE=Release \
  -DFLIPDIST_REUSE_PROFILE_COUNTERS=ON
cmake --build build-profile -j

FLIPDIST_PROFILE=1 FLIPDIST_PROFILE_ABORT_MS=2200 \
./build-profile/flipdist --case random --n 27 --seed 5 \
  --count 1 --max-k 81 --bfs-cap 1 \
  > results/profile_n27_seed5_reuse_counters_2200.txt
```

Run the same profile command for the hard seeds being studied. The current retained summary is `benchmarks/random_n26_27_profile_reuse_counters_summary.csv`. These counters show how much work is repeated and how much is new branch growth inside `TreeDistS/S.empty()`. The counters are compile-time optional so default benchmark builds do not carry extra hot-path instrumentation.

Focused `n=27..30` 90% feasibility checks:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 27 --n-max 30 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n27_30_s0_20_t2_m3.csv

python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 27 --n-max 30 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 10 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n27_30_s0_20_t10_m3.csv

FLIPDIST_DEBUG_ENABLE_PARTITION_CACHE=1 FLIPDIST_EMPTY_S_ORDER=conflict \
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 27 --n-max 30 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n27_30_s0_20_t2_cache_conflict.csv

FLIPDIST_DEBUG_ENABLE_PARTITION_SPLIT_CACHE=1 \
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 27 --n-max 30 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 2 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n27_30_s0_20_t2_splitcache.csv
```

The focused feasibility summary is retained in `benchmarks/random_n27_30_90_feasibility_summary.csv`.

## Java Parity

Compile the Java oracle:

```bash
mkdir -p oracle/java/out
javac -cp oracle/java/lib/acm.jar -d oracle/java/out oracle/java/src/*.java
```

Feasible oracle checks:

```bash
python3 tools/run_flipdist_java_parity.py \
  --case random --n 12 --seed 0 --count 5 \
  --cpp-binary ./build/flipdist --max-k 36 \
  --java-out oracle/java/out --java-lib oracle/java/lib/acm.jar \
  --output results/parity_flipdist_java_n12_seeds0_5.csv

python3 tools/run_flipdist_java_parity.py \
  --case random --n 13 --seed 0 --count 5 \
  --cpp-binary ./build/flipdist --max-k 39 \
  --java-out oracle/java/out --java-lib oracle/java/lib/acm.jar \
  --output results/parity_flipdist_java_n13_seeds0_5.csv
```

Brute-force parity uses the same oracle:

```bash
python3 tools/run_bruteforce_java_parity.py \
  --case random --n 12 --seed 0 --count 5 \
  --cpp-binary ./build/bf_bst \
  --java-out oracle/java/out --java-lib oracle/java/lib/acm.jar \
  --output results/parity_bruteforce_java_n12_seed0.csv
```

## AStarFlipDistance

AStarFlipDistance is optional. Clone and build it outside the tracked source surface, for example under ignored `third_party/`:

```bash
mkdir -p third_party
# clone/build AStarFlipDistance here with its own upstream instructions
```

Then run shared-convex comparison with an explicit binary:

```bash
python3 tools/run_shared_convex_flipdist_vs_astar_sweep.py \
  --case random --n-min 22 --n-max 25 \
  --seed-min 0 --seed-max 100 \
  --case-timeout-sec 10 \
  --flipdist-binary ./build/flipdist \
  --astar-binary /path/to/A_star_for_flipdistance \
  --output results/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10.csv \
  --summary-output results/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv
```

Plot a retained comparison:

```bash
python3 tools/plot_flipdist_vs_astar_compare.py \
  --input results/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10.csv \
  --output results/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10.svg \
  --summary-output results/shared_convex_flipdist_vs_astar_n22_25_seeds0_100_timeout10_summary.csv
```

A local no-Gurobi AStar build can be used for `simple` and `combined` comparisons when Gurobi is unavailable, but the checkout and patch stay under ignored `third_party/`. The current local comparable summary is retained in `benchmarks/shared_convex_local_flipdist_vs_astar_n22_30_summary.csv`.

The local summary was produced from two bounded sweeps:

```bash
python3 tools/run_shared_convex_flipdist_vs_astar_sweep.py \
  --case random --n-min 22 --n-max 25 \
  --seed-min 0 --seed-max 20 \
  --case-timeout-sec 10 \
  --flipdist-max-k 90 \
  --astar-binary third_party/AStarFlipDistance/build-nogurobi/A_star_for_flipdistance \
  --astar-algos simple,combined \
  --shared-root results/shared_convex_bench \
  --output results/shared_convex_flipdist_vs_astar_n22_25_s0_20.csv \
  --summary-output results/shared_convex_flipdist_vs_astar_n22_25_s0_20_summary.csv

python3 tools/run_shared_convex_flipdist_vs_astar_sweep.py \
  --case random --n-min 26 --n-max 30 \
  --seed-min 0 --seed-max 10 \
  --case-timeout-sec 10 \
  --flipdist-max-k 100 \
  --astar-binary third_party/AStarFlipDistance/build-nogurobi/A_star_for_flipdistance \
  --astar-algos simple,combined \
  --shared-root results/shared_convex_bench \
  --output results/shared_convex_flipdist_vs_astar_n26_30_s0_10.csv \
  --summary-output results/shared_convex_flipdist_vs_astar_n26_30_s0_10_summary.csv
```

When timeout subsets differ, prefer paired-row medians over raw solved-case medians.

## Metric Notes

Random sweep rows are directed. A single generated seed normally produces `a->b` and `b->a` rows. Directed coverage counts both rows. Pair coverage counts the seed only when both directions return `ok`.

`status=ok` means the exact solver found a distance within the configured `max_k`. `status=timeout` means the Python harness stopped the run at the configured wall-clock cap. In the full hard-limit sweep, the cap is applied to the `flipdist` process for one seed; if the process does not return both JSON rows before the cap, the harness records both directions as timeout. Timing near a fixed threshold such as 2s can vary slightly between runs and machines.
