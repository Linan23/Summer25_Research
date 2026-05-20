# Benchmarks

Run all commands from the repository root after building with CMake or `./setup.sh`.

## Random Solver Sweep

Sanity sweep:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 23 --n-max 25 \
  --seed-min 0 --seed-max 20 \
  --timeout-sec 4 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n23_25_seeds0_20_t4_m3.csv
```

Full retained benchmark command:

```bash
python3 tools/sweep_flipdist_limits.py \
  --case random --n-min 23 --n-max 25 \
  --seed-min 0 --seed-max 100 \
  --timeout-sec 2.5 --max-k-mult 3 \
  --cpp-binary ./build/flipdist \
  --output results/random_n23_25_seeds0_100_t2p5_m3.csv
```

The retained headline summary is copied to `benchmarks/random_n23_25_seeds0_100_t2p5_m3.csv`.

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

## Metric Notes

Random sweep rows are directed. A single generated instance produces `a->b` and `b->a` rows. First-direction coverage counts only `a->b`; directed coverage counts both rows.

`status=ok` means the exact solver found a distance within the configured `max_k`. `status=timeout` means the Python harness stopped the process for that directed case. Timing near a fixed threshold such as 2s can vary slightly between runs and machines.
