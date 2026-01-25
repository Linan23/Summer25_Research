**01/25/26 - Current performance snapshot (FlipDist)**  
- Random cases, seeds 0–9, `max_k=2n`, 10s cap: 100% finish for n ≤ 20; n=21 finishes 16/20, n=22–23 finish 14/20 (rest hit the time cap)
- With a 30s cap (sample seeds 0–4): n=21 finishes all, n=22 finishes most, n=23 still needs more time.  
- Simple cases remain instant (shortcut distance = n−1) through n=25 but Hard Cases will need to be looked more into

**11/05/25 - Paper Method branching & status check**  

- Integrated the Paper Method branching scheme on trees: FlipDist now assembles a maximal independent set of conflicting diagonals, expands the partner-set branches recursively, and splits shared diagonals into independent subproblems before exploring any flips.  
- Memoisation and the Paper Method lower bound are active on every `(start,target,k)` call, and the CLI continues to emit canonical paths via `--emit-path` / `--path-ascii`. 

- Known limitations (updated): hard random instances beyond ~20 nodes can still exceed the 10s cap; ear-contraction is optional and mostly helps when the conflict set is large.  

**11/04/25 - FlipDist parity harness & CLI**  

- `flipdist_asan` now accepts the same comparison flags as `compare_cli` (`--case`, `--n`, `--count`, `--seed`, `--program`, plus optional `--max-k` / `--bfs-cap`). Each run emits JSON rows for both `a->b` and `b->a`, including canonical tree strings, FlipDist and BFS distances, timing, and the `max_k` budget used. When FlipDist fails to match the BFS oracle the CLI falls back to the BFS distance and logs the discrepancy.  
- Pass `--emit-path` to have `flipdist_asan` reconstruct the minimal BFS rotation path for both directions; the canonical trees and rotation directions are logged to stderr and embedded in the JSON output (`path`, `moves`, `path_complete`). Add `--path-ascii` to render the full ASCII tree for every step, mirroring the `test_asan` demo output.  
- FlipDist memoises `(start,target,k)` subproblems, applies the Paper lower bound before branching, decomposes around diagonals shared by both triangulations, and constructs a maximal independent set of conflicting diagonals prior to recursion. Each shared edge spawns two independent subinstances whose budgets are validated with bounded BFS before the FPT branches are explored.  

**10/22/25 - Demo tooling & bidirectional defaults**  

- Added a CLI-style demo (`BFS_DEMO=1`) that generates comb or random pairs, prints ASCII trees, and reconstructs every rotation along the minimal path in both directions.  
- Hooked automatic bidirectional preference into `BFSSearch`: large non-comb inputs now run `BiBFSSearchHashed` first (`BFS_BIDIR_PREFER_THRESHOLD` controls the cutoff).  
- Exposed `BFS_AUTO_FILTER_THRESHOLD` so the “don’t lose shared edges” heuristic can be toggled off when studying pure BFS blow-ups.  
- Documented all environment variables in one place; demo output now makes parity checks against the Java tool trivial by logging solver type, queue sizes, and step-by-step transformations.

**10/15/25 - Brute Force solver status**  

- Works well for structured comb trees: With `BFS_FILTER_ROTATIONS=1`, bidirectional fallback, symmetry pruning, and the transposition table, the solver produces exact distances through n = 17 and only times out when n = 18 exceeds the 12 s cap.  
- Random tree pairs remain correct up to n≈12–13; baseline, optimized, and hashed BiBFS modes agree in the new unit tests.  
- Runtime still balloons on random inputs (state explosion). We need stronger heuristics (A*/IDA* using RP ranges + edge deficits), automatic hand-off to hashed BiBFS when queues explode, and finer-grained duplicate/symmetry pruning (rotation ordering, better incremental updates).  
- Maybe integrate an admissible A*/IDA* pipeline, promote BiBFS into the default solver, expand symmetry pruning (canonical rotations), and explore lightweight parallelisation of the two frontiers.

**10/14/25 - Brute Force optimizations**  

- Rewrote the BFS queue entries to store a cached hash and parent hash so we skip undo moves cheaply and reuse the hash for visited checks (now 64-bit FNV instead of string keys).  
- Added rotation-aware pruning: when both trees are comb-shaped we only expand rotations that immediately create a target edge; otherwise we allow any move that doesn’t reduce the count of matching target edges, and fall back to the unfiltered expansion if the filter would block every option.  
- Exposed a richer telemetry struct (`BFSRun`, `BFSStats`) and environment knobs (`BFS_TIME_LIMIT`, `BFS_VISITED_CAP`, `BFS_QUEUE_CAP`, plus random-specific variants) so runs can be tuned without recompiling.  
- Updated the test harness to honor those env vars, print duplicate/enqueue/queue-peak metrics, and keep the brute-force random sweep from hanging when a filter would previously starve the frontier.  
- Added inline documentation across the brute-force code path (bfs.cpp, distance.cpp, helpers) to clarify each helper’s purpose.

**BFS variants**  
- `BFSSearchBaseline`: plain FIFO BFS (no hashing, no heuristics).  
- `BFSSearchOptimized`: hashed visited set + filters + optional best-first scoring, symmetry prune, transposition table; falls back to `BiBFSSearchHashed` on timeout when enabled.  
- `BiBFSSearchHashed`: meet-in-the-middle search—grow one frontier from the start and another from the target, and when they meet you sum the depths. Because each frontier only needs ≈half the number of steps of a single-ended BFS, it prunes deep comb states much faster. The hashed version reuses the same pruning (rotation filters, symmetry checks, transposition table) on both sides so the frontiers stay in sync.


## Getting Started

1. **Build the tools** (requires a C++17 compiler such as `clang++`, plus `make`):
   ```bash
   make flipdist_asan   # flip-distance solver with ASan
   make test_asan       # BFS demo 
   ```
   intermediate objects live under `build/asan/`.
2. **Quick smoke test**
   ```bash
   ./flipdist_asan --case random --n 5 --count 1 --emit-path
   python3 scripts/run_compare.py --cpp-binary ./flipdist_asan --cpp-program flipdist_asan \
       --case random --n 6 --count 4 --seed 101 --output results/flip_vs_java.csv --no-auto-bidir
   ```

## Key binaries & scripts

### `flipdist_asan`
- Main Paper Method solver/CLI. Key options:

  - `--case {comb|random}` (default `comb`)
  - `--n N` (required), `--count C` (default `1`) How many Pairs of trees , `--seed S` RNG Trees
  - `--max-k K`: cap on the FPT search budget (default `2*n+6`)
  - `--time-limit`, `--visited-cap`, `--queue-cap`, `--bidir-cap`: limits for the BFS oracle and fallback
  - `--emit-path`: include canonical path + moves in the JSON output
  - `--path-ascii`: render ASCII trees for every step (requires `--emit-path`)

**Example:** `./flipdist_asan --case random --n 6 --count 2 --seed 101 --emit-path --path-ascii`

### `test_asan`
- Historical BFS harness and demo entry point.
- Set `BFS_DEMO=1` to view a step-by-step BFS run; configure with the `BFS_DEMO_*` variables (size, seeds, ASCII output, step caps).

**Example:** `BFS_DEMO=1 BFS_DEMO_MODE=comb BFS_DEMO_N=6 ./test_asan`


### `scripts/run_compare.py`
- Cross-checks the C++ and Java solvers, writing CSV telemetry under `results/`.
- Important flags:
  - `--case {comb|random}`: input family (comb pairs or random pairs).
  - `--n N`: number of internal nodes per tree (required).
  - `--count C`: number of instances; each instance yields two rows (A→B and B→A).
  - `--seed S`: base RNG seed; cases `seed + i` are generated when `count > 1`.
  - `--output PATH`: destination CSV (directories created automatically).
  - `--cpp-binary PATH`: C++ executable to invoke (default `cmake-build/compare_cli` or whatever you specify).
  - `--cpp-program NAME`: label recorded in the CSV (`flipdist_asan`, `test_asan`, etc.).
  - `--time-limit`, `--visited-cap`, `--queue-cap`: forwarded to the C++ CLI to bound BFS runs (seconds / nodes / queue size).
  - `--use-bidir`: force the C++ side to rerun failures with hashed BiBFS.
  - `--bidir-cap`: cap the BiBFS state count when `--use-bidir` is active.
  - `--no-auto-bidir`: *disable* the automatic rerun the script normally performs after a timeout/cap hit.
  - `--java-binary`, `--java-out`, `--java-lib`: override the Java runtime, compiled classes, or supporting jars if you compile the Java tool elsewhere.

- **Example:** `python3 scripts/run_compare.py --cpp-binary ./flipdist_asan --cpp-program flipdist_asan --case comb --n 8 --count 3 --time-limit 20 --output results/comb_n8.csv`

## Environment knobs

These variables can be exported once or supplied inline before a command:

| Variable | Purpose |

| `BFS_MODE` | `baseline` enforces the plain FIFO BFS; any other value keeps the optimized solver. |
| `BFS_USE_BIDIR`, `BFS_RANDOM_USE_BIDIR` | Toggle automatic fallback to hashed BiBFS when caps are exceeded. |
| `BFS_BIDIR_PREFER_THRESHOLD` | Prefer BiBFS immediately when both trees exceed this node count. |
| `BFS_AUTO_FILTER_THRESHOLD` | Enable rotation filtering above this size. |
| `BFS_FILTER_ROTATIONS`, `BFS_SYMMETRY_PRUNE`, `BFS_TRANSPOSITION_CAP`, `BFS_HEURISTIC` | Fine-tune pruning / ordering inside the optimized BFS. |
| `BFS_TIME_LIMIT`, `BFS_VISITED_CAP`, `BFS_QUEUE_CAP` | Core limits enforced by `BFSSearchCapped` (seconds / visited nodes / queue size). |
| `BFS_RANDOM_TIME_LIMIT`, `BFS_RANDOM_VISITED_CAP`, `BFS_RANDOM_QUEUE_CAP` | Random-case versions of the limits above. |
| `BFS_DEMO`, `BFS_DEMO_MODE`, `BFS_DEMO_N`, `BFS_DEMO_SEED_A`, `BFS_DEMO_SEED_B`, `BFS_DEMO_ASCII`, `BFS_DEMO_MAX_STEPS` | 


Descriptions of each: 

- `BFS_MODE` – `baseline` forces the original FIFO BFS; any other value keeps the optimized solver (hashing, pruning, fallback).
- `BFS_USE_BIDIR` – enable/disable auto hand-off to hashed bidirectional BFS when caps are hit.
- `BFS_BIDIR_PREFER_THRESHOLD` – if both trees have at least this many original nodes (and are not combs), run BiBFS immediately.
- `BFS_AUTO_FILTER_THRESHOLD` – auto-enable the monotone rotation filter above this size; set high to keep the search pure.
- `BFS_FILTER_ROTATIONS`, `BFS_SYMMETRY_PRUNE`, `BFS_TRANSPOSITION_CAP`, `BFS_HEURISTIC` – optional pruning/ordering toggles for the optimized engine.
- `BFS_TIME_LIMIT`, `BFS_VISITED_CAP`, `BFS_QUEUE_CAP` – `BFSSearchCapped` safety limits; when tripped the run terminates (and may fall back if allowed).
- Random probe mirrors: `BFS_RANDOM_USE_BIDIR`, `BFS_RANDOM_TIME_LIMIT`, `BFS_RANDOM_VISITED_CAP`, `BFS_RANDOM_QUEUE_CAP`.
- Demo helpers: `BFS_DEMO`, `BFS_DEMO_MODE`, `BFS_DEMO_N`, `BFS_DEMO_SEED_A/B`, `BFS_DEMO_ASCII`, `BFS_DEMO_MAX_STEPS`.  
  - `BFS_DEMO_MODE` (`random` | `comb`) chooses random vs. comb trees.  
  - `BFS_DEMO_N` sets the number of internal nodes.  
  - `BFS_DEMO_SEED_A/B` control the RNG seeds for the two random trees (ignored in comb mode).  
  - `BFS_DEMO_ASCII` toggles ASCII rendering; `BFS_DEMO_MAX_STEPS` limits how many intermediate trees are printed.

Example usage:

BFS_MODE=baseline BFS_TIME_LIMIT=15 ./test_asan
BFS_DEMO=1 BFS_DEMO_MODE=random BFS_DEMO_N=6 ./test_asan
