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


**Runtime knobs** (set via environment variables):
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

**BFS variants**  
- `BFSSearchBaseline`: plain FIFO BFS (no hashing, no heuristics).  
- `BFSSearchOptimized`: hashed visited set + filters + optional best-first scoring, symmetry prune, transposition table; falls back to `BiBFSSearchHashed` on timeout when enabled.  
- `BiBFSSearchHashed`: meet-in-the-middle search—grow one frontier from the start and another from the target, and when they meet you sum the depths. Because each frontier only needs ≈half the number of steps of a single-ended BFS, it prunes deep comb states much faster. The hashed version reuses the same pruning (rotation filters, symmetry checks, transposition table) on both sides so the frontiers stay in sync.

**Comparison script (`scripts/run_compare.py`)**  
- Runs both C++ (`test_asan`) and Java (TriangulationFlip) solvers on the same seeds, capturing solver, timing, and queue stats.  
- Default flow: launch pure BFS; if it times out or hits caps, automatically re-run the same instance with `BiBFSSearchHashed` to recover the distance.  
- `--no-auto-bidir` disables the replay, so the CSV shows raw timeouts (`distance_cpp=-1`).  
- Accepts `--case random|comb`, `--n`, `--count`, `--seed`, `--output`, plus the env vars above for finer control (e.g., to force pure BFS or raise time/space caps).  
- CSV columns include `solver_cpp/solver_java` (bfs/bidir), `status_cpp/status_java` (ok/timeout/cap), and the canonical tree strings so mismatches are reproducible.
