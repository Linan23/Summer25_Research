**10/14/25 - Brute Force optimizations**  
- Rewrote the BFS queue entries to store a cached hash and parent hash so we skip undo moves cheaply and reuse the hash for visited checks (now 64-bit FNV instead of string keys).  
- Added rotation-aware pruning: when both trees are comb-shaped we only expand rotations that immediately create a target edge; otherwise we allow any move that doesn’t reduce the count of matching target edges, and fall back to the unfiltered expansion if the filter would block every option.  
- Exposed a richer telemetry struct (`BFSRun`, `BFSStats`) and environment knobs (`MY_BFS_TIME_LIMIT`, `MY_BFS_VISITED_CAP`, `MY_BFS_QUEUE_CAP`, plus random-specific variants) so runs can be tuned without recompiling.  
- Updated the test harness to honor those env vars, print duplicate/enqueue/queue-peak metrics, and keep the brute-force random sweep from hanging when a filter would previously starve the frontier.  
- Added inline documentation across the brute-force code path (bfs.cpp, distance.cpp, helpers) to clarify each helper’s purpose.

**10/15/25 - Brute Force solver status**  
- Works well for structured comb trees: With `MY_BFS_FILTER_ROTATIONS=1`, bidirectional fallback, symmetry pruning, and the transposition table, the solver produces exact distances through n = 17 and only times out when n = 18 exceeds the 12 s cap.  
- Random tree pairs remain correct up to n≈12–13; baseline, optimized, and hashed BiBFS modes agree in the new unit tests.  
- Runtime still balloons on random inputs (state explosion). We need stronger heuristics (A*/IDA* using RP ranges + edge deficits), automatic hand-off to hashed BiBFS when queues explode, and finer-grained duplicate/symmetry pruning (rotation ordering, better incremental updates).  
- Maybe integrate an admissible A*/IDA* pipeline, promote BiBFS into the default solver, expand symmetry pruning (canonical rotations), and explore lightweight parallelisation of the two frontiers.




**Runtime knobs** (set via environment variables):
- `MY_BFS_MODE`: `baseline` forces the original FIFO BFS; anything else uses the optimized solver (hashes, rotation filters, bidir fallback).
- `MY_BFS_FILTER_ROTATIONS=1`: enable comb-specific pruning (only allow rotations that add target edges, with a fallback for dead fronts).
- `MY_BFS_TIME_LIMIT`: per-run time cap used by `BFSSearchCapped` (seconds).
- `MY_BFS_VISITED_CAP` / `MY_BFS_QUEUE_CAP`: hard limits on the visited set / queue size; hitting them terminates the search with a `CAP` flag.
- `MY_BFS_USE_BIDIR=1`: after the single-ended BFS times out, automatically try the hashed bidirectional search before giving up.
- `MY_BFS_TRANSPOSITION_CAP`: size of the per-hash table that stores the best depth seen; bypassing worse duplicates avoids re-expansion.
- `MY_BFS_SYMMETRY_PRUNE=1`: skip rotations that would only mirror identical subtrees, keeping search on a canonical path.
- Other flags: `MY_BFS_RANDOM_USE_BIDIR`, `MY_BFS_RANDOM_*` variants mirror the main switches for the random benchmarks.

**BFS variants**  
- `BFSSearchBaseline`: plain FIFO BFS (no hashing, no heuristics).  
- `BFSSearchOptimized`: hashed visited set + filters + optional best-first scoring, symmetry prune, transposition table; falls back to `BiBFSSearchHashed` on timeout when enabled.  
- `BiBFSSearchHashed`: meet-in-the-middle search—grow one frontier from the start and another from the target, and when they meet you sum the depths. Because each frontier only needs ≈half the number of steps of a single-ended BFS, it prunes deep comb states much faster. The hashed version reuses the same pruning (rotation filters, symmetry checks, transposition table) on both sides so the frontiers stay in sync.
