# FlipDist Browser Visualizer

This tool opens a separate browser window for generating and inspecting
FlipDist cases. It shows Tree A, Tree B, the selected directed row, and a
BST-style animation when the browser can reconstruct a small rotation path.

The browser asks the local `serve.py` process to run the existing `flipdist`
binary. The exact distance still comes from `flipdist`; browser-side BFS is
used only to find a display path for small cases.

For larger solved cases, the server may reconstruct a witness path by repeatedly
calling the existing solver on one-rotation neighbors. This is a visualizer-only
certificate path; it does not add a public `flipdist` CLI mode and does not
change the solver's search logic.

## Build First

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Open The Visualizer

```bash
python3 tools/visualizer/serve.py
```

The server prints the local URL and opens the browser by default. Use
`--no-open` when running in a terminal-only environment.

If the solver binary is somewhere else, pass it explicitly:

```bash
python3 tools/visualizer/serve.py --solver /path/to/flipdist
```

## Generate Cases

Use the controls at the top of the page:

- `n`: number of tree nodes, from 5 to 60.
- `Case type`: `Easy case` uses the solver's `simple` generator; `Hard case`
  uses the solver's random generator.
- `Seed`: used for hard/random cases.
- `max_k`: optional; leave blank to use the Li-Xia/STT diameter budget
  `2n - 6` for this repository's tree-node `n`.

Click `Generate`. The visualizer runs one solver case and renders both directed
rows, `a->b` and `b->a`.

The server receives normal `flipdist` JSON-lines rows with these fields:

- `case_type`
- `n`
- `seed`
- `direction`
- `distance`
- `status`
- `time_ms`
- `tree_a`
- `tree_b`
- `max_k`

## Animation Limits

Path animation is shortest-path only. When the server can certify witness steps,
the page animates those steps. Otherwise the browser automatically tries to
reconstruct a path whose length matches the solver distance for `n <= 14`.
Hard cases can still reach the state cap.

When a case is too large, not solved, or reaches the cap, the visualizer still
shows Tree A, Tree B, directed status, distance, and runtime.
