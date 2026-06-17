# Third-Party Dependencies

`third_party/` is reserved for optional local checkouts and builds. Contents are ignored by git by default so external code does not become part of the tracked research artifact by accident.

AStarFlipDistance is optional. Keep it as an external dependency and pass its binary explicitly to comparison tools:

```bash
python3 tools/run_shared_convex_flipdist_vs_astar_sweep.py \
  --astar-binary /path/to/A_star_for_flipdistance
```

Do not vendor external source or build outputs into the tracked repository unless the project policy changes.
