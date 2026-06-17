#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tests/benchmark_slice.sh [--build-dir DIR] [--output PATH]

Run a small random FlipDist benchmark slice for path and harness checks.
EOF
}

build_dir="build"
output="results/refactor_benchmark_slice_random_n23_seed0_2.csv"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      build_dir="${2:?missing value for --build-dir}"
      shift 2
      ;;
    --output)
      output="${2:?missing value for --output}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

cd "$(dirname "$0")/.."

if [[ "${build_dir}" = /* ]]; then
  flipdist="${build_dir}/flipdist"
else
  flipdist="./${build_dir}/flipdist"
fi

python3 tools/sweep_flipdist_limits.py \
  --case random \
  --n-min 23 --n-max 23 \
  --seed-min 0 --seed-max 2 \
  --timeout-sec 4 \
  --max-k-mult 3 \
  --cpp-binary "${flipdist}" \
  --bfs-cap 1 \
  --output "${output}"
