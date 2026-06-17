#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tests/smoke.sh [--build-dir DIR]

Run the fast C++ smoke checks from the repository root.
EOF
}

build_dir="build"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      build_dir="${2:?missing value for --build-dir}"
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
  bf_bst="${build_dir}/bf_bst"
  flipdist="${build_dir}/flipdist"
else
  bf_bst="./${build_dir}/bf_bst"
  flipdist="./${build_dir}/flipdist"
fi

"${bf_bst}"
"${flipdist}" --case random --n 12 --seed 0 --count 1 --max-k 30 --bfs-cap 1
