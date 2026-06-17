#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: tests/java_parity.sh [--build-dir DIR] [--output PATH]

Run Java oracle parity for random n=12..13, seeds 0..5.
EOF
}

build_dir="build"
output="results/parity_flipdist_java_n12_13_seeds0_5.csv"
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

mkdir -p oracle/java/out
javac -cp oracle/java/lib/acm.jar -d oracle/java/out oracle/java/src/*.java

python3 tools/run_flipdist_java_parity_sweep.py \
  --case random \
  --n-min 12 --n-max 13 \
  --seed-min 0 --seed-max 5 \
  --count 1 \
  --cpp-binary "${flipdist}" \
  --bfs-cap 1 \
  --java-out oracle/java/out \
  --java-lib oracle/java/lib/acm.jar \
  --output "${output}"
