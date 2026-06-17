#!/usr/bin/env python3
"""Check that documented command-line scripts expose --help."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate --help for maintained FlipDist scripts.")
    parser.add_argument(
        "--include-archive",
        action="store_true",
        help="Also check tools/research_archive/*.py.",
    )
    return parser.parse_args()


def main() -> int:
    cfg = parse_args()
    scripts = [ROOT / "scripts" / "setup_dev.py", *sorted((ROOT / "tools").glob("*.py"))]
    if cfg.include_archive:
        scripts.extend(sorted((ROOT / "tools" / "research_archive").glob("*.py")))

    failures: list[Path] = []
    for script in scripts:
        proc = subprocess.run(
            [sys.executable, str(script.relative_to(ROOT)), "--help"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        status = "ok" if proc.returncode == 0 else f"failed rc={proc.returncode}"
        print(f"{status}: {script.relative_to(ROOT)}")
        if proc.returncode != 0:
            failures.append(script.relative_to(ROOT))

    if failures:
        print("scripts without working --help:", file=sys.stderr)
        for script in failures:
            print(f"  {script}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
