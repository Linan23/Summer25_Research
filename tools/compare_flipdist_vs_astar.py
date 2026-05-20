#!/usr/bin/env python3
import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def parse_int_list(s: str) -> List[int]:
    out: List[int] = []
    for part in s.split(','):
        part = part.strip()
        if not part:
            continue
        out.append(int(part))
    return out


def quantile(values: List[float], q: float) -> float:
    if not values:
        return float('nan')
    if len(values) == 1:
        return float(values[0])
    values = sorted(values)
    pos = (len(values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    frac = pos - lo
    return values[lo] * (1.0 - frac) + values[hi] * frac


def load_flipdist_jsonl(paths: List[Path]) -> Dict[Tuple[int, int], Dict[str, object]]:
    rows: Dict[Tuple[int, int], Dict[str, object]] = {}
    for path in paths:
        with path.open('r', encoding='utf-8') as fh:
            for line in fh:
                line = line.strip()
                if not line.startswith('{'):
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError:
                    continue
                n = int(obj.get('n', -1))
                seed = int(obj.get('seed', -1))
                direction = obj.get('direction', '')
                if n < 0 or seed < 0:
                    continue
                key = (n, seed)
                rec = rows.setdefault(key, {
                    'n': n,
                    'seed': seed,
                    'flipdist_a2b_ms': None,
                    'flipdist_b2a_ms': None,
                    'flipdist_a2b_dist': None,
                    'flipdist_b2a_dist': None,
                    'flipdist_a2b_status': None,
                    'flipdist_b2a_status': None,
                })
                if direction == 'a->b':
                    rec['flipdist_a2b_ms'] = float(obj.get('time_ms', float('nan')))
                    rec['flipdist_a2b_dist'] = int(obj.get('distance', -1))
                    rec['flipdist_a2b_status'] = str(obj.get('status', ''))
                elif direction == 'b->a':
                    rec['flipdist_b2a_ms'] = float(obj.get('time_ms', float('nan')))
                    rec['flipdist_b2a_dist'] = int(obj.get('distance', -1))
                    rec['flipdist_b2a_status'] = str(obj.get('status', ''))

    for rec in rows.values():
        a = rec.get('flipdist_a2b_ms')
        b = rec.get('flipdist_b2a_ms')
        if isinstance(a, float) and isinstance(b, float):
            rec['flipdist_max_ms'] = max(a, b)
            rec['flipdist_mean_ms'] = (a + b) / 2.0
        else:
            rec['flipdist_max_ms'] = None
            rec['flipdist_mean_ms'] = None
        sa = rec.get('flipdist_a2b_status')
        sb = rec.get('flipdist_b2a_status')
        rec['flipdist_ok_both'] = (sa == 'ok' and sb == 'ok')
    return rows


def load_astar_set_stats(astar_root: Path, n: int, seed: int, variant: str) -> Dict[str, Optional[float]]:
    p = astar_root / f'n_{n}' / f'set_{seed}' / f'results_{variant}.csv'
    if not p.exists():
        return {
            'exists': False,
            'pairs': None,
            'runtime_mean_ms': None,
            'runtime_median_ms': None,
            'runtime_p95_ms': None,
            'flip_mean': None,
            'flip_median': None,
            'timeout_rows': None,
        }

    runtimes: List[float] = []
    flips: List[float] = []
    timeout_rows = 0
    with p.open('r', encoding='utf-8') as fh:
        header_seen = False
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if not header_seen:
                header_seen = True
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                flip = float(parts[2])
                runtime = float(parts[4])
            except ValueError:
                continue
            flips.append(flip)
            runtimes.append(runtime)
            if int(flip) < 0:
                timeout_rows += 1

    if not runtimes:
        return {
            'exists': True,
            'pairs': 0,
            'runtime_mean_ms': None,
            'runtime_median_ms': None,
            'runtime_p95_ms': None,
            'flip_mean': None,
            'flip_median': None,
            'timeout_rows': 0,
        }

    return {
        'exists': True,
        'pairs': len(runtimes),
        'runtime_mean_ms': statistics.fmean(runtimes),
        'runtime_median_ms': statistics.median(runtimes),
        'runtime_p95_ms': quantile(runtimes, 0.95),
        'flip_mean': statistics.fmean(flips),
        'flip_median': statistics.median(flips),
        'timeout_rows': timeout_rows,
    }


def safe_ratio(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    if b == 0:
        return None
    return a / b


def main() -> int:
    ap = argparse.ArgumentParser(description='Compare flipdist timing vs bundled A* timing by n/seed')
    ap.add_argument('--flipdist-jsonl', nargs='+', required=True,
                    help='One or more flipdist JSONL output files')
    ap.add_argument('--astar-root', default='third_party/AStarFlipDistance/data/random_experiments_paper',
                    help='Path to A* random_experiments_paper root')
    ap.add_argument('--n-values', required=True, help='Comma-separated n list, e.g. 15,20,25')
    ap.add_argument('--seeds', required=True, help='Comma-separated seed list, e.g. 0,1,2,3,4')
    ap.add_argument('--output', required=True, help='Output CSV path')
    args = ap.parse_args()

    flip_paths = [Path(x) for x in args.flipdist_jsonl]
    n_values = parse_int_list(args.n_values)
    seeds = parse_int_list(args.seeds)
    astar_root = Path(args.astar_root)

    flip_rows = load_flipdist_jsonl(flip_paths)

    out_rows: List[Dict[str, object]] = []
    for n in n_values:
        for seed in seeds:
            rec = dict(flip_rows.get((n, seed), {'n': n, 'seed': seed}))

            a_simple = load_astar_set_stats(astar_root, n, seed, 'simple')
            a_combined = load_astar_set_stats(astar_root, n, seed, 'combined')

            rec['astar_simple_pairs'] = a_simple['pairs']
            rec['astar_simple_mean_ms'] = a_simple['runtime_mean_ms']
            rec['astar_simple_median_ms'] = a_simple['runtime_median_ms']
            rec['astar_simple_p95_ms'] = a_simple['runtime_p95_ms']
            rec['astar_simple_timeout_rows'] = a_simple['timeout_rows']

            rec['astar_combined_pairs'] = a_combined['pairs']
            rec['astar_combined_mean_ms'] = a_combined['runtime_mean_ms']
            rec['astar_combined_median_ms'] = a_combined['runtime_median_ms']
            rec['astar_combined_p95_ms'] = a_combined['runtime_p95_ms']
            rec['astar_combined_timeout_rows'] = a_combined['timeout_rows']

            rec['ratio_flipdist_max_vs_astar_simple_median'] = safe_ratio(
                rec.get('flipdist_max_ms'), a_simple['runtime_median_ms'])
            rec['ratio_flipdist_max_vs_astar_combined_median'] = safe_ratio(
                rec.get('flipdist_max_ms'), a_combined['runtime_median_ms'])

            out_rows.append(rec)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        'n', 'seed',
        'flipdist_a2b_ms', 'flipdist_b2a_ms', 'flipdist_max_ms', 'flipdist_mean_ms',
        'flipdist_a2b_dist', 'flipdist_b2a_dist', 'flipdist_a2b_status', 'flipdist_b2a_status', 'flipdist_ok_both',
        'astar_simple_pairs', 'astar_simple_mean_ms', 'astar_simple_median_ms', 'astar_simple_p95_ms', 'astar_simple_timeout_rows',
        'astar_combined_pairs', 'astar_combined_mean_ms', 'astar_combined_median_ms', 'astar_combined_p95_ms', 'astar_combined_timeout_rows',
        'ratio_flipdist_max_vs_astar_simple_median',
        'ratio_flipdist_max_vs_astar_combined_median',
    ]

    with out_path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in out_rows:
            w.writerow({k: row.get(k) for k in fieldnames})

    print(f'Wrote {out_path} ({len(out_rows)} rows)')

    by_n: Dict[int, List[Dict[str, object]]] = {}
    for row in out_rows:
        by_n.setdefault(int(row['n']), []).append(row)

    print('\nSummary by n (median over seeds):')
    for n in sorted(by_n.keys()):
        rows = by_n[n]
        fmax = [float(r['flipdist_max_ms']) for r in rows if r.get('flipdist_max_ms') is not None]
        asim = [float(r['astar_simple_median_ms']) for r in rows if r.get('astar_simple_median_ms') is not None]
        acomb = [float(r['astar_combined_median_ms']) for r in rows if r.get('astar_combined_median_ms') is not None]
        ok_both = sum(1 for r in rows if r.get('flipdist_ok_both'))
        print(
            f"n={n}: flipdist_med_max_ms={statistics.median(fmax) if fmax else 'NA'} "
            f"astar_simple_med_ms={statistics.median(asim) if asim else 'NA'} "
            f"astar_combined_med_ms={statistics.median(acomb) if acomb else 'NA'} "
            f"flipdist_ok_both={ok_both}/{len(rows)}"
        )

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
