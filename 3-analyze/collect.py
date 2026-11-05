#!/usr/bin/env python
"""Measure every calculation in the group and write data/summary.json + the CSVs."""

from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida import load_profile          # noqa: E402

from aiida_ld1_scripts import collect as collect_module   # noqa: E402
from aiida_ld1_scripts import recipe as recipe_module      # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recipe', default=None)
    parser.add_argument('--group', default=recipe_module.DEFAULT_GROUP)
    parser.add_argument('--out', default=None)
    arguments = parser.parse_args()

    load_profile()
    summary = collect_module.collect(arguments.group, arguments.recipe)
    out = collect_module.write(summary, arguments.out)

    generation = summary['generation']
    print(f"generation pk {generation['pk']}: "
          f"{generation['number_of_states']} states, Etotps = {generation['etotps']} Ry")
    derivatives = generation.get('log_derivatives') or {}
    if derivatives:
        print(f"  poles coincide: {derivatives.get('poles_all_coincide')}, "
              f"pooled median |dL| = {derivatives.get('pooled_median'):.3e}")
    transferability = summary.get('transferability') or {}
    if transferability.get('rows'):
        print(f"transferability: worst |Delta E| = "
              f"{transferability['max_abs_mry']:.1f} mRy")
    scan = summary.get('rlderiv_scan') or []
    if scan:
        best = summary['best_rlderiv']
        print(f"rlderiv scan: {len(scan)} points, best = {best['rlderiv']:g} "
              f"(pooled median {best['pooled_median']:.3e})")
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
