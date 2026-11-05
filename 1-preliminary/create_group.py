#!/usr/bin/env python
"""Create the results group and check that this profile can run ld1.x.

Does not submit anything: it prints the plan derived from the recipe, so a typo is
caught before seven calculations are in the graph.

    python 1-preliminary/create_group.py --code ld1.x@localhost
"""

from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida import load_profile, orm          # noqa: E402

from aiida_ld1_scripts import recipe as recipe_module   # noqa: E402
from aiida_ld1_scripts import workflow                  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recipe', default=None)
    parser.add_argument('--group', default=recipe_module.DEFAULT_GROUP)
    parser.add_argument('--code', default='ld1.x@localhost')
    arguments = parser.parse_args()

    load_profile()
    recipe = recipe_module.load(arguments.recipe)
    code = workflow.get_code(arguments.code)
    group = workflow.get_group(arguments.group)

    scan = recipe.get('rlderiv_scan') or {}
    print(f'recipe      : {recipe["name"]}')
    print(f'  description: {recipe.get("description", "")}')
    print(f'  states     : {len(recipe["generation"]["states"])} '
          f'({", ".join(str(state[0]) for state in recipe["generation"]["states"])})')
    print(f'  file_pseudopw: {recipe["generation"]["inputp"]["file_pseudopw"]}')
    print(f'  transferability: {len(recipe["configts"]["test"]["configts"])} configurations')
    print(f'  rlderiv scan: {scan.get("values")}')
    print(f'code        : {code.label}@{code.computer.label} ({code.pk})')
    print(f'group       : {group.label} ({group.pk}, {len(group.nodes)} nodes)')
    print()
    print('planned calculations: 1 generation + '
          f'{len(scan.get("values") or [])} rlderiv points + 1 transferability = '
          f'{2 + len(scan.get("values") or [])}')
    print()
    print('next: python 2-submit/00_submit_all.py --wait')


if __name__ == '__main__':
    main()
