#!/usr/bin/env python
"""Run the whole pipeline: generation -> rlderiv scan -> transferability.

With ``--wait`` everything runs in process (no daemon needed, which is what makes this
usable on a laptop); without it the calculations are submitted and the script returns.
Analysis afterwards is ``3-analyze/``.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida import load_profile          # noqa: E402

from aiida_ld1_scripts import recipe as recipe_module   # noqa: E402
from aiida_ld1_scripts import workflow                  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recipe', default=None)
    parser.add_argument('--group', default=recipe_module.DEFAULT_GROUP)
    parser.add_argument('--code', default='ld1.x@localhost')
    parser.add_argument('--wait', action='store_true')
    parser.add_argument('--reuse', action='store_true',
                        help='skip anything already finished in the group')
    arguments = parser.parse_args()

    load_profile()
    recipe = recipe_module.load(arguments.recipe)

    generation = workflow.submit_generation(recipe, arguments.code, arguments.group,
                                            wait=arguments.wait, reuse=arguments.reuse)
    workflow.submit_rlderiv_scan(recipe, arguments.code, arguments.group,
                                 wait=arguments.wait, reuse=arguments.reuse)
    workflow.submit_transferability(recipe, arguments.code, arguments.group,
                                    wait=arguments.wait, reuse=arguments.reuse,
                                    generation_node=generation)
    print()
    print(f'all submitted into group {arguments.group!r}')
    print('next: python 3-analyze/collect.py && python 3-analyze/figures.py '
          '&& python 3-analyze/report.py')


if __name__ == '__main__':
    main()
