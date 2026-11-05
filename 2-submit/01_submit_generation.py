#!/usr/bin/env python
"""Submit the generation run (iswitch = 3) for the recipe's pseudopotential."""

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
    parser.add_argument('--wait', action='store_true',
                        help='run in process instead of submitting to the daemon')
    parser.add_argument('--reuse', action='store_true',
                        help='reuse a finished generation already in the group')
    arguments = parser.parse_args()

    load_profile()
    recipe = recipe_module.load(arguments.recipe)
    node = workflow.submit_generation(recipe, arguments.code, arguments.group,
                                      wait=arguments.wait, reuse=arguments.reuse)
    if node.outputs.get('pseudopotential'):
        print(f"  pseudopotential: pk {node.outputs['pseudopotential'].pk}")


if __name__ == '__main__':
    main()
