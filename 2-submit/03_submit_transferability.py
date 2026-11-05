#!/usr/bin/env python
"""Submit the iswitch = 2 transferability test of the generated pseudopotential.

The ``pseudo_in`` of the run is the generation's ``pseudopotential`` **output node**, so
the test is linked in the provenance graph to the exact file it validates.
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
    parser.add_argument('--reuse', action='store_true')
    arguments = parser.parse_args()

    load_profile()
    recipe = recipe_module.load(arguments.recipe)
    workflow.submit_transferability(recipe, arguments.code, arguments.group,
                                    wait=arguments.wait, reuse=arguments.reuse)


if __name__ == '__main__':
    main()
