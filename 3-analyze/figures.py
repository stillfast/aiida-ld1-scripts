#!/usr/bin/env python
"""Draw the three figures into data/."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida import load_profile          # noqa: E402

from aiida_ld1_scripts import collect as collect_module   # noqa: E402
from aiida_ld1_scripts import report as report_module      # noqa: E402
from aiida_ld1_scripts import recipe as recipe_module      # noqa: E402
from aiida_ld1_scripts import workflow                     # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', default=recipe_module.DEFAULT_GROUP)
    parser.add_argument('--out', default=None)
    arguments = parser.parse_args()

    load_profile()
    out = pathlib.Path(arguments.out) if arguments.out else recipe_module.DATA_DIR
    data_dir = pathlib.Path(arguments.out) if arguments.out else recipe_module.DATA_DIR
    summary = json.loads((data_dir / 'summary.json').read_text(encoding='utf-8'))

    # the figures need the ld1.dlog / ld1.wfc files, i.e. the nodes themselves
    nodes = {'generation': collect_module.newest(arguments.group, 'generation')}
    for path in report_module.make_figures(summary, nodes, out):
        print(f'wrote {path}')


if __name__ == '__main__':
    main()
