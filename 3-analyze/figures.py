#!/usr/bin/env python
"""Draw the figures into data/, and optionally into the manuscript's figs/ directory.

    python 3-analyze/figures.py
    python 3-analyze/figures.py --paper ../../uranium-ppgen/ai_v1/6p25/figs
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida import load_profile          # noqa: E402

from aiida_ld1_scripts import collect as collect_module   # noqa: E402
from aiida_ld1_scripts import report as report_module      # noqa: E402
from aiida_ld1_scripts import recipe as recipe_module      # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', default=recipe_module.DEFAULT_GROUP)
    parser.add_argument('--out', default=None)
    parser.add_argument('--paper', default=None,
                        help="also write fig3_1_radial.* and fig3_3_logderiv.* here")
    arguments = parser.parse_args()

    load_profile()
    out = pathlib.Path(arguments.out) if arguments.out else recipe_module.DATA_DIR
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))

    # the figures need the ld1.x files themselves, i.e. the retrieval of the generation run
    node = collect_module.newest(arguments.group, 'generation')
    with tempfile.TemporaryDirectory(prefix='ld1-figures-') as temporary:
        directory = collect_module.dump_retrieved(node, pathlib.Path(temporary))
        for path in report_module.make_figures(summary, directory, out):
            print(f'wrote {path}')

        if arguments.paper:
            rcuts = {label: entry['rcut'] for label, entry
                     in (summary['generation'].get('wavefunctions') or {}).items()}
            written = report_module.export_paper_figures(
                directory, arguments.paper, rcuts=rcuts,
                window=summary['criteria']['window'],
                threshold=summary['criteria']['pole_threshold'])
            for base in written:
                print(f'wrote {base}.{{png,svg,pdf}}')


if __name__ == '__main__':
    main()
