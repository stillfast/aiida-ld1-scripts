#!/usr/bin/env python
"""Write data/REPORT.md from data/summary.json."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from aiida_ld1_scripts import report as report_module      # noqa: E402
from aiida_ld1_scripts import recipe as recipe_module      # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default=None)
    arguments = parser.parse_args()

    out = pathlib.Path(arguments.out) if arguments.out else recipe_module.DATA_DIR
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    print(f'wrote {report_module.write_report(summary, out)}')


if __name__ == '__main__':
    main()
