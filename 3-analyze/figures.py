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

from aiida_ld1_scripts import logderiv_scan                 # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', default=None)
    parser.add_argument('--out', default=None)
    parser.add_argument('--paper', default=None,
                        help="also write fig3_1_radial.* and fig3_3_logderiv.* here")
    parser.add_argument('--scan-only', action='store_true',
                        help='redraw S1 without AiiDA; default input: tests/data/rlderiv_scan')
    parser.add_argument('--scan-root', type=pathlib.Path,
                        help='original rlderiv_4p0 ... rlderiv_8p0 directories for S1')
    parser.add_argument('--jump-threshold', type=float, default=5.0,
                        help='display-only connector break; not an accuracy/pole threshold')
    arguments = parser.parse_args()
    if not (0 < arguments.jump_threshold < float('inf')):
        parser.error('--jump-threshold must be positive and finite.')

    out = pathlib.Path(arguments.out) if arguments.out else ROOT / 'data'
    if arguments.scan_only:
        scan_root = arguments.scan_root or ROOT / 'tests/data/rlderiv_scan'
        logderiv_scan.make_figure(scan_root, out, arguments.jump_threshold)
        print(f'wrote {out}/fig-S1-logderiv-radius-scan.{{png,pdf,svg}}')
        return

    from aiida import load_profile
    from aiida_ld1_scripts import collect as collect_module
    from aiida_ld1_scripts import report as report_module
    from aiida_ld1_scripts import recipe as recipe_module
    arguments.group = arguments.group or recipe_module.DEFAULT_GROUP
    load_profile()
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))

    # the figures need the ld1.x files themselves, i.e. the retrieval of the generation run
    node = collect_module.newest(arguments.group, 'generation')
    with tempfile.TemporaryDirectory(prefix='ld1-figures-') as temporary:
        directory = collect_module.dump_retrieved(node, pathlib.Path(temporary))
        for path in report_module.make_figures(summary, directory, out):
            print(f'wrote {path}')

        # Use the same finished scan nodes as the collected summary, never the
        # shipped fixtures when rendering a live calculation group.
        scan_root = arguments.scan_root
        if scan_root is None and summary.get('rlderiv_scan'):
            nodes = {node.pk: node for node in
                     collect_module._finished(arguments.group, 'rlderiv_scan')}
            scan_root = pathlib.Path(temporary) / 'scan'
            seen = set()
            for entry in summary['rlderiv_scan']:
                radius = float(entry['rlderiv'])
                if radius in seen:
                    raise ValueError(f'duplicate scan radius {radius}; recollect one run per radius')
                seen.add(radius)
                node = nodes[entry['pk']]
                actual = float(node.inputs.parameters.get_dict()['input']['rlderiv'])
                if actual != radius:
                    raise ValueError('scan node radius differs from summary; recollect the group')
                name = f"rlderiv_{radius:.1f}".replace('.', 'p')
                collect_module.dump_retrieved(node, scan_root / name)
            if seen != set(logderiv_scan.RADII):
                raise ValueError('S1 requires the reference 4,5,6,7,8 bohr scan; use --scan-root')
        if scan_root is not None:
            logderiv_scan.make_figure(scan_root, out, arguments.jump_threshold)
            print(f'wrote {out}/fig-S1-logderiv-radius-scan.{{png,pdf,svg}}')

        if arguments.paper:
            rcuts = {label: entry['rcut'] for label, entry
                     in (summary['generation'].get('wavefunctions') or {}).items()}
            valence = {str(state['label'])
                       for state in summary['generation'].get('states') or []} or None
            written = report_module.export_paper_figures(
                directory, arguments.paper, rcuts=rcuts, valence=valence,
                window=summary['criteria']['window'],
                threshold=summary['criteria']['pole_threshold'])
            for base in written:
                print(f'wrote {base}.{{png,svg,pdf}}')


if __name__ == '__main__':
    main()
