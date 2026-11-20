"""S1 must preserve the real scan, including samples beside broken connectors."""
import csv
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from aiida_ld1_scripts import logderiv_scan as scan

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data/rlderiv_scan'


def test_entry_point_preserves_all_reference_values_without_aiida(tmp_path):
    # This CLI must not import AiiDA: CI installs only plotting/analysis dependencies.
    bootstrap = """
import importlib.abc, runpy, sys
class NoWorkflowImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('aiida', 'yaml'):
            raise RuntimeError('plot-only path imported ' + fullname)
sys.meta_path.insert(0, NoWorkflowImports())
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name='__main__')
"""
    result = subprocess.run(
        [sys.executable, '-c', bootstrap, str(ROOT / '3-analyze/figures.py'), '--scan-only',
         '--out', str(tmp_path)], capture_output=True, text=True, check=True)
    assert 'fig-S1-logderiv-radius-scan' in result.stdout
    metadata = json.loads((tmp_path / 'fig-S1-metadata.json').read_text())
    assert metadata['source_value_count'] == 12040
    assert metadata['excluded_samples'] == 0
    assert metadata['display_ylimit'] == [-10000, 10000]
    for item in metadata['input_files']:
        assert hashlib.sha256((DATA / item['relative_path']).read_bytes()).hexdigest() == item['sha256']
    with (tmp_path / 'fig-S1-source.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 6020
    for radius in scan.RADII:
        ae = scan.read_dlog(DATA / f'rlderiv_{radius}p0/ld1.dlog')
        ps = scan.read_dlog(DATA / f'rlderiv_{radius}p0/ld1ps.dlog')
        for channel in range(4):
            subset = [r for r in rows if int(r['radius_bohr']) == radius and int(r['l']) == channel]
            np.testing.assert_array_equal([float(r['energy_Ry']) for r in subset], ae[:, 0])
            np.testing.assert_array_equal([float(r['L_AE']) for r in subset], ae[:, channel + 1])
            np.testing.assert_array_equal([float(r['L_PS']) for r in subset], ps[:, channel + 1])
    for suffix in ('.png', '.pdf', '.svg', '-preview.png'):
        assert (tmp_path / f'fig-S1-logderiv-radius-scan{suffix}').stat().st_size > 1000
    alignment = json.loads((tmp_path / 'fig-S1-alignment.json').read_text())
    assert alignment['verdict'] == 'PASS'
    assert len(alignment['plot_rectangles_pt']) == 20


def test_display_breaks_preserve_both_endpoints():
    energy = np.arange(4.)
    values = np.array([10., -12., 1., -2.])
    mask = scan.break_mask(values, 5.)
    np.testing.assert_array_equal(mask, [False, True, False, False])
    x, y = scan.separated_line(energy, values, mask)
    np.testing.assert_array_equal(x[np.isfinite(x)], energy)
    np.testing.assert_array_equal(y[np.isfinite(y)], values)
    assert np.isnan(y[1]) and y[0] == 10 and y[2] == -12


def test_live_entry_uses_selected_group_scan_nodes(tmp_path, monkeypatch):
    """Exercise live dispatch with fake nodes; real AiiDA/QE is not needed in CI."""
    from aiida_ld1_scripts import report

    summary = {'rlderiv_scan': [{'pk': r, 'rlderiv': r} for r in scan.RADII]}
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    aiida = ModuleType('aiida')
    aiida.load_profile = lambda: None
    collect = ModuleType('aiida_ld1_scripts.collect')
    collect.newest = lambda group, role: SimpleNamespace(pk=99)
    groups = []

    def finished(group, role):
        groups.append((group, role))
        return [SimpleNamespace(pk=r, inputs=SimpleNamespace(parameters=SimpleNamespace(
            get_dict=lambda r=r: {'input': {'rlderiv': r}}))) for r in scan.RADII]

    def dump(node, target):
        target.mkdir(parents=True, exist_ok=True)
        (target / 'node-pk').write_text(str(node.pk))
        return target

    def render(data_root, output, threshold):
        assert output == tmp_path
        assert threshold == 5
        for radius in scan.RADII:
            assert (data_root / f'rlderiv_{radius}p0/node-pk').read_text() == str(radius)

    collect._finished, collect.dump_retrieved = finished, dump
    monkeypatch.setitem(sys.modules, 'aiida', aiida)
    monkeypatch.setitem(sys.modules, 'aiida_ld1_scripts.collect', collect)
    monkeypatch.setattr(report, 'make_figures', lambda *args: [])
    monkeypatch.setattr(scan, 'make_figure', render)
    monkeypatch.setattr(sys, 'argv', ['figures.py', '--group', 'current-group', '--out', str(tmp_path)])
    runpy.run_path(str(ROOT / '3-analyze/figures.py'), run_name='__main__')
    assert groups == [('current-group', 'rlderiv_scan')]


@pytest.mark.parametrize('contents', [
    '0 1 2\n1 2 3\n',  # missing channels
    '0 1 2 3 4\n0 2 3 4 5\n',  # repeated energy
    '0 1 2 3 4\n1 nan 3 4 5\n',
])
def test_invalid_raw_input_is_rejected(tmp_path, contents):
    path = tmp_path / 'ld1.dlog'
    path.write_text(contents)
    with pytest.raises(ValueError):
        scan.read_dlog(path)


@pytest.mark.parametrize('threshold', [0, -1, float('nan'), float('inf')])
def test_public_api_rejects_invalid_display_threshold(tmp_path, threshold):
    with pytest.raises(ValueError, match='positive and finite'):
        scan.make_figure(DATA, tmp_path, threshold)
