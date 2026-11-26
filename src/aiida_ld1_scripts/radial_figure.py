"""Replot AE/PS reference wavefunctions without an AiiDA profile or QE runtime."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import zipfile

from . import analysis, recipe, report


def extract_generation(archive: Path, directory: Path) -> dict:
    """Read the unique generation result from a SQLite-format AiiDA export.

    Select by the recipe's generation label and retrieved link, not a local PK.
    Multiple cached generation runs are allowed only when their wavefunction hashes
    are identical. Never select the first of several scientifically different runs.
    """
    with zipfile.ZipFile(archive) as source:
        metadata = json.loads(source.read('metadata.json'))
        if metadata.get('key_format') != 'sha256':
            raise ValueError('Expected a SHA256, SQLite-format AiiDA archive.')
        database = directory / 'archive.sqlite3'
        database.write_bytes(source.read('db.sqlite3'))
        expected = recipe.load()['name'] + ' generation'
        with sqlite3.connect(f'{database.as_uri()}?mode=ro', uri=True) as connection:
            rows = connection.execute(
                "SELECT p.uuid, r.repository_metadata "
                "FROM db_dbnode p JOIN db_dblink l ON l.input_id=p.id "
                "JOIN db_dbnode r ON r.id=l.output_id "
                "WHERE l.label='retrieved' AND p.label=? ORDER BY p.uuid", (expected,)
            ).fetchall()
        if not rows:
            raise ValueError(f'No retrieved generation result with label {expected!r}.')
        sources = [(uuid, json.loads(meta)['o']) for uuid, meta in rows]
        signatures = {tuple(objects[name]['k'] for name in ('ld1.wfc', 'ld1ps.wfc'))
                      for _, objects in sources}
        if len(signatures) != 1:
            raise ValueError('Generation runs have different wavefunctions; select an unambiguous archive.')
        hashes = {}
        for name in ('ld1.wfc', 'ld1ps.wfc'):
            key = sources[0][1][name]['k']
            contents = source.read('repo/' + key)
            if hashlib.sha256(contents).hexdigest() != key:
                raise ValueError(f'Archive object checksum failed for {name}.')
            (directory / name).write_bytes(contents)
            hashes[name] = key
    return {'generation_uuids': [uuid for uuid, _ in sources], 'input_sha256': hashes}


def make_figure(archive: str | Path, out: str | Path) -> dict:
    """Export raw signed u(r), with the existing phase convention and all line samples.

    Only markers are spaced (34 target positions), never the continuous curve.
    No smoothing, interpolation, amplitude rescaling, or uncertainty bars are used.
    """
    archive, out = Path(archive).resolve(), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    # The per-orbital cutoffs are collected from the same generation record.
    summary_path = archive.parent / 'summary.json'
    summary = json.loads(summary_path.read_text(encoding='utf8'))
    rcuts = {label: values['rcut'] for label, values in
             summary['generation']['wavefunctions'].items()}
    plt = report._matplotlib()
    plt.rcParams['savefig.dpi'] = 600
    with tempfile.TemporaryDirectory(prefix='ld1-radial-') as name:
        directory = Path(name)
        provenance = extract_generation(archive, directory)
        labels, grid, ae = analysis.read_wfc(directory / 'ld1.wfc')
        ps_labels, ps_grid, ps, _ = analysis.read_wfc_halves(directory / 'ld1ps.wfc')
        import numpy as np
        if not np.array_equal(grid, ps_grid):
            raise ValueError('AE and PS radial grids must match.')
        report.figure_radial(plt, directory, out / 'fig_radial_logx',
                             rcuts=rcuts, xscale='log', xlim=(1e-2, 1e2),
                             pseudo_line=True, formats=('png', 'pdf', 'svg'))
        report.figure_radial(plt, directory, out / 'fig_radial',
                             rcuts=rcuts, pseudo_line=True, formats=('png',))
        provenance.update({
            'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
            'summary_sha256': hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            'orbitals': ps_labels, 'grid_points_per_curve': len(grid),
            'display': {'xscale': 'log', 'xlim_bohr': [0.01, 100],
                        'ordinate': 'signed u(r) = r R(r), original amplitudes',
                        'AE': 'solid', 'PS': 'dashed + hollow circles',
                        'cutoff': 'vertical dotted',
                        'marker_target_count': 34, 'continuous_line': 'all visible samples',
                        'dpi': 600},
            'uncertainty': 'Not applicable: deterministic reference wavefunctions.',
        })
    (out / 'fig_radial_logx-metadata.json').write_text(
        json.dumps(provenance, indent=2) + '\n', encoding='utf8')
    return provenance
