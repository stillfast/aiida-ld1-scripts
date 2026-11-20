#!/usr/bin/env python3
"""Replot Figure S1 from the original ld1.x logarithmic-derivative files.

Run: python 3-analyze/figures.py --scan-only
Or:  python -m aiida_ld1_scripts.logderiv_scan --data-root tests/data/rlderiv_scan
Requires only NumPy and Matplotlib; no AiiDA or Quantum ESPRESSO installation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import platform

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
RADII = (4, 5, 6, 7, 8)
CHANNELS = ('s', 'p', 'd', 'f')
AE_STYLE = dict(color='#303030', lw=0.95, ls='-', zorder=2)
PS_STYLE = dict(color='#D55E00', lw=0.65, ls=(0, (3, 2)), zorder=3)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_dlog(path: Path) -> np.ndarray:
    """Columns: energy (Ry), l=0, l=1, l=2, l=3. Keep the QE values unchanged."""
    if not path.is_file():
        raise FileNotFoundError(f'Missing original data: {path}')
    content = path.read_text(encoding='utf-8').replace('D', 'E').replace('d', 'e')
    data = np.loadtxt(io.StringIO(content), comments='#', ndmin=2)
    if data.shape[1] != 5:
        raise ValueError(f'{path}: expected 5 columns (E, s, p, d, f), found {data.shape[1]}')
    if len(data) < 2 or not np.isfinite(data).all():
        raise ValueError(f'{path}: at least two finite rows are required.')
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f'{path}: energy must increase strictly; no sorting is applied.')
    return data


def break_mask(values: np.ndarray, threshold: float) -> np.ndarray:
    """Break only a displayed connector at a large opposite-sign jump.

    This heuristic is a plotting parameter, not an accuracy threshold or a pole
    locator. Both endpoint samples remain drawn and exported. No value is removed.
    """
    before = np.zeros(len(values), dtype=bool)
    before[1:] = ((np.signbit(values[:-1]) != np.signbit(values[1:]))
                  & (np.abs(values[:-1]) >= threshold)
                  & (np.abs(values[1:]) >= threshold))
    return before


def separated_line(energy: np.ndarray, values: np.ndarray, breaks: np.ndarray) -> tuple:
    """Insert NaN between original samples; do not interpolate or discard them."""
    indices = np.flatnonzero(breaks)
    return np.insert(energy, indices, np.nan), np.insert(values, indices, np.nan)


def tick_label(value: float, _position: float) -> str:
    if value == 0:
        return '0'
    exponent = int(round(np.log10(abs(value))))
    return '$'+('-' if value < 0 else '')+r'10^{'+str(exponent)+'}$'


def local_alignment_check(fig, json_out, tolerance_pt=1.5, **_kwargs):
    """Self-contained check of this fixed 5x4 rectangular layout in PDF points.

    Checks equal plot widths/heights, row/column edges and repeated gutters.
    The optional external audit adds further typography/label diagnostics.
    """
    fig.canvas.draw()
    width, height = fig.get_size_inches()*72
    boxes = np.array([[a.get_position().x0*width, a.get_position().y0*height,
                       a.get_position().x1*width, a.get_position().y1*height] for a in fig.axes])
    grid = boxes.reshape(len(RADII), len(CHANNELS), 4)
    deviations = {
        'width': float(np.ptp(boxes[:, 2]-boxes[:, 0])),
        'height': float(np.ptp(boxes[:, 3]-boxes[:, 1])),
        'row_bottom': float(np.max(np.ptp(grid[:, :, 1], axis=1))),
        'row_top': float(np.max(np.ptp(grid[:, :, 3], axis=1))),
        'column_left': float(np.max(np.ptp(grid[:, :, 0], axis=0))),
        'column_right': float(np.max(np.ptp(grid[:, :, 2], axis=0))),
        'horizontal_gutter': float(np.ptp(grid[:, 1:, 0]-grid[:, :-1, 2])),
        'vertical_gutter': float(np.ptp(grid[:-1, :, 1]-grid[1:, :, 3])),
    }
    report = {'verdict': 'PASS' if max(deviations.values()) <= tolerance_pt else 'FAIL',
              'check': 'uniform 5x4 plot rectangles', 'tolerance_pt': tolerance_pt,
              'max_deviations_pt': deviations, 'plot_rectangles_pt': boxes.tolist()}
    Path(json_out).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    if report['verdict'] != 'PASS':
        raise RuntimeError('Panel rectangles failed alignment: '+str(deviations))
    return report


try:
    # Optional, for authoring QA; the complete plotting script works without it.
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    require_matplotlib_panel_alignment = local_alignment_check


def make_figure(data_root: Path, output: Path, jump_threshold: float = 5.0) -> dict:
    """Export the reference five-radius S1 grid; retain every raw QE sample.

    ``data_root`` contains rlderiv_4p0 through rlderiv_8p0, each with an AE/PS pair.
    The display-only jump threshold is independent of analysis.py's pole settings.
    """
    data_root, output = Path(data_root), Path(output)
    if not np.isfinite(jump_threshold) or jump_threshold <= 0:
        raise ValueError('jump_threshold must be positive and finite.')
    with plt.rc_context({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                         'font.size': 8, 'axes.linewidth': 0.6,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8.5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'}):
        return _make_figure(data_root, output, jump_threshold)


def _make_figure(data_root: Path, output: Path, jump_threshold: float) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    datasets, input_files = {}, []
    for radius in RADII:
        folder = data_root/f'rlderiv_{radius}p0'
        ae, ps = [read_dlog(folder/name) for name in ('ld1.dlog', 'ld1ps.dlog')]
        if ae.shape != ps.shape or not np.allclose(ae[:, 0], ps[:, 0], atol=1e-12, rtol=0):
            raise ValueError(f'r={radius}: AE and PS require identical energy grids.')
        if ae[0, 0] != -3 or ae[-1, 0] != 3:
            raise ValueError(f'r={radius}: this Figure S1 expects the full -3 to 3 Ry window.')
        datasets[radius] = ae, ps
        for name in ('ld1.dlog', 'ld1ps.dlog'):
            path = folder/name
            input_files.append({'relative_path': f'rlderiv_{radius}p0/{name}',
                                'sha256': sha256(path), 'rows': len(ae)})
    max_value = max(float(np.abs(a[:, 1:]).max()) for pair in datasets.values() for a in pair)
    ylimit = 10.0**np.ceil(np.log10(max(max_value, 10.0)))
    fig, axes = plt.subplots(5, 4, figsize=(178/25.4, 210/25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.125, right=0.985, bottom=0.075, top=0.905, wspace=0.16, hspace=0.40)
    breaks_summary = []
    csv_path = output/'fig-S1-source.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['radius_bohr', 'l', 'energy_Ry', 'L_AE', 'L_PS',
                         'break_connector_before_AE', 'break_connector_before_PS'])
        for row, radius in enumerate(RADII):
            ae, ps = datasets[radius]
            energy = ae[:, 0]
            for col, channel in enumerate(CHANNELS):
                ax = axes[row, col]
                va, vp = ae[:, col+1], ps[:, col+1]
                ba, bp = break_mask(va, jump_threshold), break_mask(vp, jump_threshold)
                xa, ya = separated_line(energy, va, ba)
                xp, yp = separated_line(energy, vp, bp)
                ax.plot(xa, ya, **AE_STYLE)
                ax.plot(xp, yp, **PS_STYLE)
                # Show ALL source samples, including both sides of every broken connector.
                ax.plot(energy, vp, ls='none', marker='o', ms=1.5, markerfacecolor='none',
                        markeredgecolor=PS_STYLE['color'], markeredgewidth=0.55, zorder=4)
                ax.plot(energy, va, ls='none', marker='.', ms=0.9, color=AE_STYLE['color'], zorder=2)
                ax.set_yscale('symlog', linthresh=1.0, linscale=0.6)
                ax.set(xlim=(-3, 3), ylim=(-ylimit, ylimit))
                ax.set_xticks([-3, 0, 3])
                ax.set_yticks([-ylimit, -ylimit/100, 0, ylimit/100, ylimit])
                ax.yaxis.set_major_formatter(FuncFormatter(tick_label))
                ax.minorticks_off()
                ax.tick_params(direction='in', length=2.5, width=0.6, pad=2,
                               labelleft=col == 0, labelbottom=row == 4)
                ax.tick_params(axis='x', pad=6)
                ax.axhline(0, color='#D5D5D5', lw=0.4, zorder=0)
                letter = chr(ord('a')+row*4+col)
                ax.text(0.0, 1.15, letter, transform=ax.transAxes, fontweight='bold',
                        fontsize=8, ha='left', va='bottom', clip_on=False)
                if col == 0:
                    ax.text(0.16, 1.15, f'$r_d={radius}$ bohr', transform=ax.transAxes,
                            fontsize=8, ha='left', va='bottom', clip_on=False)
                breaks_summary.append({'radius_bohr': radius, 'l': col,
                                       'AE_break_indices': np.flatnonzero(ba).tolist(),
                                       'PS_break_indices': np.flatnonzero(bp).tolist()})
                for index, e in enumerate(energy):
                    writer.writerow([radius, col, e, va[index], vp[index], int(ba[index]), int(bp[index])])
    for col, channel in enumerate(CHANNELS):
        box = axes[0, col].get_position()
        fig.text((box.x0+box.x1)/2, 0.955, f'${channel}$  ($l={col}$)',
                 fontsize=9, ha='center', va='center')
    handles = [Line2D([], [], color=AE_STYLE['color'], lw=1.0, label='AE'),
               Line2D([], [], color=PS_STYLE['color'], lw=0.7, ls='--', marker='o',
                      ms=3, markerfacecolor='none', markeredgewidth=0.7, label='PS')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.57, 0.998),
               ncol=2, frameon=False, handlelength=3, columnspacing=2)
    fig.text(0.55, 0.025, 'Energy (Ry)', ha='center', fontsize=9)
    fig.text(0.022, 0.50, r'$L_l(E,r_d)$ (QE output)', rotation=90,
             rotation_mode='anchor', ha='center', va='center', fontsize=9)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=output/'fig-S1-alignment.json',
                                      tolerance_pt=1.5, gutter_tolerance_pt=1.5,
                                      require_panel_labels=True, strict=True)
    stem = output/'fig-S1-logderiv-radius-scan'
    fig.savefig(str(stem)+'.png', dpi=600, facecolor='white')
    fig.savefig(str(stem)+'.pdf', facecolor='white')
    fig.savefig(str(stem)+'.svg', facecolor='white')
    fig.savefig(str(stem)+'-preview.png', dpi=300, facecolor='white')
    plt.close(fig)
    report = {'repository': 'https://github.com/stillfast/aiida-ld1-scripts', 'schema_version': 1,
              'input_files': input_files, 'radii_bohr': list(RADII), 'angular_momenta': [0, 1, 2, 3],
              'source_value_count': sum(len(a)*4+len(b)*4 for a, b in datasets.values()),
              'energy_window_Ry': [-3, 3], 'symlog': {'linthresh': 1.0, 'linscale': 0.6},
              'display_ylimit': [-ylimit, ylimit], 'jump_connector_threshold': jump_threshold,
              'broken_connectors': breaks_summary, 'excluded_samples': 0,
              'processing': 'unchanged QE columns; no smoothing, fitting, resampling or ranking; heuristic line breaks retain all points',
              'source_csv_sha256': sha256(csv_path), 'python': platform.python_version(),
              'numpy': np.__version__, 'matplotlib': matplotlib.__version__}
    (output/'fig-S1-metadata.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=HERE.parents[1]/'tests/data/rlderiv_scan')
    parser.add_argument('--output', type=Path, default=HERE.parents[1]/'data')
    parser.add_argument('--jump-threshold', type=float, default=5.0,
                        help='Display-only: break opposite-sign connectors with both |L| >= this value.')
    args = parser.parse_args()
    if not np.isfinite(args.jump_threshold) or args.jump_threshold <= 0:
        parser.error('--jump-threshold must be positive and finite.')
    report = make_figure(args.data_root, args.output, args.jump_threshold)
    print(json.dumps({'source_values': report['source_value_count'], 'excluded_samples': 0,
                      'output': str(args.output.resolve())}, indent=2))


if __name__ == '__main__':
    main()
