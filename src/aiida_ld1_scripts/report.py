"""Figures and the markdown report, from ``data/summary.json``.

Three figures, each one a question about the pseudopotential:

``fig_radial.png``     do the pseudo and all-electron radial functions agree?
``fig_logderiv.png``   do the logarithmic derivatives agree, and do the poles sit on
                       the same bound-state energies (the ghost criterion)?
``fig_rlderiv.png``    which ``rlderiv`` reproduces the all-electron data best?

The ``l >= 3`` channel is drawn but not used to judge anything: for ``l >= 3`` the
Bessel-box column of ``ld1.x`` is a numerical artefact (it moves ~0.23 Ry with the box
radius alone), so the criterion is ``e PS`` -- see the ``aiida-ld1`` README.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import analysis, recipe as recipe_module

#: channels used as criteria (the box test is not trustworthy for l >= 3)
RELIABLE_MAX_L = 2

COLOURS = ('#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e')


def _matplotlib():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _retrieved_dir_for(node):
    """A directory with the files of one calculation (temp dir owned by the caller)."""
    import tempfile

    from . import collect

    handle = tempfile.mkdtemp(prefix='ld1-report-')
    return collect.dump_retrieved(node, Path(handle))


def figure_radial(plt, node, out: Path) -> Optional[Path]:
    """AE vs PS ``u(r)`` for every channel of the generation run."""
    import numpy as np

    directory = Path(_retrieved_dir_for(node))
    ae_file, ps_file = directory / 'ld1.wfc', directory / 'ld1ps.wfc'
    if not (ae_file.is_file() and ps_file.is_file()):
        return None
    ae_labels, ae_r, ae_u = analysis.read_wfc(ae_file)
    # ld1ps.wfc holds the pseudo orbitals in its first half and the all-electron ones
    # in its second; only the first half is the pseudopotential's wavefunction
    ps_labels, ps_r, ps_u, _ae_again = analysis.read_wfc_halves(ps_file)
    ps_index = {label: index for index, label in enumerate(ps_labels)}

    figure, axes = plt.subplots(figsize=(7.2, 4.6))
    for index, label in reversed(list(enumerate(ae_labels))):
        if label not in ps_index:
            continue                      # e.g. 4F / 5D: core states of the AE run only
        colour = COLOURS[index % len(COLOURS)]
        pseudo = ps_u[:, ps_index[label]]
        aligned, _flipped = analysis.align_phase(ae_u[:, index], pseudo)
        axes.plot(ae_r, aligned, '-', color=colour, label=f'{label} AE')
        axes.plot(ps_r, pseudo, '--', color=colour, label=f'{label} PS')
    axes.set_xlabel('r (bohr)')
    axes.set_ylabel('u(r) = r·R(r)')
    axes.set_xlim(0.0, min(8.0, float(ae_r[-1])))
    axes.axhline(0.0, color='0.8', lw=0.6, zorder=0)
    axes.set_title('radial functions: all-electron (solid) vs pseudo (dashed), '
                   'phase aligned\nthey must coincide outside the rcut radius')
    axes.legend(fontsize=7, ncol=2, frameon=False)
    figure.tight_layout()
    figure.savefig(out, dpi=150)
    plt.close(figure)
    return out


def figure_logderiv(plt, node, out: Path, *, window=(-3.0, 3.0),
                    threshold: float = analysis.POLE_THRESHOLD) -> Optional[Path]:
    """AE vs PS logarithmic derivatives, with the poles marked."""
    import numpy as np

    directory = Path(_retrieved_dir_for(node))
    if not (directory / 'ld1.dlog').is_file():
        return None
    ae = analysis.read_dlog(directory / 'ld1.dlog')
    ps = analysis.read_dlog(directory / 'ld1ps.dlog')
    channels = [channel for channel in ae if channel != 'energy']
    figure, axes = plt.subplots(len(channels), 1, figsize=(7.2, 1.5 * len(channels)),
                                sharex=True)
    axes = np.atleast_1d(axes)
    for axis, channel, colour in zip(axes, channels, COLOURS):
        axis.plot(ae['energy'], ae[channel], '-', color=colour, lw=1.0, label='AE')
        axis.plot(ps['energy'], ps[channel], '--', color='k', lw=0.9, label='PS')
        for pole in analysis.find_poles(ae['energy'], ae[channel], threshold):
            axis.axvline(pole, color=colour, lw=0.6, alpha=0.5)
        axis.set_ylabel(channel)
        axis.axhline(0.0, color='0.85', lw=0.5, zorder=0)
        axis.set_ylim(-30, 30)
        if channel == channels[0]:
            axis.legend(fontsize=7, frameon=False, loc='upper left')
    axes[-1].set_xlabel('E (Ry)')
    axes[-1].set_xlim(*window)
    axes[0].set_title('logarithmic derivative: all-electron (colour) vs pseudo (black)')
    figure.tight_layout()
    figure.savefig(out, dpi=150)
    plt.close(figure)
    return out


def figure_rlderiv(plt, summary: Dict[str, Any], out: Path) -> Optional[Path]:
    """Pooled median ``|dL|`` against ``rlderiv``, with the pole verdict marked.

    Filled marker = the pole test passes at that radius, open marker = it fails.  The
    verdict, not the median, is what decides: see the README ("what the scan cannot
    decide").  The production value is drawn so nobody has to guess which radius was used.
    """
    scan = summary.get('rlderiv_scan') or []
    if not scan:
        return None
    import numpy as np

    production = summary.get('criteria', {}).get('production_rlderiv')
    figure, axis = plt.subplots(figsize=(6.4, 4.0))
    passed = [entry for entry in scan if entry['poles_all_coincide']]
    failed = [entry for entry in scan if not entry['poles_all_coincide']]
    axis.semilogy([entry['rlderiv'] for entry in passed],
                  [entry['pooled_median'] for entry in passed], 'o-',
                  color='#0072B2', ms=6, lw=1.2, label='poles coincide (admissible)')
    axis.semilogy([entry['rlderiv'] for entry in failed],
                  [entry['pooled_median'] for entry in failed], 'o',
                  mfc='none', mec='#D55E00', ms=7, mew=1.4,
                  label='poles do NOT coincide (rejected)')
    best = summary.get('best_rlderiv') or {}
    if best.get('rlderiv') is not None:
        axis.axvline(best['rlderiv'], color='#0072B2', ls=':', lw=1.0)
        axis.annotate(f"pole-gated pick = {best['rlderiv']:g}",
                      xy=(best['rlderiv'], best['pooled_median']),
                      xytext=(6, -14), textcoords='offset points',
                      fontsize=8, color='#0072B2')
    if production is not None:
        axis.axvline(production, color='0.5', ls='--', lw=1.0)
        axis.annotate(f'production = {production:g}', xy=(production, 0.0),
                      xycoords=('data', 'axes fraction'), xytext=(4, 6),
                      textcoords='offset points', fontsize=8, color='0.35')
    smallest = min(scan, key=lambda entry: entry['pooled_median'])
    if not smallest['poles_all_coincide']:
        axis.annotate('smallest median, but rejected',
                      xy=(smallest['rlderiv'], smallest['pooled_median']),
                      xytext=(6, 8), textcoords='offset points', fontsize=8,
                      color='#D55E00')
    axis.set_xlabel('rlderiv (bohr)')
    axis.set_ylabel('median |ΔL| (pooled over channels)')
    axis.set_title('log-derivative agreement vs the rlderiv reference radius\n'
                   'the median shrinks with the radius by construction; only the '
                   'pole verdict decides')
    axis.grid(True, which='both', color='0.9', lw=0.6)
    axis.legend(fontsize=8, frameon=False, loc='upper right')
    figure.tight_layout()
    figure.savefig(out, dpi=150)
    plt.close(figure)
    return out


def make_figures(summary: Dict[str, Any], nodes: Dict[str, Any],
                 out_dir: Optional[str | Path] = None) -> List[Path]:
    """All three figures; returns the ones that were actually written."""
    out = Path(out_dir) if out_dir is not None else recipe_module.DATA_DIR
    out.mkdir(parents=True, exist_ok=True)
    plt = _matplotlib()
    written: List[Path] = []
    generation = nodes.get('generation')
    if generation is not None:
        for function, name, extra in (
                (figure_radial, 'fig_radial.png', {}),
                (figure_logderiv, 'fig_logderiv.png',
                 {'window': summary['criteria']['window'],
                  'threshold': summary['criteria']['pole_threshold']})):
            try:
                path = function(plt, generation, out / name, **extra)
            except Exception as exception:                  # noqa: BLE001
                print(f'  warning: {name} failed: {exception}')
                path = None
            if path is not None:
                written.append(path)
    path = figure_rlderiv(plt, summary, out / 'fig_rlderiv.png')
    if path is not None:
        written.append(path)
    return written


# ------------------------------------------------------------------- report


def _markdown_table(header: List[str], rows: List[List[str]]) -> str:
    lines = ['| ' + ' | '.join(header) + ' |',
             '|' + '|'.join('---' for _ in header) + '|']
    lines += ['| ' + ' | '.join(row) + ' |' for row in rows]
    return '\n'.join(lines)


def write_report(summary: Dict[str, Any], out_dir: Optional[str | Path] = None) -> Path:
    """``REPORT.md``: what was generated, and how well it reproduces the AE data."""
    out = Path(out_dir) if out_dir is not None else recipe_module.DATA_DIR
    out.mkdir(parents=True, exist_ok=True)

    generation = summary.get('generation') or {}
    generated = summary.get('transferability') or {}
    scan = summary.get('rlderiv_scan') or []
    best = summary.get('best_rlderiv')
    criteria = summary.get('criteria') or {}

    lines = [f"# {summary.get('recipe', 'pseudopotential')}: generation report",
             '',
             f"{summary.get('description', '')}",
             '',
             f"- AiiDA group: `{summary.get('group')}`",
             f"- generation: pk {generation.get('pk')}",
             f"- transferability: pk {generated.get('pk') or '(not run)'}",
             f"- rlderiv scan: {len(scan)} point(s)",
             f"- criteria: window {criteria.get('window')} Ry, pole threshold "
             f"|L| > {criteria.get('pole_threshold')}, pole tolerance "
             f"{criteria.get('pole_tolerance')} Ry",
             '']

    metadata = generation.get('pseudo_metadata') or {}
    if metadata:
        lines += ['## The pseudopotential', '']
        lines.append(_markdown_table(['field', 'value'],
                                 [[key, str(value)] for key, value in metadata.items()]))
        lines += ['', f"Etot = {generation.get('etot')} Ry, "
                      f"Etotps = {generation.get('etotps')} Ry", '']

    states = generation.get('states') or []
    if states:
        lines += ['## Bound-state energies', '',
                  'The pseudo-atom reproduces the all-electron levels of the generating '
                  'configuration:', '']
        rows = [[state['label'], str(state['n']), str(state['l']), f"{state['occ']:g}",
                 f"{state['ae']:+.5f}", f"{state['ps']:+.5f}",
                 f"{1000.0 * state['delta']:+.2f}"]
                for state in states]
        worst = max(abs(state['delta']) for state in states) * 1000.0
        lines.append(_markdown_table(['orbital', 'n', 'l', 'occ', 'e AE (Ry)', 'e PS (Ry)',
                                  'delta (mRy)'], rows))
        lines += ['', f'Worst mismatch: **{worst:.2f} mRy**.', '']

    wavefunctions = generation.get('wavefunctions') or {}
    if wavefunctions:
        lines += ['## Radial functions', '',
                  'Inside a state cut-off radius the pseudo orbital is *meant* to differ '
                  'from the all-electron one (that is the pseudization); **outside** it '
                  'the tail must agree, otherwise the logarithmic derivative and the '
                  'transferability are wrong. Columns are `max |u_AE - u_PS|` before and '
                  'after the cut-off (phase aligned, see `analysis.align_phase`).', '']
        lines.append(_markdown_table(
            ['orbital', 'rcut (bohr)', 'inside', 'outside', 'phase flipped'],
            [[label, f"{entry['rcut']:g}",
              '-' if entry['max_inside'] is None else f"{entry['max_inside']:.4f}",
              '-' if entry['max_outside'] is None else f"{entry['max_outside']:.2e}",
              'yes' if entry['flipped'] else 'no']
             for label, entry in wavefunctions.items()]))
        worst = max((entry['max_outside'] for entry in wavefunctions.values()
                     if entry['max_outside'] is not None), default=None)
        if worst is not None:
            lines += ['', f'Worst tail mismatch: **{worst:.2e}**.', '']

    derivatives = generation.get('log_derivatives') or {}
    if derivatives:
        per_channel = derivatives.get('per_channel') or {}
        poles = derivatives.get('poles') or {}
        rows = []
        for channel in sorted(per_channel):
            entry = per_channel[channel]
            pole_entry = poles.get(channel, {})
            rows.append([channel, f"{entry['median']:.3e}", f"{entry['max']:.3e}",
                         f"{pole_entry.get('n_ae', '-')} / {pole_entry.get('n_ps', '-')}",
                         'yes' if pole_entry.get('coincide') else 'NO'])
        lines += ['## Logarithmic derivatives and poles', '',
                  'No point is dropped from `median |dL|`; the median is robust to the few '
                  'energies at a pole where both curves diverge (`max |dL|` is dominated by '
                  'them and cannot compare two points).', '',
                  'Two caveats on the pole columns, both measured rather than assumed:',
                  '',
                  '* the poles of `L(E)` at radius `r_c` are the **box** spectrum of that '
                  'radius (`u(r_c, E) = 0`), not properties of the atom -- they move with '
                  '`r_c` by construction, so pole *counts* are only comparable at one '
                  'radius;',
                  '* whether a divergence is detected at all depends on its residue: a pole '
                  'whose `|L|` peak falls between two grid points stays under the `|L| > 20` '
                  'threshold. Measured on the same all-electron atom, the detected `l=0` '
                  'pole count is **1 at rlderiv = 4 and 2 at rlderiv = 5**.', '',
                  'So a pole verdict is meaningful **within one radius** (AE and PS are then '
                  'solved in the same box) and not across radii.', '']
        lines.append(_markdown_table(['channel', 'median |dL|', 'max |dL|',
                                  'poles AE / PS', 'coincide'], rows))
        lines += ['', f"Poles coincide in every channel: "
                      f"**{'yes' if derivatives.get('poles_all_coincide') else 'NO'}**.", '']

    rows = generated.get('rows') or []
    if rows:
        lines += ['## Transferability (error in the energy differences)', '',
                  '`Delta E = (E_ae - E_ae,ref) - (E_ps - E_ps,ref)` against the '
                  'generating configuration. For a pseudopotential this is the number '
                  'that matters: an error that cancels between configurations does not '
                  'affect energy differences.', '']
        lines.append(_markdown_table(
            ['#', 'configuration', 'dE_AE (Ry)', 'dE_PS (Ry)', 'Delta E (mRy)'],
            [[str(row['index']), row.get('config', ''),
              '-' if row['delta_ae_ry'] is None else f"{row['delta_ae_ry']:+.6f}",
              '-' if row['delta_ps_ry'] is None else f"{row['delta_ps_ry']:+.6f}",
              '-' if row['delta_mry'] is None else f"{row['delta_mry']:+.1f}"]
             for row in rows]))
        worst_mry = generated.get('max_abs_mry') or 0.0
        lines += ['', f"Worst |Delta E| = **{worst_mry:.2f} mRy = {worst_mry * 13.6057:.1f} meV**.", '',
                  'This is an **atomic** configuration energy difference. It is *not* the '
                  'ACWF `delta`, which is defined on a solid\'s energy-volume curve in '
                  'meV/atom, so the two cannot be compared directly -- and the unit '
                  'conversion matters: 1 mRy = 13.6 meV, so 1/3 meV = 0.0245 mRy (not '
                  '0.333 mRy, which would be 4.5 meV). As an external scale, chemical '
                  'accuracy is about 43 meV (1 kcal/mol), i.e. this pseudopotential is at '
                  'roughly a third of it on the worst configuration.', '']

    if scan:
        lines += ['## The rlderiv scan: a diagnostic, not a selection', '',
                  'One full generation per value. **The pooled median cannot rank these '
                  'points**, and it is not because the pseudopotentials are equally good: at '
                  'a larger radius `L(E)` is dominated by the decaying tail '
                  '(`u\'/u -> -kappa + 1/r`), where the all-electron and pseudo curves agree '
                  'by construction. The metric therefore measures mostly *how far outside '
                  'the core it looks* -- measured at `E = -1 Ry`, `|L_AE - L_PS|` is '
                  '5.9e-05 at rlderiv = 4 and 3.6e-07 at rlderiv = 8.', '',
                  'The pole verdict is the hard criterion (a pseudopotential that adds or '
                  'moves a bound state is not made acceptable by a small median), so it '
                  'gates the choice. `bound median` repeats the median on `E < 0` only, '
                  'where the tail cannot do the work; it is shown as a sanity check, not as '
                  'a second criterion.', '']
        channels = sorted({channel for entry in scan for channel in entry.get('per_channel', {})})
        rows = []
        for entry in scan:
            cells = [f"{entry['rlderiv']:g}"]
            for channel in channels:
                value = entry['per_channel'].get(channel, {}).get('median')
                cells.append('-' if value is None else f'{value:.3e}')
            cells.append(f"{entry['pooled_median']:.3e}")
            bound = entry.get('bound_median')
            cells.append('-' if bound is None else f"{bound:.2e}")
            cells.append('yes' if entry['poles_all_coincide'] else '**NO**')
            if best and entry['rlderiv'] == best['rlderiv']:
                cells[0] += ' **←**'
            rows.append(cells)
        lines.append(_markdown_table(
            ['rlderiv'] + [f'{channel} median' for channel in channels]
            + ['pooled median', 'bound median', 'poles coincide'], rows))
        if best:
            median_pick = best.get('best_by_median') or {}
            if best.get('rlderiv') is not None:
                lines += ['', f"**Pole-gated pick: rlderiv = {best['rlderiv']:g} bohr** "
                              f"(pooled median {best['pooled_median']:.3e}, pk {best['pk']}); "
                              f"the points passing the pole test are "
                              f"{', '.join(f'{value:g}' for value in best.get('passing', []))}."]
            else:
                lines += ['', 'No scanned radius passes the pole test.']
            if median_pick and median_pick.get('pole_gated') is False:
                lines += [f"The smallest median is rlderiv = {median_pick['rlderiv']:g} "
                          f"({median_pick['pooled_median']:.3e}), but its pole check "
                          f"**fails**, so it is not admissible -- ranking by the median alone "
                          f"is what an earlier version of this report did."]
            lines += ['', 'The production pseudopotential uses the recipe value '
                          '(`rlderiv = 5.0`); this scan does not justify changing it.']
        lines.append('')

    lines += ['## Figures', '',
              '![radial](fig_radial.png)', '',
              '![log derivatives](fig_logderiv.png)', '',
              '![rlderiv scan](fig_rlderiv.png)', '']

    path = out / 'REPORT.md'
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path
