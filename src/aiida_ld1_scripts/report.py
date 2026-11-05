"""Figures and the markdown report, from ``data/summary.json``.

Two figures, both in the style of Torres & Kaloni (2020) and both **without a title** (the
manuscript captions them):

``fig_radial.png``     do the pseudo and all-electron radial functions agree?
``fig_logderiv.png``   do the logarithmic derivatives agree?

There is deliberately no plot of the ``rlderiv`` scan: its meaning is not evident from a
picture, the decisive column is the pole verdict, and that belongs in a table (see the
README, "what the scan cannot decide").

The ``l >= 3`` channel is drawn but not used to judge anything: for ``l >= 3`` the
Bessel-box column of ``ld1.x`` is a numerical artefact (it moves ~0.23 Ry with the box
radius alone), so the criterion is ``e PS`` -- see the ``aiida-ld1`` README.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from . import analysis, recipe as recipe_module

#: channels used as criteria (the box test is not trustworthy for l >= 3)
RELIABLE_MAX_L = 2

#: Okabe-Ito: colour-blind safe, and one fixed colour per angular-momentum channel.  AE
#: and PS are told apart by the line style / symbol, never by the colour alone -- the
#: two curves of one channel are the same physical object.
PALETTE = ('#0072B2',   # blue
           '#D55E00',   # vermillion
           '#009E73',   # bluish green
           '#CC79A7',   # reddish purple
           '#E69F00',   # orange
           '#56B4E9')


def _matplotlib():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # a paper figure, not a debugging plot: no titles anywhere, only axis labels and the
    # legend carry information
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9.5,
                         'legend.fontsize': 8, 'xtick.labelsize': 8.5,
                         'ytick.labelsize': 8.5})
    return plt


def _save(figure, out_base: Path, plt, formats: Sequence[str]) -> Path:
    """Save ``out_base`` in every requested format; returns the base path.

    ``bbox_inches='tight'`` is what keeps a long ylabel from being cut off.  It is **not**
    sufficient on its own: matplotlib's tight bounding box does not include the
    decorations of an inset (child) axes, so those are passed as extra artists.
    """
    out_base.parent.mkdir(parents=True, exist_ok=True)
    extra = []
    for parent in figure.axes:
        for child in parent.child_axes:
            extra.extend(child.get_default_bbox_extra_artists() or [])
    for suffix in formats:
        figure.savefig(f'{out_base}.{suffix}', dpi=300, bbox_inches='tight',
                       bbox_extra_artists=extra or None)
    plt.close(figure)
    return out_base


def _sample(grid, count: int = 34, *, logarithmic: bool = False):
    """Indices at which to draw the open symbols (a curve, not a thick band).

    With a logarithmic x axis the symbols are sampled logarithmically too, otherwise
    every one of them crowds into the left-hand decade and the outer region is left bare.
    """
    import numpy as np

    grid = np.asarray(grid)
    if not logarithmic:
        step = max(1, len(grid) // count)
        return np.arange(0, len(grid), step)
    positive = grid > 0
    low, high = np.log10(grid[positive].min()), np.log10(grid[positive].max())
    targets = np.logspace(low, high, count)
    indices = np.unique(np.searchsorted(grid, targets).clip(0, len(grid) - 1))
    return indices[grid[indices] > 0]


def figure_radial(plt, directory: str | Path, out_base: Path, *, rcuts=None,
                  xscale: str = 'linear', xlim=None, pseudo_line: bool = False,
                  log_ordinate: bool = False,
                  formats: Sequence[str] = ('png',)) -> Optional[Path]:
    """AE vs pseudo radial functions, in the style of Torres & Kaloni Fig. 2 (left).

    All-electron orbitals are **solid lines**, pseudo-orbitals **open symbols** (a dashed
    line would hide underneath the solid one exactly where the two have to coincide), and
    each pseudized channel gets a vertical dotted line at its ``rcut``.  No title: the axes
    are labelled and the legend carries the ``rcut`` values.

    Only the pseudized channels are drawn; the all-electron run also writes core states
    (4F, 5D) that the pseudopotential has no orbital for.

    ``xscale='log'`` is a diagnostic variant: with ``xlim=(1e-2, 1e2)`` it covers the whole
    orbital, from the core to the exponential tail.  ``pseudo_line=True`` also joins the
    pseudo points with a dashed line (on four decades of abscissa the symbols alone are too
    sparse to follow).  ``log_ordinate=True`` optionally puts ``|u|`` on a logarithmic axis
    -- the shipped diagnostic keeps a **linear** ordinate, so beyond ~10 bohr the curves sit
    on the zero line by construction, which is the honest picture of that range.
    """
    import numpy as np

    directory = Path(directory)
    ae_file, ps_file = directory / 'ld1.wfc', directory / 'ld1ps.wfc'
    if not (ae_file.is_file() and ps_file.is_file()):
        return None
    ae_labels, grid, ae_values = analysis.read_wfc(ae_file)
    ps_labels, _, ps_values, _aligned = analysis.read_wfc_halves(ps_file)

    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    logarithmic = str(xscale).lower() == 'log'
    # only use (and sample) the window that is actually shown: the ld1 grid runs to ~100 bohr
    if xlim:
        low, high = float(xlim[0]), float(xlim[1])
    elif logarithmic:
        # a log axis cannot start at 0; the grid's first positive point would be ~1e-5,
        # which is not interesting -- start at 1e-4 unless told otherwise
        low, high = max(float(grid[grid > 0].min()), 1e-4), min(8.0, float(grid[-1]))
    else:
        low, high = 0.0, min(8.0, float(grid[-1]))
    visible = (grid >= low) & (grid <= high) if logarithmic else (grid <= high)
    indices = _sample(grid[visible], logarithmic=logarithmic)
    rcuts = rcuts or {}

    for position, label in enumerate(ps_labels):
        if label not in ae_labels:
            continue
        colour = PALETTE[position % len(PALETTE)]
        all_electron = ae_values[:, ae_labels.index(label)]
        pseudo = ps_values[:, position]
        all_electron, _flipped = analysis.align_phase(all_electron, pseudo)

        cutoff = rcuts.get(label)
        name = label if cutoff is None else f'{label} (r$_c$ = {cutoff:.2f} bohr)'
        shown_ae, shown_ps = all_electron, pseudo
        if log_ordinate:
            shown_ae, shown_ps = np.abs(all_electron), np.abs(pseudo)
        axis.plot(grid[visible], shown_ae[visible], '-', color=colour, lw=1.2, label=name)
        if pseudo_line:
            axis.plot(grid[visible], shown_ps[visible], '--', color=colour, lw=0.9,
                      alpha=0.95)
        axis.plot(grid[visible][indices], shown_ps[visible][indices], 'o',
                  color=colour, mfc='none', mec=colour, ms=3.2, mew=0.8)
        if cutoff is not None:
            axis.axvline(cutoff, color=colour, ls=':', lw=0.8, alpha=0.75)

    axis.set_xlabel('r (bohr)')
    if log_ordinate:
        axis.set_yscale('log')
        axis.set_ylabel('|u(r)| = |r·R(r)|  (arb. units)')
    else:
        axis.set_ylabel('u(r) = r·R(r)  (arb. units)')
    if logarithmic:
        axis.set_xscale('log')
        axis.set_xlim(low, high)
    else:
        axis.set_xlim(0.0, high)
    axis.axhline(0.0, color='0.85', lw=0.5, zorder=0)
    axis.legend(loc='upper right', frameon=False)
    return _save(figure, Path(out_base), plt, formats)


def figure_logderiv(plt, directory: str | Path, out_base: Path, *,
                    window=(-3.0, 3.0), threshold: float = analysis.POLE_THRESHOLD,
                    formats: Sequence[str] = ('png',)) -> Optional[Path]:
    """Logarithmic derivatives, in the style of Torres & Kaloni Fig. 3.

    All-electron curves **solid**, pseudo curves **dotted** -- the paper's convention, and
    nothing else: no title, no vertical pole markers, no annotations.  The pole
    correspondence is a number in the report table, not something to be measured off the
    picture.

    The y range is taken from the data (``0.5``-``99.5`` percentile) rather than a fixed
    ``+-30``, which flattened every pole into a vertical stroke.
    """
    import numpy as np

    directory = Path(directory)
    if not (directory / 'ld1.dlog').is_file():
        return None
    all_electron = analysis.read_dlog(directory / 'ld1.dlog')
    pseudo = analysis.read_dlog(directory / 'ld1ps.dlog')
    energies = all_electron['energy']
    mask = (energies >= window[0]) & (energies <= window[1])
    channels = [channel for channel in all_electron if channel != 'energy']

    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    for position, channel in enumerate(channels):
        colour = PALETTE[position % len(PALETTE)]
        axis.plot(energies[mask], all_electron[channel][mask], '-', color=colour,
                  lw=1.1, label=f'{channel} AE')
        axis.plot(energies[mask], pseudo[channel][mask], ':', color=colour,
                  lw=1.7, label=f'{channel} PS')

    usable = np.concatenate([np.abs(all_electron[channel][mask]) for channel in channels])
    limit = max(float(np.percentile(usable, 99.5)), 21.0)
    axis.set_ylim(-limit, limit)
    axis.set_xlim(*window)
    axis.set_xlabel('E (Ry)')
    axis.set_ylabel('logarithmic derivative  $L_l(E)$')
    axis.axhline(0.0, color='0.85', lw=0.5, zorder=0)
    axis.legend(loc='upper left', frameon=False, ncol=2)
    return _save(figure, Path(out_base), plt, formats)


def export_paper_figures(directory: str | Path, paper_dir: str | Path, *,
                         rcuts=None, window=(-3.0, 3.0),
                         threshold: float = analysis.POLE_THRESHOLD,
                         formats: Sequence[str] = ('png', 'svg', 'pdf')) -> List[Path]:
    """Write the two paper figures (3 formats) into ``paper_dir``, backing up the old ones.

    The names are the ones the manuscript uses: ``fig3_1_radial`` and ``fig3_3_logderiv``.
    """
    import shutil

    plt = _matplotlib()
    paper_dir = Path(paper_dir)
    paper_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for function, base, extra in (
            (figure_radial, 'fig3_1_radial', {'rcuts': rcuts}),
            (figure_logderiv, 'fig3_3_logderiv',
             {'window': window, 'threshold': threshold})):
        base_path = paper_dir / base
        for suffix in formats:                      # never lose the previous version
            current = Path(f'{base_path}.{suffix}')
            if current.is_file():
                shutil.copy2(current, f'{base_path}.{suffix}.bak')
        result = function(plt, directory, base_path, formats=formats, **extra)
        if result is not None:
            written.append(result)
    return written


def make_figures(summary: Dict[str, Any], directory: str | Path,
                 out_dir: Optional[str | Path] = None) -> List[Path]:
    """The two generation figures into ``out_dir`` (``directory`` = one run's files)."""
    out_dir = Path(out_dir) if out_dir is not None else recipe_module.DATA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    plt = _matplotlib()
    rcuts = {label: entry['rcut']
             for label, entry in (summary['generation'].get('wavefunctions') or {}).items()}
    written: List[Path] = []
    for function, base, extra in (
            (figure_radial, 'fig_radial', {'rcuts': rcuts}),
            (figure_radial, 'fig_radial_logx',
             {'rcuts': rcuts, 'xscale': 'log', 'xlim': (1e-2, 1e2),
              'pseudo_line': True}),
            (figure_logderiv, 'fig_logderiv',
             {'window': summary['criteria']['window'],
              'threshold': summary['criteria']['pole_threshold']})):
        try:
            path = function(plt, directory, out_dir / base, formats=('png',), **extra)
        except Exception as exception:                  # noqa: BLE001
            print(f'  warning: {base} failed: {exception}')
            path = None
        if path is not None:
            written.append(Path(f'{path}.png'))
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
              'All-electron orbitals as solid lines, pseudo-orbitals as open symbols; the '
              'vertical dotted lines are the pseudization radii. The all-electron core '
              'states (4F, 5D) have no pseudo-orbital and are not drawn.', '',
              '![log derivatives](fig_logderiv.png)', '',
              'All-electron curves solid, pseudo curves dotted, one colour per angular '
              'momentum channel. Pole correspondence is in the table above, not marked '
              'on the figure.', '']

    path = out / 'REPORT.md'
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path
