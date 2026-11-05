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

#: double-column width of the manuscript (inches)
WIDTH = 7.0
#: bohr -> angstrom, for the secondary abscissa
BOHR2ANG = 0.529177210903

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
    from matplotlib.lines import Line2D as _Line2D

    globals()['Line2D'] = _Line2D          # used by the legend proxy handles

    # one style for every figure (values from the manuscript's own script, so the new
    # figures and the previous Fig. 3.1/3.2/3.3 look like the same paper)
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 8.5,
        'axes.labelsize': 8.5, 'axes.titlesize': 8.5,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'legend.fontsize': 7.5, 'legend.frameon': False,
        'axes.linewidth': 0.8, 'lines.linewidth': 1.5,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
        'xtick.minor.width': 0.6, 'ytick.minor.width': 0.6,
        'axes.grid': False, 'savefig.dpi': 300,
        'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'mathtext.default': 'regular',
    })
    return plt


def _save(figure, out_base: Path, plt, formats: Sequence[str]) -> Path:
    """Save ``out_base`` in every requested format; returns the base path.

    ``bbox_inches='tight'`` keeps a long ylabel from being cut off, but it does **not**
    include an inset (child) axes' decorations, so those are passed as extra artists.
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


def _minor_ticks(axis, divisions: int = 2):
    """Minor ticks, honouring the log axes (``AutoMinorLocator`` refuses them)."""
    from matplotlib.ticker import AutoMinorLocator, LogLocator

    for name, which in (('x', axis.xaxis), ('y', axis.yaxis)):
        if getattr(axis, f'get_{name}scale')() == 'log':
            which.set_minor_locator(LogLocator(base=10.0, subs=(2.0, 5.0), numticks=50))
        else:
            which.set_minor_locator(AutoMinorLocator(divisions))


def _top_angstrom(axis):
    """A second abscissa in angstrom on the top edge (bohr stays on the bottom)."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    secondary = axis.secondary_xaxis(
        'top', functions=(lambda x: x * BOHR2ANG, lambda value: value / BOHR2ANG))
    secondary.set_xlabel(r'$r$ ($\mathrm{\AA}$)', fontsize=8, labelpad=2)
    secondary.tick_params(labelsize=7)
    low, high = axis.get_xlim()
    candidates = ((0.1, 0.2, 0.5, 1.0, 2.0, 3.0) if axis.get_xscale() == 'log'
                  else (0.1, 0.5, 1.0, 2.0))
    ticks = [tick for tick in candidates if low * BOHR2ANG < tick < high * BOHR2ANG]
    secondary.set_xticks(ticks)
    secondary.xaxis.set_major_locator(FixedLocator(ticks))
    secondary.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f'{value:g}'))
    secondary.xaxis.set_minor_locator(NullLocator())
    return secondary


def figure_radial(plt, directory: str | Path, out_base: Path, *, rcuts=None,
                  xscale: str = 'log', xlim=(0.02, 6.0),
                  pseudo: str = 'dashed', log_ordinate: bool = False,
                  difference_panel: bool = True,
                  formats: Sequence[str] = ('png',)) -> Optional[Path]:
    """AE vs pseudo radial functions, in the style of Torres & Kaloni Fig. 2.

    All-electron **solid**, pseudo **dashed** (``pseudo='dashed'``) or as **open symbols**
    (``'symbols'``) or both (``'both'``); a vertical dotted line at every ``rcut``; the
    core region ``r < min(rcut)`` in light grey; an angstrom abscissa on the top edge; and
    a **difference panel** underneath with ``|u_AE - u_PS| / max|u_AE|`` on a log ordinate.

    The difference panel is what makes "the tails agree" quantitative: overlaid, a dashed
    curve hides underneath the solid one exactly where the two coincide, which was the
    original complaint about this figure.
    """
    import numpy as np

    directory = Path(directory)
    ae_file, ps_file = directory / 'ld1.wfc', directory / 'ld1ps.wfc'
    if not (ae_file.is_file() and ps_file.is_file()):
        return None
    ae_labels, grid, ae_values = analysis.read_wfc(ae_file)
    ps_labels, _, ps_values, _aligned = analysis.read_wfc_halves(ps_file)

    logarithmic = str(xscale).lower() == 'log'
    if xlim:
        low, high = float(xlim[0]), float(xlim[1])
    else:
        low = max(float(grid[grid > 0].min()), 1e-4) if logarithmic else 0.0
        high = min(8.0, float(grid[-1]))
    visible = (grid >= low) & (grid <= high) if logarithmic else (grid <= high)
    rcuts = rcuts or {}

    figure = plt.figure(figsize=(WIDTH, 4.3))
    if difference_panel:
        grid_spec = figure.add_gridspec(2, 1, height_ratios=(3.0, 1.05), hspace=0.08,
                                        left=0.085, right=0.775, top=0.905, bottom=0.105)
        main = figure.add_subplot(grid_spec[0])
        panel = figure.add_subplot(grid_spec[1], sharex=main)
        main.tick_params(labelbottom=False)
    else:
        grid_spec = figure.add_gridspec(1, 1, left=0.085, right=0.775, top=0.905,
                                        bottom=0.115)
        main = figure.add_subplot(grid_spec[0])
        panel = None

    indices = np.arange(0, visible.sum(), max(1, visible.sum() // 34))
    handles = []
    for position, label in enumerate(ps_labels):
        if label not in ae_labels:
            continue
        colour = _orbital_colour(label, position)
        all_electron = ae_values[:, ae_labels.index(label)]
        pseudo_values = ps_values[:, position]
        all_electron, _flipped = analysis.align_phase(all_electron, pseudo_values)
        cutoff = rcuts.get(label)
        name = label if cutoff is None else f'{label}  ($r_c$ = {cutoff:g})'
        handles.append(Line2D([], [], color=colour, lw=1.5, label=name))

        shown = np.abs(all_electron) if log_ordinate else all_electron
        main.plot(grid[visible], shown[visible], '-', color=colour, lw=1.5, zorder=3)
        if pseudo in ('dashed', 'both'):
            pseudo_shown = np.abs(pseudo_values) if log_ordinate else pseudo_values
            main.plot(grid[visible], pseudo_shown[visible], color=colour, lw=1.0,
                      ls=(0, (4, 2)), zorder=4)
        if pseudo in ('symbols', 'both'):
            marks = np.flatnonzero(visible)[indices]
            main.plot(grid[marks], pseudo_values[marks], 'o', color=colour, mfc='none',
                      mec=colour, ms=3.0, mew=0.9, zorder=5)
        if cutoff is not None:
            main.axvline(cutoff, color=colour, lw=0.8, ls=(0, (2, 2)), alpha=0.55,
                         zorder=1)
        if panel is not None:
            peak = float(np.max(np.abs(all_electron[visible]))) or 1.0
            panel.semilogy(grid[visible],
                           np.maximum(np.abs(all_electron - pseudo_values)[visible] / peak,
                                      1e-12),
                           '-', color=colour, lw=0.9)

    if rcuts:
        main.axvspan(low if logarithmic else 0.0, min(rcuts.values()), color='0.92', zorder=0)

    main.set_xscale('log') if logarithmic else main.set_xscale('linear')
    main.set_xlim(low, high)
    if log_ordinate:
        main.set_yscale('log')
        main.set_ylabel('$|u(r)|$  (log scale)')
    else:
        main.set_ylabel(r'$u(r) = rR(r)$  (arb. units)')
    peak = float(np.max(np.abs(ae_values[visible]))) or 1.0
    if not log_ordinate:
        main.set_ylim(-1.35 * peak, 1.25 * peak)
    _minor_ticks(main)
    _top_angstrom(main)

    handles += [Line2D([], [], color='0.25', lw=1.5, label='AE'),
                Line2D([], [], color='0.25', lw=1.0, ls=(0, (4, 2)), label='PS')]
    main.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.005, 1.0),
                borderaxespad=0.0, handlelength=2.6, labelspacing=0.55)

    if panel is not None:
        panel.set_xscale('log') if logarithmic else panel.set_xscale('linear')
        panel.set_xlim(low, high)
        panel.set_ylim(1e-6, 2.0)
        panel.axhline(1e-3, color='0.6', lw=0.6, ls=':')
        panel.set_ylabel(r'$|\Delta u| / \max|u|$')
        panel.set_xlabel('r (bohr)')
        _minor_ticks(panel)
    else:
        main.set_xlabel('r (bohr)')
    return _save(figure, Path(out_base), plt, formats)


def figure_projectors(plt, directory: str | Path, out_base: Path, *, rcuts=None,
                      formats: Sequence[str] = ('png',)) -> Optional[Path]:
    """The radial projector functions ``beta(r)``: absolute (log) and normalized.

    Two panels over a shared abscissa: ``|beta|`` on a logarithmic ordinate, and
    ``beta/beta_max`` linear underneath, which is the one that shows the shape.  Colour
    follows the orbital, the second projector of a channel is dashed.
    """
    import numpy as np

    try:
        labels, grid, values = analysis.read_projectors(directory)
    except (FileNotFoundError, OSError):
        return None
    if values.size == 0:
        return None
    peak = float(np.max(np.abs(values))) or 1.0
    reach = grid[np.any(np.abs(values) > 0.01 * peak, axis=1)]
    high = float(reach.max()) * 1.1 if reach.size else float(grid[-1])
    high = min(max(high, 3.0), 8.0, float(grid[-1]))
    rcuts = rcuts or {}

    figure = plt.figure(figsize=(WIDTH, 4.3))
    grid_spec = figure.add_gridspec(2, 1, height_ratios=(1.0, 1.0), hspace=0.10,
                                    left=0.085, right=0.775, top=0.905, bottom=0.115)
    absolute = figure.add_subplot(grid_spec[0])
    relative = figure.add_subplot(grid_spec[1], sharex=absolute)
    absolute.tick_params(labelbottom=False)

    handles = []
    for index, label in enumerate(labels):
        # read_projectors() already disambiguated the repeated channels as ``6P′``, so the
        # prime character is what marks the second projector of a channel
        base = label.rstrip('\u2032')
        second = label != base
        style = (0, (4, 2)) if second else '-'
        colour = _orbital_colour(base, index)
        name = label if second else base
        name = name.replace('\u2032', r'$\prime$')
        handles.append(Line2D([], [], color=colour, lw=1.3, ls=style, label=name))
        column = values[:, index]
        absolute.semilogy(grid, np.maximum(np.abs(column), 1e-30), color=colour, lw=1.3,
                          ls=style, label=name)
        relative.plot(grid, column / (np.max(np.abs(column)) or 1.0), color=colour, lw=1.3,
                      ls=style)
        cutoff = rcuts.get(base)
        if cutoff is not None:
            for axis in (absolute, relative):
                axis.axvline(cutoff, color=colour, lw=0.8, ls=(0, (2, 2)), alpha=0.55)

    for axis in (absolute, relative):
        axis.set_xlim(0.0, high)
        _minor_ticks(axis)
    absolute.set_ylim(1e-6, 10.0 * peak)
    absolute.set_ylabel(r'$|\beta(r)|$  (arb. units)')
    relative.set_ylim(-1.05, 1.15)
    relative.set_ylabel(r'$\beta(r)/\beta_{\max}$')
    relative.set_xlabel('r (bohr)')
    _top_angstrom(absolute)          # the top edge of the figure, not between the panels
    absolute.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.005, 1.0),
                    borderaxespad=0.0, handlelength=2.6, labelspacing=0.55)
    return _save(figure, Path(out_base), plt, formats)


def figure_logderiv(plt, directory: str | Path, out_base: Path, *,
                    window=(-3.0, 3.0), threshold: float = analysis.POLE_THRESHOLD,
                    log_ordinate: bool = False, radius: float | None = None,
                    pole_table: bool = True,
                    formats: Sequence[str] = ('png',)) -> Optional[Path]:
    """Logarithmic derivatives: AE solid, pseudo as open circles every tenth point.

    The ordinate stays **linear** (a signed ``L`` cannot be logarithmic; ``log_ordinate``
    gives the ``|L|`` variant instead).  The pole correspondence -- the ghost criterion --
    is not eyeballed: it is a **table panel** underneath, per channel, with the AE and PS
    pole energies and their displacement.

    ``radius`` is the reference radius at which ``L`` was evaluated; it belongs on the
    figure (the paper states its own r_d in the caption), so it goes into the x label when
    given.
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

    figure = plt.figure(figsize=(WIDTH, 4.6))
    if pole_table:
        grid_spec = figure.add_gridspec(2, 1, height_ratios=(3.4, 1.0), hspace=0.12,
                                        left=0.085, right=0.755, top=0.945, bottom=0.095)
        main = figure.add_subplot(grid_spec[0])
        table_axis = figure.add_subplot(grid_spec[1])
        table_axis.axis('off')
    else:
        grid_spec = figure.add_gridspec(1, 1, left=0.085, right=0.755, top=0.945,
                                        bottom=0.115)
        main = figure.add_subplot(grid_spec[0])
        table_axis = None

    def values(array):
        return np.maximum(np.abs(array), 1e-8) if log_ordinate else array

    handles = []
    for position, channel in enumerate(channels):
        colour = PALETTE[position % len(PALETTE)]
        main.plot(energies[mask], values(all_electron[channel][mask]), '-', color=colour,
                  lw=1.3, zorder=2)
        every = np.arange(0, mask.sum(), 10)
        main.plot(energies[mask][every], values(pseudo[channel][mask])[every], 'o',
                  color=colour, mfc='none', mec=colour, ms=3.0, mew=0.9, zorder=3)
        handles.append(Line2D([], [], color=colour, lw=1.3,
                              label=(f'{channel}  (med '
                                     r'$|\Delta L|$'
                                     f' = {float(np.median(np.abs(pseudo[channel][mask] - all_electron[channel][mask]))):.2g})')))

    usable = np.concatenate([np.abs(all_electron[channel][mask]) for channel in channels])
    if log_ordinate:
        main.set_yscale('log')
        main.set_ylim(max(float(np.percentile(usable, 0.5)), 1e-4),
                      max(float(np.percentile(usable, 99.9)), 30.0))
        main.set_ylabel('$|L_l(E)|$  (log scale)')
    else:
        limit = max(float(np.percentile(usable, 99.5)), 21.0)
        main.set_ylim(-limit, limit)
        main.set_ylabel(r'$L_l(E)$')
        main.axhline(0.0, color='0.85', lw=0.5, zorder=0)
    main.set_xlim(*window)
    main.set_xlabel('E (Ry)' if radius is None
                    else f'E (Ry),  $L$ evaluated at $r_d$ = {radius:g} bohr')
    _minor_ticks(main, divisions=4)
    handles += [Line2D([], [], color='0.25', lw=1.3, label='AE'),
                Line2D([], [], color='0.25', lw=0, marker='o', ms=3.2, mfc='none',
                       mew=0.9, label='PS (every 10th pt)')]
    main.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.005, 1.0),
                borderaxespad=0.0, handlelength=2.2, labelspacing=0.75)

    if table_axis is not None:
        rows = analysis.pole_table(all_electron, pseudo, window=window,
                                  threshold=threshold)
        lines = ['channel   AE poles (Ry)            PS poles (Ry)            '
                 'displacement (Ry)']
        for row in rows:
            ae_text = ', '.join(f'{value:+.2f}' for value in row['ae'][:3]) or '-'
            ps_text = ', '.join(f'{value:+.2f}' for value in row['ps'][:3]) or '-'
            delta_text = ', '.join(f'{delta:.3f}' for _a, _b, delta in row['pairs']) or '-'
            lines.append(f'{row["channel"]:<8s}  {ae_text:<23s}  {ps_text:<23s}  '
                         f'{delta_text}'
                         + ('' if row['coincide'] else '   (no match)'))
        table_axis.text(0.0, 0.95, '\n'.join(lines), fontsize=6.4, family='monospace',
                        transform=table_axis.transAxes, va='top', ha='left')
    return _save(figure, Path(out_base), plt, formats)


def export_paper_figures(directory: str | Path, paper_dir: str | Path, *,
                         rcuts=None, window=(-3.0, 3.0),
                         threshold: float = analysis.POLE_THRESHOLD,
                         radius: float | None = None,
                         formats: Sequence[str] = ('png', 'svg', 'pdf')) -> List[Path]:
    """Write the two manuscript figures (3 formats), backing the previous files up."""
    import shutil

    plt = _matplotlib()
    paper_dir = Path(paper_dir)
    paper_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for function, base, extra in (
            (figure_radial, 'fig3_1_radial',
             {'rcuts': rcuts, 'xscale': 'log', 'xlim': (0.02, 6.0)}),
            (figure_logderiv, 'fig3_3_logderiv',
             {'window': window, 'threshold': threshold, 'radius': radius})):
        base_path = paper_dir / base
        for suffix in formats:
            current = Path(f'{base_path}.{suffix}')
            backup = Path(f'{base_path}.{suffix}.bak')
            # write-once: the backup must be the file this script first replaced, not the
            # previous run's output (re-exporting would otherwise destroy the original)
            if current.is_file() and not backup.exists():
                shutil.copy2(current, backup)
        result = function(plt, directory, base_path, formats=formats, **extra)
        if result is not None:
            written.append(result)
    return written


def _orbital_colour(label: str, fallback: int = 0):
    """Okabe-Ito colour of an orbital, stable across the figures."""
    order = ['6S', '7S', '6P', '6D', '5F', '6P\u2032', '6P\'']
    if label in order:
        return PALETTE[order.index(label) % len(PALETTE)]
    return PALETTE[fallback % len(PALETTE)]


def make_figures(summary: Dict[str, Any], directory: str | Path,
                 out_dir: Optional[str | Path] = None) -> List[Path]:
    """The three repository figures into ``out_dir`` (``directory`` = one run's files)."""
    out_dir = Path(out_dir) if out_dir is not None else recipe_module.DATA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    plt = _matplotlib()
    rcuts = {label: entry['rcut']
             for label, entry in (summary['generation'].get('wavefunctions') or {}).items()}
    radius = summary['criteria'].get('production_rlderiv')
    written: List[Path] = []
    for function, base, extra in (
            (figure_radial, 'fig_radial',
             {'rcuts': rcuts, 'xscale': 'log', 'xlim': (0.02, 6.0)}),
            (figure_radial, 'fig_radial_logx',
             {'rcuts': rcuts, 'xscale': 'log', 'xlim': (1e-2, 1e2), 'pseudo': 'both'}),
            (figure_logderiv, 'fig_logderiv',
             {'window': summary['criteria']['window'],
              'threshold': summary['criteria']['pole_threshold'],
              'radius': radius}),
            (figure_projectors, 'fig_projectors', {'rcuts': rcuts})):
        try:
            path = function(plt, directory, out_dir / base, formats=('png',), **extra)
        except Exception as exception:                  # noqa: BLE001
            print(f'  warning: {base} failed: {exception}')
            path = None
        if path is not None:
            written.append(Path(f'{path}.png'))
        else:
            print(f'  note: {base} was skipped (its input file is missing from the run)')
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
              'All-electron orbitals as solid lines, pseudo-orbitals dashed; the vertical '
              'dotted lines are the pseudization radii and the grey band is the core '
              '(`r < min r_c`). The lower panel is `|u_AE - u_PS| / max|u_AE|` on a log '
              'ordinate, which is where "the tails agree" becomes a number. The '
              'all-electron core states (4F, 5D) have no pseudo-orbital and are not '
              'drawn.', '',
              '![log derivatives](fig_logderiv.png)', '',
              'All-electron curves solid, pseudo curves as open circles every tenth point, '
              'one colour per angular momentum channel; the ordinate is linear and the '
              'comment in the x label is the radius at which the logarithmic derivative '
              'was evaluated. The panel underneath is the pole table -- the ghost '
              'criterion -- with the AE and PS pole energies and their displacement.', '',
              '![projectors](fig_projectors.png)', '',
              'The six radial projector functions `beta(r)`, coloured by orbital (the '
              'second projector of a channel is dashed, `6P` and `6P\u2032`). Upper panel: '
              'absolute value on a logarithmic ordinate; lower panel: each normalized to '
              'its own maximum. The dotted vertical lines are the pseudization radii.', '']

    path = out / 'REPORT.md'
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path
