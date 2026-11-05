"""The figures must follow the paper's conventions and must not clip their labels.

Conventions (Torres & Kaloni 2020, Figs. 2-3): all-electron = **solid line**, pseudo =
**open symbols** (radial) or **dotted line** (log derivative); a vertical dotted line at
every pseudization radius; the pole correspondence shown explicitly.  The labels are
checked against the figure's *tight* bounding box, which is what ``bbox_inches='tight'``
uses when saving -- an earlier version of the rlderiv figure cut its ylabel and title.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from aiida_ld1_scripts import report as report_module            # noqa: E402
from aiida_ld1_scripts import analysis                           # noqa: E402


def write_run(directory: Path) -> Path:
    """A tiny run directory: 3 all-electron columns, 2 pseudized channels, a pole."""
    directory.mkdir(parents=True, exist_ok=True)
    r = np.linspace(0.01, 8.0, 200)
    core = np.exp(-2.0 * r)               # an all-electron core state: no pseudo partner
    s = np.exp(-r) * (1.0 - np.exp(-6 * r))
    f = r ** 3 * np.exp(-r)
    # ld1.wfc: AE, including one channel the pseudopotential has no orbital for
    (directory / 'ld1.wfc').write_text(
        '#            r                  6S                  4F                  5F\n'
        + ''.join(f'{x:20.12f}{a:20.12f}{c:20.12f}{b:20.12f}\n'
                  for x, a, b, c in zip(r, s, f, core)), encoding='utf-8')
    # ld1ps.wfc: PS half then phase-aligned AE half (the all-electron one is negated,
    # which is what align_phase has to detect)
    (directory / 'ld1ps.wfc').write_text(
        '#            r                  6S                  5F'
        '                  6S                  5F\n'
        + ''.join(f'{x:20.12f}{a * 0.8 if x < 1.5 else a:20.12f}'
                  f'{b * 0.8 if x < 1.5 else b:20.12f}{-a:20.12f}{-b:20.12f}\n'
                  for x, a, b in zip(r, s, f)), encoding='utf-8')
    # dlog: a pole for l=0 near -0.3 Ry (sign flip with |L| > 20), smooth otherwise
    energy = np.round(np.arange(-3.0, 3.01, 0.1), 10)
    l0 = np.where(energy < -0.3, 60.0, -60.0) + 0.01 * energy
    l1 = 0.1 * energy
    for name, shift in (('ld1.dlog', 0.0), ('ld1ps.dlog', 0.003)):
        lines = [f'{e + shift:20.12f}{a:20.12f}{b:20.12f}'
                 for e, a, b in zip(energy, l0, l1)]
        (directory / name).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return directory


@pytest.fixture()
def captured(monkeypatch, tmp_path):
    """Render a figure and hand it over instead of closing it."""
    from matplotlib.text import Text

    figures = {}

    def fake_save(figure, out_base, plt, formats):
        figures[Path(out_base).name] = figure
        return out_base

    monkeypatch.setattr(report_module, '_save', fake_save)
    plt = report_module._matplotlib()
    directory = write_run(tmp_path / 'run')

    def render(which, **kwargs):
        figures.clear()
        function = getattr(report_module, which)
        function(plt, directory, tmp_path / which, formats=('png',), **kwargs)
        figure = figures[which]
        renderer = figure.canvas.get_renderer()
        figure.canvas.draw()
        return figure, renderer, Text

    return render, plt


def outside_tight_bbox(figure, renderer, Text=None, tolerance: float = 0.06):
    """Decorations outside the box ``bbox_inches='tight'`` would keep (inches).

    Built from the **current** tight bounding boxes of the axes (and of their inset child
    axes) plus the figure's own texts, not from ``findobj(Text)``: matplotlib keeps tick
    labels of previous limits in the artist tree with stale positions, so walking every
    ``Text`` reports labels that are not drawn at all.

    Note the unit trap: ``Figure.get_tightbbox`` returns **inches** while
    ``Axes.get_tightbbox`` returns **pixels**, so the axes boxes have to be converted.

    ``tolerance`` (0.06 in) absorbs the half-a-label-height overhang of a tick label on
    the frame; a cut-off ylabel or a two-line title misses by 0.2 in or more, which is the
    case this guards (it was cut in an earlier version of ``fig_rlderiv``).
    """
    tight = figure.get_tightbbox(renderer)
    to_inches = figure.dpi_scale_trans.inverted()
    candidates = []
    for axes in figure.axes:
        box = axes.get_tightbbox(renderer)
        if box is not None:                        # pixels
            candidates.append(box.transformed(to_inches))
        for child in axes.child_axes:              # the inset is not in figure.axes
            box = child.get_tightbbox(renderer)
            if box is not None:
                candidates.append(box.transformed(to_inches))
    for text in figure.texts:                       # the "not drawn" footnote
        if text.get_visible() and text.get_text().strip():
            candidates.append(text.get_window_extent(renderer)
                              .transformed(to_inches))
    return [box for box in candidates
            if (box.x0 < tight.x0 - tolerance or box.y0 < tight.y0 - tolerance
                or box.x1 > tight.x1 + tolerance or box.y1 > tight.y1 + tolerance)]


def test_radial_uses_solid_lines_and_open_symbols(captured):
    """Paper Fig. 2: AE solid, pseudo open symbols, one rcut line per pseudized channel."""
    render, _plt = captured
    figure, renderer, Text = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    axes = figure.axes[0]
    markers = [line.get_marker() for line in axes.get_lines()]
    open_symbol_series = [line for line in axes.get_lines() if line.get_marker() == 'o']
    assert len(open_symbol_series) == 2                       # 6S and 5F only
    assert all(line.get_markerfacecolor() in ('none', 'None', (0, 0, 0, 0))
               for line in open_symbol_series)
    vertical = [line for line in axes.get_lines()
                if len(set(line.get_xdata())) == 1 and line.get_linestyle() == ':']
    assert len(vertical) == 2                                 # one per rcut
    # the all-electron core state 4F has no pseudo partner and must not be drawn
    labels = [text.get_text() for text in axes.get_legend().get_texts()]
    assert all('4F' not in label for label in labels)
    assert outside_tight_bbox(figure, renderer, Text) == []


def test_the_figures_carry_no_title(captured):
    """The figures are captioned by the manuscript, not by their own title."""
    render, _plt = captured
    for which, extra in (('figure_radial', {'rcuts': {'6S': 1.5, '5F': 1.3}}),
                         ('figure_logderiv', {'window': (-3.0, 3.0), 'threshold': 20.0})):
        figure, _renderer, _Text = render(which, **extra)
        assert figure.axes[0].get_title() == '', which
        assert not figure.texts, which          # no figure-level titles either


def test_radial_has_no_difference_inset(captured):
    """The |u_AE - u_PS| inset was removed: the report table carries those numbers."""
    render, _plt = captured
    figure, _renderer, _Text = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    assert figure.axes[0].child_axes == []


def test_logderiv_is_curves_only(captured):
    """Paper Fig. 3: AE solid, PS dotted, and nothing else -- no pole markers."""
    render, _plt = captured
    figure, renderer, Text = render('figure_logderiv', window=(-3.0, 3.0), threshold=20.0)
    axes = figure.axes[0]
    dotted = [line for line in axes.get_lines() if line.get_linestyle() == ':']
    assert len(dotted) == 2                    # one per channel (l=0, l=1)
    vertical = [line for line in axes.get_lines()
                if len(set(line.get_xdata())) == 1]
    assert vertical == []                      # no vertical pole lines
    assert [line for line in axes.get_lines() if line.get_marker() == '^'] == []
    assert len(axes.texts) == 0                # no pole-energy annotations
    assert outside_tight_bbox(figure, renderer, Text) == []


def test_export_paper_figures_backs_up_the_previous_files(tmp_path, monkeypatch):
    """Overwriting the manuscript figures must not lose the previous version."""
    directory = write_run(tmp_path / 'run')
    paper = tmp_path / 'figs'
    paper.mkdir()
    (paper / 'fig3_1_radial.png').write_text('old radial', encoding='utf-8')
    (paper / 'fig3_3_logderiv.png').write_text('old logderiv', encoding='utf-8')

    plt = report_module._matplotlib()
    monkeypatch.setattr(report_module, '_matplotlib', lambda: plt)
    written = report_module.export_paper_figures(
        directory, paper, rcuts={'6S': 1.5, '5F': 1.3}, formats=('png',))
    assert len(written) == 2
    for name in ('fig3_1_radial', 'fig3_3_logderiv'):
        assert Path(f'{paper}/{name}.png.bak').read_text(encoding='utf-8').startswith('old')
        assert Path(f'{paper}/{name}.png').stat().st_size > 1000


def test_save_uses_the_tight_bounding_box(tmp_path, monkeypatch):
    """``bbox_inches='tight'`` is what guarantees a long ylabel cannot be cut off."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    seen = {}
    figure, axis = plt.subplots()
    axis.set_ylabel('a very long ylabel that would be clipped by a fixed canvas')
    original = figure.savefig

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(figure, 'savefig', spy)
    report_module._save(figure, tmp_path / 'f', plt, ('png',))
    assert seen.get('bbox_inches') == 'tight'
    assert seen.get('dpi', 0) >= 300


def test_the_log_x_variant_only_samples_the_visible_range(captured):
    """The diagnostic log-x radial figure, and no symbols outside the axes."""
    render, _plt = captured
    figure, _renderer, _Text = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3},
                                      xscale='log')
    axes = figure.axes[0]
    assert axes.get_xscale() == 'log'
    low, high = axes.get_xlim()
    assert low > 0.0 and high <= 8.0
    symbols = [line for line in axes.get_lines() if line.get_marker() == 'o']
    assert len(symbols) == 2
    for line in symbols:
        xdata = line.get_xdata()
        assert xdata.min() >= low and xdata.max() <= high


def test_the_diagnostic_variant_is_log_log_over_four_decades(captured):
    """``fig_radial_logx``: x from 1e-2 to 1e2, |u| logarithmic, PS dashed + symbols."""
    render, _plt = captured
    figure, _renderer, _Text = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3},
                                      xscale='log', xlim=(1e-2, 1e2),
                                      pseudo_line=True, log_ordinate=True)
    axes = figure.axes[0]
    assert axes.get_xscale() == 'log' and axes.get_yscale() == 'log'
    assert axes.get_xlim() == pytest.approx((0.01, 100.0))
    dashed = [line for line in axes.get_lines() if line.get_linestyle() == '--']
    symbols = [line for line in axes.get_lines() if line.get_marker() == 'o']
    assert len(dashed) == 2 and len(symbols) == 2
    # the symbols must carry the channel colour, not matplotlib's default cycle
    solid = [line for line in axes.get_lines()
             if line.get_linestyle() == '-' and line.get_color() != '0.85']
    assert [line.get_color() for line in symbols] == [line.get_color() for line in solid]
    assert axes.get_ylabel().startswith('|u(r)|')
    assert axes.get_title() == '' and len(axes.child_axes) == 0
