"""The figures must follow the house style and must not clip their labels.

House style (the one the manuscript's own figure script uses): 7.0 in double-column
width, 8.5 pt base, minor ticks, 0.8 pt axes, AE **solid** lw 1.5 and PS **dashed**
lw 1.0, a vertical dotted line per pseudization radius, the core region shaded, an
angstrom abscissa on the top edge, legends **outside** the axes, and no title.

The log derivative is linear in L (a signed quantity cannot be logarithmic) and carries a
**pole table** panel; the projectors are shown absolute (log) and normalized.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use('Agg')

from aiida_ld1_scripts import report as report_module            # noqa: E402


def write_run(directory: Path) -> Path:
    """A tiny run directory: 3 all-electron columns, 2 pseudized channels, a pole."""
    directory.mkdir(parents=True, exist_ok=True)
    r = np.linspace(0.01, 8.0, 200)
    core = np.exp(-2.0 * r)
    s = np.exp(-r) * (1.0 - np.exp(-6 * r))
    f = r ** 3 * np.exp(-r)
    (directory / 'ld1.wfc').write_text(
        '#            r                  6S                  4F                  5F\n'
        + ''.join(f'{x:20.12f}{a:20.12f}{c:20.12f}{b:20.12f}\n'
                  for x, a, b, c in zip(r, s, f, core)), encoding='utf-8')
    (directory / 'ld1ps.wfc').write_text(
        '#            r                  6S                  5F'
        '                  6S                  5F\n'
        + ''.join(f'{x:20.12f}{a * 0.8 if x < 1.5 else a:20.12f}'
                  f'{b * 0.8 if x < 1.5 else b:20.12f}{-a:20.12f}{-b:20.12f}\n'
                  for x, a, b in zip(r, s, f)), encoding='utf-8')

    beta_labels = ['6S', '7S', '6P', '6D', '5F', '6P']
    beta = np.stack([np.exp(-(1.0 + 0.3 * k) * r) * (1.0 + 0.1 * k)
                     for k in range(len(beta_labels))], axis=1)
    (directory / 'U.beta').write_text(
        '#            r' + ''.join(f'{label:>20s}' for label in beta_labels) + '\n'
        + ''.join(f'{x:20.12f}' + ''.join(f'{v:20.12f}' for v in row) + '\n'
                  for x, row in zip(r, beta)), encoding='utf-8')

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
        getattr(report_module, which)(plt, directory, tmp_path / which, formats=('png',),
                                      **kwargs)
        figure = figures[which]
        renderer = figure.canvas.get_renderer()
        figure.canvas.draw()
        return figure, renderer, Text

    return render, plt


def outside_tight_bbox(figure, renderer, Text=None, tolerance: float = 0.06):
    """Decorations outside the box ``bbox_inches='tight'`` would keep (inches)."""
    tight = figure.get_tightbbox(renderer)
    to_inches = figure.dpi_scale_trans.inverted()
    candidates = []
    for axes in figure.axes:
        box = axes.get_tightbbox(renderer)
        if box is not None:
            candidates.append(box.transformed(to_inches))
        for child in axes.child_axes:
            box = child.get_tightbbox(renderer)
            if box is not None:
                candidates.append(box.transformed(to_inches))
    for text in figure.texts:
        if text.get_visible() and text.get_text().strip():
            candidates.append(text.get_window_extent(renderer).transformed(to_inches))
    return [box for box in candidates
            if (box.x0 < tight.x0 - tolerance or box.y0 < tight.y0 - tolerance
                or box.x1 > tight.x1 + tolerance or box.y1 > tight.y1 + tolerance)]


def curves(axis):
    """The data curves of a panel (a reference line has only two distinct abscissae)."""
    return [line for line in axis.get_lines() if len(set(line.get_xdata())) > 2]


def is_solid(line) -> bool:
    return line.get_linestyle() in ('-', 'solid')


def is_dashed(line) -> bool:
    """``ls=(0, (4, 2))`` is a tuple, not the ``'--'`` shorthand matplotlib also accepts."""
    style = line.get_linestyle()
    return isinstance(style, tuple) or style in ('--', 'dashed', '-.')


# ---------------------------------------------------------------- house style

def test_the_style_is_the_manuscript_one(captured):
    render, plt = captured
    render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    params = matplotlib.rcParams
    assert params['font.size'] == pytest.approx(8.5)
    assert params['axes.linewidth'] == pytest.approx(0.8)
    assert params['xtick.direction'] == 'out'
    assert params['svg.fonttype'] == 'none'
    assert params['pdf.fonttype'] == 42


def test_the_figures_carry_no_title(captured):
    render, _plt = captured
    for which, extra in (('figure_radial', {'rcuts': {'6S': 1.5, '5F': 1.3}}),
                         ('figure_projectors', {'rcuts': {'6S': 1.5}}),
                         ('figure_logderiv', {'window': (-3.0, 3.0)})):
        figure, _renderer, _Text = render(which, **extra)
        for axis in figure.axes:
            assert axis.get_title() == '', which
        assert not figure.texts, which


def test_every_figure_is_seven_inches_wide(captured):
    render, _plt = captured
    for which, extra in (('figure_radial', {'rcuts': {'6S': 1.5}}),
                         ('figure_projectors', {'rcuts': {'6S': 1.5}}),
                         ('figure_logderiv', {})):
        figure, _r, _T = render(which, **extra)
        assert figure.get_size_inches()[0] == pytest.approx(7.0), which
        assert outside_tight_bbox(figure, _r) == []


# ------------------------------------------------------------------- radial

def test_radial_is_a_two_panel_figure_with_the_difference_underneath(captured):
    render, _plt = captured
    figure, renderer, Text = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    assert len(figure.axes) == 2
    main, panel = figure.axes
    assert main.get_xscale() == 'log'
    assert main.get_yscale() == 'linear'
    assert panel.get_yscale() == 'log'            # |du|/max|u|
    assert panel.get_xlabel() == 'r (bohr)'
    assert main.get_xlabel() == ''                # shared abscissa
    assert outside_tight_bbox(figure, renderer, Text) == []


def test_radial_uses_solid_all_electron_and_dashed_pseudo(captured):
    render, _plt = captured
    figure, _r, _T = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    main = figure.axes[0]
    solid = [line for line in curves(main) if is_solid(line)]
    dashed = [line for line in curves(main) if is_dashed(line)]
    assert len(solid) == 2 and len(dashed) == 2
    for line in solid:
        assert line.get_linewidth() == pytest.approx(1.5)
    for line in dashed:
        assert line.get_linewidth() == pytest.approx(1.0)
    vertical = [line for line in main.get_lines() if len(set(line.get_xdata())) == 1]
    assert len(vertical) == 2                     # one per pseudization radius


def test_radial_has_an_angstrom_axis_and_a_shaded_core(captured):
    render, _plt = captured
    figure, _r, _T = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    main = figure.axes[0]
    secondary = main.child_axes
    assert len(secondary) == 1
    assert 'AA' in secondary[0].get_xlabel() or 'mathrm' in secondary[0].get_xlabel()
    assert secondary[0].get_xticks().size >= 1
    assert main.patches, 'the core region should be shaded'


def test_radial_skips_the_all_electron_core_states(captured):
    render, _plt = captured
    figure, _r, _T = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3})
    labels = [text.get_text() for text in figure.axes[0].get_legend().get_texts()]
    assert all('4F' not in label for label in labels)
    assert any('6S' in label for label in labels)


def test_the_diagnostic_variant_spans_four_decades(captured):
    render, _plt = captured
    figure, _r, _T = render('figure_radial', rcuts={'6S': 1.5, '5F': 1.3},
                            xscale='log', xlim=(1e-2, 1e2), pseudo='both')
    main = figure.axes[0]
    assert main.get_xlim() == pytest.approx((0.01, 100.0))
    assert [line for line in curves(main) if line.get_marker() == 'o']


# ---------------------------------------------------------------- projectors

def test_projectors_have_an_absolute_and_a_normalized_panel(captured):
    render, _plt = captured
    figure, renderer, Text = render('figure_projectors', rcuts={'6S': 1.5})
    assert len(figure.axes) == 2
    absolute, relative = figure.axes
    assert absolute.get_yscale() == 'log'
    assert relative.get_yscale() == 'linear'
    assert len(curves(absolute)) == 6
    # only one channel carries two projectors (6P and the unoccupied 6P'), so exactly one
    # curve is dashed -- the other four channels have a single projector each
    assert sum(1 for line in curves(absolute) if is_dashed(line)) == 1
    labels = [text.get_text() for text in absolute.get_legend().get_texts()]
    assert len(labels) == 6
    assert any('6P' in label for label in labels)
    assert outside_tight_bbox(figure, renderer, Text) == []


# --------------------------------------------------------------- logderiv

def test_logderiv_is_linear_with_open_circles_and_a_pole_table(captured):
    render, _plt = captured
    figure, renderer, Text = render('figure_logderiv', window=(-3.0, 3.0),
                                    radius=5.0)
    assert len(figure.axes) == 2
    main, table = figure.axes
    assert main.get_yscale() == 'linear'
    assert main.get_xscale() == 'linear'
    circles = [line for line in main.get_lines() if line.get_marker() == 'o']
    assert len(circles) == 2                      # one per channel, open circles
    for line in circles:
        assert line.get_markerfacecolor() in ('none', 'None', (0, 0, 0, 0))
    assert len(curves(main)) == 4                 # 2 AE solid + 2 PS circle series
    # the evaluation radius belongs on the figure
    assert 'r_d' in main.get_xlabel() or 'r_d' in main.get_xlabel().replace('$', '')
    # the pole table panel
    assert table.axison is False
    text = table.texts[0].get_text()
    assert 'AE poles' in text and 'PS poles' in text
    assert outside_tight_bbox(figure, renderer, Text) == []


def test_logderiv_can_still_be_logarithmic(captured):
    render, _plt = captured
    figure, _r, _T = render('figure_logderiv', window=(-3.0, 3.0), log_ordinate=True)
    main = figure.axes[0]
    assert main.get_yscale() == 'log'
    assert main.get_ylabel().startswith('$|L_l(E)|$')


# ------------------------------------------------------------------ pipeline

def test_the_shipped_configuration_is_what_we_think_it_is(monkeypatch, tmp_path):
    """Guard the arguments make_figures passes, not a convenient hand-made call."""
    directory = write_run(tmp_path / 'run')
    summary = {
        'generation': {'wavefunctions': {'6S': {'rcut': 1.5}, '5F': {'rcut': 1.3}}},
        'criteria': {'window': (-3.0, 3.0), 'pole_threshold': 20.0,
                     'production_rlderiv': 5.0},
    }
    captured = {}

    def fake_save(figure, out_base, plt, formats):
        captured[Path(out_base).name] = figure
        return out_base

    monkeypatch.setattr(report_module, '_save', fake_save)
    written = report_module.make_figures(summary, directory, tmp_path / 'data')

    assert sorted(captured) == ['fig_logderiv', 'fig_projectors', 'fig_radial',
                                'fig_radial_logx']
    assert sorted(Path(path).stem for path in written) == sorted(captured)
    radial = captured['fig_radial'].axes[0]
    assert (radial.get_xscale(), radial.get_yscale()) == ('log', 'linear')
    assert radial.get_xlim() == pytest.approx((0.02, 6.0))
    assert len(captured['fig_radial'].axes) == 2          # the difference panel
    logx = captured['fig_radial_logx'].axes[0]
    assert (logx.get_xscale(), logx.get_yscale()) == ('log', 'linear')
    assert logx.get_xlim() == pytest.approx((0.01, 100.0))
    derivative = captured['fig_logderiv'].axes[0]
    assert derivative.get_yscale() == 'linear'
    assert derivative.get_xscale() == 'linear'
    projectors = captured['fig_projectors']
    assert projectors.axes[0].get_yscale() == 'log'
    assert len(projectors.axes) == 2


def test_export_paper_figures_backs_up_and_stays_linear(tmp_path, monkeypatch):
    directory = write_run(tmp_path / 'run')
    paper = tmp_path / 'figs'
    paper.mkdir()
    (paper / 'fig3_1_radial.png').write_text('old radial', encoding='utf-8')
    (paper / 'fig3_3_logderiv.png').write_text('old logderiv', encoding='utf-8')

    captured = {}
    real = report_module._save

    def fake_save(figure, out_base, plt, formats):
        captured[Path(out_base).name] = figure
        return real(figure, out_base, plt, formats)

    monkeypatch.setattr(report_module, '_save', fake_save)
    written = report_module.export_paper_figures(
        directory, paper, rcuts={'6S': 1.5, '5F': 1.3}, radius=5.0, formats=('png',))
    assert len(written) == 2
    for name in ('fig3_1_radial', 'fig3_3_logderiv'):
        assert Path(f'{paper}/{name}.png.bak').read_text(encoding='utf-8').startswith('old')
        assert Path(f'{paper}/{name}.png').stat().st_size > 1000
    assert captured['fig3_3_logderiv'].axes[0].get_yscale() == 'linear'
    assert captured['fig3_1_radial'].axes[0].get_xscale() == 'log'

    # re-exporting must NOT overwrite the backup with the previous output
    report_module.export_paper_figures(directory, paper, rcuts={'6S': 1.5},
                                       formats=('png',))
    assert Path(f'{paper}/fig3_1_radial.png.bak').read_text(encoding='utf-8') == 'old radial'
