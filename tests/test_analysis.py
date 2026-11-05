"""The analysis must reproduce the numbers of the reference ``pp-gen`` report.

The fixtures in ``tests/data/rlderiv_scan`` are the ten ``.dlog`` files of the five
``rlderiv`` runs of the 6p25 work (``uranium-ppgen/ai_v1/6p25/test_rlderiv``); the
expected values are the rows of that work's ``report.md``.  Agreement here is what
"reproduced with the new pipeline" means for the analysis side -- the AiiDA side is
covered by running the pipeline itself.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from aiida_ld1_scripts import analysis

DATA = Path(__file__).parent / 'data' / 'rlderiv_scan'

#: rows of the reference report: per-channel median |dL| and the pooled median
REFERENCE = {
    4.0: (0.00088539, 0.00014845, 0.002797, 0.0038712, 0.00087752),
    5.0: (0.0023266, 8.1214e-05, 0.0024538, 0.0010934, 0.0010307),
    6.0: (0.0006382, 0.00010648, 0.00055293, 0.00041187, 0.00056483),
    7.0: (0.00028473, 0.00028562, 0.00020432, 0.00017836, 0.00023586),
    8.0: (0.00016667, 0.00014074, 0.00011317, 7.6425e-05, 0.00013758),
}


def directory(value: float) -> Path:
    return DATA / f'rlderiv_{int(value)}p0'


@pytest.mark.parametrize('value', sorted(REFERENCE))
def test_pooled_median_matches_the_reference_report(value):
    measured = analysis.compare_run(directory(value), window=(-3.0, 3.0),
                                    pole_tolerance=0.04)
    expected = REFERENCE[value]
    assert measured['pooled_median'] == pytest.approx(expected[4], rel=1e-4)
    for channel in range(4):
        assert measured['per_channel'][f'l={channel}']['median'] == \
            pytest.approx(expected[channel], rel=1e-4)


def test_the_scan_ranks_rlderiv_8_best():
    """The reference report picks 8.0 a.u. (0.00013758, the smallest)."""
    pooled = {}
    for value in REFERENCE:
        pooled[value] = analysis.compare_run(directory(value), window=(-3.0, 3.0),
                                             pole_tolerance=0.04)['pooled_median']
    assert min(pooled, key=pooled.get) == 8.0
    # the criterion is not monotone in rlderiv -- 5.0 is worse than 4.0, which is why the
    # scan has to measure every point instead of extrapolating
    assert pooled[5.0] > pooled[4.0] > pooled[6.0] > pooled[7.0] > pooled[8.0]


def test_poles_and_their_verdicts_match_the_reference_table():
    """l=0 and l=1 coincide; l=2 and l=3 are displaced by more than the tolerance."""
    measured = analysis.compare_run(directory(8.0), window=(-3.0, 3.0),
                                    pole_tolerance=0.04)
    poles = measured['poles']
    for channel in ('l=0', 'l=1'):
        assert poles[channel]['coincide'] is True, channel
    for channel in ('l=2', 'l=3'):
        assert poles[channel]['coincide'] is False, channel
        # the reference lists a displaced pole: l=2 +2.649 vs +2.693, l=3 +2.392 vs +2.345
        assert poles[channel]['extra_in_ps'], channel
    assert measured['poles_all_coincide'] is False


def test_pole_positions_are_interpolated_not_quantised():
    """The reference table has l=0 at +0.430 and +1.535 Ry, not on the 0.02 grid."""
    ae = analysis.read_dlog(directory(8.0) / 'ld1.dlog')
    poles = analysis.find_poles(ae['energy'], ae['l=0'], threshold=20.0)
    assert len(poles) == 2
    assert poles[0] == pytest.approx(0.430, abs=2e-3)
    assert poles[1] == pytest.approx(1.535, abs=2e-3)
    # off-grid: a 0.02 Ry grid cannot place a pole at 0.430 exactly
    assert abs((poles[0] / 0.02) - round(poles[0] / 0.02)) > 1e-3


def test_the_median_is_used_because_the_mean_is_not_robust():
    ae = analysis.read_dlog(directory(4.0) / 'ld1.dlog')
    ps = analysis.read_dlog(directory(4.0) / 'ld1ps.dlog')
    per_channel = analysis.log_derivative_error(ae, ps, window=(-3.0, 3.0))
    for channel, entry in per_channel.items():
        assert entry['max'] > 10 * entry['median'], channel      # poles blow up the max
        assert entry['mean'] > entry['median'], channel


# ---------------------------------------------------------------- transferability

def test_transferability_delta_e_sign_and_values():
    """Same sign as ld1.x's own Delta E column; reference 6p25 values in mRy."""
    reference = {'etot': -55958.130614, 'etotps': -181.732843}
    shifts = [(0.0, 0.0), (0.396703, 0.397216), (0.492878, 0.493174),
              (0.145115, 0.145599), (0.604082, 0.605156)]
    sections = [{'etot': reference['etot'] + a, 'etotps': reference['etotps'] + b,
                 'states': []} for a, b in shifts]
    rows = analysis.transferability(sections)
    assert len(rows) == 5
    assert rows[0]['delta_ry'] == pytest.approx(0.0)
    assert [round(row['delta_mry'], 1) for row in rows[1:]] == \
        [-0.5, -0.3, -0.5, -1.1]          # -0.512, -0.296, -0.483, -1.074, rounded


def test_transferability_of_an_empty_run_is_empty():
    assert analysis.transferability([]) == []


# ---------------------------------------------------------------------- formats


def test_read_wfc_returns_labels_and_values():
    path = DATA / 'synthetic.wfc'
    path.write_text('#            r                  5F                  6S\n'
                    '      0.001    0.0     1.0\n'
                    '      0.002    0.5     0.5\n', encoding='utf-8')
    try:
        labels, r, values = analysis.read_wfc(path)
        assert labels == ['5F', '6S']
        assert list(r) == pytest.approx([0.001, 0.002])
        assert values.shape == (2, 2)
    finally:
        path.unlink()


# ------------------------------------------------------------------ wavefunctions

def test_split_wfc_halves_detects_the_repeated_header():
    """ld1ps.wfc holds PS then AE, so the header repeats itself."""
    labels = ['5F', '6S', '6P', '6D', '7S']
    pseudo, aligned = analysis.split_wfc_halves(labels + labels)
    assert pseudo == labels and aligned == labels
    # a normal ld1.wfc (no repeat) stays whole
    assert analysis.split_wfc_halves(labels) == (labels, [])


def test_align_phase_flips_an_opposite_orbital():
    r = np.linspace(0.1, 5.0, 50)
    u = np.exp(-r)
    same, flipped = analysis.align_phase(u, u)
    assert flipped is False
    assert same is u or np.allclose(same, u)
    opposite, flipped = analysis.align_phase(u, -u)
    assert flipped is True
    assert np.allclose(opposite, -u)


def test_wavefunction_agreement_is_small_outside_the_cutoff(tmp_path):
    """PS and AE must differ inside rcut (by design) and agree outside it."""
    r = np.linspace(0.01, 6.0, 200)
    ae = np.exp(-r)                              # the "true" orbital
    ps = np.where(r < 1.5, 0.7 * np.exp(-r), np.exp(-r))   # altered only inside rcut
    header = '#            r                  6S\n'
    # the all-electron run wrote the orbital with the opposite overall sign, which is
    # what align_phase has to detect
    (tmp_path / 'ld1.wfc').write_text(
        header + ''.join(f'{x:20.12f}{-y:20.12f}\n' for x, y in zip(r, ae)), encoding='utf-8')
    # the pseudo file repeats the header: PS first, then the phase-aligned AE
    (tmp_path / 'ld1ps.wfc').write_text(
        '#            r                  6S                  6S\n'
        + ''.join(f'{x:20.12f}{p:20.12f}{a:20.12f}\n' for x, p, a in zip(r, ps, ae)),
        encoding='utf-8')
    measured = analysis.wavefunction_agreement(tmp_path / 'ld1.wfc',
                                               tmp_path / 'ld1ps.wfc',
                                               cutoffs={'6S': 1.5})
    entry = measured['6S']
    assert entry['rcut'] == 1.5
    assert entry['flipped'] is True              # ld1.wfc phase vs pseudo
    assert entry['max_inside'] > 0.1             # differs where it is allowed to
    assert entry['max_outside'] < 1e-9           # agrees in the tail


def test_wavefunction_agreement_skips_states_missing_from_the_pseudo_file(tmp_path):
    """ld1.wfc can hold core states (4F, 5D) that the pseudopotential has no orbital for."""
    r = np.linspace(0.01, 4.0, 20)
    (tmp_path / 'ld1.wfc').write_text(
        '#            r                  4F                  6S\n'
        + ''.join(f'{x:18.10f}{y:18.10f}{y:18.10f}\n' for x, y in zip(r, np.exp(-r))),
        encoding='utf-8')
    (tmp_path / 'ld1ps.wfc').write_text(
        '#            r                  6S\n'
        + ''.join(f'{x:18.10f}{y:18.10f}\n' for x, y in zip(r, np.exp(-r))),
        encoding='utf-8')
    measured = analysis.wavefunction_agreement(tmp_path / 'ld1.wfc',
                                               tmp_path / 'ld1ps.wfc',
                                               cutoffs={'6S': 1.5})
    assert set(measured) == {'6S'}
