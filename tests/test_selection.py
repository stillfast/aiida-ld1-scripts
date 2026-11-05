"""The rlderiv selection rule: the pole verdict gates the median."""

from __future__ import annotations

from aiida_ld1_scripts.collect import _select


def entry(rlderiv, median, coincide, bound=1e-6):
    return {'rlderiv': rlderiv, 'pooled_median': median,
            'poles_all_coincide': coincide, 'pk': int(rlderiv * 10),
            'bound_median': bound, 'per_channel': {}, 'poles': {}}


def test_a_smaller_median_does_not_beat_a_failed_pole_check():
    """The case that was wrong: r=8 has the smallest median but fails the poles."""
    summary = _select([entry(4, 8.8e-4, False), entry(5, 1.0e-3, True),
                       entry(7, 2.4e-4, True), entry(8, 1.4e-4, False)])
    assert summary['rlderiv'] == 7.0
    assert summary['pk'] == 70
    assert summary['passing'] == [5.0, 7.0]
    assert summary['best_by_median']['rlderiv'] == 8.0
    assert summary['best_by_median']['pole_gated'] is False


def test_no_point_passing_means_no_pick():
    summary = _select([entry(4, 8.8e-4, False), entry(8, 1.4e-4, False)])
    assert summary['rlderiv'] is None
    assert summary['pole_gated'] is False
    assert summary['best_by_median']['rlderiv'] == 8.0


def test_empty_scan_selects_nothing():
    assert _select([]) is None
