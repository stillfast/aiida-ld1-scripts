"""Reading the ``ld1.x`` artifacts and measuring how well the pseudopotential matches.

These are plain functions over the files ``ld1.x`` writes -- ``ld1.dlog``,
``ld1ps.dlog``, ``ld1.wfc``, ``ld1ps.wfc`` -- so they can be tested without AiiDA and
reused on an old run directory.

The ``.dlog`` file is the logarithmic derivative ``L_l(E) = d ln u_l / d ln r`` at the
reference radius, computed for each angular momentum on the energy grid of ``&input``
(``eminld``, ``emaxld``, ``deld``)::

         -3.000000000000      1.618190149386      1.645891606330   ...
         (energy / Ry)        L(l=0)              L(l=1)           ...

Column ``i + 1`` is ``l = i``: the columns are indexed by angular momentum, not by the
orbital card.  A channel with no projector shows the *local potential* alone, and the
curves diverge at the energies of the bound states ("poles").

Two numbers are used to judge a set of parameters:

``median |dL|``
    the median of ``|L_PS - L_AE|`` over the window -- no point is skipped, but the
    median is what makes it robust to the few diverging points at a pole.  Lower is
    better.  The ranking uses the median over **all channels concatenated**.
``poles``
    the bound-state energies of each channel, located by interpolating the sign change
    of ``L``; AE and PS must have the same number, each within ``tolerance``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: |L| above which the logarithmic derivative is at a pole
POLE_THRESHOLD = 20.0


def read_dlog(path: str | Path) -> Dict[str, np.ndarray]:
    """``{'energy': E, 'l=0': L, ...}`` from one ``.dlog`` file."""
    data = np.loadtxt(path)
    if data.ndim == 1:                       # a single grid point
        data = data[None, :]
    columns = {'energy': data[:, 0]}
    for index in range(1, data.shape[1]):
        columns[f'l={index - 1}'] = data[:, index]
    return columns


def split_wfc_halves(labels: Sequence[str]) -> Tuple[List[str], List[str]]:
    """Split an ``ld1ps.wfc`` header into its ``(pseudo, all-electron)`` halves.

    ``write_resultsps.f90`` writes both sets of orbitals in one file: the first half
    holds the **pseudo** orbitals and the second half the matching **all-electron**
    ones (phase-aligned).  The header therefore repeats itself::

        # r    5F  6S  6P  6D  7S   5F  6S  6P  6D  7S

    Treating the repeat as an ordinary duplicate and building a ``{label: column}``
    dict silently keeps the *last* copy -- i.e. the all-electron orbital where the
    pseudo one is wanted, which is how a first version of the figure ended up
    comparing AE with AE.
    """
    labels = list(labels)
    half = len(labels) // 2
    if labels and len(labels) % 2 == 0 and labels[:half] == labels[half:]:
        return labels[:half], labels[half:]
    return labels, []


def read_wfc_halves(path: str | Path):
    """``(labels, r, pseudo, all_electron)``; ``all_electron`` is ``None`` for ld1.wfc."""
    labels, grid, values = read_wfc(path)
    pseudo_labels, all_labels = split_wfc_halves(labels)
    if not all_labels:
        return pseudo_labels, grid, values, None
    half = len(pseudo_labels)
    return pseudo_labels, grid, values[:, :half], values[:, half:]


def read_wfc(path: str | Path) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """``(labels, r, u)`` from one ``.wfc`` file (the first column is ``r``)."""
    labels: List[str] = []
    rows = []
    with open(path, encoding='utf-8', errors='replace') as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                labels = stripped.lstrip('#').split()[1:]     # drop the 'r' column
                continue
            rows.append([float(value) for value in stripped.split()])
    data = np.array(rows)
    return labels, data[:, 0], data[:, 1:]


def align_phase(ae: np.ndarray, ps: np.ndarray) -> Tuple[np.ndarray, bool]:
    """Return the all-electron array in the phase of the pseudo one.

    ``u(r)`` has no defined sign, and the all-electron orbital written by the AE run
    (``ld1.wfc``) need not share the phase of the pseudo orbital.  QE itself solves
    this for the *second* half of ``ld1ps.wfc`` (``write_resultsps.f90`` writes the
    phase-aligned all-electron orbital), but comparing with ``ld1.wfc`` means doing it
    here: if ``|AE + PS| < |AE - PS|`` the two are in opposite phase.
    """
    if np.max(np.abs(ae + ps)) < np.max(np.abs(ae - ps)):
        return -ae, True
    return ae, False


def wavefunction_agreement(ae_path: str | Path, ps_path: str | Path,
                           cutoffs: Dict[str, float] | None = None,
                           default_cutoff: float = 2.0) -> Dict[str, Dict[str, object]]:
    """How well the pseudo orbitals reproduce the all-electron ones.

    Inside a state's cut-off radius the two are *meant* to differ (that is what the
    pseudization does); **outside** it they must agree, up to the phase.  So the
    interesting number is ``max |u_AE - u_PS|`` for ``r > r_c`` -- a large value there
    means the tail of the pseudo orbital is wrong, which is what spoils a
    logarithmic derivative and the transferability.

    ``cutoffs`` maps a state label to its ``rcut`` (from the recipe's state cards).
    """
    ae_labels, ae_r, ae_values = read_wfc(ae_path)
    ps_labels, ps_r, ps_values, _aligned = read_wfc_halves(ps_path)
    r = ae_r
    result: Dict[str, Dict[str, object]] = {}
    for index, label in enumerate(ps_labels):
        if label not in ae_labels:
            continue
        ae = ae_values[:, ae_labels.index(label)]
        ps = ps_values[:, index]
        ae, flipped = align_phase(ae, ps)
        cutoff = float((cutoffs or {}).get(label, default_cutoff))
        inside = r < cutoff
        outside = ~inside & (r <= 10.0)         # the tail that matters
        result[label] = {
            'rcut': cutoff,
            'flipped': flipped,
            'max_inside': float(np.max(np.abs(ae[inside] - ps[inside]))) if inside.any() else None,
            'max_outside': float(np.max(np.abs(ae[outside] - ps[outside]))) if outside.any() else None,
        }
    return result


#: the spectroscopic letters, ``s p d f g h i k`` -> l = 0..7
LETTERS = 'spdfghik'


def state_l(label: str) -> int:
    """Angular momentum of a state label: ``6S`` -> 0, ``5F`` -> 3, ``6P'`` -> 1.

    The trailing letters of a label may carry a prime to disambiguate a repeated channel
    (``6P`` / ``6P'``); only the first alphabetic character matters.  ``-1`` means "no
    letter", which never matches a channel.
    """
    letters = ''.join(character for character in str(label) if character.isalpha())
    if not letters:
        return -1
    try:
        return LETTERS.index(letters[0].lower())
    except ValueError:
        return -1


# --------------------------------------------------------------------- poles


def find_poles(energies: np.ndarray, values: np.ndarray,
               threshold: float = POLE_THRESHOLD) -> List[float]:
    """Energies where the logarithmic derivative diverges.

    A pole is a sign flip between adjacent points with ``|L| > threshold`` on at least
    one side; ``20`` rather than ``1e3`` because a 0.02 Ry grid often just misses the
    peak.  When *both* sides exceed the threshold the position comes from the
    divergence model ``L ~ c/(E - E_p)``, i.e.
    ``E_p = E_k + dE·|L_k+1| / (|L_k| + |L_k+1|)``; otherwise the cell midpoint is used.
    Without this the two poles of AE and PS would both snap to the midpoint and a real
    displacement would be reported as zero.
    """
    poles: List[float] = []
    magnitudes = np.abs(values)
    for index in range(len(values) - 1):
        if values[index] * values[index + 1] >= 0.0:
            continue
        here, there = magnitudes[index], magnitudes[index + 1]
        if max(here, there) <= threshold:
            continue
        if min(here, there) > threshold:
            span = float(energies[index + 1]) - float(energies[index])
            poles.append(float(energies[index]) + span * there / (here + there))
        else:
            poles.append(0.5 * (float(energies[index]) + float(energies[index + 1])))
    return poles


def find_zeros(energies: np.ndarray, values: np.ndarray) -> List[float]:
    """Energies where ``L = 0`` (a dip on a logarithmic axis), interpolated."""
    zeros: List[float] = []
    for index in range(1, len(energies)):
        left, right = values[index - 1], values[index]
        if left == 0.0:
            zeros.append(float(energies[index - 1]))
        elif np.sign(left) != np.sign(right):
            span = right - left
            zeros.append(float(energies[index - 1] - left / span * (energies[index] - energies[index - 1])))
    return zeros


def match_poles(left: Sequence[float], right: Sequence[float],
                tolerance: float) -> Dict[str, object]:
    """Do two pole lists describe the same states?"""
    matched, unmatched = [], list(right)
    for position in left:
        closest = min(unmatched, key=lambda other: abs(other - position), default=None)
        if closest is not None and abs(closest - position) <= tolerance:
            matched.append((position, closest))
            unmatched.remove(closest)
    return {
        'n_ae': len(left),
        'n_ps': len(right),
        'matched': len(matched),
        'coincide': len(matched) == len(left) == len(right),
        'missing_in_ps': [position for position in left
                          if not any(abs(position - pair[0]) < 1e-12 for pair in matched)],
        'extra_in_ps': unmatched,
    }


# ------------------------------------------------------------------- metrics


def _differences(ae: Dict[str, np.ndarray], ps: Dict[str, np.ndarray],
                 window: Tuple[float, float]):
    """``(channel, |L_PS - L_AE|)`` inside the window, per channel."""
    energies = ae['energy']
    mask = (energies >= window[0]) & (energies <= window[1])
    for channel in ae:
        if channel == 'energy':
            continue
        yield channel, np.abs(ps[channel][mask] - ae[channel][mask])


def log_derivative_error(ae: Dict[str, np.ndarray], ps: Dict[str, np.ndarray],
                         window: Tuple[float, float] = (-3.0, 3.0)
                         ) -> Dict[str, Dict[str, float]]:
    """Per channel, the statistics of ``|L_PS - L_AE|`` inside the window.

    **No point is excluded.**  Near a pole both curves diverge, so a few points carry
    huge differences; the *median* is what ranks parameter sets because it is robust to
    them, while the *mean* is not -- do not use ``mean`` to compare points.

    These are exactly the numbers of ``pp-gen``'s ``report.md`` (``--scope scanned``),
    which is the analysis this pipeline reproduces.
    """
    result: Dict[str, Dict[str, float]] = {}
    for channel, difference in _differences(ae, ps, window):
        if difference.size == 0:
            continue
        result[channel] = {
            'median': float(np.median(difference)),
            'mean': float(np.mean(difference)),
            'min': float(np.min(difference)),
            'max': float(np.max(difference)),
            'n_points': int(difference.size),
        }
    return result


def pooled_median(ae: Dict[str, np.ndarray], ps: Dict[str, np.ndarray],
                  window: Tuple[float, float] = (-3.0, 3.0)) -> float:
    """The ranking criterion: median of ``|dL|`` over **all channels concatenated**.

    Pooling the raw differences (rather than taking the median of the per-channel
    medians) is what ``pp-gen``'s ``--scope all`` does, and it is the number the
    ``rlderiv`` scan is ranked by.
    """
    arrays = [difference for _channel, difference in _differences(ae, ps, window)
              if difference.size]
    if not arrays:
        return float('nan')
    return float(np.median(np.concatenate(arrays)))


def compare_run(directory: str | Path, window: Tuple[float, float] = (-3.0, 3.0),
                pole_threshold: float = POLE_THRESHOLD,
                pole_tolerance: float = 0.04) -> Dict[str, object]:
    """Everything measured on one ``ld1.x`` run directory."""
    directory = Path(directory)
    ae = read_dlog(directory / 'ld1.dlog')
    ps = read_dlog(directory / 'ld1ps.dlog')
    per_channel = log_derivative_error(ae, ps, window=window)

    # the pole comparison is done inside the same window as the statistics: outside it
    # the two files may not even share a grid, and a feature there is not what the scan
    # is judging
    energies = ae['energy']
    mask = (energies >= window[0]) & (energies <= window[1])
    poles: Dict[str, Dict[str, object]] = {}
    for channel in ae:
        if channel == 'energy':
            continue
        poles[channel] = match_poles(
            find_poles(energies[mask], ae[channel][mask], pole_threshold),
            find_poles(ps['energy'][mask], ps[channel][mask], pole_threshold),
            pole_tolerance)
    return {
        'directory': str(directory),
        'per_channel': per_channel,
        'pooled_median': pooled_median(ae, ps, window=window),
        'poles': poles,
        'poles_all_coincide': all(entry['coincide'] for entry in poles.values()),
    }


def select_rlderiv(scanned: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Pick an ``rlderiv`` in **two tiers**: pole coincidence first, median second.

    Ranking by the pooled median alone was wrong, and visibly so: rlderiv = 8 has the
    smallest median (1.38e-4) while its own pole check *fails* -- exactly the criterion
    the report calls the ghost criterion.  A pseudopotential that adds or moves a bound
    state is not made acceptable by a small median, so the pole verdict gates the choice.

    ``best_by_median`` is kept in the summary so the difference between the two rules is
    visible instead of silently resolved.
    """
    if not scanned:
        return None
    by_median = min(scanned, key=lambda entry: entry['pooled_median'])
    passing = [entry for entry in scanned if entry['poles_all_coincide']]
    chosen = min(passing, key=lambda entry: entry['pooled_median']) if passing else None
    return {
        'rlderiv': chosen['rlderiv'] if chosen else None,
        'pooled_median': chosen['pooled_median'] if chosen else None,
        'pk': chosen['pk'] if chosen else None,
        'rule': 'poles must coincide, then smallest pooled median',
        'pole_gated': bool(passing),
        'passing': [entry['rlderiv'] for entry in passing],
        'best_by_median': {
            'rlderiv': by_median['rlderiv'],
            'pooled_median': by_median['pooled_median'],
            'pole_gated': by_median['poles_all_coincide'],
        },
    }


def read_projectors(directory: str | Path
                    ) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """The radial projector functions ``beta(r)`` of a run.

    ``&inputp.file_beta`` gets the same layout as ``ld1.wfc``: a ``#`` header naming one
    column per state, then ``r`` and one column per projector.  The state labels repeat
    when a channel carries two projectors -- this recipe has ``6P`` twice, the second one
    with ``occ = 0`` -- so the repeats are disambiguated as ``6P``, ``6P'``, ... rather
    than left as two identical legend entries.
    """
    candidates = sorted(Path(directory).glob('*.beta'))
    if not candidates:
        raise FileNotFoundError(f'no *.beta file in {directory}')
    labels, grid, values = read_wfc(candidates[0])

    counted: Dict[str, int] = {}
    unique: List[str] = []
    for label in labels:
        counted[label] = counted.get(label, 0) + 1
        unique.append(label if counted[label] == 1 else f'{label}{chr(0x2032) * (counted[label] - 1)}')
    return unique, grid, values


def pole_table(ae: Dict[str, np.ndarray], ps: Dict[str, np.ndarray],
               window: Tuple[float, float] = (-3.0, 3.0),
               threshold: float = POLE_THRESHOLD,
               tolerance: float = 0.04) -> List[Dict[str, object]]:
    """Per channel: the AE and PS pole energies and the displacement of each pair.

    This is the table the log-derivative figure shows and the report prints, so that the
    ghost criterion is a number rather than something read off a curve.
    """
    energies = ae['energy']
    mask = (energies >= window[0]) & (energies <= window[1])
    rows: List[Dict[str, object]] = []
    for channel in ae:
        if channel == 'energy':
            continue
        left = find_poles(energies[mask], ae[channel][mask], threshold)
        right = find_poles(ps['energy'][mask], ps[channel][mask], threshold)
        matched = match_poles(left, right, tolerance)
        pairs = []
        remaining = list(right)
        for position in left:
            nearest = min(remaining, key=lambda other: abs(other - position), default=None)
            if nearest is not None and abs(nearest - position) <= tolerance:
                pairs.append((position, nearest, abs(nearest - position)))
                remaining.remove(nearest)
        rows.append({
            'channel': channel,
            'ae': left,
            'ps': right,
            'pairs': pairs,
            'coincide': matched['coincide'],
            'max_delta': max((delta for _a, _b, delta in pairs), default=0.0),
        })
    return rows


# --------------------------------------------------------- transferability test


def transferability(sections: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    """Delta E for each test configuration relative to the first (the reference).

    ``Delta E = (E_ae - E_ae,ref) - (E_ps - E_ps,ref)``, i.e. the error the
    pseudopotential makes in the energy *difference* between two configurations --
    the quantity that matters for transferability.  ``delta_mry`` is the same in mRy.

    The sign is the one ``ld1.x`` prints in its own ``dEtot_ps ... Delta E=`` column, so
    a value computed here can be checked against ``ld1.test`` directly (and matches the
    ``delta`` column of the older ``pp-gen`` reports).
    """
    if not sections:
        return []
    reference = sections[0]
    rows: List[Dict[str, object]] = []
    for index, section in enumerate(sections, 1):
        delta_ae = (section['etot'] - reference['etot']
                    if section['etot'] is not None and reference['etot'] is not None else None)
        delta_ps = (section['etotps'] - reference['etotps']
                    if section['etotps'] is not None and reference['etotps'] is not None else None)
        delta = (delta_ae - delta_ps
                 if delta_ae is not None and delta_ps is not None else None)
        rows.append({
            'index': index,
            'delta_ae_ry': delta_ae,
            'delta_ps_ry': delta_ps,
            'delta_ry': delta,
            'delta_mry': None if delta is None else 1000.0 * delta,
        })
    return rows


def states_table(sections: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    """The AE/PS level table of the first section, one row per channel."""
    if not sections:
        return []
    return [dict(state) for state in sections[0]['states']]


def pseudo_metadata(path: str | Path) -> Dict[str, object]:
    """A few ``<PP_HEADER>`` fields of the generated UPF (for the report)."""
    import re

    text = Path(path).read_text(encoding='utf-8', errors='replace')
    match = re.search(r'<PP_HEADER(.*?)/?>', text, re.DOTALL)
    if not match:
        return {}
    body = match.group(1)
    fields: Dict[str, object] = {}
    for key in ('element', 'pseudo_type', 'relativistic', 'functional', 'z_valence'):
        found = re.search(rf'{key}\s*=\s*"([^"]*)"', body)
        if found:
            fields[key] = found.group(1).strip()
    return fields
