"""Turning the AiiDA group into plain data files.

Everything downstream (the report, the figures) reads ``data/summary.json`` and the
CSVs this module writes, so the analysis can be re-run without touching the database
and the numbers in the repository can be checked against the archive.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from aiida import orm

from . import analysis, recipe as recipe_module

#: files copied out of a ``retrieved`` folder for the analysis
#: ``ld1.x`` names the beta file after ``&inputp.file_beta`` (this recipe: ``U.beta``), so
#: it is picked up by suffix rather than by the literal name it is given here
RETRIEVED_FILES = ('ld1.dlog', 'ld1ps.dlog', 'ld1.wfc', 'ld1ps.wfc', 'ld1.test', 'ld1.pwe')


def _finished(group_label: str, role: str) -> List[Any]:
    from . import workflow

    return workflow.already_submitted(group_label, role)


def newest(group_label: str, role: str):
    nodes = _finished(group_label, role)
    if not nodes:
        raise SystemExit(f'no finished {role!r} calculation in group {group_label!r}')
    return nodes[-1]


def dump_retrieved(node, target: Path) -> Path:
    """Copy the analysis-relevant files of one calculation into ``target``."""
    if 'retrieved' not in node.outputs:
        raise SystemExit(f'{node.pk} has no retrieved folder')
    retrieved = node.outputs['retrieved']
    target.mkdir(parents=True, exist_ok=True)
    available = set(retrieved.list_object_names())
    names = [name for name in RETRIEVED_FILES if name in available]
    names += [name for name in sorted(available) if name.endswith('.beta')]
    for name in dict.fromkeys(names):
        (target / name).write_bytes(retrieved.get_object_content(name, mode='rb'))
    return target


def settings_of(recipe: Dict[str, Any]) -> Dict[str, Any]:
    return recipe_module.analysis_settings(recipe)


def collect(group_label: str = recipe_module.DEFAULT_GROUP,
            recipe_path: Optional[str] = None) -> Dict[str, Any]:
    """Measure every role in the group; returns the ``summary`` dictionary."""
    recipe = recipe_module.load(recipe_path)
    settings = settings_of(recipe)
    orm.load_group(group_label)          # fail early if the group does not exist

    summary: Dict[str, Any] = {
        'recipe': recipe['name'],
        'description': recipe.get('description', ''),
        'group': group_label,
        'criteria': settings,
        'generation': {},
        'transferability': {},
        'rlderiv_scan': [],
    }

    with tempfile.TemporaryDirectory(prefix='ld1-collect-') as temporary:
        temporary = Path(temporary)

        # ---------------------------------------------------------- generation
        node = newest(group_label, 'generation')
        directory = dump_retrieved(node, temporary / 'generation')
        parameters = node.outputs.output_parameters.get_dict()
        pseudo = node.outputs['pseudopotential'] if 'pseudopotential' in node.outputs else None
        pseudo_path = temporary / 'generation' / 'U.pbe-n-nc.UPF'
        if pseudo is not None:
            pseudo_path.write_bytes(pseudo.get_content(mode='rb'))
        summary['generation'] = {
            'pk': node.pk,
            'etot': parameters.get('etot'),
            'etotps': parameters.get('etotps'),
            'number_of_states': parameters.get('number_of_states'),
            'states': analysis.states_table(parameters.get('test_sections') or []),
            'pseudo_metadata': analysis.pseudo_metadata(pseudo_path) if pseudo_path.is_file() else {},
            'log_derivatives': (analysis.compare_run(
                directory, window=settings['window'],
                pole_threshold=settings['pole_threshold'],
                pole_tolerance=settings['pole_tolerance'])
                if (directory / 'ld1.dlog').is_file() else {}),
            'bessel_cutoffs': sorted((parameters.get('bessel_test') or {}), key=float),
            'wavefunctions': (analysis.wavefunction_agreement(
                directory / 'ld1.wfc', directory / 'ld1ps.wfc',
                cutoffs={str(state[0]): float(state[5])
                         for state in recipe['generation']['states']})
                if (directory / 'ld1.wfc').is_file()
                and (directory / 'ld1ps.wfc').is_file() else {}),
        }

        # ------------------------------------------------------ transferability
        nodes = _finished(group_label, 'configts')
        if nodes:
            node = nodes[-1]
            parameters = node.outputs.output_parameters.get_dict()
            sections = parameters.get('test_sections') or []
            rows = analysis.transferability(sections)
            configurations = recipe['configts']['test'].get('configts') or []
            for row in rows:
                if row['index'] - 1 < len(configurations):
                    row['config'] = configurations[row['index'] - 1]
            summary['transferability'] = {
                'pk': node.pk,
                'number_of_configurations': parameters.get('number_of_test_configurations'),
                'rows': rows,
                'max_abs_mry': max((abs(row['delta_mry']) for row in rows
                                    if row['delta_mry'] is not None), default=None),
            }

        # -------------------------------------------------------- rlderiv scan
        nodes = _finished(group_label, 'rlderiv_scan')
        scanned = []
        for node in nodes:
            value = float(node.inputs.parameters.get_dict()['input']['rlderiv'])
            directory = dump_retrieved(node, temporary / f'rlderiv={value:g}')
            measured = analysis.compare_run(
                directory, window=settings['window'],
                pole_threshold=settings['pole_threshold'],
                pole_tolerance=settings['pole_tolerance'])
            # the same median restricted to the *bound* region: unlike the whole window it
            # cannot be improved by moving the probe into the decaying tail
            ae = analysis.read_dlog(directory / 'ld1.dlog')
            ps = analysis.read_dlog(directory / 'ld1ps.dlog')
            scanned.append({
                'rlderiv': value,
                'pk': node.pk,
                'per_channel': measured['per_channel'],
                'pooled_median': measured['pooled_median'],
                'bound_median': analysis.pooled_median(ae, ps, window=(-3.0, -0.05)),
                'poles_all_coincide': measured['poles_all_coincide'],
                'poles': measured['poles'],
            })
        scanned.sort(key=lambda entry: entry['rlderiv'])
        summary['rlderiv_scan'] = scanned
        summary['best_rlderiv'] = analysis.select_rlderiv(scanned)

    return summary


# ------------------------------------------------------------------- writing


def write(summary: Dict[str, Any], out_dir: Optional[str | Path] = None) -> Path:
    """Write ``summary.json`` plus one CSV per measured quantity."""
    out = Path(out_dir) if out_dir is not None else recipe_module.DATA_DIR
    out.mkdir(parents=True, exist_ok=True)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=False) + '\n',
                                      encoding='utf-8')

    rows = ['label,n,l,occ,e_ae_ry,e_ps_ry,delta_ry,delta_mry']
    for state in summary['generation'].get('states') or []:
        rows.append(f"{state['label']},{state['n']},{state['l']},{state['occ']},"
                    f"{state['ae']},{state['ps']},{state['delta']},"
                    f"{1000.0 * state['delta']:.3f}")
    (out / 'states.csv').write_text('\n'.join(rows) + '\n', encoding='utf-8')

    rows = ['index,config,dE_ae_ry,dE_ps_ry,delta_ry,delta_mry']
    for entry in summary['transferability'].get('rows') or []:
        def cell(value):
            return '' if value is None else (f'{value:.6f}' if isinstance(value, float) else str(value))
        rows.append(','.join([str(entry['index']), repr(entry.get('config', '')),
                              cell(entry.get('delta_ae_ry')), cell(entry.get('delta_ps_ry')),
                              cell(entry.get('delta_ry')), cell(entry.get('delta_mry'))]))
    (out / 'transferability.csv').write_text('\n'.join(rows) + '\n', encoding='utf-8')

    channels = sorted({channel for entry in summary.get('rlderiv_scan') or []
                       for channel in (entry.get('per_channel') or {})})
    rows = ['rlderiv,pk,' + ','.join(f'{channel}_median' for channel in channels)
            + ',pooled_median,poles_coincide']
    for entry in summary.get('rlderiv_scan') or []:
        per_channel = entry.get('per_channel') or {}
        cells = [('' if channel not in per_channel else f"{per_channel[channel]['median']:.6e}")
                 for channel in channels]
        rows.append(f"{entry['rlderiv']:g},{entry['pk']}," + ','.join(cells)
                    + f",{entry['pooled_median']:.6e},{entry['poles_all_coincide']}")
    (out / 'rlderiv_scan.csv').write_text('\n'.join(rows) + '\n', encoding='utf-8')
    return out
