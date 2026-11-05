"""Loading the recipe that drives every stage.

One YAML file is the single source of truth: the generation run, the transferability
test and the ``rlderiv`` scan are all derived from it, so a number is never copied into
a script.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, List

import yaml

#: repository root (this file lives in ``src/aiida_ld1_scripts/``)
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RECIPE = REPO_ROOT / '0-recipe' / 'uranium-14ve-6p25.yml'
DATA_DIR = REPO_ROOT / 'data'
DEFAULT_GROUP = 'aiida-ld1-6p25'


def load(path: str | Path | None = None) -> Dict[str, Any]:
    """Read the recipe (default: the shipped 6p25 uranium one)."""
    recipe_path = Path(path) if path is not None else DEFAULT_RECIPE
    with open(recipe_path, encoding='utf-8') as handle:
        recipe = yaml.safe_load(handle)
    for key in ('name', 'generation'):
        if key not in recipe:
            raise ValueError(f'{recipe_path} has no {key!r} section')
    return recipe


def generation_parameters(recipe: Dict[str, Any]) -> Dict[str, Any]:
    """The ``parameters`` input of the generation ``CalcJob``."""
    generation = copy.deepcopy(recipe['generation'])
    return generation


def rlderiv_variants(recipe: Dict[str, Any]) -> List[Dict[str, Any]]:
    """One ``parameters`` dict per scanned ``rlderiv`` value."""
    scan = recipe.get('rlderiv_scan') or {}
    values = scan.get('values') or []
    variants = []
    for value in values:
        generation = copy.deepcopy(recipe['generation'])
        generation['input']['rlderiv'] = float(value)
        variants.append({'rlderiv': float(value), 'parameters': generation})
    return variants


def configts_parameters(recipe: Dict[str, Any]) -> Dict[str, Any]:
    """The ``parameters`` input of the ``iswitch = 2`` transferability run."""
    return copy.deepcopy(recipe['configts'])


def analysis_settings(recipe: Dict[str, Any]) -> Dict[str, Any]:
    """Window and pole criteria used by :mod:`aiida_ld1_scripts.analysis`."""
    scan = recipe.get('rlderiv_scan') or {}
    config = recipe.get('generation', {}).get('input', {})
    return {
        'window': tuple(scan.get('window', (-3.0, 3.0))),
        'pole_threshold': float(scan.get('pole_threshold', 20.0)),
        'pole_tolerance': float(scan.get('pole_tolerance', 0.04)),
        'grid_step': float(config.get('deld', 0.02)),
        'production_rlderiv': float(config.get('rlderiv', 0.0)),
    }
