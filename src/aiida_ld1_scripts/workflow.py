"""Submitting the recipe: generation, the ``rlderiv`` scan and the transferability test.

Every computation is a :class:`aiida_ld1.calculations.ld1.Ld1Calculation` and everything
lands in one AiiDA group, so the provenance of a pseudopotential is one query away:

    from aiida import orm
    group = orm.load_group('aiida-ld1-6p25')
    generation = next(node for node in group.nodes
                      if node.process_label == 'Ld1Calculation'
                      and node.base.attributes.get('ld1_role') == 'generation')

The two entry points are :func:`submit_generation` (also used by the ``rlderiv`` scan,
which is just the generation recipe with one number changed) and
:func:`submit_transferability`, whose ``pseudo_in`` is the *output node* of the
generation run -- that link is the whole reason to keep this in AiiDA.
"""

from __future__ import annotations

import sys
from typing import Any, Dict, Iterable, List, Optional

from aiida import orm
from aiida.engine import run_get_node, submit
from aiida.plugins import CalculationFactory

from . import recipe as recipe_module

PLUGIN = 'ld1.ld1'
#: roles are recorded as **subgroups** (``<group>.generation``, ``<group>.rlderiv_scan``,
#: ``<group>.configts``), not as node attributes: a finished ``CalcJob`` is sealed and
#: its attributes are immutable, while a group accepts a node at any time.
ROLES = ('generation', 'rlderiv_scan', 'configts')


def get_code(label: str) -> orm.AbstractCode:
    try:
        return orm.load_code(label)
    except Exception as exception:                          # noqa: BLE001
        raise SystemExit(
            f'cannot load the code {label!r}: {exception}\n'
            f'Create it once with (ld1.x is serial):\n'
            f'  verdi code create core.code.installed --label {label.split("@")[0]} \\\n'
            f'      --computer {label.split("@")[-1]} \\\n'
            f'      --filepath-executable /path/to/q-e/bin/ld1.x \\\n'
            f'      --default-calc-job-plugin {PLUGIN} --no-with-mpi') from exception


def get_group(label: str) -> orm.Group:
    """The results group, created on first use."""
    try:
        return orm.load_group(label)
    except Exception:                                       # noqa: BLE001
        return orm.Group(label).store()


def role_group(group_label: str, role: str) -> orm.Group:
    """The subgroup holding one role (``generation``, ``rlderiv_scan``, ``configts``)."""
    return get_group(f'{group_label}.{role}')


def nodes_with_role(group_label: str, role: str) -> List[orm.ProcessNode]:
    """Finished calculations of one role, oldest first."""
    try:
        group = orm.load_group(f'{group_label}.{role}')
    except Exception:                                       # noqa: BLE001
        return []
    nodes = [node for node in group.nodes
             if getattr(node, 'process_label', None) == 'Ld1Calculation']
    return sorted(nodes, key=lambda node: node.pk)


def build(parameters: Dict[str, Any], code: orm.AbstractCode, *,
          label: str, description: str = '',
          pseudo_in: Optional[orm.SinglefileData] = None,
          wallclock: int = 600):
    """A ready ``ProcessBuilder`` for one ``ld1.x`` run."""
    calculator = CalculationFactory(PLUGIN)
    builder = calculator.get_builder()
    builder.code = code
    builder.parameters = orm.Dict(parameters)
    if pseudo_in is not None:
        builder.pseudo_in = pseudo_in
    # a serial atomic-physics code: no MPI, and the direct scheduler still wants two
    # of the three resource keys
    builder.metadata.options.resources = {'num_machines': 1,
                                          'num_mpiprocs_per_machine': 1}
    builder.metadata.options.max_wallclock_seconds = wallclock
    builder.metadata.label = label
    if description:
        builder.metadata.description = description
    return builder


def _launch(builder, *, role: str, wait: bool, group_label: str):
    """Submit (or run in process when ``wait``) and file the node under its role."""
    if wait:
        _results, node = run_get_node(builder)
    else:
        node = submit(builder)
    get_group(group_label).add_nodes(node)
    role_group(group_label, role).add_nodes(node)
    return node


def already_submitted(group_label: str, role: str) -> List[orm.ProcessNode]:
    """Finished calculations of this role, oldest last (for idempotent re-runs).

    A node can only be asked about its role through its subgroup, so this is the inverse
    of :func:`role_group`.
    """
    try:
        group = orm.load_group(f'{group_label}.{role}')
    except Exception:                                       # noqa: BLE001
        return []
    nodes = [node for node in group.nodes
             if getattr(node, 'process_label', None) == 'Ld1Calculation'
             and node.is_finished_ok]
    return sorted(nodes, key=lambda node: node.pk)


def _summarise(node) -> str:
    state = node.process_state
    exit_status = node.exit_status
    return f'pk={node.pk} {state} exit={exit_status}'


def submit_generation(recipe: Dict[str, Any], code_label: str, group_label: str,
                      *, wait: bool = False, reuse: bool = False):
    """The generation run (``role = generation``)."""
    if reuse:
        previous = already_submitted(group_label, 'generation')
        if previous:
            print(f'reusing generation {_summarise(previous[-1])}')
            return previous[-1]
    builder = build(recipe_module.generation_parameters(recipe), get_code(code_label),
                    label=f"{recipe['name']} generation",
                    description='ld1.x generation of the full pseudopotential',
                    )
    node = _launch(builder, role='generation', wait=wait, group_label=group_label)
    print(f'generation: {_summarise(node)}')
    return node


def submit_rlderiv_scan(recipe: Dict[str, Any], code_label: str, group_label: str,
                        *, wait: bool = False, reuse: bool = False
                        ) -> List[orm.ProcessNode]:
    """One generation per ``rlderiv`` value (``role = rlderiv_scan``)."""
    nodes: List[orm.ProcessNode] = []
    for variant in recipe_module.rlderiv_variants(recipe):
        value = variant['rlderiv']
        role = 'rlderiv_scan'
        if reuse:
            previous = [node for node in already_submitted(group_label, role)
                        if abs(float(node.inputs.parameters.get_dict()
                                     ['input']['rlderiv']) - value) < 1e-9]
            if previous:
                print(f'rlderiv = {value}: reusing {_summarise(previous[-1])}')
                nodes.append(previous[-1])
                continue
        builder = build(variant['parameters'], get_code(code_label),
                        label=f"{recipe['name']} rlderiv={value:g}",
                        description=f'generation with rlderiv = {value:g} bohr',
                        )
        node = _launch(builder, role=role, wait=wait, group_label=group_label)
        print(f'rlderiv = {value:g}: {_summarise(node)}')
        nodes.append(node)
    return nodes


def submit_transferability(recipe: Dict[str, Any], code_label: str, group_label: str,
                           *, wait: bool = False, reuse: bool = False,
                           generation_node=None):
    """The ``iswitch = 2`` test of the generated pseudopotential (``role = configts``).

    ``pseudo_in`` is the generation's ``pseudopotential`` output node, so the test is
    linked to the exact file it validates.
    """
    if generation_node is None:
        previous = already_submitted(group_label, 'generation')
        if not previous:
            raise SystemExit('no finished generation run in the group: submit 2-submit/01 first')
        generation_node = previous[-1]

    if 'pseudopotential' not in generation_node.outputs:
        raise SystemExit(f'generation {generation_node.pk} has no pseudopotential output')
    pseudo = generation_node.outputs['pseudopotential']

    if reuse:
        previous = already_submitted(group_label, 'configts')
        if previous:
            print(f'reusing transferability run {_summarise(previous[-1])}')
            return previous[-1]

    builder = build(recipe_module.configts_parameters(recipe), get_code(code_label),
                    label=f"{recipe['name']} transferability (nconf=5)",
                    description='iswitch = 2 test of the generated pseudopotential '
                                f'against {len(recipe["configts"]["test"]["configts"])} '
                                'configurations',
                    pseudo_in=pseudo)
    node = _launch(builder, role='configts', wait=wait, group_label=group_label)
    print(f'transferability: {_summarise(node)}')
    return node
