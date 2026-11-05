# aiida-ld1-scripts

Scripts that **generate and verify a pseudopotential** with
[`aiida-ld1`](https://github.com/stillfast/aiida-ld1), storing the whole chain in AiiDA.

This is the reproduction, in AiiDA, of the "6p25" uranium pseudopotential work of
`uranium-ppgen/ai_v1/6p25`: one generation run, a scan of the logarithmic-derivative
reference radius `rlderiv`, and a five-configuration transferability test — plus the
analysis of the three, and an AiiDA archive that carries all of it.

It follows the layout of
[`acwf-verification-scripts`](https://github.com/aiidateam/aiida-common-workflows):
one numbered folder per stage, a single input file, and an archive exported at the end.

## The result

The recipe is a 14-valence-electron uranium pseudopotential (5f³6s²6p⁶6d¹7s²),
scalar-relativistic, norm-conserving, PBE, Troullier–Martins, with a sixth (unoccupied)
6P projector pinned at −0.25 Ry. Reproduced from AiiDA:

| quantity | value |
|---|---|
| Etot / Etotps of the generating configuration | −55958.130614 / **−181.732842** Ry |
| worst bound-state error (PS vs AE) | **0.08 mRy** (5F) |
| worst radial-function tail mismatch (`r > rcut`) | **9.4 × 10⁻⁴** |
| logarithmic derivatives, median \|ΔL\| | 1.03 × 10⁻³ |
| poles of AE and PS | coincide in every channel |
| transferability, worst \|ΔE\| over 5 configurations | **1.08 mRy** |
| best `rlderiv` | **8.0 bohr** (pooled median 1.38 × 10⁻⁴, the best of 5) |

The full report is [data/REPORT.md](data/REPORT.md), generated from the AiiDA group.

**Does it reproduce the earlier work?** Yes, number for number:

* the five `rlderiv` rows (20 per-channel medians and 5 pooled medians) match the old
  `report.md` to `1e-4` relative — `tests/test_analysis.py` asserts this;
* the pole verdicts match (`l=0`, `l=1` coincide, `l=2`, `l=3` are displaced by more than
  the 0.04 Ry tolerance at `rlderiv = 8`);
* Etotps reproduces to the last printed digit (−181.732842 vs −181.732843), and the
  transferability ΔE to 1 µRy (the old numbers were computed from the six decimals
  printed in `ld1.test`, these from full precision).

## Content of the folders

```
0-recipe/     the single source of truth: 0-recipe/uranium-14ve-6p25.yml
1-preliminary/  create the group, check the code, print the plan (submits nothing)
2-submit/     submit the generation, the rlderiv scan and the transferability test
3-analyze/    collect -> CSVs + summary.json, figures, REPORT.md
4-export/     verdi archive create -> data/aiida-ld1-6p25.aiida
data/         the committed results (REPORT.md, CSVs, figures, the archive)
src/          the code the scripts share (recipe loading, analysis, report)
```

Everything is derived from the recipe. Change a number there, re-run `2-submit/` and the
whole chain is regenerated with new provenance — no number is copied into a script.

## Starting

```bash
pip install -e .                      # needs aiida-ld1 >= 0.2 as well
verdi status                          # a profile must be configured
```

Register the code once (ld1.x is serial, so `--no-with-mpi`):

```bash
verdi code create core.code.installed \
    --label ld1.x --computer localhost \
    --filepath-executable /path/to/q-e/bin/ld1.x \
    --default-calc-job-plugin ld1.ld1 --no-with-mpi
```

Then:

```bash
python 1-preliminary/create_group.py --code ld1.x@localhost   # what will be submitted
python 2-submit/00_submit_all.py --wait                        # 7 calculations, seconds each
python 3-analyze/collect.py                                    # data/*.csv + summary.json
python 3-analyze/figures.py                                    # data/fig_*.png
python 3-analyze/report.py                                     # data/REPORT.md
bash   4-export/export_group.sh                                # data/aiida-ld1-6p25.aiida
```

`--wait` runs in process, so no daemon is needed. Without it the calculations are
submitted and `verdi process list` shows them. Each stage also has its own script
(`2-submit/01_…`, `02_…`, `03_…`) and `--reuse` skips anything already finished.

## What is stored, and where

| AiiDA entity | content |
|---|---|
| group `aiida-ld1-6p25` | every calculation |
| group `aiida-ld1-6p25.generation` | the run that produced the pseudopotential |
| group `aiida-ld1-6p25.rlderiv_scan` | one run per `rlderiv` value |
| group `aiida-ld1-6p25.configts` | the `iswitch = 2` transferability test |

The roles are **subgroups** rather than node attributes on purpose: a finished `CalcJob`
is sealed and its attributes are immutable, while a group accepts a node at any time.

The transferability test's `pseudo_in` is the generation's `pseudopotential` **output
node**, so the test is linked in the graph to the exact file it validates — that link is
the reason to keep this in AiiDA rather than in a directory of files.

The archive can be re-imported anywhere:

```bash
verdi archive import data/aiida-ld1-6p25.aiida
verdi group list
```

## How the numbers are defined

* **ΔE (transferability)** = `(E_AE − E_AE,ref) − (E_PS − E_PS,ref)`, with the generating
  configuration as the reference: the error the pseudopotential makes in the energy
  *difference* between two configurations, which is the quantity that matters. The sign
  is the one `ld1.x` prints in its own `Delta E=` column. A pseudopotential with a large
  absolute-energy error but small ΔE is fine for energy differences (the ACWF screening
  scale is 1/3 meV/atom = 0.333 mRy).
* **median |ΔL|**: median of `|L_PS − L_AE|` at the `rlderiv` radius over the energy
  window, **pooled over all channels** for ranking. No point is dropped: near a pole both
  curves diverge, which is why the *median* is the criterion and the *mean* is not.
* **poles**: sign change of `L(E)` with `|L| > 20`, position from the divergence model
  `L ≈ c/(E−E_P)`. AE and PS must have the same number, each within
  `pole_tolerance + grid step` (0.04 Ry). This is the ghost criterion — not the Bessel
  box test, whose `l ≥ 3` column is an artefact of the diagonalisation (see the
  `aiida-ld1` README).
* **phase matching**: `u(r)` has no defined sign, and `ld1.wfc` (all-electron run) need
  not share the phase of the pseudo orbital, so `analysis.align_phase` flips it when
  `|AE + PS| < |AE − PS|`. Note also that `ld1ps.wfc` contains the pseudo orbitals in
  its **first** half and the phase-aligned all-electron orbitals in its second — building
  a `{label: column}` dict over that file silently keeps the *all-electron* copy.

## Tests

```bash
pytest        # 12 tests, no database needed
```

The fixtures in `tests/data/` are the ten `.dlog` files of the reference `rlderiv` scan,
so the analysis is checked against the numbers of the original report rather than against
itself.

## License

MIT.
