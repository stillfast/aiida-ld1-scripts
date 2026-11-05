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
| `rlderiv` | diagnostic only — the pole-gated pick is **7 bohr**, the production value stays 5 (see below) |

The full report is [data/REPORT.md](data/REPORT.md), generated from the AiiDA group.

**Does it reproduce the earlier work?** Yes, number for number:

* the five `rlderiv` rows (20 per-channel medians and 5 pooled medians) match the old
  `report.md` to `1e-4` relative — `tests/test_analysis.py` asserts this. Note that this
  reproduces the *numbers*, not the old conclusion: see "what the scan cannot decide";
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
  absolute-energy error but small ΔE is fine for energy differences. **This is an atomic
  quantity and is not the ACWF `delta`**, which is defined on a solid's energy–volume
  curve in meV/atom. The conversion matters: `1 mRy = 13.6 meV`, so our worst value
  (1.07 mRy = **14.6 meV**) is about a third of chemical accuracy (~43 meV), while
  `1/3 meV` would be 0.0245 mRy, not 0.333 mRy.
* **median |ΔL|**: median of `|L_PS − L_AE|` at the `rlderiv` radius over the energy
  window, pooled over all channels. No point is dropped: near a pole both curves diverge,
  which is why the *median* is used and the *mean* is not. It is **not** used to rank
  `rlderiv` — see below.
* **poles**: sign change of `L(E)` with `|L| > 20`, position from the divergence model
  `L ≈ c/(E−E_P)`. AE and PS must have the same number, each within
  `pole_tolerance + grid step` (0.04 Ry). This is the ghost criterion — not the Bessel
  box test, whose `l ≥ 3` column is an artefact of the diagonalisation (see the
  `aiida-ld1` README). It is only meaningful **within one radius**: the poles of `L` at
  `r_c` are the *box* spectrum of that radius, so they move with `r_c` by construction.
* **phase matching**: `u(r)` has no defined sign, and `ld1.wfc` (all-electron run) need
  not share the phase of the pseudo orbital, so `analysis.align_phase` flips it when
  `|AE + PS| < |AE − PS|`. Note also that `ld1ps.wfc` contains the pseudo orbitals in
  its **first** half and the phase-aligned all-electron orbitals in its second — building
  a `{label: column}` dict over that file silently keeps the *all-electron* copy.

## What the `rlderiv` scan cannot decide

The scan is a **diagnostic, not a selection**, for two measured reasons:

1. **The pooled median rewards looking outside the core.** At larger `r_c`, `L(E)` is
   dominated by the decaying tail (`u'/u → −κ + 1/r`), where all-electron and pseudo
   curves agree by construction. At `E = −1 Ry`, `|L_AE − L_PS|` is 5.9 × 10⁻⁵ at
   `rlderiv = 4` and 3.6 × 10⁻⁷ at `rlderiv = 8` — so "the median falls with `r`" says
   nothing about the pseudopotential.
2. **The pole criterion is not stationary in `r`.** The detected `l = 0` pole count of the
   *same* all-electron atom is 1 at `rlderiv = 4` and 2 at `rlderiv = 5`: the poles are the
   box spectrum of radius `r_c`, and whether a divergence survives the `|L| > 20` test
   depends on its residue. Comparing "how many poles coincide" across radii compares
   different tests.

So the pole verdict gates the choice **within** each point, and no point is allowed to win
on the median alone. That matters: the smallest median belongs to `rlderiv = 8`
(1.38 × 10⁻⁴) but its pole check **fails** (`l = 2` displaced by 0.043 Ry, `l = 3` by
0.047 Ry against a 0.04 Ry tolerance), while the points that pass the pole test are 5 and
7. The gated pick is therefore **7** (median 2.36 × 10⁻⁴), which the bound-region median
(`E < 0`, where the tail cannot help) independently ranks best as well (2.54 × 10⁻⁶).

**The production pseudopotential stays at `rlderiv = 5`**: the scan does not justify
changing it, and the choice of `r` remains a diagnostic setting. In a paper the scan
belongs in the SI as a 5-row table — plotting it hides the column that decides.

## Tests

```bash
pytest        # 16 tests, no database needed
```

The fixtures in `tests/data/` are the ten `.dlog` files of the reference `rlderiv` scan,
so the analysis is checked against the numbers of the original report rather than against
itself.

CI runs these tests only: the analysis layer needs nothing but numpy. The submission
layer (`2-submit/`) needs an AiiDA profile and a real `ld1.x`, so it is verified by
running the pipeline — which is how every number in this README was produced.

## The archive

`data/aiida-ld1-6p25.aiida` (2.1 MB) holds 14 group links, the 7 calculations with their
inputs and outputs, and **five generated UPFs** (one per generation run: the main one and
the five `rlderiv` points share four distinct files, the configts run writes none), plus
the retrieved `ld1.x` output including `ld1.dlog`, `ld1.test` and the wavefunctions.

## License

MIT.
