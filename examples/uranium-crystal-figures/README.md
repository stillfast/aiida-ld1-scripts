# Uranium crystal figures used in the manuscript

```bash
python -m pip install numpy matplotlib
python examples/uranium-crystal-figures/plot_crystal_figures.py --data /path/to/structures.json
# Optional input/output locations:
python examples/uranium-crystal-figures/plot_crystal_figures.py --data structures.json --out figures
```

This example is independent of the ld1 generation and analysis pipeline. It needs
no AiiDA, ASE, QE or database. Crystal coordinates originate in the separate
bcc-uranium-verification repository; their source commit and UUIDs are preserved
in the locally supplied JSON. This repository includes code and documentation only;
structure coordinates and generated images are supplied/stored locally. No elastic
data are inferred from these geometry drawings.

The plotter generates four manuscript figures in the selected output directory:

- `alpha-cell-origin`: one U atom placed at (0,0,0); all atoms use filled circles.
- `alpha-bc-2y`: independent fractional bc projection showing 2y. All alpha atoms
  are in the same 4c Wyckoff orbit. Orange rings only highlight dimension endpoints.
- `beta-composite`: full cell + z=0 layer + III-site ac projection, with one shared
  I--V site-colour legend. Parameter labels/leaders stay inside their own views.
- `gamma-cell`: ideal bcc conventional cell, with no independent internal parameter.

`beta-cell-sites`, `beta-ab-z0`, `beta-ac-III` are additional independent views,
with local legends/direct site labels, for possible separate placement.

All images omit embedded titles, captions, phase headings and panel letters.
Exports include 600 dpi PNG, editable PDF/SVG, 300 dpi previews and metadata.
Captions belong to the manuscript. Sphere radii are schematic; dotted coordinate
guides are not bonds. Fractional-coordinate diagrams do not use isotropic Angstrom
scales. Cell drawings contain boundary-equivalent images of the 4/30/2 unique atoms.

For alpha, the origin shift preserves all periodic pair distances. The paired atom
becomes (0,2y,1/2); 2y is the b fractional component, corresponding to length 2yb,
not the complete U--U distance. For beta, (x,y,z) -> (x,-y,z) mod 1 aligns the
representative coordinates with Qiu Table 1; all pair distances are checked.

Coordinate definitions/layout references: Beeler et al. (2013), Fig. 1 / Table 2;
Qiu et al. (2023), Fig. 1 / Tables 1 and 4. Numerical values are this work's LCAO
relaxed structures, not the reference papers' calculations.

All seven final PNGs redraw identically without optional authoring tools. The beta
composite has a render-time alignment check (1.5 pt tolerance); independent views
have one axis. Local authoring checks found PDF glyph sizes of at least 5.6 pt and 0 FAIL/0 WARN
in collision checks. Generated images and QA files are not committed here.

## Local input schema

The JSON is an object with a `structures` mapping containing `alpha`, `beta`,
`gamma`. Each structure supplies `fractional_raw` (N by 3), `cell_lengths_A` (3),
`cell_volume_A3`, `uuid` and per-atom `labels`. The conventional cells have 4, 30
and 2 atoms. Use this workflow's local archived structure JSON or prepare the
same schema from your own calculations; the source coordinates are not bundled.
`repository` and `commit` may be supplied as provenance metadata. The supplied
site ordering is validated against the representative coordinate conventions
before drawing. Output directories and local inputs are ignored by Git.
