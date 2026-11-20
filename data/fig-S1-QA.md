# Figure S1 verification (2026-11-20)

Question: how do AE/PS logarithmic-derivative curves compare across the five
fixed-recipe diagnostic radii? A quantitative 5 × 4 grid separates the radius
and angular-momentum dependence. This is a visual diagnostic, not a radius ranking.

- Python/Matplotlib; 178 × 210 mm; PNG 600 dpi, PDF/SVG editable text.
- All ten original files, all 12,040 values retained; CSV round-trip tested against
  the committed raw fixtures. Input SHA256 hashes are in the metadata JSON.
- Signed symlog, linear |L| ≤ 1; common ±10,000 range retains all finite samples.
- Only connecting lines are broken at large opposite-sign jumps. Neither endpoint
  is removed; no interpolation, smoothing, fit or statistical summary is applied.
- Authoring alignment audit: PASS, 20 panels, 1.5 pt tolerance, 0 failures/warnings.
- PDF font audit: minimum 5.6 pt, no glyph runs below 5 pt.
- PDF collision audit: PASS, 0 failures and 0 warnings.
- Final preview inspected panel by panel: shared axes, row radii, channel labels,
  legend and endpoint samples are visible without clipped labels.
- Static source audit: 0 failures. Two warnings reviewed: TIFF is unnecessary for
  these requested PNG/PDF/SVG exports; logarithm inputs are guarded (tick zero is
  handled separately and the axis amplitude is at least 10).
- Plot-only CLI tested with AiiDA and YAML imports blocked; live-group dispatch
  tested with fake nodes. No new AiiDA/QE calculation was performed for this sync.

`fig-S1-alignment.json`, `fig-S1-font-audit.json` and
`fig-S1-collision-audit.json` preserve the rendered checks. External authoring
auditors are optional; the distributed plotter includes its own rectangle check.
