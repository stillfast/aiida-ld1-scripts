# AE/PS radial comparison: figure contract and QA

- Role: compare the five occupied reference-state AE/PS radial functions outside
  the pseudization region, while showing the smoothing inside it.
- Input: the committed AiiDA archive; wavefunction object hashes are verified
  against its SHA256 repository keys. No new calculation or modified source data.
- Encoding: orbital colour + independent AE/PS legend; AE solid, PS dashed with
  hollow circles; dotted vertical cutoffs. Signed u(r), unscaled amplitudes,
  original phase alignment, logarithmic 0.01–100 bohr horizontal axis.
- Continuous curves: every visible radial sample. Only circle positions are spaced
  to avoid marker clutter. No smoothing, interpolation or error bars.
- Exports: 600-dpi PNG, PDF with TrueType fonts, SVG with editable text.
- Static source preflight: 18 pass, 3 advisory warnings, no failures.
  The warnings concern a generic TIFF submission preference, dynamically supplied
  DPI, and journal-specific physical width. This Markdown review uses PNG at
  its retained Obsidian width; PDF/SVG are also provided, PNG is actually 600 dpi,
  and the existing 157.5 mm source width is preserved. No target-journal
  submission size has been assumed.
- PDF glyph audit: minimum 5.46 pt; no text below 5 pt.
- Rendered collision audit: 0 failures and 0 warnings.
- Alignment: NOT APPLICABLE, one plot area; geometry recorded in alignment JSON.
- Visual inspection: method legend occupies the unused upper-left area; orbital
  legend is upper-right. Both are clear of curves and cutoff lines; axis labels
  and limits are visible. Overlapping outer AE/PS curves are intentional and
  distinguished by PS markers.
- Tests: the full repository suite passes (45 tests), including hollow marker,
  full-grid dashed curve and explicit AE/PS legend checks.

Reproduce with: python 3-analyze/figures.py --radial-only
