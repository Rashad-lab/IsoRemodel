# IsoRemodel 0.2.1

- Add explicit paired/unpaired design selection and optional named or rolling
  leave-one-subject/sample-out sensitivity analyses, with design-specific
  validation and inference.
- Add the study BioProject, RNA run download guidance and bioRxiv citation.
  Present public documentation without internal manuscript-revision language.

- Identify the paired/unpaired extension as 0.2.1, distinct from the 0.2.0
  plotting release. Update release metadata, output names and configuration paths.

# IsoRemodel 0.2.0

- Retain the confirmed DTU RNG contract: L'Ecuyer-CMRG, Inversion normal
  generation and Rejection sampling, with the original cohort seeds.
  The study-machine refit reproduced the original gene tests, FDR calls
  and stageR results at both reporting thresholds.
- Add a complete replot-only entry point that combines the completed early
  RNA results with the verified isolated DTU refit. It writes a separate
  output tree, figure audit, logs, progress records and review bundle.
- Adopt compact DejaVu Sans, 16-point titles and axis titles, 14-point body
  labels, 600-dpi raster exports, and PDF/SVG/PNG/JPEG output.
- Replace long intersection labels with aligned count bars and an exact
  membership matrix. Retain complete counts and meaningful component panels.
- Allocate enrichment row spacing from the number of wrapped label lines.
  Use integer count legends, visible minimum point sizes and fraction axes
  starting at zero. Preserve nominal-only pathway captions.
- Repair filtering-flow endpoint labels, technical-QC legends/sample labels,
  heatmap guides and individual gene-example labels.
- Record the agreed statistical and plotting decisions in the release guide.
  A final review of the complete study-machine figure output is still needed.

# IsoRemodel 0.1.2

- Set the DTU RNG algorithm explicitly to L'Ecuyer-CMRG, with Inversion normal
  generation and Rejection sampling. Keep the original cohort seeds.
- Restore the caller's random-number settings and seed after DTU, including
  failed fits. Include the effective algorithm and seed in DTU cache identity
  and run records.
- Add a targeted DTU refit that reads completed import and remodeling caches,
  checks ordered counts and model design before fitting, and writes to a new
  timestamped output directory.
- Compare common precision, final gene and isoform test results, FDR calls, and
  stageR results at FDR 0.05 and 0.10 with the original analysis and the previous
  IsoRemodel run. Retain a compact review bundle on success or failure.
- Defer all plot-style and layout changes. Figure export is disabled in the
  targeted confirmation run.

# IsoRemodel 0.1.1

- Introduced the IsoRemodel package name, banner, run logging and progress
  records.
