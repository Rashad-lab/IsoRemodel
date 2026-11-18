# 0.2.1 - Paired/unpaired analysis and exclusion sensitivity

- Release the independent-sample analysis branch separately from 0.2.0.
- Support paired and unpaired designs with no automatic exclusions, selected
  exclusions, or complete rolling leave-one-subject/sample-out analyses.
- Add independent prevalence, QC, bootstrap/permutation remodeling and
  condition-only native DTU, with resumable checkpoints and matching reports.
- Add BioProject/RNA-run guidance, study barcodes, upstream preparation steps
  and the bioRxiv citation to the public README.
- Record version 0.2.1 in package/citation metadata, banners and release output
  names. Analytical methods and validated settings are unchanged by this
  version-record correction; the 0.2.0 layout contract remains compatible.

# 0.2.0 - Confirmed DTU and approved figure layouts

- Retain the confirmed explicit DTU RNG algorithm, original cohort seeds and
  caller-state restoration. Record the exact study-machine agreement in the
  updated reproducibility handoff.
- Add a replot-only update combining completed early RNA results with the
  verified isolated DTU refit, without repeating statistical analyses.
- Use compact DejaVu Sans, 16-point titles/axis titles, 14-point labels and
  legends, 600-dpi rasters, PDF/SVG/PNG/JPEG and meaningful separate panels.
- Replace long overlap labels with exact count bars and membership matrices.
- Space enrichment rows by wrapped label height; use integer count legends,
  visible minimum points and fraction axes starting at zero.
- Repair endpoint labels in filtering-flow composites and general QC,
  heatmap and gene-example layouts. Preserve nominal-only pathway captions.
- Isolate new analyses in `outputs_v0.2.0`; keep the existing results and
  pathway reference cache available. Full study-output visual review remains
  required after the combined export.

# 0.1.2 — Explicit DTU RNG and targeted confirmation refit

- Pin DTU to the paper's RNG algorithm and retain its cohort seeds.
- Restore caller RNG state and include effective settings in DTU cache keys.
- Read completed import/remodeling caches without rerunning upstream work.
- Check saved and newly prepared ordered model inputs before fitting.
- Write fresh DTU results to a separate timestamped directory, with plots disabled.
- Compare model estimates, raw tests, FDR and stageR results at both thresholds.
- Preserve unavailable tests and require exact discovery membership agreement.
- Save comparison evidence, logs and a review ZIP on success or failure.
- Record the controlled RNG check and deferred plotting work in the DTU handoff.

# 0.1.1 — IsoRemodel name and run progress

- Renamed the public package and delivery directory to IsoRemodel.
- Preserved the existing R entry points and scientific settings.
- Added a versioned startup banner and timestamped severity labels.
- Retained console messages, warnings and errors in each run's log.
- Added progress history and latest-progress tables with stage position and
  elapsed time, plus updates at configured analysis checkpoints.
- Kept progress separate from a runtime estimate; checkpoint messages do not
  imply fixed-interval heartbeats.
- Preserved historical source filenames in provenance records.

# 0.1.0 — Initial reusable R workflow

- Converted the source R Markdown analyses into ordinary R stage functions.
- Replaced fixed sample/condition assumptions with CSV metadata and YAML settings.
- Preserved the original twelve Pigeon filenames and patient-2 sensitivity preset.
- Added selected/all/none one-patient exclusions without silently relaxing support.
- Preserved all raw full-length counts before canonical annotation is attached.
- Kept legacy descriptive complexity denominators separate from the revised core matrix.
- Removed ORF, NMD, CDS, protein and DNA prerequisites.
- Added file-content/seed/code-aware resume, explicit unavailable states and strict completion.
- Added actual-input synthetic fixtures and meaningful import, QC, gain/loss,
  remodeling, DTU and runner regression checks.
- Centralized manuscript typography/colors and PDF/SVG/PNG/JPEG exports.
- Added separate panels, saved plot objects and a statistics-free replot command.

This is a refactor for independent execution. It is not a certification that
study results have already been reproduced under a newly installed environment.
Record the PC environment and compare with the publication outputs before a
reproduction release. The exact-version DRIMSeq/stageR branch requires native
validation on the analysis machine.
