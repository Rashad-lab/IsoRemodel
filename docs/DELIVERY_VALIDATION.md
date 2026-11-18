# Delivery validation

## Current release: IsoRemodel 0.2.1

The paired/unpaired extension is released as 0.2.1 to distinguish it from
the original 0.2.0 plotting package. This version-record correction updates
metadata, release/output names and current documentation; analytical methods
and settings are unchanged. The unchanged display layout contract remains
0.2.0, and completed 0.2.0 results remain supported for replotting.

The complete extension validation is recorded in
[VALIDATION_v0.2.0_PAIRED_UNPAIRED.md](VALIDATION_v0.2.0_PAIRED_UNPAIRED.md).
Those checks ran before this numbering correction: all 16 bounded test files,
real DRIMSeq 1.38.0/stageR 1.32.0, and both rolling-exclusion smoke designs passed.
The historical logs retain their original version labels. Focused checks for
the 0.2.1 metadata, logging, configuration and run records are retained in
`validation_logs/v0.2.1_version_record_checks.log`.

## Historical plotting release: IsoRemodel 0.2.0

The v0.2.0 release incorporates the DTU RNG correction confirmed on the study
PC and the plotting layouts approved in the focused review. The existing
study refit, rather than a new local production fit, supplies the exact
DTU/stageR reproduction evidence recorded in DTU_REPRODUCIBILITY_HANDOFF.md.

New checks executed for this release used native R 4.3.3 with ggplot2 4.0.3 and
patchwork 1.3.2. The supported analysis environment remains R 4.5.x with the
audited DRIMSeq 1.38.0 adapter; this test interpreter does not redefine that.

- Nine bounded test files passed: configuration, DTU reference/matching
  utilities, targeted-refit comparisons, RNG state/contract, logging, plot
  layouts, complete replot integration, remodeling enrichment rules, and
  remodeling mathematics/classes. Actual DRIMSeq 1.38.0 and stageR execution
  subtests explicitly skipped because those packages were unavailable.
- The plotting-only integration exported 74 fixture figures covering both
  cohorts and FDR thresholds, including corrected DTU reports, meaningful
  standalone panels and saved-ORA views. It rejected incorrect RNG/completion
  provenance, excluded an intentionally invalid old DTU plot, preserved source
  hashes and restored caller RNG/logging options. Review PDF volumes and ZIP
  paths passed checks. No statistical stage was invoked.
- Native PDF, SVG, PNG and JPEG exporters passed at 600 dpi, including separate
  facets and composite components. Raster dimensions/DPI and PDF structure
  were verified; editable SVG export also succeeded.
- New enrichment/UpSet/flow layouts were rendered and visually inspected using
  the uploaded small study plotting tables. The Hallmark example starts at 0%,
  has integer count breaks and a visible minimum point size.
- All five vendored numerical-engine files are byte-identical to v0.1.2.
  AST checks confirmed that non-reporting DTU functions remain unchanged; the
  run wrapper changes only the recorded workflow version and removes an unused
  optional plotting dependency. Six other analytical R files are byte-identical.
  Scientific YAML settings are identical after excluding output/reporting.
- All R source/Rd files parsed and R CMD build produced a source package.
  This is not a claim of a complete native R CMD check or package installation.

The attempted legacy import/gain-loss/QC/runner tests could not execute in
this test runtime because data.table, Matrix and/or vegan were unavailable.
Those analytical sources were unchanged and their earlier validation remains
historical evidence below; these tests should still run in the study environment.
The full study data and corrected DTU model objects were not recomputed here.

After running replot_isoremodel() on the study PC, inspect the complete output
collection and share IsoRemodel_v0.2.0_replot_review.zip. The pilot and fixture
checks do not certify every full-study plot as free of clipping.

## Historical checks for version 0.1.2

That release was **IsoRemodel 0.1.2**. Native R 4.3.3 was used for
bounded source-mode regression checks. This is a test environment; the study
confirmation refit must run on the existing Windows/R 4.5.0 environment.

- All R source files and documentation Rd files parsed.
- RNG tests passed for success/error restoration, present/absent caller seed,
  exact reproduction of the original sampling sequence, invalid settings, and
  an effective RNG contract independent of the caller's algorithm.
- Targeted-refit tests passed for complete keyed comparisons, missing/extra
  identities, duplicate keys, numeric tolerances, unavailable versus P=1,
  exact count/sample/design ordering, and production cache validation.
- The source-only refit loader loaded without starting analysis and preserved
  caller RNG state and the code-root option.
- AST comparison confirmed that the numerical adapter, filters, reference
  ordering, stageR implementation and plotting functions were unchanged.
- Full DRIMSeq refitting, real stageR execution, enrichment, and the complete
  production review-bundle path were not executed here. Their confirmation
  uses the study PC and the returned `DTU_refit_review.zip`.

The earlier validation below is retained as historical release information.
It does not claim that all earlier numerical or visual tests were repeated
for this focused RNG change.

## Historical version 0.1.1 delivery

The public package name was **IsoRemodel**, version **0.1.1**. That update retained
the scientific workflow and adds a startup banner, timestamped console/run logs
and stage/checkpoint progress records. The analytical validation below describes
the checks actually executed; a new package name or progress message does not
extend the scope of those checks.

This package was assembled from the author-supplied RNA scripts and revision
modules. The supplied study YAML uses the publication thresholds and twelve
original Pigeon filenames. Raw patient data are not bundled.

## Checks for version 0.1.1

All R files parsed under the renamed package. All nine bounded test files
passed, including the expanded runner and new banner/logging tests. These cover
branded output, cache reuse, warning/error logs, progress history, explicit
stage indices and restoration of caller options. The PDF/SVG/PNG/JPEG typography
fixture also passed using the IsoRemodel resource paths. Native DRIMSeq/stageR
tests remain explicitly skipped in this environment.

## Checks executed during initial preparation

An actual R interpreter compiled to WebAssembly (R 4.6.0, WebR 0.6.0) was used
for syntax and bounded runtime checks. This is a validation environment, not
the paper's Windows/R 4.5.0 environment and not a dependency of this repository.

- Every delivered R source and executable test file parsed successfully.
- Eight bounded test files passed: configuration/metadata, DTU reference and
  matching/ORA utilities, gain/loss, import, QC, remodeling mathematics/classes,
  runner/checkpoint/failure handling, and remodeling foreground/background rules.
- The runner tests checked actual import, content-based resume, font-only cache
  reuse, failure audit persistence, and recovery after restoring raw input.
- A six-pair raw-Pigeon synthetic workflow executed through import, gain/loss,
  technical QC and remodeling for the complete cohort and patient-2 exclusion.
  The test used reduced bootstrap/null/depth draws and completed without a
  runtime error. Its calibration gates correctly retained audit-only status.
- An additional bounded run enabled variance modeling, independent-refinement
  selection, GPD diagnostics and all 16 combinations of two isoform counts,
  two depths, two biological variances and null/alternative truth. It completed
  without a runtime error. These tiny simulations do not establish calibration.
  No refinement candidate qualified under that low-resolution covariance setup;
  execution of a nonempty production refinement family remains a PC check.
- The plotting fixture exported PDF, outlined SVG, editable SVG, PNG and JPEG,
  plus independent facet and patchwork components. The PNG/JPEG dimensions were
  verified against the requested 600 pixels/inch. SVG text/path structure and
  PDF page structure were checked, and raster previews were visually reviewed.
  Long composite titles were wrapped to their individual panel width.

## Branches requiring native analysis-machine validation

The exact DRIMSeq 1.38.0 adapter and stageR numerical execution could not run
in the WebAssembly environment. Their tests explicitly reported skips rather
than being represented as executed. Run `scripts/run_tests.R` and
`scripts/run_smoke_test.R --include-dtu` in the paper's native R environment.

Production enrichment with AnnotationDbi/org.Hs.eg.db/clusterProfiler and actual
MSigDB resources requires the native packages and the captured resource files.
Identifier/foreground and hypergeometric utility checks do not substitute for
running those database-backed branches.

The full original study data were not available for a new production run here.
No claim of exact biological-result reproduction, numerical equivalence across
R versions, or a complete native `R CMD check` is made. Use the comparison list
in VALIDATION.md before labeling a GitHub release as a reproduced paper run.

## Deliberate differences from older descriptive code

The revised core matrix preserves all raw counts before canonical annotation
is attached. The earlier winning-gene-only descriptive complexity calculation
is retained separately. ORF/protein outcomes are excluded, so the technical
screen uses the actual RNA-only metric family and recomputes its FDR rather
than inheriting an old hard-coded number of properties. Missing optional
annotations remain unavailable. No patient-level technical comparison is
presented as an additional biological replicate.
