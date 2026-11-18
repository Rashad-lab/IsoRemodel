# IsoRemodel v0.2.1 release guide

## Paired and unpaired analyses

Version 0.2.1 includes the complete input-to-result workflow for paired subjects
and independent samples. Start with `config/analysis.yaml` for paired data or
`config/analysis_unpaired.yaml` for independent samples. Edit input/reference
paths, comparison labels and metadata before running:

```r
source("load.R")
run_rna_analysis("config/analysis_unpaired.yaml", validate_only = TRUE)
run_rna_analysis("config/analysis_unpaired.yaml")
```

`sensitivity.mode` selects no automatic omissions (`none`), separately named
subject/sample omissions (`selected`), or full rolling leave-one-out (`all`).
Named joint-exclusion cohorts are also supported. Every cohort retains the
configured sample/support minima. Full rolling runs repeat cohort inference
and native DTU fitting; they are opt-in and can resume compatible checkpoints.

Unpaired mode uses independent-group prevalence/QC, within-group bootstrap,
label-permutation composition tests and condition-only DTU. Its results do not
receive paired architecture classes or patient-direction support. Saved-result
replotting checks design compatibility and uses design-specific reports.

The README includes BioProject PRJNA1492232, RNA run SRR39530479, the twelve cDNA
barcodes, a short upstream Pigeon preparation workflow and the bioRxiv citation.
See `docs/METHODS.md` and `docs/VALIDATION.md` for method and test details.

## Run the update on the completed study

Extract `IsoRemodel_v0.2.1` under `F:/IsoRemodel`. In the existing R environment:

```r
source("F:/IsoRemodel/IsoRemodel_v0.2.1/run_replot.R")
result <- replot_isoremodel()
```

This is the recommended update path for the completed study. It reads early
RNA results from `F:/IsoRemodel/outputs` and DTU results from the confirmed
refit at `F:/IsoRemodel/DTU_refit/20261026T043514_45632/M06_dtu`. It writes a
new timestamped folder under `F:/IsoRemodel/replot_v0.2.1`.

The new folder includes individual figures, `figure_index.tsv`,
`run_summary.tsv`, `reporting_settings.yaml`, source hashes, `run.log`,
`progress.tsv`, `progress_latest.tsv` and review PDF volumes under
`review_volumes/`. Upload `IsoRemodel_v0.2.1_replot_review.zip` for the full
visual audit. The ZIP retains review material and small plotted-data tables,
not model objects or every large raster export.

Use an empty destination outside both source trees when supplying
`output_dir`. Source result folders are not overwritten. Plotting-only
execution does not refit DTU, run stageR, resample, or recompute enrichment
tests. Missing required corrected DTU evidence is an error, not permission to
fall back to the old DTU output tree.

## Statistical decisions retained

The DTU implementation explicitly selects `L'Ecuyer-CMRG` with `Inversion`
normal generation and `Rejection` sampling. It retains study master seed
260605 and the effective cohort seeds 261605 and 262605. The caller's RNG
state is restored, and RNG settings are included in cache keys and records.

The completed isolated refit reproduced the original gene/feature tests,
FDR values, unavailable results, model estimates and stageR results exactly.
The retained study counts are:

| Cohort | FDR threshold | DTU genes | Confirmed isoforms |
|---|---:|---:|---:|
| All six pairs | <0.05 | 497 | 593 |
| All six pairs | <0.10 | 681 | 865 |
| Patient 2 excluded | <0.05 | 377 | 452 |
| Patient 2 excluded | <0.10 | 571 | 714 |

Both reporting thresholds remain included. All six pairs are the primary
RNA cohort; patient-2 exclusion is a sensitivity analysis. These count values
are records of the validated outputs, not constraints used to create results.
Enrichment was outside the exact final-result comparison; existing corrected
DTU enrichment tables remain the plotting source.

## Approved plotting decisions

| Item | Decision |
|---|---|
| Typeface and sizes | Compact DejaVu Sans; titles/axis titles 16 pt; labels, subtitles and legends 14 pt |
| Colors | Retain the configured article palettes and their meanings |
| Formats | PDF, SVG, PNG and JPEG; raster exports 600 dpi |
| Composite figures | Also export meaningful individual panels; omit blank/spacer panels |
| Labels | Concise scientific language; FDR in reader-facing text; no module numbers or "exploratory" labels |
| Overlap plots | Aligned exact-intersection count bars and membership dots |
| Enrichment | Allocate vertical row spacing from wrapped term height; reserve label width |
| Enrichment count legends | Integer breaks only; preserve values and mapped count meaning |
| Hallmark refinements | Start fraction axis at 0%; minimum point size 3 for visible four-gene points |
| Filtering flow | Reserve sufficient room for endpoint counts in the composite |
| QC, heatmaps and examples | Apply sample-label wrapping, compact legends, readable guides and facet labels |

`reporting.maximum_upset_intersections` defaults to 40. If a plot displays
only the largest intersections, the caption states that restriction and all
exact counts remain in the accompanying data. Set sizes are checked against
the reconstructed exact intersections. This display cap does not change
membership or analytical results.

`reporting.enrichment_label_width` defaults to 40 characters and
`reporting.minimum_enrichment_point_size` to 3. Height adapts to the actual
wrapped lines. Retain a suitably large panel or separate enrichment figure
instead of shrinking long text until it is unreadable. Fractions remain in
their valid range. Nominal-only findings keep the nominal-P caption; no
pathway is presented as an FDR discovery merely because it is displayed.

## New analyses

For a new full analysis, edit the inputs and references in the YAML, then use
the normal entry point:

```r
source("F:/IsoRemodel/IsoRemodel_v0.2.1/load.R")
run_rna_analysis("F:/IsoRemodel/IsoRemodel_v0.2.1/config/analysis.yaml")
```

The supplied v0.2.1 YAML writes new analyses to
`F:/IsoRemodel/outputs_v0.2.1`, separate from the existing completed outputs.
The previous frozen pathway-reference cache remains configured for reuse.
No raw data, reference files or fitted models are distributed in this release.

## Validation and remaining review

The earlier four-example study-machine plotting pilot completed without
warnings or errors at 600 dpi. Source hashes and saved plotted values remained
unchanged. Its exact-intersection totals reproduced 497, 28, 377 and 27, and
the Hallmark count legend displayed 4 and 5. These checks support the selected
layout rules; they do not establish that all plots in the complete collection
are free of clipping or overlap.

The release validation record describes bounded code/rendering checks. The
original full production workflow has not been rerun in the delivery runtime.
After the update is exported on the study PC, inspect every figure family,
both cohorts, both FDR thresholds and separate panels. Recheck crowded labels,
legend ordering, very small point symbols and dense intersections at the
intended manuscript size. The analytical results remain locked during this
review.
