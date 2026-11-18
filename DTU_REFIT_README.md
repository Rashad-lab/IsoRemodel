# Historical DTU confirmation procedure

**The study confirmation is complete.** The corrected refit reproduced the
original DTU and stageR results exactly. This procedure is retained to explain
or deliberately repeat that comparison; it is not a required v0.2.0 update step.
For the completed study, use `run_replot.R` and `replot_isoremodel()` as described
in the main README. See `DTU_REPRODUCIBILITY_HANDOFF.md` for the locked results.

The historical v0.1.2 correction fixed the random-number settings used by DTU
and provided the targeted refit below. It did not rerun the complete RNA
workflow.

## Run on the study PC

Extract `IsoRemodel_v0.1.2.zip` into `F:/IsoRemodel`. Keep the existing
`IsoRemodel_v0.1.1` folder and its study configuration in place. In a fresh R
session, run:

```r
source("F:/IsoRemodel/IsoRemodel_v0.1.2/run_DTU_refit.R")
result <- refit_dtu_rng()
```

Sourcing the file loads the new code but does not start the analysis. The
second line starts it. Installation of the R package is optional for this
source-mode run. Keep the existing analysis packages; do not upgrade them
before this comparison. The refit checks for DRIMSeq 1.38.0 and stageR 1.32.0.

The defaults are the actual study locations:

| Setting | Path |
|---|---|
| Existing study YAML | `F:/IsoRemodel/IsoRemodel_v0.1.1/config/analysis.yaml` |
| Completed IsoRemodel outputs | `F:/IsoRemodel/outputs` |
| Original DTU results | `F:/PacBio revision analysis/outputs/M06_DTU_robustness` |
| New refit parent directory | `F:/IsoRemodel/DTU_refit` |

If the study YAML is elsewhere, pass its actual path explicitly:

```r
result <- refit_dtu_rng(
  config_file = "F:/IsoRemodel/IsoRemodel_v0.1.1/config/analysis.yaml",
  previous_output = "F:/IsoRemodel/outputs",
  original_dtu = "F:/PacBio revision analysis/outputs/M06_DTU_robustness",
  output_root = "F:/IsoRemodel/DTU_refit"
)
```

The repository's generic YAML remains available for new full analyses. This
targeted run deliberately reads the existing study YAML so it retains the
actual local annotation and reference paths.

## What runs

1. Read the completed production import and remodeling caches. No raw Pigeon
   files are imported and no upstream statistical analysis is repeated.
2. Verify that the saved original and previous IsoRemodel models used the
   same ordered counts, samples, gene identities, and paired design.
3. Prepare the new DTU families using the study settings and check them
   against the original models before either new fit starts.
4. Refit DRIMSeq for all six pairs and for the patient-2-excluded cohort.
5. Run stageR at both FDR 0.05 and 0.10, and regenerate DTU-dependent tables
   and the configured enrichment analyses.
6. Compare the final results with both saved runs and create a compact review
   ZIP.

The required cached upstream files are:

```text
F:/IsoRemodel/outputs/M01_import/cache/stage_result.rds
F:/IsoRemodel/outputs/M04_remodeling/cache/stage_result.rds
```

For each cohort, the original and previous IsoRemodel fitted caches are read
from their `cache/drimseq_fit_<cohort>.rds` files. The original files are never
used as replacement fitted results: each new DTU model is fitted afresh.

The old output directories and their completion markers stay in place.
Frozen pathway resources are copied into the new run before use. They are
not refreshed from a newer database. Full fitted objects are retained in the
new run directory and excluded from the review ZIP.

## Random-number correction

```yaml
dtu:
  rng:
    kind: "L'Ecuyer-CMRG"
    normal_kind: Inversion
    sample_kind: Rejection
```

The existing DTU master seed is 260605. It produces cohort seeds 261605
(all six pairs) and 262605 (patient 2 excluded), exactly as before. The
effective settings are recorded in the cache identity and QC tables. The
caller’s RNG settings and seed are restored after success or error.

No filtering thresholds, statistical model, numerical adapter, FDR procedure,
or plot styling was changed to make the results agree.

## What to send back

Each run creates a new timestamped directory:

```text
F:/IsoRemodel/DTU_refit/<run_id>/
  M06_dtu/                  # new models, tables, and QC
  references/              # local copies of frozen resources
  review/                  # compact comparison evidence and logs
  DTU_refit_review.zip      # upload this file
  refit_status.rds
```

Upload **`DTU_refit_review.zip`**. It includes the input checks, effective
configuration, package versions, source hashes, comparison summaries, gene
and isoform comparisons, changed discovery lists, and the run log.

The final status is `original_results_reproduced` only when the DTU completion
checks pass and the original model results and discovery sets pass the
comparisons. Differences from the previous IsoRemodel run are expected and
are reported separately. Read the comparison tables for the numerical
tolerances and the largest observed differences.

If a fit or comparison fails after a run directory has been created, completed
new outputs are retained. The default strict mode ends with an error after
writing the review bundle; upload that bundle rather than repeating the full
workflow. Missing dependencies or an invalid configuration can stop before a
run directory is created.

## Plotting and future runs

Figure export is disabled for this confirmation run. No plotting changes are
included. Once final DTU agreement is reviewed, the retained results can be
used for the separate figure-layout update.

Each call starts a fresh isolated confirmation run. There is no implicit
resume of a failed confirmation run. Keep its outputs for review before
starting another call.
