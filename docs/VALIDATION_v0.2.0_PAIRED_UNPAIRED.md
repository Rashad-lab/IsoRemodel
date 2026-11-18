# Paired and unpaired extension validation

Date: 2026-11-02. Test-run label: 0.2.0. This extension is released as 0.2.1;
the original validation logs retain their test-run labels. The subsequent
version-record correction changes no analytical methods or settings.

All 16 bounded test files passed after the independent-sample extension and
final audit fixes. The native DTU tests executed with DRIMSeq 1.38.0 and
stageR 1.32.0, without skipping those branches. The full log is
`validation_logs/v0.2.0_paired_unpaired_tests.log`.

The source-mode validation runtime used R 4.3.3 on Linux. This checks the
implemented source workflow and exact DTU adapter, not installation of the
package in the supported paper environment (R 4.5.x/Bioconductor 3.22).
The full paper dataset was not rerun for this update.

## Native independent-model fixture

A bounded dispersed-count fixture contains 12 genes and 36 isoforms. Both
the full 3-versus-4 model and sample-excluded 3-versus-3 model produced finite
tests for all 12 genes. Checks inspect the actual fitted null-model slot
(intercept only), full-design rank, sample identities, conserved gene/isoform
families, full-family BH correction and stageR at alpha 0.05 and 0.10.
Missing feature evidence is conservatively retained as P=1 in its frozen family.

## Raw-input rolling-exclusion smoke runs

Synthetic Pigeon files were imported and processed through all stages, including
native DTU. Resampling counts were deliberately reduced for bounded execution.

| Design | Main plus omission cohorts | DTU finite genes per cohort | Unexplained numerical failures |
|---|---:|---:|---:|
| Paired, six subject pairs | 7 | 4 of 4 | 0 |
| Unpaired, twelve independent samples | 13 | 3 of 4 | 0 |

All five DTU completion gates passed in every cohort. The unpaired fixture has
one condition-specific-zero PTEN boundary per cohort, explicitly audited and
conservatively retained as P=1. StageR preserves the four-gene family at both
alpha levels, and its gene-screen FDR agrees with the full-family correction.
Compact design, RNG, fitting and stageR audits are retained under
`validation_logs/v0.2.0_native_smoke/`.

The paired whole-workflow smoke has `audit_only` status because its reduced
remodeling simulations do not meet production calibration gates. The unpaired
whole-workflow smoke has `diagnostic_complete` status. Neither smoke is a
production result or writes a production completion claim.

## Reporting and reproducibility

Saved-result replot integration exported 74 paired and 23 independent figures,
completed review PDFs/ZIPs and verified unchanged source files. Mismatched
configuration/result designs are rejected. Independent figures do not use
paired architecture classes or patient-support reports.

Additional checks cover correct exclusion units and minimum group sizes,
exact balanced-group permutation resolution, configured bootstrap confidence,
the shared JSD-or-fraction effect rule, caller RNG restoration, rarefaction count
conservation, matching cohort-cache reuse and interrupted gene-block recovery.
Independent cache keys include scientific settings, source helpers, RNG contract
and runtime/package versions.

The invalid-count runner test deliberately emits an import error and verifies
failure retention and recovery. Synthetic rarefaction and numeric-reshape
notices in the logs are expected fixture diagnostics.
