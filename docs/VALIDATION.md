# Validation and result comparison

The completed paired/unpaired extension checks (released as v0.2.1) are recorded in
[VALIDATION_v0.2.0_PAIRED_UNPAIRED.md](VALIDATION_v0.2.0_PAIRED_UNPAIRED.md).

The test runner sources `tests/test_*.R` and records failures without starting
the full study analysis. Run it from the repository using:

```sh
Rscript scripts/run_tests.R
```

Then run the synthetic smoke test:

```sh
Rscript scripts/run_smoke_test.R
Rscript scripts/run_smoke_test.R --unpaired --leave-one-out
```

Add `--include-dtu` to exercise the exact DRIMSeq 1.38.0/stageR branch, and
`--plots` to include figure export. The default synthetic smoke run excludes
DTU. Bounded test files also print explicit skips if the exact DRIMSeq version
or stageR is unavailable; a passed file with a skipped branch does not mean
that branch was tested.

The synthetic dataset tests program behavior and data accounting. Reduced
simulation counts used in the smoke configuration do not satisfy production
calibration requirements and must not be interpreted biologically.

The independent-sample smoke mode assigns distinct synthetic subject IDs and
uses condition-only inference. `--leave-one-out` adds every subject omission
in paired mode or every sample omission in unpaired mode. Combine these flags
with `--include-dtu` to check the native model on each cohort when the required
packages are installed.

## Input tests

- Classification/junction keys are unique and valid.
- Integer FL values are finite and nonnegative.
- Local-end variants sharing a chain aggregate without losing FL support.
- Core counts retain all local rows before canonical annotations are joined.
- A chain's competing gene labels are resolved deterministically and audited.
- No-junction identifiers remain sample-specific.
- Missing optional annotation stays missing.
- Reordering metadata does not change the identity/count results.
- Duplicated libraries, incomplete pairs, missing core columns and fractional
  counts produce useful errors.

## Gain/loss tests

- The complete count matrix defines rarefaction depth.
- Presence requires both FL and within-gene fraction thresholds.
- A gene unavailable in one sample is not called an absent isoform.
- All mapped isoforms contribute to the gene denominator.
- Opposite-direction events veto the cohort call as specified.
- A robust event is a full-depth cohort event recovered in the required
  fraction of all requested rarefactions.
- Every exclusion uses the primary-cohort depth target.

## Remodeling tests

- Base-2 JSD identities, exact-zero handling and metric bounds.
- The effect gate is JSD or max absolute delta IF.
- Dominance, secondary-remodeling and class boundary examples.
- F, U and NQ represent different states.
- Paired bootstrap retains each patient's two specimens together.
- The vectorized and reference bootstrap routines agree on bounded examples.
- Covariance/profile and null-inference status remain visible.
- Independent refinement never pools selection and refinement draws.
- Observed-only leave-one-out is not labelled a bootstrap or model refit.
- A changed seed invalidates the stochastic cache, whereas a font change does
  not alter numerical analysis identity.

## DTU tests

- The paired design is full rank and tests the intended condition coefficient.
- Filtering and feature identifiers are reproducible.
- Unavailable scheduled gene and feature tests remain in their frozen families.
- Audited boundary cases are distinguished from unexplained numerical failures.
- The adapter checks DRIMSeq 1.38.0 before modifying local helper copies.
- StageR includes the full gene/feature families at each requested alpha.
- Patient directional support is distinct from a biological-replicate DTU fit.
- Optional ORF fields are not required or fabricated.

## Independent-sample and exclusion tests

- Paired cohorts omit complete subject pairs; unpaired cohorts omit only the
  specified sample. Wrong-design exclusion fields are rejected.
- None, selected and rolling exclusions produce the requested main/omission
  cohorts. Every retained condition meets the minimum biological sample count.
- Independent DTU uses a full-rank condition-only design and an intercept-only
  null; subject identifiers do not affect the model or mean-fraction effects.
- Independent remodeling uses within-group bootstrap and independent-label
  permutation, preserves caller RNG state, honors configured interval confidence
  and reports the finite permutation resolution.
- Raw-Pigeon integration exercises the full independent cohort and all eight
  omissions from an eight-sample fixture.
- Matching completed cohorts and interrupted gene blocks resume correctly.
  Incompatible scientific settings, design or runtime versions invalidate caches.
- Replotting renders independent effects and usage without paired classes or
  support reports, rejects mismatched designs and preserves source files.

## Study-data regression

On the analysis machine, compare the refactor with the completed paper outputs:

1. Per-sample original and aggregated FL totals; local/chain/gene ID sets.
2. Canonical annotation and eligible gain/loss universe.
3. Full-depth, common-depth and final robust gain/loss tables.
4. Per-gene quantifiability, metrics, class, bootstrap flags and support.
5. Cohort-test family, numerical availability and empirical/refined P/FDR.
6. DTU gene/feature family, numerical status and alpha-specific stageR results.
7. Structural-category denominators, matching controls and pathway backgrounds.
8. Patient-exclusion cohorts and their independent filtering/refitting outputs.
9. Interrupted versus uninterrupted execution for the same input/configuration.

Exact agreement is required for identifiers, categories and integer counts.
Specify an appropriate numerical tolerance for continuous estimates and record
the R/package/OS environment. Reuse the paper's deterministic seed namespaces
where exact simulation continuity is intended. Do not select a tolerance after
seeing a discrepancy or substitute a smoke result for this comparison.

No code here claims new proteomic validation, RNA integrity measurements,
independent biological replication or an experimental confirmation of predicted
protein consequences. Those are outside this workflow.
