# DTU reproducibility decisions - 26 October 2026

## Confirmed outcome

The isolated study-machine refit reproduced the original DTU and stageR
results exactly. Comparison covered gene and isoform P values, FDR values,
model estimates, unavailable results, stageR gene/feature values and the
identities of significant genes and isoforms. Agreement was not judged from
discovery counts alone. All 96 input checks and five completion gates passed.

Use the corrected DTU results here for the combined plotting update:

```text
F:/IsoRemodel/DTU_refit/20261026T043514_45632/M06_dtu
```

| Cohort | FDR threshold | DTU genes | Confirmed isoforms |
|---|---:|---:|---:|
| All six pairs | <0.05 | 497 | 593 |
| All six pairs | <0.10 | 681 | 865 |
| Patient 2 excluded | <0.05 | 377 | 452 |
| Patient 2 excluded | <0.10 | 571 | 714 |

These values summarize the validated study results; they are not inserted as
targets into model fitting. The manuscript's primary RNA cohort remains all
six patient pairs. Patient-2 exclusion remains a sensitivity analysis.

The exact comparison covered DTU and stageR. Pathway-enrichment tables were
not directly compared in the confirmation bundle; plotting reads the existing
corrected-run enrichment tables without recomputing their tests.

## Locked RNG contract

- RNG kind: `L'Ecuyer-CMRG`.
- Normal generation: `Inversion`.
- Sampling: `Rejection`.
- Study DTU master seed: 260605; effective cohort seeds: 261605 and 262605.
- Construct the serial backend before selecting the scoped RNG and seed,
  retaining the audited order of operations.
- Restore the caller's RNG kind and `.Random.seed` after success or error,
  including when no `.Random.seed` existed before the call.
- Include the effective RNG algorithm and cohort seeds in DTU cache identity
  and saved run records. Do not inherit the ambient RNG algorithm.
- Retain the DRIMSeq 1.38.0 numerical adapter, filtering, paired design,
  FDR calculation, stageR implementation and numerical-failure accounting.

## Explanation established by the controlled check

The original gain/loss module selected `L'Ecuyer-CMRG` globally. Its later DTU
module set a seed without selecting an algorithm, and inherited that setting.
IsoRemodel restored the caller's earlier RNG state, so the first public run's
DTU step instead used Mersenne-Twister. A seed alone does not define the same
sampling sequence across algorithms.

DRIMSeq selects a random subset of genes for initial common-precision
estimation. Changing only the RNG algorithm reproduced the two saved
common-precision estimates from identical ordered inputs:

| Cohort | L'Ecuyer-CMRG: original | Mersenne-Twister: earlier IsoRemodel |
|---|---:|---:|
| All six pairs | 624.492832892777 | 357.346436290016 |
| Patient 2 excluded | 55.7280900008412 | 93.5032784264444 |

The subsequent corrected refit established agreement in the final tests and
stageR results as well. No change to counts, testing families, measured
effects, filtering definitions or multiplicity procedures was needed.

## v0.2.0 update

The corrected RNG behavior remains part of the package's DTU implementation.
For the completed study, `replot_isoremodel()` reads the already verified DTU
results and the completed early RNA results into a new export folder. It does
not repeat the refit, upstream analysis, resampling, stageR or enrichment
tests. Prior results and caches remain available.

The approved plotting decisions are described in `docs/RELEASE_v0.2.1.md`.
The four-example pilot succeeded with unchanged source data; full-collection
visual review remains required after exporting the update.

ORF consequences, DNA analyses and ASE/ASTU remain outside this package.
