# IsoRemodel

An R workflow for paired and unpaired long-read RNA isoform analysis from
Pigeon classification and junction files. IsoRemodel harmonizes transcript
identities across samples, summarizes RNA quality and isoform complexity,
measures condition-associated changes in isoform composition and detection,
and performs differential transcript usage (DTU) with DRIMSeq and stageR.

Version 0.2.1 supports optional named exclusions and complete leave-one-out
sensitivity analyses. Pairing is an explicit study-design choice. Paired
analyses preserve subject matching; unpaired analyses compare independent
biological samples without constructing artificial pairs.

## Quick start

1. Prepare filtered Pigeon classification and junction files for each library.
2. Install dependencies in R 4.5.x and edit the configuration and sample sheet.
3. Validate the inputs and run the selected analysis.

From the extracted or cloned repository:

```r
source("scripts/install_dependencies.R")
check_rna_dependencies()
# Run when dependencies need to be installed:
# install_rna_dependencies()

source("load.R")
run_rna_analysis("config/analysis.yaml", validate_only = TRUE)
result <- run_rna_analysis("config/analysis.yaml")
```

`config/analysis.yaml` contains the paired astrocytoma study settings and original
input filenames. Adjust paths before running; raw data and reference genomes
are not included in the software repository. For another study, edit
`comparison`, the metadata and applicable sample-support requirements.

The command-line equivalent is:

```sh
Rscript scripts/run_analysis.R --config=config/analysis.yaml --validate-only
Rscript scripts/run_analysis.R --config=config/analysis.yaml
```

## Study data and citation

The associated study profiled matched primary and recurrent astrocytomas from
six patients using PacBio HiFi DNA and Kinnex full-length RNA sequencing.
Raw sequencing data are deposited under
[NCBI BioProject PRJNA1492232](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1492232).

For the RNA analysis in this repository, select
[SRA run SRR39530479](https://www.ncbi.nlm.nih.gov/sra/?term=SRR39530479),
labelled **Merged long read RNA sequencing runs**. Its submitted record lists
**three source files**. These are the RNA sequencing inputs to demultiplex into
the twelve libraries below; the run is not a separate run for each patient.
The other blood and primary/recurrent tumor records labelled long-read DNA
are not the RNA inputs for IsoRemodel. Use the original PacBio BAM files for
upstream processing so that the required read metadata are retained.

The supplied SRA submission screenshot marks the records **To be released**
with a scheduled date of **2027-08-31**. This describes that submission record,
not a guarantee of present public availability. Check NCBI for the current
release status; downloads may be unavailable until the data are released.

For now, cite the associated preprint:

Rashad S, Ando D, Yamashita S, Shimoda Y, Kanamori M, Endo H, Niizuma K.
**Long read multi-omics sequencing reveals DNA-to-RNA evolutionary remodeling
trajectories of recurrent astrocytoma.** bioRxiv. 2026;2026.07.12.738108.
Posted July 14, 2026. [doi:10.64898/2026.07.12.738108](https://doi.org/10.64898/2026.07.12.738108).

The software citation is provided in `CITATION.cff`. A final article citation
and an archived software-release DOI can be added when available.

## Software and installation

The supported analysis environment is R 4.5.x with Bioconductor 3.22 and the
bundled adapter for **DRIMSeq 1.38.0**. The adapter checks that exact DRIMSeq
version because it uses version-specific numerical helpers. The installer
retains packages already at the requested versions and does not automatically
update unrelated packages.

```r
source("scripts/install_dependencies.R")
install_rna_dependencies()
```

Restart R after installing packages. Archived source packages may require
Rtools on Windows or compilation tools and system libraries on Linux/macOS.
To deliberately update compatible dependencies while retaining the required
DRIMSeq version:

```r
install_rna_dependencies(upgrade = TRUE, pin_known_versions = FALSE)
```

The principal dependencies are data.table, dplyr, tidyr, readr, Matrix, vegan,
yaml, digest, DRIMSeq, stageR, BiocParallel, ggplot2, patchwork, scales, ragg,
svglite, systemfonts, sysfonts and showtext. Optional annotation/enrichment
branches use AnnotationDbi, org.Hs.eg.db, GO.db, clusterProfiler and msigdbr.
The installer contains the complete dependency specification.

Installing the package is optional for source-mode use. After dependency
installation, `R CMD INSTALL .` installs the package for use with
`library(IsoRemodel)`. The `load.R` route is convenient for a repository checkout.

## Inputs

### Preparing Pigeon input files

IsoRemodel starts from filtered Pigeon outputs. The upstream bulk Iso-Seq
workflow that produces these files is:

1. **Remove primers and assign samples.** Start from PacBio HiFi reads, or
   segmented reads for the corresponding MAS-Seq/Kinnex library preparation.
   Use `lima` to remove cDNA primers and demultiplex the sample barcodes.
2. **Refine full-length reads.** Use `isoseq refine` to remove poly(A) tails
   and artificial concatemers, producing full-length non-concatemer (FLNC)
   reads. Combine sequencing movies belonging to the same sample for
   downstream clustering, keeping biological samples separate.
3. **Cluster and align.** Use `isoseq cluster2` to generate transcript
   consensus sequences, then align them to the reference genome with `pbmm2`
   using the Iso-Seq preset.
4. **Collapse and prepare.** Use `isoseq collapse` to obtain nonredundant
   transcript models. For bulk data, supply the FLNC BAM to obtain correct
   `flnc_count.txt` counts. Use `pigeon prepare`
   to prepare the collapsed GFF and reference annotation files.
5. **Classify and filter.** Use `pigeon classify --fl flnc_count.txt`
   to annotate transcript structures and gene assignments. Use
   `pigeon filter` to remove potential artifacts. The resulting filtered
   classification and matching junction files are the main IsoRemodel inputs.

For complete commands, reference preparation and software-specific options,
follow the official [PacBio Iso-Seq and Pigeon workflow documentation](https://isoseq.how/).
The exact study input filenames are listed below and in `config/samples.csv`.

### Study cDNA barcodes

The twelve primary and recurrent tumor libraries use the following assignments.
Retain these sample identities when combining technical sequencing runs.

| Patient | Sample | cDNA barcode |
|---|---|---|
| ptn1 | ptn1_primary | bc01 |
| ptn1 | ptn1_recurrent | bc02 |
| ptn2 | ptn2_primary | bc03 |
| ptn2 | ptn2_recurrent | bc04 |
| ptn3 | ptn3_primary | bc05 |
| ptn3 | ptn3_recurrent | bc06 |
| ptn4 | ptn4_primary | bc07 |
| ptn4 | ptn4_recurrent | bc08 |
| ptn5 | ptn5_primary | bc09 |
| ptn5 | ptn5_recurrent | bc10 |
| ptn6 | ptn6_primary | bc11 |
| ptn6 | ptn6_recurrent | bc12 |

### Sample metadata

`config/samples.csv` contains the twelve study libraries.
`config/samples_template.csv` is a starting template for a paired study.
For independent samples, use `config/analysis_unpaired.yaml` with
`config/samples_unpaired_template.csv`. There must be one row per included
biological library.

| Field | Meaning |
|---|---|
| sample_id | Unique library ID |
| subject_id | Required subject ID for paired specimens; optional in unpaired mode, where blank/missing IDs become sample IDs and supplied IDs must identify independent subjects |
| condition | One of the two labels defined in the YAML comparison |
| classification_file | Filtered Pigeon classification file |
| junctions_file | Matching Pigeon junction file |
| fl_count_column | Integer full-length count column, such as `FL.sample` |
| pigeon_summary_file | Optional Pigeon summary file |
| filter_reasons_file | Optional Pigeon filtering-reasons file |
| barcode | cDNA barcode annotation |
| include | Whether to include the library |

The study input filenames follow this pattern, shown for `ptn1_primary`:

```text
ptn1_primary_classification.filtered_lite_classification.txt
ptn1_primary_classification.filtered_lite_junctions.txt
ptn1_primary_classification.filtered.summary.txt
ptn1_primary_classification.filtered_lite_reasons.txt
```

Classification and junction files must describe compatible transcript IDs and
reference coordinates. Counts must be finite, nonnegative integers. Set the
actual Pigeon count-column name in the sample sheet: `FL.sample` is the study
setting, while other runs may produce `FL.<sample>`. TPM and `fl_assoc` are not
automatic substitutes for the configured integer FL counts. Leave an optional
file path blank when that file is unavailable.

Relative sample-file paths resolve from `inputs.root_directory`; paths in YAML
resolve from the YAML file's directory. Absolute paths are also accepted.
Reference genome/annotation versions must remain consistent with the upstream
processing. The study preset specifies GRCh38 and a GENCODE v49 GTF for its
annotation-dependent analyses. Preserve the reference/database versions and
exported pathway membership snapshots for reproducibility.

## Choose paired or unpaired analysis

Use the design that matches how the biological samples were collected.

| Setting | Paired mode | Unpaired mode |
|---|---|---|
| `comparison.paired` | `true` | `false` |
| Biological unit | Matched subject | Independent sample |
| Metadata requirement | One library in each condition per subject | Independent subject IDs and at least the configured minimum samples per condition |
| DTU design | Subject plus condition | Condition only |
| Composition inference | Paired resampling and the paired study model | Independent-group bootstrap and condition-label permutation |
| Exclusion unit | Entire subject pair | One independent sample |

For paired data:

```yaml
comparison:
  reference: primary
  comparison: recurrent
  paired: true
  incomplete_pairs: error
  technical_replicates: error
```

For unpaired data, use `config/analysis_unpaired.yaml`, set `paired: false`
and use the two appropriate condition labels. The sample sheet can omit
`subject_id`; missing values are filled from `sample_id`. If supplied, subject
IDs must be unique across the independent samples. Unpaired mode does not
model repeated measures, batch covariates or technical replicates. Repeated
specimens from one subject must not be declared independent simply to obtain
more replicates. After an exclusion, each condition must still meet the
configured minimum sample count (default three). Consequently, a complete
leave-one-sample-out analysis needs at least four samples in each condition
when that minimum remains three.

The paired remodeling classes A–F and within-patient directional support have
paired definitions. Unpaired results instead report independent-group effects
and their own method/class labels; they must not be interpreted as matched
patient trajectories. Likewise, unpaired detection/prevalence differences are
condition-associated changes, not observed within-patient gain/loss events.
See [docs/METHODS.md](docs/METHODS.md) for the statistical definitions.

## Exclusion sensitivity and leave one out

Every run includes the main cohort. Sensitivity analyses are optional:

| `sensitivity.mode` | Paired mode | Unpaired mode |
|---|---|---|
| `none` | Main cohort only | Main cohort only |
| `selected` | Omit each ID in `sensitivity.patients` separately | Omit each ID in `sensitivity.samples` separately |
| `all` | Omit every included subject once | Omit every included sample once |

Paired example:

```yaml
sensitivity:
  mode: selected
  patients: [ptn2]
  analysis: full
```

Unpaired example, using IDs from that study's sample sheet:

```yaml
sensitivity:
  mode: selected
  samples: [sample_03, sample_07]
  analysis: full
```

Use `mode: all` for the rolling leave-one-out analysis. Each omission is a
separate sensitivity cohort; this does not exclude all listed units together
or enumerate every possible subset. For a joint exclusion, define a named
cohort with `exclude_subjects` in paired mode or `exclude_samples` in unpaired
mode. The first cohort must contain all included units. Explicit named cohorts
are always included; `mode` controls the additional automatic omissions. Keep
only the main cohort definition when using `mode: none` for a run with no exclusions.

Full leave-one-out can be computationally expensive. With N subjects or
independent samples it adds N cohort analyses to the main analysis. Each cohort
repeats the applicable resampling/inference and DTU filtering/fitting. Leave it
disabled unless needed; start with named exclusions when possible. Compatible
completed work can resume from caches. Exclusions reuse the main cohort's
common-depth target and retained-isoform universe. Cache identities include
the design, included samples, settings and random seeds.

Invalid exclusions fail validation; support requirements and group-size minima
are not silently reduced. In paired mode, omitting a subject removes both
members. In unpaired mode, omitting a sample removes only that sample.

The optional paired observed-only influence check
(`sensitivity.observed_leave_one_patient_out.enabled`) is separate from full
leave-one-out. It does not repeat bootstrap or cohort inference and is not
available as a substitute for full unpaired sensitivity analysis.

## Run and validate

```sh
Rscript scripts/run_tests.R
Rscript scripts/run_smoke_test.R
Rscript scripts/run_smoke_test.R --unpaired --leave-one-out
Rscript scripts/run_analysis.R --config=config/analysis.yaml --validate-only
Rscript scripts/run_analysis.R --config=config/analysis.yaml
```

The smoke test uses synthetic files and reduced settings; it does not reproduce
the full study or establish production calibration. Add `--include-dtu` for
its DTU branch when the required native packages are installed. To run selected
stages, use `--stages=import,qc`; to deliberately bypass compatible caches use
`--force`. Force can repeat expensive work.

The paired study settings use 100 common-depth rarefactions, 5,000 paired
bootstrap resamples, 1,000 split-library repetitions and 20,000 reference null
draws, with optional independent refinement. DTU uses full-depth counts and
reports FDR thresholds 0.05 and 0.10. These paired production settings and
class definitions do not automatically define the unpaired inference.

## Outputs and resume

The YAML `output.directory` controls the output destination. Each stage has
separate tables, figures, audit records and reusable objects. The directory
names are identifiers; they do not change with the comparison labels.

| Folder | Contents |
|---|---|
| `M01_import` | Local-to-harmonized transcript map, counts, annotations and identity/count-conservation audits |
| `M02_gain_loss` | Paired gain/loss or unpaired condition-level detection/prevalence and depth-control results |
| `M03_qc` | Depth, detection, structural composition and design-appropriate technical summaries |
| `M04_remodeling` | Gene-level composition effects, method-specific inference and sensitivity outputs |
| `M06_dtu` | Frozen tested families, model-fit audits, gene/isoform P values, FDR and stageR results |
| `runs/<run_id>` | Resolved settings, included samples/cohorts, versions, logs and progress |

Missing annotation and unavailable tests retain their reasons and denominators.
They are not converted into zero effects or tested negatives. ORF prediction,
protein consequences, DNA integration and ASE/ASTU are outside this package.
Enrichment outputs depend on the selected design and available annotation
resources; inspect branch status records rather than assuming that every
optional analysis ran.

The startup banner identifies the package version. `run.log`, `progress.tsv`
and `progress_latest.tsv` record work completed and active stages. Progress
percentages count completed stages rather than predict elapsed/remaining time.
Keep the run records when reporting a problem. A completion marker is written
only after the applicable completion checks pass.

Re-running the same configuration reuses compatible completed work. A changed
design, sample membership, threshold, seed or numerical contract invalidates
the affected cache. Figure styling is separate from statistical calculations.

## Figures and replotting

Default exports use DejaVu Sans with 16-point titles/axis titles and 14-point
labels, subtitles and legends. PDF, SVG, PNG and JPEG are exported, with
600-dpi rasters. Meaningful individual panels accompany composites.

Overlap plots use exact membership dots and count bars. Enrichment plots wrap
labels with per-term spacing, use integer gene-count legends and visible
minimum point sizes, and start fraction axes at zero. Nominal-only pathway
findings retain that label. Configuration controls include
`reporting.maximum_upset_intersections`, `reporting.enrichment_label_width`
and `reporting.minimum_enrichment_point_size`.

Saved results can be replotted without repeating statistical analyses. For a
repository checkout, load `run_replot.R`, then call `replot_isoremodel()` with
the source and destination paths appropriate to the saved run. Use a new
output folder outside the source trees. The result includes figure indices,
review PDF volumes and a compact review ZIP. The source results are preserved.
Inspect the complete exported collection at the intended display size.

## Reproducibility and citation records

After a successful analysis:

```r
source("scripts/capture_environment.R")
capture_rna_environment("environment", project_root = ".")
```

This captures installed versions, session/OS/library information and, when
requested, a real `renv.lock`. Also retain input/reference checksums and the
upstream PacBio tool versions. A successful installation or small smoke test
does not establish reproduction of the study's results or validate a new
biological dataset. Compare eligibility, count accounting, effects, P/FDR,
stageR and sensitivity results when assessing reproducibility.

The reference paired study fit was checked for agreement of gene/isoform
P values, FDR calls and stageR results. New unpaired analyses use a different
biological design and inference, and require their own validation. See
[docs/DELIVERY_VALIDATION.md](docs/DELIVERY_VALIDATION.md) for the checks actually
executed for this release and [docs/VALIDATION.md](docs/VALIDATION.md) for the
analysis-machine validation procedure.

## License

GPL version 3 or later. See `LICENSE` and `NOTICE`. The numerical adapter
contains adaptations of upstream DRIMSeq helpers and retains the required
version guard. Reference genome, annotation and pathway resources have their
own provenance and terms. For GitHub distribution, include source, templates,
tests, documentation and tested environment records, but not raw patient files
or large output caches.
