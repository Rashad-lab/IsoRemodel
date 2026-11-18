# Reference resources and citations

## Annotation

The study uses the configured GRCh38/GENCODE v49 annotation. The annotation
file is not bundled. Record its source/release and content checksum. Gene-label
mapping must preserve original labels and expose unresolved/ambiguous maps.

The core importer uses the coordinates already supplied by Pigeon. It does not
rerun alignment, transcript collapsing or Pigeon classification. Record those
upstream software/reference versions in the public data-processing methods.

## Pathway sets

The original remodeling follow-up and DTU follow-up use distinct GO resources:
the former uses annotation-database GO membership and the latter uses the
MSigDB GO biological-process collection. Retain that distinction for study
reproduction. Hallmark and KEGG collections are also used where configured.

Export and preserve the term-to-gene membership and term-name tables actually
used, together with package/database versions and ID-mapping provenance.
Updates to current remote collections are a changed reference, not merely a
new plotting run. Human-specific annotation is the study preset; other species
need their own explicit mapping and pathway resources and validation.

The eligible universe is analysis-specific. Standard ORA must not replace the
finite-tested/evaluable foreground's background with every annotated gene.
Matched enrichment uses deterministic unique controls and records matching
coverage and standardized mean differences. A single matched set is not a
permutation-based biological null.

## Software references

Use citation() from the installed versions and preserve the generated records.
Principal method references include:

- Nowicka M, Robinson MD. DRIMSeq: a Dirichlet-multinomial framework for
  multivariate count outcomes in genomics. F1000Research. 2016;5:1356.
  doi:10.12688/f1000research.8900.2
- Van den Berge K, Soneson C, Robinson MD, Clement L. stageR: a general
  stage-wise method for controlling the gene-level false discovery rate in
  differential expression and differential transcript usage. Genome Biology.
  2017;18:151. doi:10.1186/s13059-017-1271-0

Primary software documentation reviewed for this release:

- https://isoseq.how/classification/pigeon-output.html
- https://bioconductor.org/packages/DRIMSeq/
- https://bioconductor.org/install/
- https://bioconductor.org/about/release-announcements/
- https://rstudio.github.io/renv/reference/snapshot.html

The final article citation and an archived software-release DOI can be added
to CITATION.cff when available. The bioRxiv DOI below identifies the preprint,
not an archived software release.


## Associated study and sequencing data

Rashad S, Ando D, Yamashita S, Shimoda Y, Kanamori M, Endo H, Niizuma K.
Long read multi-omics sequencing reveals DNA-to-RNA evolutionary remodeling
trajectories of recurrent astrocytoma. bioRxiv. 2026;2026.07.12.738108.
Posted July 14, 2026. https://doi.org/10.64898/2026.07.12.738108

- BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1492232
- RNA run: https://www.ncbi.nlm.nih.gov/sra/?term=SRR39530479
- Upstream PacBio Iso-Seq/Pigeon workflow: https://isoseq.how/

The supplied SRA submission record labels SRR39530479 as merged long-read RNA
sequencing runs and lists three source files. It marks the submission as to be
released on2027-08-31; check NCBI for current access. DNA records in the same
BioProject are not inputs to the RNA-only workflow. The paired study's
sample-to-cDNA-barcode assignments are supplied in the README and sample sheet.
