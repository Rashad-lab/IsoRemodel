\name{refit_dtu_rng}
\alias{refit_dtu_rng}
\title{Confirm Publication DTU Results with Explicit RNG Settings}
\usage{
refit_dtu_rng(
  config_file = "F:/IsoRemodel/IsoRemodel_v0.1.1/config/analysis.yaml",
  previous_output = "F:/IsoRemodel/outputs",
  original_dtu = "F:/PacBio revision analysis/outputs/M06_DTU_robustness",
  output_root = "F:/IsoRemodel/DTU_refit", strict = TRUE
)
}
\arguments{
\item{config_file}{Existing study YAML. Raw Pigeon files are not reopened.}
\item{previous_output}{Completed IsoRemodel output root with production import and remodeling stage caches.}
\item{original_dtu}{Original publication DTU output directory containing fitted-model caches.}
\item{output_root}{Separate parent for a newly timestamped confirmation run.}
\item{strict}{Stop after retaining the review ZIP if refitting fails or original final results are not reproduced.}
}
\description{
Loads completed import and remodeling handoffs, validates ordered model inputs,
and refits only differential transcript usage using the publication RNG
algorithm and cohort seeds. Tests and stageR are compared with the original
and previous IsoRemodel caches at both FDR thresholds. Figure exports are
suppressed and existing result directories remain read only.
}
\value{Invisibly, a list containing the run status, new run directory, and review ZIP path.}
\details{
This focused confirmation requires production settings, DRIMSeq 1.38.0,
stageR 1.32.0, master seed 260605, L'Ecuyer-CMRG with Inversion and Rejection,
and both the all-six and patient-2-excluded cohorts. Raw unavailable tests are
compared separately from conservative P values. No upstream analyses run.
}
