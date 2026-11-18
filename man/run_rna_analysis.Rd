\name{run_rna_analysis}
\alias{run_rna_analysis}
\title{Run the configured paired or unpaired RNA workflow}
\usage{
run_rna_analysis(config_file = "config/analysis.yaml", stages = NULL,
  validate_only = FALSE, force = FALSE)
}
\arguments{
\item{config_file}{Path to the YAML analysis configuration. Relative YAML paths
are interpreted from the configuration directory.}
\item{stages}{Optional character vector of analysis stages. The supported stage
names are \code{import}, \code{qc}, \code{gain_loss}, \code{remodeling} and
\code{dtu}. The runner resolves the required upstream inputs. The default
\code{NULL} uses the enabled configuration stages.}
\item{validate_only}{If \code{TRUE}, validate configuration, metadata and input
paths without running the analyses.}
\item{force}{If \code{TRUE}, recompute requested work rather than accepting a
compatible completed stage. This may repeat expensive stochastic analyses.}
}
\description{
Imports Pigeon files, harmonizes splice chains and runs configured quality
checks, depth-controlled detection changes, isoform remodeling and differential
transcript usage. The YAML selects a paired or independent-sample design and
optional subject or sample exclusions, including full leave-one-out runs. ORF
prediction and protein-consequence analyses are outside this workflow.
}
\value{The pipeline run record, returned invisibly. Analysis tables, figures,
completion audits and reusable objects are written under the configured output
directory. See the repository README for stage-specific outputs.}
\details{
The paper settings can require substantial computation. Run the supplied
bounded tests and synthetic smoke test before beginning a new full-data run.
The numerical DTU adapter requires DRIMSeq version 1.38.0.

Paired designs require one sample per condition per subject. Unpaired designs
require independent biological samples and use condition-only DTU, within-group
bootstrap and condition-label permutation. Paired architecture classes are not
assigned to unpaired results. See \code{config/analysis_unpaired.yaml} for an
independent-sample template. All exclusion cohorts must retain the configured
minimum support; these thresholds are not reduced automatically.

The startup banner identifies IsoRemodel and its version. Timestamped run
messages are retained in \code{runs/<run_id>/run.log}; the same directory
contains \code{progress.tsv} and \code{progress_latest.tsv}. Progress records
describe stages and elapsed time, not an estimate of remaining runtime.
}
\examples{
\dontrun{
run_rna_analysis("config/analysis.yaml", validate_only = TRUE)
run_rna_analysis("config/analysis.yaml", stages = c("import", "qc"))
run_rna_analysis("config/analysis.yaml")
}
}
