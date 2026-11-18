\name{replot_rna_figures}
\alias{replot_rna_figures}
\title{Export saved RNA figures with current reporting settings}
\usage{
replot_rna_figures(config_file = "config/analysis.yaml", source_dir = NULL,
  output_dir = NULL)
}
\arguments{
\item{config_file}{Path to the YAML configuration containing reporting settings.}
\item{source_dir}{Existing output directory containing saved \code{*.plot.rds}
files. The default is the configured analysis output directory.}
\item{output_dir}{Empty destination outside the source directory tree. The
default is a timestamped sibling of the source directory.}
}
\description{Re-exports saved plot objects using the common typography and
graphics devices. This function does not rerun statistical analyses. For trees
containing DTU results, use \code{replot_isoremodel()} to verify and use the
corrected saved DTU results.}
\value{Invisibly, an audit data frame identifying exported figures and any
failures. New files and their audit are written to the destination.}
\details{Successful exports are retained if another figure fails. Review
\code{replot_audit.tsv} if the operation reports incomplete export. The original
outputs are not overwritten.}
\examples{
\dontrun{
replot_rna_figures("config/analysis.yaml", output_dir = "RNA_restyled")
}
}
