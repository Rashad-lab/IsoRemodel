\name{replot_isoremodel}
\alias{replot_isoremodel}
\title{Export RNA figures from completed paired or unpaired results}
\usage{
replot_isoremodel(
  config_file = NULL,
  source_dir = "F:/IsoRemodel/outputs",
  dtu_source_dir = "F:/IsoRemodel/DTU_refit/20261026T043514_45632/M06_dtu",
  output_dir = NULL, formats = NULL, dpi = NULL, strict = TRUE
)
}
\arguments{
\item{config_file}{YAML configuration providing reporting settings. By default,
use the configuration in the loaded repository or installed package.}
\item{source_dir}{Completed early RNA result tree containing saved plot objects.}
\item{dtu_source_dir}{Verified isolated paired DTU refit or a completed native
v0.2.0-or-newer DTU result directory. The design must match the configuration. Old DTU
figures in \code{source_dir} are not substituted for this source.}
\item{output_dir}{New empty destination outside both source trees. The default
is a new timestamped folder under \code{replot_v0.2.1} beside \code{source_dir}.}
\item{formats}{Optional figure-format override. The configured defaults are
PDF, SVG, PNG and JPEG.}
\item{dpi}{Optional raster-resolution override; the configured default is 600.}
\item{strict}{Stop after retaining review outputs when required exports or
integrity checks fail.}
}
\description{Re-exports saved early RNA plots and reconstructs DTU graphics
from the verified completed DTU results and existing enrichment tables.
The operation does not refit models, rerun stageR, resample or recompute
enrichment tests.}
\details{Exports use compact DejaVu Sans, approved count/membership overlap
layouts, wrapped enrichment rows, integer count legends, visible minimum
enrichment point sizes, fraction axes starting at zero and repaired filtering
counts. Meaningful individual panels are retained alongside composites.
Independent-sample results use independent effect/usage reports and do not
receive paired architecture or patient-support reports.
The source results are not overwritten. The new output includes figure
indices, logs, progress records, source hashes, review PDF volumes and a review
ZIP. A final visual audit of the complete study output remains required.}
\value{Invisibly, a list identifying the new output and review files.}
\examples{
\dontrun{
result <- replot_isoremodel()
}
}
