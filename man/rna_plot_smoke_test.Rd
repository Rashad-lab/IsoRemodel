\name{rna_plot_smoke_test}
\alias{rna_plot_smoke_test}
\title{Check RNA figure fonts and export devices}
\usage{
rna_plot_smoke_test(output_root, cfg = list())
}
\arguments{
\item{output_root}{Directory for the small synthetic figure exports.}
\item{cfg}{Optional configuration list, including reporting settings. The
default uses the shared style defaults.}
}
\description{Creates a small synthetic paired plot and separate facet panels to
check the bundled font and PDF, SVG, PNG and JPEG export. It does not read
patient data or fit statistical models.}
\value{The output directory, returned invisibly.}
\examples{
\dontrun{
rna_plot_smoke_test("RNA_font_check")
}
}
