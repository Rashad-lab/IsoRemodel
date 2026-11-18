# Standard startup message for an installed R package.
.onAttach <- function(libname, pkgname) {
  rna_banner(startup = TRUE)
}
