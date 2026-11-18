# Rscript scripts/refit_DTU.R [--config=FILE] [--previous-output=DIR]
#   [--original-dtu=DIR] [--output-root=DIR]
args <- commandArgs(trailingOnly = TRUE)
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(entry)) stop("Use Rscript for scripts/refit_DTU.R.", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", entry[[1L]]), winslash = "/", mustWork = TRUE)))
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/refit_DTU.R [--config=FILE] [--previous-output=DIR] [--original-dtu=DIR] [--output-root=DIR]\n")
} else {
  if (any(!grepl("^--(config|previous-output|original-dtu|output-root)=", args))) stop("Unknown command-line argument.")
  get <- function(key) {
    x <- grep(paste0("^--", key, "="), args, value = TRUE)
    if (length(x) > 1L) stop("Repeated argument: ", key)
    if (length(x)) sub(paste0("^--", key, "="), "", x) else NULL
  }
  values <- list(config_file = get("config"), previous_output = get("previous-output"), original_dtu = get("original-dtu"), output_root = get("output-root"))
  values <- values[!vapply(values, is.null, logical(1))]
  source(file.path(root, "run_DTU_refit.R"), local = TRUE)
  do.call(refit_dtu_rng, values)
}
