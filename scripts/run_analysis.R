# Run from any working directory using Rscript /path/to/scripts/run_analysis.R.
args <- commandArgs(trailingOnly = TRUE)
script_file <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_file)) stop("Use Rscript for this command; source load.R for interactive use.", call. = FALSE)
project <- dirname(dirname(normalizePath(sub("^--file=", "", script_file[[1L]]), winslash = "/", mustWork = TRUE)))
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/run_analysis.R [--config=FILE] [--stages=import,gain_loss,qc,remodeling,dtu] [--validate-only] [--force]\n",
      "Without --stages, all enabled stages in YAML run. Relative config paths use your working directory.\n",
      "Default config: <repository>/config/analysis.yaml. --force ignores statistical stage caches.\n", sep = "")
} else {
  known <- grepl("^--(config|stages)=", args) | args %in% c("--validate-only", "--force")
  if (any(!known)) stop("Unknown command-line argument: ", paste(args[!known], collapse = ", "), call. = FALSE)
  value <- function(flag, default = NULL) {
    hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (length(hit) > 1L) stop("Repeated --", flag, " argument.", call. = FALSE)
    if (length(hit)) sub(paste0("^--", flag, "="), "", hit[[1L]]) else default
  }
  source(file.path(project, "load.R"), local = TRUE)
  config <- value("config", file.path(project, "config", "analysis.yaml"))
  selected <- value("stages")
  if (!is.null(selected)) selected <- trimws(strsplit(selected, ",", fixed = TRUE)[[1L]])
  run_rna_analysis(config_file = config, stages = selected,
    validate_only = "--validate-only" %in% args, force = "--force" %in% args)
}
