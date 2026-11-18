# Exports previously saved plots with current YAML reporting settings.
# It never calls the statistical analysis engines or modifies source figures.
args <- commandArgs(trailingOnly = TRUE)
script_file <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_file)) stop("Use Rscript for this command; source load.R for interactive use.", call. = FALSE)
project <- dirname(dirname(normalizePath(sub("^--file=", "", script_file[[1L]]), winslash = "/", mustWork = TRUE)))
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/replot.R [--config=FILE] [--source=DIR] [--output=NEW_DIR]\n",
      "Default source: YAML output directory. Default output: a new timestamped sibling directory.\n",
      "All *.plot.rds objects are exported with the current reporting settings; statistics are unchanged.\n", sep = "")
} else {
  known <- grepl("^--(config|source|output)=", args)
  if (any(!known)) stop("Unknown command-line argument: ", paste(args[!known], collapse = ", "), call. = FALSE)
  value <- function(flag, default = NULL) {
    hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
    if (length(hit) > 1L) stop("Repeated --", flag, " argument.", call. = FALSE)
    if (length(hit)) sub(paste0("^--", flag, "="), "", hit[[1L]]) else default
  }
  source(file.path(project, "load.R"), local = TRUE)
  replot_rna_figures(config_file = value("config", file.path(project, "config", "analysis.yaml")),
    source_dir = value("source"), output_dir = value("output"))
}
