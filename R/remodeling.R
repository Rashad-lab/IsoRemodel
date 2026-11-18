#' Gene-level isoform remodeling from harmonized integer full-length counts.
#'
#' The vendored engine retains Module 04 v0.6.0 metric definitions, paired
#' bootstrap, complementary split null, additive finite-count REML model,
#' independent refinement, GPD diagnostics, patient LRT and truth-known oracles.
#' Counts are never transformed before descriptive fraction metrics.
run_remodeling <- function(inputs, cfg, cohorts, stage_dir) {
  required <- c("data.table","digest","dplyr","forcats","ggplot2","patchwork",
                "purrr","readr","scales","stringr","tibble","tidyr")
  rna_require(required)
  root <- cfg$.code_root %||% getOption("IsoRemodel.code_root") %||% getwd()
  engine_dir <- file.path(root, "inst", "engines")
  if (!file.exists(file.path(engine_dir, "remodeling_engine.R"))) {
    installed <- system.file("engines", package = "IsoRemodel")
    if (nzchar(installed)) engine_dir <- installed
  }
  if (!file.exists(file.path(engine_dir, "remodeling_engine.R"))) {
    stop("Cannot locate the bundled remodeling engine. Set cfg$.code_root to the repository root.")
  }
  source_files <- c(wrapper=file.path(root,"R","remodeling.R"),
    configuration=file.path(engine_dir,"remodeling_config.R"),
    statistics=file.path(engine_dir,"remodeling_functions.R"),
    orchestration=file.path(engine_dir,"remodeling_engine.R"))
  source_sha <- vapply(source_files,function(path) {
    if (file.exists(path)) digest::digest(file=path,algo="sha256") else NA_character_
  },character(1))
  # Installed packages may contain byte-compiled R code rather than the source wrapper.
  if (is.na(source_sha[["wrapper"]])) source_sha[["wrapper"]] <-
    digest::digest(body(sys.function()),algo="sha256")
  cfg$.remodeling_code_sha256 <- source_sha
  # Namespace imports are local; running this stage does not attach packages or
  # overwrite helpers in the caller's global environment.
  env <- new.env(parent = environment())
  for (pkg in required) {
    ns <- asNamespace(pkg)
    exports <- getNamespaceExports(pkg)
    for (nm in exports) assign(nm, getExportedValue(pkg,nm), envir=env)
  }
  sys.source(file.path(engine_dir,"remodeling_functions.R"), envir=env)
  sys.source(file.path(engine_dir,"remodeling_config.R"), envir=env)
  env$flatten_config_for_record <- function(x, prefix="") {
    rows <- list()
    for (nm in names(x)) {
      key <- if(nzchar(prefix)) paste(prefix,nm,sep=".") else nm
      value <- x[[nm]]
      if (is.list(value)) rows[[length(rows)+1L]] <- env$flatten_config_for_record(value,key)
      else rows[[length(rows)+1L]] <- data.frame(parameter=key,value=paste(value,collapse=","))
    }
    dplyr::bind_rows(rows)
  }
  sys.source(file.path(engine_dir,"remodeling_engine.R"), envir=env)
  result <- env$run_remodeling_engine(inputs,cfg,cohorts,stage_dir)
  result$enrichment <- run_remodeling_enrichment(inputs, result, cfg, stage_dir)
  result$complete <- isTRUE(result$complete) && isTRUE(result$enrichment$complete)
  result
}

# Internal function loader used by the small mathematical regression tests.
.rna_remodeling_test_environment <- function(code_root, cfg=list()) {
  env <- new.env(parent=environment())
  for (pkg in c("dplyr","purrr","tibble","tidyr")) {
    if (!requireNamespace(pkg,quietly=TRUE)) stop("Tests need ",pkg)
    for(nm in getNamespaceExports(pkg)) assign(nm,getExportedValue(pkg,nm),envir=env)
  }
  sys.source(file.path(code_root,"inst","engines","remodeling_functions.R"),envir=env)
  sys.source(file.path(code_root,"inst","engines","remodeling_config.R"),envir=env)
  env$cfg <- env$rna_remodeling_engine_config(cfg)
  env
}
