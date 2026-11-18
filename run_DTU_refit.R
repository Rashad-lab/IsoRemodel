# Source this file; no analysis runs until refit_dtu_rng() is called.
.isoremodel_refit_entry <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
if (is.null(.isoremodel_refit_entry)) stop("Use source() with the full path to run_DTU_refit.R.", call. = FALSE)
refit_dtu_rng <- local({
  root <- dirname(normalizePath(.isoremodel_refit_entry, winslash = "/", mustWork = TRUE))
  env <- new.env(parent = globalenv())
  previous <- getOption("IsoRemodel.code_root")
  on.exit(options(IsoRemodel.code_root = previous), add = TRUE)
  sys.source(file.path(root, "R", "load.R"), envir = env)
  env$load_rna_code(root, envir = env, startup = FALSE)
  function(...) {
    previous <- getOption("IsoRemodel.code_root")
    on.exit(options(IsoRemodel.code_root = previous), add = TRUE)
    options(IsoRemodel.code_root = root)
    env$refit_dtu_rng(...)
  }
})
rm(.isoremodel_refit_entry)
message("IsoRemodel DTU confirmation loaded. Call refit_dtu_rng() to start the targeted refit.")
