# Bounded end-to-end checks for provenance, statistical cache identity, and
# retained failures. Only the real import stage runs; no fitting or rendering.
test_rna_runner <- function() {
  rna_require(c("data.table", "digest", "readr", "yaml"))
  root <- .repo_root
  work <- tempfile("rna_runner_test_")
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixture <- make_rna_fixture(file.path(work, "pigeon"), n_pairs = 6L)
  cfg <- yaml::read_yaml(file.path(root, "config", "analysis.yaml"), eval.expr = FALSE)
  cfg$inputs$sample_metadata <- fixture$metadata_file
  cfg$inputs$root_directory <- dirname(fixture$metadata_file)
  cfg$output$directory <- file.path(work, "workflow")
  cfg$references <- list(annotation_gtf = NULL, gene_id_map = NULL,
    pathway_term2gene = NULL, pathway_term_names = NULL)
  cfg$enrichment$enabled <- FALSE
  cfg$enrichment$reference_cache_directory <- file.path(work, "references")
  cfg$execution$profile <- "production"
  cfg$execution$resume <- TRUE
  cfg$execution$force_refit <- FALSE
  cfg$execution$strict_completion <- TRUE
  cfg$reporting$enabled <- FALSE
  cfg$stages <- as.list(setNames(rep(TRUE, 5L), c("import", "gain_loss", "qc", "remodeling", "dtu")))
  yaml_file <- file.path(work, "analysis.yaml")
  yaml::write_yaml(cfg, yaml_file)
  previous <- getOption("IsoRemodel.code_root")
  options(IsoRemodel.code_root = root)
  on.exit(options(IsoRemodel.code_root = previous), add = TRUE)

  validation <- run_rna_analysis(yaml_file, stages = "import", validate_only = TRUE)
  stopifnot(isTRUE(validation$valid), !dir.exists(cfg$output$directory))
  # The manifest helper also accepts a data.table without treating column names
  # as a data.table expression. Import itself must leave caller metadata intact.
  parsed <- read_rna_config(yaml_file, validate_inputs = FALSE)
  manifest_dt <- data.table::as.data.table(fixture$samples)
  files <- rna_run_inputs_manifest(parsed, manifest_dt)
  stopifnot(nrow(files) >= 25L, all(nchar(files$sha256) == 64L))
  before <- serialize(fixture$samples, NULL)
  import_pigeon(parsed, fixture$samples, file.path(work, "direct_import"))
  stopifnot(identical(before, serialize(fixture$samples, NULL)))

  progress_options <- c("IsoRemodel.log_file", "IsoRemodel.progress_file",
    "IsoRemodel.progress_latest_file", "IsoRemodel.run_id", "IsoRemodel.run_start")
  options_before <- stats::setNames(lapply(progress_options, getOption), progress_options)
  first <- run_rna_analysis(yaml_file, stages = "import")
  options_after <- stats::setNames(lapply(progress_options, getOption), progress_options)
  stopifnot(identical(options_before, options_after),
    identical(first$run_record$package_name, "IsoRemodel"),
    identical(first$run_record$package_version, rna_package_version()))
  first_progress_path <- file.path(first$run_record$run_directory, "progress.tsv")
  first_latest_path <- file.path(first$run_record$run_directory, "progress_latest.tsv")
  stopifnot(file.exists(first_progress_path), file.exists(first_latest_path))
  first_progress <- utils::read.delim(first_progress_path, stringsAsFactors = FALSE)
  first_latest <- utils::read.delim(first_latest_path, stringsAsFactors = FALSE)
  required_progress <- c("timestamp_utc", "run_id", "stage", "stage_index", "stage_total",
    "status", "elapsed_seconds", "detail")
  stopifnot(all(required_progress %in% names(first_progress)), nrow(first_progress) >= 3L,
    all(first_progress$run_id == first$run_record$run_id),
    all(grepl("Z$", first_progress$timestamp_utc)),
    all(!is.na(as.POSIXct(first_progress$timestamp_utc, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
    all(first_progress$elapsed_seconds >= 0), nrow(first_latest) == 1L,
    identical(first_latest$status, "partial_complete"))
  if ("stages_percent" %in% names(first_progress))
    stopifnot(all(first_progress$stages_percent[first_progress$status == "starting"] < 100))
  first_log <- readLines(file.path(first$run_record$run_directory, "run.log"), warn = FALSE)
  stopifnot(any(grepl(paste0("IsoRemodel v", rna_package_version()), first_log, fixed = TRUE)))
  stopifnot(identical(first$run_record$status, "partial_complete"),
    !isTRUE(first$run_record$completed),
    !isTRUE(first$run_record$stage_records$import$reused),
    !file.exists(file.path(cfg$output$directory, "COMPLETE")))
  first_key <- first$run_record$stage_records$import$key
  second <- run_rna_analysis(yaml_file, stages = "import")
  stopifnot(isTRUE(second$run_record$stage_records$import$reused),
    identical(second$run_record$stage_records$import$key, first_key),
    !identical(first$run_record$run_directory, second$run_record$run_directory))
  reused_progress <- utils::read.delim(file.path(second$run_record$run_directory, "progress.tsv"), stringsAsFactors = FALSE)
  stopifnot(any(reused_progress$status == "reused"),
    all(reused_progress$stage_index[reused_progress$status == "reused" & reused_progress$stage == "import"] == 1L),
    all(reused_progress$run_id == second$run_record$run_id),
    identical(utils::tail(reused_progress$status, 1L), "partial_complete"))
  cfg$reporting$font_family <- "a font that does not need installation for import"
  cfg$reporting$font_size_pt <- 22
  cfg$reporting$dpi <- 144
  yaml::write_yaml(cfg, yaml_file)
  third <- run_rna_analysis(yaml_file, stages = "import")
  stopifnot(isTRUE(third$run_record$stage_records$import$reused),
    identical(third$run_record$stage_records$import$key, first_key))

  # A stage warning must be logged as WARNING while retaining normal R warning
  # behavior. The outer suppression only keeps this intentional fixture quiet.
  engine_env <- environment(run_rna_analysis)
  original_import <- get("import_pigeon", envir = engine_env, inherits = FALSE)
  warning_import <- function(...) {
    warning("runner warning audit fixture", call. = FALSE)
    original_import(...)
  }
  assign("import_pigeon", warning_import, envir = engine_env)
  on.exit(assign("import_pigeon", original_import, envir = engine_env), add = TRUE)
  warned <- suppressWarnings(run_rna_analysis(yaml_file, stages = "import", force = TRUE))
  assign("import_pigeon", original_import, envir = engine_env)
  warning_log <- readLines(file.path(warned$run_record$run_directory, "run.log"), warn = FALSE)
  stopifnot(any(grepl("[WARNING]", warning_log, fixed = TRUE) &
    grepl("runner warning audit fixture", warning_log, fixed = TRUE)))

  source_file <- fixture$samples$classification_file[[1L]]
  original <- readBin(source_file, "raw", n = file.info(source_file)$size)
  bad <- data.table::fread(source_file, data.table = FALSE, check.names = FALSE)
  bad[["FL.sample"]][1L] <- -1L
  utils::write.table(bad, source_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  error <- tryCatch({run_rna_analysis(yaml_file, stages = "import"); NULL}, error = function(e) e)
  stopifnot(inherits(error, "error"))
  failed <- readRDS(file.path(cfg$output$directory, "latest_run_manifest.rds"))
  stopifnot(identical(failed$status, "failed"), identical(failed$failed_stage, "import"),
    !isTRUE(failed$completed), file.exists(file.path(failed$run_directory, "failure_audit.rds")),
    file.exists(file.path(failed$run_directory, "failure.txt")),
    file.exists(file.path(failed$run_directory, "stage_completion.tsv")),
    !file.exists(file.path(cfg$output$directory, "COMPLETE")))
  failed_progress <- utils::read.delim(file.path(failed$run_directory, "progress.tsv"), stringsAsFactors = FALSE)
  failed_latest <- utils::read.delim(file.path(failed$run_directory, "progress_latest.tsv"), stringsAsFactors = FALSE)
  failed_log <- readLines(file.path(failed$run_directory, "run.log"), warn = FALSE)
  stopifnot(any(failed_progress$status == "failed"), identical(failed_latest$status, "failed"),
    any(grepl("[ERROR]", failed_log, fixed = TRUE) & grepl("FAILED at import", failed_log, fixed = TRUE)))
  if ("stages_percent" %in% names(failed_latest)) stopifnot(failed_latest$stages_percent < 100)
  # Restore byte-identical raw input. The last complete cache is still valid;
  # a failed new attempt never overwrites it with a claim of completion.
  writeBin(original, source_file)
  recovered <- run_rna_analysis(yaml_file, stages = "import")
  stopifnot(isTRUE(recovered$run_record$stage_records$import$reused),
    identical(recovered$run_record$stage_records$import$key, first_key),
    identical(recovered$run_record$status, "partial_complete"))
  invisible(TRUE)
}
test_rna_runner()
