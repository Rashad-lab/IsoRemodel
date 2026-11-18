# Orchestration, content-addressed whole-stage resume, and separate plot export.
# Statistical engines retain their own fine-grained checkpoints and audit files.

rna_run_absolute <- function(path, base = getwd()) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) rna_stop("A nonempty path is required.")
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  normalizePath(if (absolute) path else file.path(base, path), winslash = "/", mustWork = FALSE)
}

rna_run_package_versions <- function(root) {
  description <- file.path(root, "DESCRIPTION")
  packages <- c("data.table", "digest", "readr", "yaml")
  if (file.exists(description)) {
    d <- read.dcf(description)
    for (field in intersect(c("Imports", "Suggests"), colnames(d))) {
      fields <- trimws(strsplit(d[1L, field], ",", fixed = TRUE)[[1L]])
      packages <- c(packages, sub("[[:space:]]*\\(.*$", "", fields))
    }
  }
  packages <- sort(unique(packages[nzchar(packages)]))
  data.frame(package = packages, version = vapply(packages, function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "Not installed")
  }, character(1)), stringsAsFactors = FALSE)
}

rna_run_file_manifest <- function(paths) {
  paths <- unique(as.character(paths[!is.na(paths) & nzchar(paths)]))
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  paths <- sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
  if (!length(paths)) return(data.frame(path = character(), bytes = numeric(), sha256 = character()))
  data.frame(path = paths, bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1)),
    stringsAsFactors = FALSE)
}

rna_run_inputs_manifest <- function(cfg, samples, references = FALSE) {
  columns <- intersect(c("classification_file", "junctions_file", "pigeon_summary_file", "filter_reasons_file"), names(samples))
  sample_frame <- as.data.frame(samples)
  paths <- c(cfg$.metadata_path, unlist(sample_frame[, columns, drop = FALSE], use.names = FALSE))
  if (references) {
    refs <- unlist(cfg$references[c("annotation_gtf", "gene_id_map", "pathway_term2gene", "pathway_term_names")], use.names = FALSE)
    refs <- refs[!is.na(refs) & nzchar(refs)]
    if (length(refs)) refs <- vapply(refs, rna_run_absolute, character(1), base = cfg$.config_dir)
    cache <- cfg$enrichment$reference_cache_directory
    if (!is.null(cache) && nzchar(cache)) {
      cache <- rna_run_absolute(cache, cfg$.config_dir)
      refs <- c(refs, list.files(cache, pattern = "frozen.*\\.(tsv|csv|rds|sqlite|db)(\\.gz)?$", full.names = TRUE, recursive = TRUE))
    }
    paths <- c(paths, refs)
  }
  rna_run_file_manifest(paths)
}

rna_run_code_manifest <- function(root) {
  files <- c(list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(root, "inst", "engines"), pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE),
    file.path(root, "DESCRIPTION"))
  out <- rna_run_file_manifest(files)
  out$path <- substring(out$path, nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L)
  out
}

rna_run_enabled_stages <- function(cfg) {
  order <- c("import", "gain_loss", "qc", "remodeling", "dtu")
  flags <- cfg$stages
  if (is.null(flags)) return(order)
  if (is.character(flags)) {
    bad <- setdiff(flags, order)
    if (length(bad)) rna_stop("Unknown configured stages: ", paste(bad, collapse = ", "))
    return(order[order %in% unique(c("import", flags))])
  }
  enabled <- vapply(order, function(stage) {
    flag <- flags[[stage]]
    if (is.list(flag)) flag <- flag$enabled
    if (is.null(flag)) TRUE else isTRUE(flag)
  }, logical(1))
  order[enabled | order == "import"]
}

rna_run_stage_complete <- function(result, stage) {
  if (identical(stage, "import")) {
    checks <- result$count_conservation
    conserved <- is.null(checks) || !nrow(checks) ||
      ("difference" %in% names(checks) && all(checks$difference == 0)) ||
      ("count_conserved" %in% names(checks) && all(checks$count_conserved))
    return(is.matrix(result$counts_matrix) && isTRUE(conserved))
  }
  if (!is.null(result$complete)) return(isTRUE(result$complete))
  if (!is.null(result$completed)) return(isTRUE(result$completed))
  FALSE
}

#' Run the configured paired or unpaired RNA workflow from raw Pigeon files.
#'
#' @param config_file Study YAML file.
#' @param stages Optional subset of import,gain_loss,qc,remodeling,dtu. Import is
#'   always available as an explicit dependency. A subset is recorded as partial.
#' @param validate_only Check configuration, metadata, and input availability.
#' @param force Ignore whole-stage caches and ask engines to refit where supported.
#' @return Invisibly, a list with run_record, inputs, and named stage results.
run_rna_analysis <- function(config_file = "config/analysis.yaml", stages = NULL,
                             validate_only = FALSE, force = FALSE) {
  rna_require(c("yaml", "digest", "readr"))
  cfg <- read_rna_config(config_file, validate_inputs = FALSE)
  samples <- as.data.frame(read_rna_metadata(cfg, validate_inputs = TRUE))
  cohorts <- rna_resolve_cohorts(samples, cfg)
  enabled <- rna_run_enabled_stages(cfg)
  if (is.null(stages)) stages <- enabled
  stages <- unique(as.character(stages))
  unknown <- setdiff(stages, c("import", "gain_loss", "qc", "remodeling", "dtu"))
  if (length(unknown)) rna_stop("Unknown requested stages: ", paste(unknown, collapse = ", "))
  if (!length(stages)) rna_stop("At least one stage must be requested.")
  stages <- c("import", "gain_loss", "qc", "remodeling", "dtu")[c("import", "gain_loss", "qc", "remodeling", "dtu") %in% unique(c("import", stages))]
  # Validate references only for requested branches. Import and RNA QC do not
  # require a local GTF or enrichment resources.
  if (any(stages %in% c("remodeling", "dtu")) && isTRUE(cfg$enrichment$enabled)) {
    reference_keys <- c("gene_id_map", "pathway_term2gene", "pathway_term_names")
    if ("dtu" %in% stages) reference_keys <- c("annotation_gtf", reference_keys)
    for (name in reference_keys) {
      path <- cfg$references[[name]]
      if (!is.null(path) && nzchar(path) && !file.exists(path))
        rna_stop("Configured reference file does not exist (references.", name, "): ", path)
    }
    if ("dtu" %in% stages && isTRUE(cfg$enrichment$matched_background) &&
        (is.null(cfg$references$annotation_gtf) || !nzchar(cfg$references$annotation_gtf)))
      rna_stop("Matched DTU enrichment requires references.annotation_gtf.")
  }
  profile <- tolower(cfg$execution$profile %||% "production")
  output <- cfg$.output_dir %||% rna_run_absolute(cfg$output$directory, cfg$.config_dir)
  root <- cfg$.code_root %||% getOption("IsoRemodel.code_root") %||% getwd()
  if (isTRUE(validate_only)) {
    rna_log("Validation passed: ", nrow(samples), " samples, ", length(unique(samples$patient)),
      if(isTRUE(cfg$comparison$paired)) " paired subjects, " else " independent samples, ", length(cohorts), " cohort(s); profile ", profile, ".")
    rna_log("Output directory: ", output)
    return(invisible(list(valid = TRUE, config = cfg, samples = samples,
      cohorts = cohorts, requested_stages = stages, output_directory = output)))
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  run_id <- paste0(gsub("[^0-9T]", "", format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC")), "_", Sys.getpid())
  run_dir <- file.path(output, "runs", run_id)
  # Subsecond clocks can still repeat on some systems. Keep every run distinct.
  suffix <- 0L
  while (dir.exists(run_dir)) {
    suffix <- suffix + 1L
    run_dir <- file.path(output, "runs", paste0(run_id, "_", suffix))
  }
  run_id <- basename(run_dir)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  option_names <- c("IsoRemodel.log_file", "IsoRemodel.progress_file",
    "IsoRemodel.progress_latest_file", "IsoRemodel.run_id", "IsoRemodel.run_start")
  old_options <- stats::setNames(lapply(option_names, getOption), option_names)
  options(IsoRemodel.log_file = file.path(run_dir, "run.log"),
    IsoRemodel.progress_file = file.path(run_dir, "progress.tsv"),
    IsoRemodel.progress_latest_file = file.path(run_dir, "progress_latest.tsv"),
    IsoRemodel.run_id = run_id, IsoRemodel.run_start = Sys.time())
  on.exit(options(old_options), add = TRUE)
  rna_banner(log = TRUE)
  rna_progress(0L, length(stages), "workflow", status = "starting", stage_index = 0L,
    detail = paste0("Preparing input fingerprints and run records; profile ", profile, "."))
  # A marker is a claim about the latest run, so retire it before starting.
  marker <- file.path(output, "COMPLETE")
  if (file.exists(marker)) {
    if (!file.copy(marker, file.path(run_dir, "previous_completion_marker.txt")))
      rna_stop("Cannot preserve the previous completion marker.")
    unlink(marker)
  }
  if (isTRUE(force)) {
    cfg$execution$force_refit <- TRUE
    cfg$dtu$force_refit <- TRUE
  }
  code_manifest <- rna_run_code_manifest(root)
  input_manifest <- rna_run_inputs_manifest(cfg, samples)
  versions <- rna_run_package_versions(root)
  scientific <- rna_scientific_config(cfg)
  sample_key_metadata <- unserialize(serialize(as.data.frame(samples), NULL))
  yaml::write_yaml(cfg[!startsWith(names(cfg), ".")], file.path(run_dir, "resolved_analysis.yaml"))
  rna_write_tsv(samples, file.path(run_dir, "resolved_samples.tsv"))
  cohort_manifest <- do.call(rbind, lapply(names(cohorts), function(id)
    data.frame(cohort_id = id, subject_id = as.character(cohorts[[id]]), stringsAsFactors = FALSE)))
  if(!isTRUE(cfg$comparison$paired))names(cohort_manifest)[2L] <- "sample_id"
  rna_write_tsv(cohort_manifest, file.path(run_dir, "resolved_cohorts.tsv"))
  cohort_design <- do.call(rbind,lapply(names(cohorts),function(id) {
    keep <- if(isTRUE(cfg$comparison$paired))samples$patient %in% cohorts[[id]] else samples$sample %in% cohorts[[id]]
    z <- samples[keep,,drop=FALSE]
    data.frame(cohort=id,design=if(isTRUE(cfg$comparison$paired))"paired"else"unpaired",
      n_samples=nrow(z),n_reference=sum(z$condition==cfg$comparison$reference),
      n_comparison=sum(z$condition==cfg$comparison$comparison),
      excluded_units=paste(setdiff(if(isTRUE(cfg$comparison$paired))unique(samples$patient)else samples$sample,cohorts[[id]]),collapse=";"))
  }))
  rna_write_tsv(cohort_design,file.path(run_dir,"cohort_design.tsv"))
  rna_write_tsv(code_manifest, file.path(run_dir, "source_code_content_sha256.tsv"))
  rna_write_tsv(input_manifest, file.path(run_dir, "input_content_sha256.tsv"))
  rna_write_tsv(versions, file.path(run_dir, "package_versions.tsv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
  record <- list(schema_version = "0.1.0", package_name = "IsoRemodel",
    package_version = rna_package_version(), run_id = run_id,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), status = "running",
    profile = profile, requested_stages = stages, enabled_stages = enabled,
    all_enabled_stages_requested = all(enabled %in% stages), completed = FALSE,
    output_directory = output, run_directory = run_dir,
    scientific_config_sha256 = rna_hash(scientific),
    input_sha256 = rna_hash(input_manifest), code_sha256 = rna_hash(code_manifest),
    R_version = R.version.string, RNGkind = RNGkind(), cohorts = cohorts,
    design = if(isTRUE(cfg$comparison$paired)) "paired" else "unpaired",
    stage_records = list())
  rna_atomic_rds(record, file.path(run_dir, "run_manifest.rds"))
  rna_log("Starting RNA analysis ", run_id, "; profile ", profile,
    "; stages ", paste(stages, collapse = ", "), ".")
  inputs <- NULL; results <- list(); keys <- list(); active_stage <- "setup"
  folders <- c(import = "M01_import", gain_loss = "M02_gain_loss", qc = "M03_qc",
               remodeling = "M04_remodeling", dtu = "M06_dtu")
  stage_key <- function(stage) {
    deps <- switch(stage, import = character(), gain_loss = "import",
      qc = c("import", "gain_loss"), remodeling = "import", dtu = c("import", "remodeling"))
    refs <- if (stage %in% c("remodeling", "dtu")) rna_run_inputs_manifest(cfg, samples, TRUE) else input_manifest
    rna_hash(list(schema = "whole-stage-cache-v1", stage = stage,
      scientific_config = scientific, samples = sample_key_metadata, cohorts = cohorts,
      dependencies = keys[intersect(deps, names(keys))],
      input_contents = refs, code_contents = code_manifest,
      R = R.version.string, package_versions = versions,
      RNGkind = if (identical(stage, "dtu")) dtu_rng_spec(cfg$dtu$rng) else RNGkind()))
  }
  persist_record <- function() {
    updated <- record
    updated$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    # This closure records the caller's current run, not a private local copy.
    record <<- updated
    rna_atomic_rds(record, file.path(run_dir, "run_manifest.rds"))
    rna_atomic_rds(record, file.path(output, "latest_run_manifest.rds"))
    rows <- lapply(record$stage_records, function(x) data.frame(
      stage = x$stage, status = x$status, complete = isTRUE(x$complete),
      reused = isTRUE(x$reused), cache_key = x$key, cache_file = x$cache_file,
      stringsAsFactors = FALSE))
    if (length(rows)) rna_write_tsv(do.call(rbind, rows), file.path(run_dir, "stage_completion.tsv"))
    invisible(record)
  }
  tryCatch({
    for (stage in stages) {
      active_stage <- stage
      finished_stage_count <- sum(vapply(record$stage_records,
        function(x) isTRUE(x$complete), logical(1)))
      rna_progress(finished_stage_count, length(stages), stage, status = "starting",
        stage_index = match(stage, stages), detail = paste0("Starting stage ", match(stage, stages), " of ", length(stages), "."))
      stage_dir <- file.path(output, folders[[stage]])
      dirs <- rna_dirs(stage_dir)
      cache_file <- file.path(dirs$cache, "stage_result.rds")
      key <- stage_key(stage)
      cached <- NULL
      if (!isTRUE(force) && !identical(cfg$execution$resume, FALSE) && file.exists(cache_file)) {
        candidate <- tryCatch(readRDS(cache_file), error = function(e) NULL)
        if (is.list(candidate) && identical(candidate$schema_version, "whole-stage-cache-v1") &&
            identical(candidate$key, key) && isTRUE(candidate$complete)) cached <- candidate
      }
      if (!is.null(cached)) {
        result <- cached$result
        rna_log("Reusing completed ", stage, " statistics. To change figure styling, use scripts/replot.R.")
        rna_progress(finished_stage_count, length(stages), stage, status = "reused",
          stage_index = match(stage, stages), detail = "Loaded matching completed statistics from the stage cache.")
        reused <- TRUE
      } else {
        rna_log("Running ", stage, " from saved inputs/raw files.")
        seed <- cfg$execution$seeds[[if (stage == "qc") "import" else stage]]
        if (!is.null(seed)) set.seed(as.integer(seed))
        result <- withCallingHandlers(switch(stage,
          import = import_pigeon(cfg, samples, stage_dir),
          gain_loss = if(isTRUE(cfg$comparison$paired)) run_gain_loss(inputs, cfg, cohorts, stage_dir) else run_unpaired_gain_loss(inputs,cfg,cohorts,stage_dir),
          qc = { inputs$gain_loss <- results$gain_loss; if(isTRUE(cfg$comparison$paired))run_qc(inputs,cfg,cohorts,stage_dir)else run_unpaired_qc(inputs,cfg,cohorts,results$gain_loss,stage_dir) },
          remodeling = if(isTRUE(cfg$comparison$paired))run_remodeling(inputs,cfg,cohorts,stage_dir)else run_unpaired_remodeling(inputs,cfg,cohorts,stage_dir),
          dtu = run_dtu(inputs, cfg, cohorts, remodeling = results$remodeling, stage_dir = stage_dir)),
          warning = function(w) rna_log("Stage ", stage, ": ", conditionMessage(w), level = "WARNING"))
        # A stage may freeze reference files on first use; record its final identity.
        key <- stage_key(stage)
        complete <- rna_run_stage_complete(result, stage)
        rna_atomic_rds(list(schema_version = "whole-stage-cache-v1", stage = stage,
          key = key, complete = complete, profile = profile,
          created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), result = result), cache_file)
        reused <- FALSE
      }
      complete <- rna_run_stage_complete(result, stage)
      record$stage_records[[stage]] <- list(stage = stage,
        status = if (complete) "complete" else "audit_only", complete = complete,
        reused = reused, key = key, cache_file = cache_file)
      keys[[stage]] <- key
      if (stage == "import") inputs <- result else results[[stage]] <- result
      persist_record()
      finished_stage_count <- sum(vapply(record$stage_records,
        function(x) isTRUE(x$complete), logical(1)))
      rna_progress(finished_stage_count, length(stages), stage,
        status = if (complete) "complete" else "audit_only", stage_index = match(stage, stages),
        detail = if (complete) "Stage completion checks passed." else "Stage retained audit outputs; completion checks did not pass.")
      if (!complete && isTRUE(cfg$execution$strict_completion))
        rna_stop("Stage ", stage, " did not pass its completion checks. Audit files were retained in ", stage_dir)
    }
    stage_pass <- all(vapply(record$stage_records, function(x) isTRUE(x$complete), logical(1)))
    record$requested_stages_complete <- stage_pass
    record$completed <- stage_pass && identical(profile, "production") && all(enabled %in% stages)
    record$status <- if (record$completed) "complete" else if (!stage_pass) "audit_only" else if (!identical(profile, "production")) "diagnostic_complete" else "partial_complete"
    persist_record()
    if (record$completed) writeLines(c(paste("Completed UTC:", record$finished_utc),
      paste("Run:", run_id), paste("Manifest:", file.path(run_dir, "run_manifest.rds"))), marker)
    rna_progress(sum(vapply(record$stage_records, function(x) isTRUE(x$complete), logical(1))),
      length(stages), "workflow", status = record$status, stage_index = length(stages),
      detail = paste0("Requested stages finished; overall run status: ", record$status, "."))
    rna_log("IsoRemodel analysis finished: ", record$status, ". Run record: ", run_dir)
    invisible(list(run_record = record, inputs = inputs, results = results))
  }, error = function(e) {
    failed_record <- record
    failed_record$status <- "failed"; failed_record$completed <- FALSE
    failed_record$failed_stage <- active_stage; failed_record$error <- conditionMessage(e)
    if (active_stage %in% names(folders)) {
      failed_record$stage_records[[active_stage]] <- list(stage = active_stage,
        status = "failed", complete = FALSE, reused = FALSE,
        key = if (exists("key", inherits = TRUE)) key else NA_character_,
        cache_file = file.path(output, folders[[active_stage]], "cache", "stage_result.rds"))
    }
    failed_record$error_call <- paste(deparse(conditionCall(e)), collapse = " ")
    record <<- failed_record
    persist_record()
    rna_atomic_rds(list(run_record = record, condition = e,
      calls = lapply(sys.calls(), deparse)), file.path(run_dir, "failure_audit.rds"))
    writeLines(c(paste("Failed stage:", active_stage), conditionMessage(e),
      "Existing stage checkpoints, tables, and audit files were retained."), file.path(run_dir, "failure.txt"))
    writeLines(capture.output(utils::sessionInfo()), file.path(run_dir, "sessionInfo_after_failure.txt"))
    rna_progress(sum(vapply(record$stage_records, function(x) isTRUE(x$complete), logical(1))),
      length(stages), active_stage, status = "failed", stage_index = match(active_stage, stages),
      detail = conditionMessage(e))
    rna_log("FAILED at ", active_stage, ": ", conditionMessage(e), "; audits: ", run_dir, level = "ERROR")
    stop(e)
  })
}

#' Export saved plot objects to a new folder without rerunning statistics.
replot_rna_figures <- function(config_file = "config/analysis.yaml", source_dir = NULL,
                               output_dir = NULL) {
  cfg <- read_rna_config(config_file, validate_inputs = FALSE)
  base_output <- cfg$.output_dir %||% rna_run_absolute(cfg$output$directory, cfg$.config_dir)
  if (is.null(source_dir)) source_dir <- base_output
  source_dir <- rna_run_absolute(source_dir)
  if (!dir.exists(source_dir)) rna_stop("Saved figure directory does not exist: ", source_dir)
  if (is.null(output_dir)) output_dir <- paste0(source_dir, "_replot_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"))
  output_dir <- rna_run_absolute(output_dir)
  same_or_nested <- identical(tolower(source_dir), tolower(output_dir)) ||
    startsWith(paste0(tolower(output_dir), "/"), paste0(tolower(source_dir), "/"))
  if (same_or_nested) rna_stop("Replot output must be a separate directory outside the original source tree.")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    rna_stop("Replot output directory is not empty: ", output_dir)
  files <- sort(list.files(source_dir, pattern = "\\.plot\\.rds$", full.names = TRUE, recursive = TRUE))
  if (!length(files)) rna_stop("No *.plot.rds files were found under ", source_dir)
  if (any(grepl("(^|/)M06_dtu(/|$)", files)))
    rna_stop("Use replot_isoremodel() for an output tree containing DTU. ",
      "It verifies the DTU RNG provenance and uses the confirmed saved results.")
  cfg$reporting$enabled <- TRUE
  cfg$reporting$save_plot_objects <- FALSE
  # The saved tree already contains individual facet exports; do not duplicate them.
  cfg$reporting$separate_panels_for_composites <- FALSE
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    relative <- substring(files[[i]], nchar(source_dir) + 2L)
    destination <- file.path(output_dir, dirname(relative))
    stem <- sub("\\.plot\\.rds$", "", basename(relative))
    rows[[i]] <- tryCatch({
      saved <- readRDS(files[[i]])
      if (!is.list(saved) || is.null(saved$plot)) rna_stop("Invalid saved plot object.")
      paths <- rna_save_plot(saved$plot, stem, destination, cfg,
        width = saved$width %||% 8, height = saved$height %||% 6)
      data.frame(source = files[[i]], output = paste(paths, collapse = ";"), status = "complete", error = "")
    }, error = function(e) data.frame(source = files[[i]], output = "", status = "failed", error = conditionMessage(e)))
  }
  audit <- do.call(rbind, rows)
  rna_write_tsv(audit, file.path(output_dir, "replot_audit.tsv"))
  yaml::write_yaml(cfg$reporting, file.path(output_dir, "reporting_settings.yaml"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  rna_log("Replot exported ", sum(audit$status == "complete"), " saved figures; ", sum(audit$status == "failed"), " failures. Output: ", output_dir)
  if (any(audit$status == "failed")) rna_stop("Replotting incomplete; review replot_audit.tsv. Successful exports were retained.")
  invisible(audit)
}
