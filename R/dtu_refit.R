# Targeted publication-reproducibility refit. No upstream statistical stage runs.
# Full fitted objects stay in the new run directory; the review ZIP is compact.
dtu_refit_compare_table <- function(new, reference, keys, numeric_fields,
                                    categorical_fields = character(),
                                    atol = 1e-12, rtol = 1e-8) {
  fields <- unique(c(keys, numeric_fields, categorical_fields))
  validate <- function(x, label) {
    if (!is.data.frame(x) || !nrow(x) || !all(fields %in% names(x)))
      stop(label, " must be a nonempty table containing: ", paste(fields, collapse = ", "), call. = FALSE)
    if (anyNA(x[keys]) || any(vapply(x[keys], function(z) any(!nzchar(as.character(z))), logical(1))) || anyDuplicated(x[keys]))
      stop(label, " has missing, blank, or duplicate comparison keys.", call. = FALSE)
    for (field in numeric_fields) if (!is.numeric(x[[field]]))
      stop(label, " comparison field is not numeric: ", field, call. = FALSE)
  }
  if (!length(keys) || anyDuplicated(keys) || any(keys %in% c(numeric_fields, categorical_fields)))
    stop("Comparison keys must be distinct from value fields.", call. = FALSE)
  if (any(!is.finite(c(atol, rtol))) || atol < 0 || rtol < 0)
    stop("Comparison tolerances must be nonnegative and finite.", call. = FALSE)
  validate(new, "New result"); validate(reference, "Reference result")
  a <- as.data.frame(new[fields]); b <- as.data.frame(reference[fields])
  a$.present_new <- TRUE; b$.present_ref <- TRUE
  z <- merge(a, b, by = keys, all = TRUE, sort = TRUE, suffixes = c("_new", "_ref"))
  present <- !is.na(z$.present_new) & !is.na(z$.present_ref)
  summaries <- list(data.frame(metric = "key_family", compared = nrow(z),
    mismatched = sum(!present), pass = all(present), max_absolute_error = NA_real_,
    max_relative_error = NA_real_, max_positive_log10_difference = NA_real_, stringsAsFactors = FALSE))
  matches <- present
  for (field in c(numeric_fields, categorical_fields)) {
    a <- z[[paste0(field, "_new")]]; b <- z[[paste0(field, "_ref")]]
    if (field %in% numeric_fields) {
      finite <- is.finite(a) & is.finite(b)
      ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & is.infinite(a) & is.infinite(b) & a == b)
      ok[finite] <- abs(a[finite] - b[finite]) <= atol + rtol * abs(b[finite])
    } else ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
    ok[is.na(ok)] <- FALSE
    ok <- present & ok
    z[[paste0(field, "_match")]] <- ok
    finite <- if (field %in% numeric_fields) present & is.finite(a) & is.finite(b) else rep(FALSE, length(a))
    positive <- finite & a > 0 & b > 0
    relative <- finite & b != 0
    summaries[[length(summaries) + 1L]] <- data.frame(metric = field, compared = sum(present),
      mismatched = sum(!ok), pass = all(ok),
      max_absolute_error = if (any(finite)) max(abs(a[finite] - b[finite])) else NA_real_,
      max_relative_error = if (any(relative)) max(abs(a[relative] - b[relative]) / abs(b[relative])) else NA_real_,
      max_positive_log10_difference = if (any(positive)) max(abs(log10(a[positive]) - log10(b[positive]))) else NA_real_,
      stringsAsFactors = FALSE)
    matches <- matches & ok
  }
  z$row_matches <- matches
  list(summary = do.call(rbind, summaries), detail = z, match = all(matches))
}

dtu_refit_validate_inputs <- function(a, b, label = "") {
  rows <- list(); add <- function(check, pass, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(comparison = label, check = check,
      pass = isTRUE(pass), detail = detail, stringsAsFactors = FALSE)
  }
  fields <- c("counts", "samples", "design", "gene_order")
  valid <- is.list(a) && is.list(b) && all(fields %in% names(a)) && all(fields %in% names(b))
  add("required_input_fields", valid)
  if (!valid) return(do.call(rbind, rows))
  x <- a$counts; y <- b$counts
  valid <- is.data.frame(x) && is.data.frame(y) && all(c("gene_id", "feature_id") %in% names(x)) && all(c("gene_id", "feature_id") %in% names(y))
  add("count_table_schema", valid)
  if (!valid) return(do.call(rbind, rows))
  add("nonempty_count_tables", nrow(x) > 0L && nrow(y) > 0L)
  add("ordered_count_columns", identical(names(x), names(y)))
  add("ordered_gene_ids", identical(as.character(x$gene_id), as.character(y$gene_id)))
  add("ordered_feature_ids", identical(as.character(x$feature_id), as.character(y$feature_id)))
  add("ordered_gene_blocks", identical(as.character(a$gene_order), as.character(b$gene_order)))
  add("valid_gene_feature_keys", !anyNA(x[c("gene_id", "feature_id")]) && !anyNA(y[c("gene_id", "feature_id")]) &&
    !anyDuplicated(x$feature_id) && !anyDuplicated(y$feature_id) && !anyDuplicated(a$gene_order) && !anyDuplicated(b$gene_order) &&
    identical(as.character(a$gene_order), unique(as.character(x$gene_id))) && identical(as.character(b$gene_order), unique(as.character(y$gene_id))))
  values <- setdiff(names(x), c("gene_id", "feature_id"))
  numeric <- length(values) > 0L && all(values %in% names(y)) && all(vapply(x[values], is.numeric, logical(1))) && all(vapply(y[values], is.numeric, logical(1)))
  add("numeric_count_columns", numeric)
  if (numeric) {
    xx <- as.matrix(x[values]); yy <- as.matrix(y[values])
    add("finite_nonnegative_integer_counts", all(is.finite(xx)) && all(is.finite(yy)) && all(xx >= 0) && all(yy >= 0) && all(xx == round(xx)) && all(yy == round(yy)))
    add("exact_ordered_counts", identical(dim(xx), dim(yy)) && identical(as.numeric(xx), as.numeric(yy)))
  }
  valid_samples <- is.data.frame(a$samples) && is.data.frame(b$samples) && all(c("sample_id", "patient", "timepoint") %in% names(a$samples)) && all(c("sample_id", "patient", "timepoint") %in% names(b$samples))
  add("sample_schema", valid_samples)
  if (valid_samples) {
    add("sample_column_alignment", identical(as.character(a$samples$sample_id), values) && identical(as.character(b$samples$sample_id), setdiff(names(y), c("gene_id", "feature_id"))))
    add("unique_nonmissing_samples", nrow(a$samples) > 0L && !anyNA(a$samples$sample_id) && !anyDuplicated(a$samples$sample_id) && !anyNA(b$samples$sample_id) && !anyDuplicated(b$samples$sample_id))
    for (field in c("sample_id", "patient", "timepoint")) {
      add(paste0("ordered_", field), identical(as.character(a$samples[[field]]), as.character(b$samples[[field]])))
      add(paste0(field, "_levels"), identical(levels(a$samples[[field]]), levels(b$samples[[field]])) && identical(is.ordered(a$samples[[field]]), is.ordered(b$samples[[field]])))
    }
  }
  valid_design <- is.matrix(a$design) && is.matrix(b$design) && is.numeric(a$design) && is.numeric(b$design)
  add("design_schema", valid_design)
  if (valid_design) {
    add("design_values", identical(dim(a$design), dim(b$design)) && identical(as.numeric(a$design), as.numeric(b$design)))
    add("design_attributes", identical(attributes(a$design), attributes(b$design)))
    add("full_rank_paired_design", all(is.finite(a$design)) && all(is.finite(b$design)) && qr(a$design)$rank == ncol(a$design) && qr(b$design)$rank == ncol(b$design) && nrow(a$design) > ncol(a$design) && nrow(b$design) > ncol(b$design))
  }
  do.call(rbind, rows)
}

dtu_refit_load_stage <- function(path, stage) {
  if (!file.exists(path)) stop("Required completed stage cache is missing: ", path, call. = FALSE)
  x <- readRDS(path)
  if (!is.list(x) || !identical(x$schema_version, "whole-stage-cache-v1") ||
      !identical(x$stage, stage) || !isTRUE(x$complete) || !identical(x$profile, "production") || !is.list(x$result))
    stop("Stage cache must be a completed production ", stage, " cache: ", path, call. = FALSE)
  x$result
}

dtu_refit_read_fit <- function(path) {
  if (!file.exists(path)) stop("Required fitted-model cache is missing: ", path, call. = FALSE)
  x <- readRDS(path)
  if (!is.list(x) || !is.list(x$result) || is.null(x$result$fitted_dm) ||
      !is.data.frame(x$result$gene) || !is.data.frame(x$result$feature))
    stop("Expected a fitted DRIMSeq cache containing gene and feature results: ", path, call. = FALSE)
  x$result
}

dtu_refit_fit_input <- function(fit) {
  dm <- fit$fitted_dm
  if (!isS4(dm) || !"counts" %in% methods::slotNames(dm)) stop("Fitted DRIMSeq object has no audited counts slot.")
  list(counts = as.data.frame(DRIMSeq::counts(dm), stringsAsFactors = FALSE),
       samples = as.data.frame(DRIMSeq::samples(dm), stringsAsFactors = FALSE),
       design = DRIMSeq::design(dm, type = "precision"),
       gene_order = as.character(names(methods::slot(dm, "counts"))))
}

dtu_refit_stageR <- function(fit, prepared, cohort, alpha, cfg) {
  e <- dtu_engine_environment(parent = environment())
  sys.source(dtu_engine_path(cfg, "dtu_core.R"), envir = e)
  e$dtu_fits <- stats::setNames(list(fit), cohort)
  e$dtu_prepared <- stats::setNames(list(prepared), cohort)
  e$cfg <- list()
  e$run_stageR_alpha(cohort, alpha)
}

dtu_refit_collect <- function(run_dir, review_dir) {
  # Only small summary tables are copied. Models and precision-grid arrays stay local.
  names <- c("dtu_rng_settings.tsv", "drimseq_fit_and_discovery_summary.tsv",
    "stageR_completion_summary.tsv", "completion_gates.tsv", "execution_backend.tsv")
  for (name in names) {
    path <- file.path(run_dir, "M06_dtu", "qc", name)
    if (file.exists(path)) file.copy(path, file.path(review_dir, paste0("refit_", name)), overwrite = TRUE)
  }
  paths <- list.files(review_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  paths <- paths[!dir.exists(paths) & basename(paths) != "FILE_SHA256SUMS.txt"]
  relative <- substring(paths, nchar(normalizePath(review_dir, winslash = "/", mustWork = TRUE)) + 2L)
  writeLines(paste(vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), relative, sep = "  "), file.path(review_dir, "FILE_SHA256SUMS.txt"))
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(review_dir)
  zip_path <- file.path(run_dir, "DTU_refit_review.zip")
  members <- list.files(".", recursive = TRUE, all.files = FALSE, full.names = FALSE)
  zip::zipr(zip_path, files = members, root = review_dir, mode = "mirror", include_directories = FALSE)
  listed <- zip::zip_list(zip_path)$filename
  if (!setequal(gsub("\\\\", "/", members), gsub("\\\\", "/", listed))) stop("Review ZIP member paths were not preserved.")
  zip_path
}

#' Refit only DTU with the publication's explicit random-number settings.
#' @export
refit_dtu_rng <- function(
    config_file = "F:/IsoRemodel/IsoRemodel_v0.1.1/config/analysis.yaml",
    previous_output = "F:/IsoRemodel/outputs",
    original_dtu = "F:/PacBio revision analysis/outputs/M06_DTU_robustness",
    output_root = "F:/IsoRemodel/DTU_refit", strict = TRUE) {
  rna_require(c("yaml", "digest", "zip", "DRIMSeq", "stageR", "BiocParallel", "MASS", "readr", "dplyr", "tidyr", "tibble", "purrr", "stringr", "scales", "forcats", "ggplot2"))
  if (!identical(as.character(utils::packageVersion("DRIMSeq")), "1.38.0") ||
      !identical(as.character(utils::packageVersion("stageR")), "1.32.0"))
    stop("This confirmation refit requires DRIMSeq 1.38.0 and stageR 1.32.0. Do not upgrade the analysis environment.", call. = FALSE)
  cfg <- read_rna_config(config_file, validate_inputs = FALSE)
  if (!identical(cfg$execution$profile, "production")) stop("The confirmation refit requires the production profile.", call. = FALSE)
  if (!identical(as.integer(cfg$execution$seeds$dtu), 260605L)) stop("The publication comparison requires execution.seeds.dtu: 260605.", call. = FALSE)
  required_rng <- list(kind = "L'Ecuyer-CMRG", normal_kind = "Inversion", sample_kind = "Rejection")
  if (!identical(dtu_rng_spec(cfg$dtu$rng), required_rng)) stop("The confirmation refit requires L'Ecuyer-CMRG/Inversion/Rejection.", call. = FALSE)
  if (!all(c(.05, .10) %in% as.numeric(cfg$multiplicity$fdr_thresholds)) || !isTRUE(cfg$dtu$run_stageR))
    stop("The confirmation refit requires both FDR 0.05/0.10 and stageR enabled.", call. = FALSE)
  normalize <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  previous_output <- normalize(previous_output); original_dtu <- normalize(original_dtu); output_root <- normalize(output_root)
  inside <- function(a, b) identical(tolower(a), tolower(b)) || startsWith(tolower(a), paste0(tolower(b), "/"))
  if (inside(output_root, previous_output) || inside(output_root, original_dtu))
    stop("output_root must be separate from the previous IsoRemodel and original DTU outputs.", call. = FALSE)
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "_", Sys.getpid())
  run_dir <- file.path(output_root, run_id)
  if (dir.exists(run_dir)) run_dir <- paste0(run_dir, "_", format(Sys.time(), "%OS6"))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  review <- file.path(run_dir, "review"); dir.create(review, recursive = TRUE)
  option_names <- c("IsoRemodel.log_file", "IsoRemodel.progress_file", "IsoRemodel.progress_latest_file", "IsoRemodel.run_id", "IsoRemodel.run_start")
  old_options <- stats::setNames(lapply(option_names, getOption), option_names)
  old_rng <- dtu_rng_snapshot()
  on.exit({options(old_options); dtu_rng_restore(old_rng)}, add = TRUE)
  options(IsoRemodel.log_file = file.path(review, "run.log"), IsoRemodel.progress_file = file.path(review, "progress.tsv"),
    IsoRemodel.progress_latest_file = file.path(review, "progress_latest.tsv"), IsoRemodel.run_id = run_id, IsoRemodel.run_start = Sys.time())
  rna_banner(log = TRUE)
  status <- list(status = "running", completed = FALSE, original_results_reproduced = FALSE,
    package_version = rna_package_version(), run_id = run_id, run_directory = run_dir,
    started_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), error = NA_character_)
  comparisons <- list(); calls <- list(); inputs_audit <- list(); precision <- list(); result <- NULL
  protected <- character(); before <- NULL; steps_passed <- 0L
  flush_comparisons <- function() {
    if (length(comparisons)) rna_write_tsv(do.call(rbind, comparisons), file.path(review, "comparison_summary.tsv"))
    if (length(precision)) rna_write_tsv(do.call(rbind, precision), file.path(review, "common_precision_comparison.tsv"))
    if (length(calls)) rna_write_tsv(do.call(rbind, calls), file.path(review, "discovery_comparison.tsv"))
  }
  tryCatch({
    rna_progress(0, 4, "preflight", status = "starting", detail = "Load completed import/remodeling caches; no upstream analyses will run.")
    paths <- c(import = file.path(previous_output, "M01_import", "cache", "stage_result.rds"),
      remodeling = file.path(previous_output, "M04_remodeling", "cache", "stage_result.rds"))
    inputs <- dtu_refit_load_stage(paths[["import"]], "import")
    remodeling <- dtu_refit_load_stage(paths[["remodeling"]], "remodeling")
    if (!all(c("counts_matrix", "sample_manifest", "isoform_annotation", "fingerprint") %in% names(inputs)) ||
        !is.matrix(inputs$counts_matrix) || !all(c("gene_table", "patient_inference") %in% names(remodeling)))
      stop("Completed upstream caches lack the required DTU handoff fields.")
    inputs <- inputs[c("counts_matrix", "sample_manifest", "isoform_annotation", "fingerprint")]
    remodeling <- remodeling[c("gene_table", "patient_inference", "fingerprint", "complete")]
    cohorts <- rna_resolve_cohorts(inputs$sample_manifest, cfg)
    if (!identical(names(cohorts), c("all_6_pairs", "excluding_ptn2_5_pairs")) ||
        length(cohorts[[1L]]) != 6L || length(cohorts[[2L]]) != 5L ||
        !identical(setdiff(cohorts[[1L]], cohorts[[2L]]), "ptn2"))
      stop("This confirmation run requires the original all-six and patient-2-excluded cohorts.")
    if (!all(names(cohorts) %in% unique(as.character(remodeling$gene_table$cohort)))) stop("Saved remodeling does not cover the requested cohorts.")
    fits <- list(original = list(), previous_IsoRemodel = list())
    for (cohort in names(cohorts)) for (source in names(fits)) {
      p <- file.path(if (source == "original") original_dtu else file.path(previous_output, "M06_dtu"), "cache", paste0("drimseq_fit_", cohort, ".rds"))
      fits[[source]][[cohort]] <- dtu_refit_read_fit(p); paths[paste(source, cohort, sep = "__")] <- p
    }
    protected <- unname(paths)
    before <- rna_run_file_manifest(protected)
    rna_write_tsv(before, file.path(review, "source_cache_manifest.tsv"))
    for (cohort in names(cohorts)) {
      check <- dtu_refit_validate_inputs(dtu_refit_fit_input(fits$original[[cohort]]), dtu_refit_fit_input(fits$previous_IsoRemodel[[cohort]]), paste(cohort, "original_vs_previous"))
      inputs_audit[[length(inputs_audit) + 1L]] <- check
    }
    rna_write_tsv(do.call(rbind, inputs_audit), file.path(review, "input_comparisons.tsv"))
    if (!all(vapply(inputs_audit, function(z) all(z$pass), logical(1)))) stop("Saved original and previous IsoRemodel inputs differ. No fit was started.")
    # Freeze the already-used pathway resources in this run, never write shared references.
    if (isTRUE(cfg$enrichment$enabled)) {
      ref_names <- c("dtu_term2gene_frozen.tsv.gz", "dtu_gene_symbol_mapping_frozen.tsv.gz", "dtu_reference_sources_frozen.rds")
      from <- file.path(cfg$enrichment$reference_cache_directory, ref_names)
      if (any(!file.exists(from))) stop("Frozen DTU pathway reference files are missing; retain the original reference cache before refitting: ", paste(from[!file.exists(from)], collapse = ", "))
      ref_dir <- file.path(run_dir, "references", "pathways"); dir.create(ref_dir, recursive = TRUE)
      if (!all(file.copy(from, file.path(ref_dir, ref_names), overwrite = FALSE))) stop("Could not copy the frozen reference resources.")
      protected <- unique(c(protected, from))
      before <- rna_run_file_manifest(protected)
      rna_write_tsv(before, file.path(review, "source_cache_manifest.tsv"))
      cfg$enrichment$reference_cache_directory <- ref_dir
      cfg$references$pathway_term2gene <- file.path(ref_dir, ref_names[[1L]])
      cfg$references$gene_id_map <- file.path(ref_dir, ref_names[[2L]])
    }
    complexity_cache <- file.path(previous_output, "M06_dtu", "cache", "reference_transcript_complexity.rds")
    if (file.exists(complexity_cache)) {
      dest <- file.path(run_dir, "M06_dtu", "cache"); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(complexity_cache, file.path(dest, basename(complexity_cache)), overwrite = FALSE)) stop("Could not copy the saved transcript-complexity cache.")
    }
    if (!file.copy(config_file, file.path(review, "source_analysis.yaml"), overwrite = FALSE)) stop("Could not record the source YAML.")
    cfg$reporting$enabled <- FALSE
    cfg$execution$force_refit <- TRUE; cfg$dtu$force_refit <- TRUE
    cfg$execution$strict_completion <- TRUE; cfg$dtu$rng <- required_rng
    cfg$output$directory <- run_dir
    yaml::write_yaml(cfg[!startsWith(names(cfg), ".")], file.path(review, "resolved_refit_analysis.yaml"))
    rna_write_tsv(rna_run_package_versions(cfg$.code_root), file.path(review, "package_versions.tsv"))
    rna_write_tsv(rna_run_code_manifest(cfg$.code_root), file.path(review, "source_code_sha256.tsv"))
    writeLines(capture.output(utils::sessionInfo()), file.path(review, "sessionInfo.txt"))
    rna_write_tsv(dtu_rng_contract(required_rng, cfg$execution$seeds$dtu, names(cohorts)), file.path(review, "requested_rng_contract.tsv"))
    cfg$.dtu_preflight <- function(prepared) {
      if (isTRUE(cfg$enrichment$enabled)) {
        mapping <- readr::read_tsv(cfg$references$gene_id_map, show_col_types = FALSE, progress = FALSE)
        required_genes <- unique(unlist(lapply(prepared, function(z) as.character(z$testing_family$gene_id))))
        if (!"gene_id" %in% names(mapping) || anyDuplicated(mapping$gene_id) || length(setdiff(required_genes, mapping$gene_id)))
          stop("Frozen gene-symbol mapping does not cover the complete prepared families; no new annotation mapping or fit was attempted.")
      }
      e <- dtu_engine_environment(parent = environment())
      sys.source(dtu_engine_path(cfg, "dtu_core.R"), envir = e)
      for (cohort in names(cohorts)) {
        p <- prepared[[cohort]]
        if (isTRUE(p$empty)) stop("The prepared DTU family is unexpectedly empty: ", cohort)
        ordered <- e$m06_reference_order(as.data.frame(p$counts), as.data.frame(p$samples))
        candidate <- list(counts = ordered$counts, samples = as.data.frame(p$samples), design = p$design, gene_order = unique(as.character(ordered$counts$gene_id)))
        check <- dtu_refit_validate_inputs(dtu_refit_fit_input(fits$original[[cohort]]), candidate, paste(cohort, "original_vs_new_prepared"))
        inputs_audit[[length(inputs_audit) + 1L]] <<- check
      }
      rna_write_tsv(do.call(rbind, inputs_audit), file.path(review, "input_comparisons.tsv"))
      if (!all(vapply(inputs_audit, function(z) all(z$pass), logical(1)))) stop("Newly prepared DTU inputs differ from the original; no model fit was started.")
      rna_log("All ordered cached and newly prepared input checks passed; starting fresh DTU fits.")
    }
    steps_passed <- 1L
    rna_progress(1, 4, "DTU_refit", status = "starting", detail = "Fresh DRIMSeq fits and stageR at FDR 0.05 and 0.10; plotting disabled.")
    result <- run_dtu(inputs, cfg, cohorts, remodeling = remodeling, stage_dir = file.path(run_dir, "M06_dtu"))
    steps_passed <- 2L
    rna_progress(2, 4, "comparison", status = "starting", detail = "Compare full families, common precision, gene/isoform tests, and stageR with both saved runs.")
    for (cohort in names(cohorts)) {
      fresh <- result$cohort_results[[cohort]]; new_fit <- fresh$fit
      for (source in names(fits)) {
        ref <- fits[[source]][[cohort]]; prefix <- paste(cohort, source, sep = "__")
        a <- as.numeric(DRIMSeq::common_precision(new_fit$fitted_dm)); b <- as.numeric(DRIMSeq::common_precision(ref$fitted_dm))
        precision[[prefix]] <- data.frame(cohort = cohort, reference = source, new = a, saved = b, delta = a - b,
          match = is.finite(a) && is.finite(b) && abs(a - b) <= 1e-8 + 1e-6 * abs(b))
        numeric_gene <- c("raw_dtu_p_value", "dtu_p_value", "dtu_bh")
        category_gene <- c("fit_status", "numerical_failure_stage", "CR_boundary_associated_P1")
        check <- dtu_refit_compare_table(as.data.frame(new_fit$gene), as.data.frame(ref$gene), "gene_id", numeric_gene, category_gene, atol = 0, rtol = 1e-8)
        check$summary$cohort <- cohort; check$summary$reference <- source; check$summary$table <- "gene_tests"
        comparisons[[paste0(prefix, "_gene")]] <- check$summary
        rna_write_tsv(check$detail, file.path(review, paste0(prefix, "_gene_comparison.tsv.gz")))
        numerical_fields <- c("precision", "likelihood_ratio", "df", "full_model_log_likelihood", "null_model_log_likelihood")
        if (!all(numerical_fields %in% names(new_fit$gene)) || !all(numerical_fields %in% names(ref$gene)))
          stop("Saved gene tests lack the fitted precision/log-likelihood fields: ", prefix)
        check <- dtu_refit_compare_table(as.data.frame(new_fit$gene), as.data.frame(ref$gene), "gene_id", numerical_fields, atol = 1e-10, rtol = 1e-8)
        check$summary$cohort <- cohort; check$summary$reference <- source; check$summary$table <- "gene_model_estimates"
        comparisons[[paste0(prefix, "_gene_model")]] <- check$summary
        rna_write_tsv(check$detail, file.path(review, paste0(prefix, "_gene_model_comparison.tsv.gz")))
        flush_comparisons()
        check <- dtu_refit_compare_table(as.data.frame(new_fit$feature), as.data.frame(ref$feature), c("gene_id", "feature_id"), "isoform_dtu_p_value", atol = 0, rtol = 1e-8)
        check$summary$cohort <- cohort; check$summary$reference <- source; check$summary$table <- "isoform_tests"
        comparisons[[paste0(prefix, "_feature")]] <- check$summary
        rna_write_tsv(check$detail, file.path(review, paste0(prefix, "_isoform_comparison.tsv.gz")))
        flush_comparisons()
        for (alpha in c(.05, .10)) {
          tag <- dtu_alpha_tag(alpha)
          new_stage <- fresh$stageR[[tag]]
          ref_stage <- dtu_refit_stageR(ref, fresh$prepared, cohort, alpha, cfg)
          if (!isTRUE(new_stage$available) || !isTRUE(ref_stage$available)) stop("stageR was unavailable in the comparison: ", prefix, "/", tag)
          for (level in c("gene", "feature")) {
            field <- paste0(if (level == "gene") "stageR_gene_adjusted_" else "stageR_isoform_adjusted_", tag)
            check <- dtu_refit_compare_table(as.data.frame(new_stage[[level]]), as.data.frame(ref_stage[[level]]),
              if (level == "gene") "gene_id" else c("gene_id", "feature_id"), field, atol = 0, rtol = 1e-8)
            check$summary$cohort <- cohort; check$summary$reference <- source; check$summary$table <- paste0("stageR_", level, "_", tag)
            comparisons[[paste(prefix, level, tag, sep = "_")]] <- check$summary
            rna_write_tsv(check$detail, file.path(review, paste0(prefix, "_stageR_", level, "_", tag, ".tsv.gz")))
            flush_comparisons()
          }
          gene_new <- new_fit$gene$gene_id[is.finite(new_fit$gene$dtu_bh) & new_fit$gene$dtu_bh < alpha]
          gene_ref <- ref$gene$gene_id[is.finite(ref$gene$dtu_bh) & ref$gene$dtu_bh < alpha]
          iso_field <- paste0("stageR_isoform_adjusted_", tag)
          tx_new <- new_stage$feature$feature_id[is.finite(new_stage$feature[[iso_field]]) & new_stage$feature[[iso_field]] < alpha]
          tx_ref <- ref_stage$feature$feature_id[is.finite(ref_stage$feature[[iso_field]]) & ref_stage$feature[[iso_field]] < alpha]
          calls[[paste(prefix, tag, sep = "_")]] <- data.frame(cohort = cohort, reference = source, FDR_threshold = alpha,
            new_DTU_genes = length(gene_new), saved_DTU_genes = length(gene_ref), shared_DTU_genes = length(intersect(gene_new, gene_ref)),
            new_only_DTU_genes = length(setdiff(gene_new, gene_ref)), saved_only_DTU_genes = length(setdiff(gene_ref, gene_new)),
            new_unavailable_gene_tests = sum(!is.finite(new_fit$gene$raw_dtu_p_value)), saved_unavailable_gene_tests = sum(!is.finite(ref$gene$raw_dtu_p_value)),
            new_confirmed_isoforms = length(tx_new), saved_confirmed_isoforms = length(tx_ref), shared_confirmed_isoforms = length(intersect(tx_new, tx_ref)),
            new_only_confirmed_isoforms = length(setdiff(tx_new, tx_ref)), saved_only_confirmed_isoforms = length(setdiff(tx_ref, tx_new)))
          changed <- data.frame(gene_id = union(setdiff(gene_new, gene_ref), setdiff(gene_ref, gene_new)))
          changed$new_significant <- changed$gene_id %in% gene_new; changed$saved_significant <- changed$gene_id %in% gene_ref
          rna_write_tsv(changed, file.path(review, paste0(prefix, "_changed_DTU_genes_", tag, ".tsv")))
          flush_comparisons()
        }
      }
    }
    summary <- do.call(rbind, comparisons); common <- do.call(rbind, precision); discovery <- do.call(rbind, calls)
    rna_write_tsv(summary, file.path(review, "comparison_summary.tsv"))
    rna_write_tsv(common, file.path(review, "common_precision_comparison.tsv"))
    rna_write_tsv(discovery, file.path(review, "discovery_comparison.tsv"))
    status$original_results_reproduced <- all(summary$pass[summary$reference == "original"]) && all(common$match[common$reference == "original"]) &&
      all(discovery$new_only_DTU_genes[discovery$reference == "original"] == 0L) && all(discovery$saved_only_DTU_genes[discovery$reference == "original"] == 0L) &&
      all(discovery$new_only_confirmed_isoforms[discovery$reference == "original"] == 0L) && all(discovery$saved_only_confirmed_isoforms[discovery$reference == "original"] == 0L)
    steps_passed <- 3L
    rna_progress(3, 4, "comparison", status = "complete", detail = paste0("Original final DTU/stageR results reproduced: ", status$original_results_reproduced))
    status$completed <- isTRUE(result$complete)
    status$status <- if (isTRUE(status$completed) && isTRUE(status$original_results_reproduced)) "original_results_reproduced" else "comparison_requires_review"
  }, error = function(e) {status$status <<- "failed"; status$error <<- conditionMessage(e); rna_log(conditionMessage(e), level = "ERROR")},
  interrupt = function(e) {status$status <<- "interrupted"; status$error <<- "Refit interrupted; completed new files were retained."; rna_log(status$error, level = "WARNING")})
  # Protect all existing numerical results. No latest-run markers are rewritten.
  if (!is.null(before)) {
    after <- rna_run_file_manifest(protected)
    unchanged <- identical(before, after)
    status$source_cache_files_unchanged <- unchanged
    rna_write_tsv(after, file.path(review, "source_cache_manifest_after.tsv"))
    if (!unchanged) {status$status <- "failed"; status$error <- paste(na.omit(c(status$error, "An existing source cache changed during the refit.")), collapse = " ")}
  }
  status$finished_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  rna_write_tsv(as.data.frame(status, stringsAsFactors = FALSE), file.path(review, "refit_status.tsv"))
  rna_atomic_rds(status, file.path(run_dir, "refit_status.rds"))
  writeLines(c(paste("Status:", status$status), paste("Original final results reproduced:", status$original_results_reproduced),
    paste("Error:", status$error), "No upstream statistical analyses or plotting updates were run.",
    "The original and previous IsoRemodel output directories were read only.",
    "Full fitted objects are retained outside this review ZIP in the new M06_dtu directory.",
    "Unavailable raw tests remain unavailable in comparisons; conservative P=1 values are compared separately.",
    "P/FDR comparisons: absolute tolerance 0, relative tolerance 1e-8; zero values must agree exactly.",
    "Gene precision, likelihood-ratio, degrees-of-freedom and log-likelihood comparisons: absolute tolerance 1e-10, relative tolerance 1e-8.",
    "Common precision comparisons: absolute tolerance 1e-8, relative tolerance 1e-6.",
    "Discovery and confirmed-isoform membership at FDR <0.05/<0.10 must agree exactly."), file.path(review, "READ_FIRST.txt"))
  rna_progress(if (isTRUE(status$completed)) 4 else steps_passed, 4, "confirmation", status = status$status,
    detail = paste0("Original final DTU/stageR results reproduced: ", status$original_results_reproduced))
  zip_path <- dtu_refit_collect(run_dir, review)
  rna_log("Review bundle: ", zip_path)
  # Rebuild once so the final log line is represented by the bundle's hashes.
  zip_path <- dtu_refit_collect(run_dir, review)
  if (isTRUE(strict) && (!isTRUE(status$completed) || !isTRUE(status$original_results_reproduced) || identical(status$status, "failed")))
    stop("DTU confirmation did not pass. Completed new outputs were retained. Upload: ", zip_path, call. = FALSE)
  invisible(list(status = status, run_directory = run_dir, review_zip = zip_path))
}
