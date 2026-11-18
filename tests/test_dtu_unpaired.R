# Independent DTU design, null, sample selection and effects are testable
# without claiming execution of the pinned native DRIMSeq numerical engine.
(function() {
  engine <- new.env(parent = globalenv())
  sys.source(file.path(.repo_root, "inst", "engines", "dtu_core.R"), engine)
  samples <- data.frame(sample_id = paste0("s", 1:7),
    patient = factor(paste0("independent", 1:7)),
    timepoint = factor(c(rep("primary", 3), rep("recurrent", 4)),
      levels = c("primary", "recurrent")))
  design <- engine$m06_design_matrix(samples, paired = FALSE)
  tested <- which(colnames(design) == "timepointrecurrent")
  null <- design[, -tested, drop = FALSE]
  stopifnot(nrow(design) == 7L, ncol(design) == 2L, qr(design)$rank == 2L,
    ncol(null) == 1L, qr(null)$rank == 1L, all(null == 1),
    !any(grepl("patient", colnames(design))))
  renamed <- samples; renamed$patient <- factor(rep("unused", 7))
  stopifnot(identical(design, engine$m06_design_matrix(renamed, paired = FALSE)))
  paired <- samples[1:6, ]; paired$patient <- factor(rep(1:3, 2))
  stopifnot(identical(engine$m06_design_matrix(paired),
    stats::model.matrix(~ patient + timepoint, paired)))
  counts <- data.frame(gene_id = c("g", "g"), feature_id = c("g|a", "g|b"),
    s1 = c(8, 0), s2 = c(5, 0), s3 = c(0, 0),
    s4 = c(5, 10), s5 = c(5, 10), s6 = c(5, 10), s7 = c(5, 10))
  ordered <- engine$m06_reference_order(counts, samples, paired = FALSE)
  again <- engine$m06_reference_order(counts, renamed, paired = FALSE)
  stopifnot(identical(ordered, again), ordered$audit$selected_reference == "g|a",
    ordered$audit$minimum_condition_reference_count == 13,
    ordered$audit$design_rank_with_gene_counts == 2L,
    ordered$audit$reference_supported_in_both_conditions,
    is.na(ordered$audit$reference_supported_in_all_pairs),
    is.na(ordered$audit$patient_pairs_without_gene_counts),
    is.na(ordered$audit$isoform_patient_pairs_with_zero_counts),
    identical(colSums(counts[, samples$sample_id]),
      colSums(ordered$counts[, samples$sample_id])))

  dependencies <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "scales", "forcats", "ggplot2")
  if (!all(vapply(dependencies, requireNamespace, logical(1), quietly = TRUE))) {
    message("SKIP independent DTU effect/selection tests: data packages unavailable")
    return(invisible(NULL))
  }
  engine <- dtu_engine_environment()
  sys.source(file.path(.repo_root, "inst", "engines", "dtu_core.R"), engine)
  engine$empty_dtu_prepared <- dtu_empty_prepared
  engine$cfg <- list(paired = FALSE)
  engine$analysis_cohorts <- list(selected = samples$sample_id[-7])
  engine$sample_manifest <- dplyr::mutate(tibble::as_tibble(samples), sample = sample_id)
  engine$counts_matrix <- matrix(0, 2, 7,
    dimnames = list(c("a", "b"), samples$sample_id))
  engine$isoform_annotation <- tibble::tibble(isoform_harmonized_id = c("a", "b"),
    associated_gene = "g", structural_category = "full-splice_match", gene_mapped = TRUE,
    cross_sample_comparable = TRUE, analyzed_for_gain_loss = TRUE)
  selected <- engine$prepare_dtu_family("selected")
  stopifnot(selected$empty, identical(as.character(selected$samples$sample_id), samples$sample_id[-7]))

  features <- tibble::tibble(gene_id = c("g", "g"), feature_id = c("g|a", "g|b"),
    isoform_harmonized_id = c("a", "b"), structural_category = "full-splice_match")
  prepared <- list(counts = tibble::as_tibble(counts), feature_map = features, samples = samples)
  fit <- list(feature = tibble::tibble(gene_id = c("g", "g"), feature_id = c("g|a", "g|b"),
    isoform_dtu_p_value = c(.01, .01)))
  stages <- list(alpha05 = list(feature = tibble::tibble(gene_id = c("g", "g"),
    feature_id = c("g|a", "g|b"), stageR_isoform_adjusted_alpha05 = c(.025, .025))))
  effects <- dtu_unpaired_effect_tables(prepared, fit, stages, c(.05, .10))
  a <- effects$feature[effects$feature$feature_id == "g|a", ]
  b <- effects$feature[effects$feature$feature_id == "g|b", ]
  stopifnot(abs(a$mean_primary_fraction - 1) < 1e-12,
    abs(a$mean_recurrent_fraction - 1/3) < 1e-12,
    abs(a$mean_delta_if + 2/3) < 1e-12,
    abs(b$mean_delta_if - 2/3) < 1e-12,
    a$evaluable_reference_samples == 2L, a$evaluable_comparison_samples == 4L,
    nrow(effects$usage) == 14L, all(is.na(effects$usage$isoform_fraction[effects$usage$sample_id == "s3"])),
    all(effects$feature$stageR_confirmed_alpha05),
    all(!effects$feature$stageR_confirmed_alpha10),
    nrow(effects$driver) == 1L, nrow(effects$paired) == 0L,
    !any(grepl("directional_support|paired_delta", names(effects$feature))),
    !any(grepl("directional_support", names(effects$driver))))
  changed <- prepared; changed$samples$patient <- renamed$patient
  stopifnot(identical(effects$feature,
    dtu_unpaired_effect_tables(changed, fit, stages, c(.05, .10))$feature))
  # The independent report path must accept no paired architecture or support
  # fields. Build each graph through ggplot rather than inspecting source text.
  plotting <- new.env(parent = environment())
  sys.source(file.path(.repo_root, "R", "dtu.R"), plotting)
  plotting$seen <- character()
  plotting$rna_save_plot <- function(plot, stem, ...) {
    ggplot2::ggplot_build(plot)
    plotting$seen <- c(plotting$seen, stem)
    invisible(NULL)
  }
  cfg <- list(comparison = list(paired = FALSE, reference = "control", comparison = "treated"),
    reporting = list(enabled = TRUE), multiplicity = list(fdr_thresholds = c(.05, .10)),
    dtu = list(maximum_gene_plots = 1L, isoforms_per_gene_example = 2L))
  result <- list(effects = effects, gene = dplyr::left_join(
    tibble::tibble(gene_id = "g", dtu_fdr = .025, raw_dtu_p_value = .01, dtu_p_value = .01),
    effects$driver, by = "gene_id"))
  plotting$dtu_unpaired_plot_outputs(result, cfg, list(), .05)
  plotting$dtu_unpaired_gene_examples(result, cfg, list())
  stopifnot(length(plotting$seen) == 2L,
    all(grepl("independent", plotting$seen)))
  message("Independent DTU design, intercept-only null, sample selection, effects and reports passed.")
  if (!requireNamespace("DRIMSeq", quietly = TRUE) ||
      as.character(utils::packageVersion("DRIMSeq")) != "1.38.0") {
    message("SKIP native independent DTU fit: exact DRIMSeq 1.38.0 unavailable")
    return(invisible(NULL))
  }
  # Exercise the real audited precision, regression and likelihood-ratio paths
  # on bounded independent samples. Positive, dispersed multinomial counts
  # avoid deliberately degenerate boundary examples already covered elsewhere.
  native_dir <- tempfile("native_independent_dtu_"); dir.create(native_dir)
  on.exit(unlink(native_dir, recursive = TRUE), add = TRUE)
  engine$cfg <- list(paired = FALSE, model_contract_version = "independent-DTU-v1",
    seed = 1357L, n_workers = 1L, rng = dtu_rng_spec(), force_refit = TRUE,
    analysis_profile = "production", pilot_max_genes = 100L,
    fdr = list(primary = .05, exploratory = .10),
    filters = list(minimum_isoforms_per_gene = 2L,
      minimum_isoform_fraction = .01, minimum_samples_isoform_fraction = 1L,
      minimum_gene_count = 20L, minimum_samples_gene_count = 3L,
      minimum_isoform_count = 1L, minimum_samples_isoform_count = 1L))
  genes <- sprintf("native%02d", 1:12)
  native_counts <- dtu_with_rng(NULL, 2468L, {
    out <- matrix(0L, length(genes) * 3L, nrow(samples))
    for (g in seq_along(genes)) for (j in seq_len(nrow(samples))) {
      expected <- if (g <= 6L && samples$timepoint[j] == "recurrent") c(.20, .35, .45) else c(.45, .35, .20)
      latent <- stats::rgamma(3L, shape = expected * 20)
      out[(3L * g - 2L):(3L * g), j] <- as.integer(stats::rmultinom(1L, 500L, latent)) + 1L
    }
    out
  })
  dimnames(native_counts) <- list(paste0(rep(genes, each = 3L), "_", rep(1:3, length(genes))), samples$sample_id)
  engine$counts_matrix <- native_counts
  engine$isoform_annotation <- tibble::tibble(isoform_harmonized_id = rownames(native_counts),
    associated_gene = rep(genes, each = 3L), structural_category = "full-splice_match",
    gene_mapped = TRUE, cross_sample_comparable = TRUE, analyzed_for_gain_loss = TRUE)
  engine$gencode_complexity <- tibble::tibble(gene_id = genes, gencode_annotated_transcripts = 3L)
  engine$analysis_cohorts <- list(independent = samples$sample_id, excluded = samples$sample_id[-7])
  engine$m06 <- rna_dirs(native_dir); engine$m06_root <- native_dir
  engine$cohort_folder <- c(independent = "independent", excluded = "excluded")
  for (cohort in names(engine$analysis_cohorts)) rna_dirs(file.path(native_dir, cohort))
  engine$write_table <- rna_write_tsv; engine$atomic_save_rds <- rna_atomic_rds
  engine$log_msg <- function(...) invisible(NULL)
  engine$make_bpparam <- function(workers) BiocParallel::SerialParam(progressbar = FALSE)
  engine$input_fingerprint <- "bounded_native_independent_fixture"
  engine$dtu_prepared <- lapply(names(engine$analysis_cohorts), engine$prepare_dtu_family)
  names(engine$dtu_prepared) <- names(engine$analysis_cohorts)
  state <- dtu_rng_snapshot()
  engine$dtu_fits <- lapply(engine$dtu_prepared, engine$fit_drimseq_model)
  stopifnot(identical(state, dtu_rng_snapshot()))
  for (cohort in names(engine$dtu_fits)) {
    prepared <- engine$dtu_prepared[[cohort]]; fit <- engine$dtu_fits[[cohort]]
    null <- methods::slot(fit$fitted_dm, "design_fit_null")
    full <- methods::slot(fit$fitted_dm, "design_fit_full")
    expected_samples <- engine$analysis_cohorts[[cohort]]
    stopifnot(identical(as.character(prepared$samples$sample_id), expected_samples),
      nrow(full) == length(expected_samples), ncol(full) == 2L,
      qr(full)$rank == 2L, ncol(null) == 1L, all(null == 1),
      identical(colnames(null), "(Intercept)"),
      setequal(fit$gene$gene_id, genes), nrow(fit$feature) == 36L,
      all(fit$gene$fit_status == "finite_test"),
      all(is.finite(fit$gene$full_model_log_likelihood)),
      all(is.finite(fit$gene$null_model_log_likelihood)),
      isTRUE(all.equal(fit$gene$dtu_bh, stats::p.adjust(fit$gene$dtu_p_value, "BH"))))
    if (requireNamespace("stageR", quietly = TRUE)) for (alpha in c(.05, .10)) {
      stage <- engine$run_stageR_alpha(cohort, alpha)
      column <- paste0("stageR_gene_adjusted_", dtu_alpha_tag(alpha))
      checked <- merge(stage$gene, fit$gene, by = "gene_id")
      stopifnot(stage$available, stage$screening_genes == 12L,
        stage$confirmation_family_genes == 12L, nrow(stage$feature) == 36L,
        setequal(stage$feature$feature_id, prepared$feature_map$feature_id),
        max(abs(checked[[column]] - checked$dtu_bh)) < 1e-8)
    }
  }
  stopifnot(!"s7" %in% DRIMSeq::samples(engine$dtu_fits$excluded$fitted_dm)$sample_id)
  if (requireNamespace("stageR", quietly = TRUE)) {
    # Missing numerical feature evidence must become conservative P=1 while
    # the gene and transcript correction families stay complete.
    conservative <- engine$dtu_fits$independent
    conservative$feature$isoform_dtu_p_value[1L] <- NA_real_
    conservative$gene$dtu_p_value[1L] <- 1
    conservative$gene$dtu_bh <- stats::p.adjust(conservative$gene$dtu_p_value, "BH")
    engine$dtu_fits$independent <- conservative
    stage <- engine$run_stageR_alpha("independent", .05)
    stopifnot(stage$conservative_feature_P1 == 1L, stage$screening_genes == 12L,
      stage$confirmation_family_genes == 12L, nrow(stage$feature) == 36L,
      setequal(stage$gene$gene_id, genes))
  }
  message("Native independent DTU 3-vs-4 and excluded 3-vs-3 fits, actual intercept-only null and complete stageR families passed.")
})()
