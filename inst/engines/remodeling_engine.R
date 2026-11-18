# Metadata-driven extraction of the complete Module 04 v0.6.0 analysis engine.
# Statistical functions are in remodeling_functions.R; no private Rmd or input checkpoint is required.
run_remodeling_engine <- function(inputs, public_cfg, cohorts, stage_dir) {
  cfg <- rna_remodeling_engine_config(public_cfg)
  params <- rna_remodeling_parameters(public_cfg)
  m04_root <- normalizePath(stage_dir, winslash = "/", mustWork = FALSE)
  m04 <- rna_dirs(stage_dir)
  m04$logs <- file.path(stage_dir, "logs")
  dir.create(m04$logs, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(m04$logs, "remodeling_run.log")
  log_msg <- function(..., level = "INFO") {
    line <- rna_log(..., level = level)
    cat(line, "\n", file = log_file, append = TRUE, sep = "")
    invisible(line)
  }
  primary_cohort_name <- names(cohorts)[[1L]]
  sample_manifest <- as_tibble(inputs$sample_manifest) %>%
    mutate(sample = as.character(sample), patient = as.character(patient),
           timepoint = as.character(timepoint)) %>% arrange(patient, timepoint)
  pair_keys <- paste(sample_manifest$patient, sample_manifest$timepoint, sep = "::")
  if (anyDuplicated(pair_keys) || anyDuplicated(sample_manifest$sample) ||
      !all(sample_manifest$timepoint %in% c("primary", "recurrent")) ||
      any(table(sample_manifest$patient) != 2L)) {
    stop("Remodeling requires one reference and one comparison sample per included subject.")
  }
  if (!length(cohorts) || is.null(names(cohorts)) || anyDuplicated(names(cohorts)) ||
      any(!nzchar(names(cohorts))) || any(!unlist(cohorts) %in% sample_manifest$patient)) {
    stop("Cohorts must be named patient-ID vectors drawn from the input metadata.")
  }
  if (!setequal(cohorts[[1L]], unique(sample_manifest$patient))) {
    stop("The primary remodeling cohort must contain all imported complete pairs; the retained-isoform mask is frozen there.")
  }
  if (cfg$thresholds$min_evaluable_pairs < 2L || cfg$bootstrap$replicates < 2L ||
      cfg$empirical_null$replicates < 1L || cfg$geometry$replicates_per_stratum < 1L ||
      cfg$gene_null$replicates < 1L || cfg$checkpoint_every_genes < 1L ||
      cfg$gene_null$covariance_initial_replicates < 2L ||
      cfg$gene_null$covariance_maximum_replicates < cfg$gene_null$covariance_initial_replicates ||
      cfg$bootstrap$quantile_type < 1L || cfg$bootstrap$quantile_type > 9L ||
      any(cfg$inference$fdr_thresholds <= 0 | cfg$inference$fdr_thresholds >= 1)) {
    stop("Invalid remodeling replicate, threshold, or covariance configuration.")
  }
  # Source function environments are scoped to this run (never .GlobalEnv).
  statistical_functions <- c("probability_vector", "entropy_bits", "jsd_base2")
  # All engine function definitions are re-bound below by the wrapper to this frame.
  engine_env <- parent.env(environment())
  for (nm in ls(engine_env, all.names = TRUE)) {
    obj <- get(nm, envir = engine_env)
    if (is.function(obj) && identical(environment(obj), engine_env) &&
        !nm %in% c("run_remodeling_engine", "rna_remodeling_parameters", "rna_remodeling_engine_config")) {
      local_obj <- obj
      environment(local_obj) <- environment()
      assign(nm, local_obj, envir = environment())
    }
  }
  scientific_cfg <- cfg
  scientific_cfg$project_root <- NULL
  scientific_cfg$strict_completion <- NULL
  scientific_cfg$resume_cached_results <- NULL
  scientific_cfg$checkpoint_every_genes <- NULL
  input_fingerprint <- digest::digest(list(
    inputs = inputs$fingerprint, counts = inputs$counts_matrix,
    annotation = inputs$isoform_annotation, manifest = sample_manifest,
    cohorts = cohorts, config = scientific_cfg,
    seed = cfg$seed, RNG = c("Mersenne-Twister", "Inversion", "Rejection"),
    engine = "remodeling-v0.6.0-public-0.1.0-seed-aware-cache",
    code_contents_sha256 = public_cfg$.remodeling_code_sha256
  ), algo = "sha256")
  cache_input_fingerprint <- input_fingerprint
  cache_thresholds <- cfg$thresholds
  cache_thresholds$patient_support <- NULL
  cache_thresholds$primary_patient_support <- NULL
  legacy_input_fingerprints_v0_5 <- character()
  legacy_thresholds_v0_5 <- cfg$thresholds
  figure_records <- list()
  plotting_enabled <- !identical(public_cfg$reporting$make_plots, FALSE) && !identical(public_cfg$reporting$enabled, FALSE)
  paper_theme <- function(base_size = NULL) rna_theme(public_cfg)
  save_publication_plot <- function(plot, filename, width = 8, height = 6, dpi = 600) {
    if (identical(public_cfg$reporting$make_plots, FALSE) ||
        identical(public_cfg$reporting$enabled, FALSE)) return(tibble())
    paths <- rna_save_plot(plot, filename, m04, public_cfg, width = width, height = height)
    if (!length(paths)) return(tibble())
    tibble(figure_name = filename, format = tolower(tools::file_ext(paths)), path = as.character(paths),
           width_inches = width, height_inches = height)
  }
  write_tsv(as_tibble(flatten_config_for_record(cfg)), file.path(m04$tables, "methods_parameters.tsv"))
  atomic_save_rds(cfg, file.path(m04$objects, "config.rds"))
  if (length(public_cfg$.remodeling_code_sha256)) write_tsv(tibble(
    component=names(public_cfg$.remodeling_code_sha256),sha256=unname(public_cfg$.remodeling_code_sha256)),
    file.path(m04$qc,"engine_source_code_sha256.tsv"))

  depth_cuts <- as.numeric(cfg$split_null_public$pooled_gene_depth_bin_lower_bounds %||% c(100,500,2000,10000))
  complexity_cuts <- as.numeric(cfg$split_null_public$isoform_complexity_breaks %||% c(2,4,8,20))
  if (is.unsorted(depth_cuts, strictly = TRUE) || is.unsorted(complexity_cuts, strictly = TRUE)) stop("Stratification cut points must be strictly increasing.")
  depth_bin_labels <- c(paste0("<", depth_cuts[[1L]]),
    if (length(depth_cuts)>1L) paste0(head(depth_cuts,-1L), "-", tail(depth_cuts,-1L)-1),
    paste0(">=", tail(depth_cuts,1L)))
  complexity_bin_labels <- c(as.character(complexity_cuts[[1L]]),
    if(length(complexity_cuts)>1L) paste0(head(complexity_cuts,-1L)+1, "-",tail(complexity_cuts,-1L)),
    paste0(">", tail(complexity_cuts,1L)))

# ---- m04-universe (source-derived computation) ----
log_msg("Remodeling step: universe.")
m02_object <- inputs
required_object_elements <- c("counts_matrix", "isoform_annotation", "sample_manifest")
if (!all(required_object_elements %in% names(m02_object))) {
  stop("The Module 02 count object lacks required elements: ",
       paste(setdiff(required_object_elements, names(m02_object)), collapse = ", "))
}

counts_matrix <- as.matrix(m02_object$counts_matrix)
isoform_annotation <- as_tibble(m02_object$isoform_annotation)

if (nrow(counts_matrix) != nrow(isoform_annotation) ||
    !identical(rownames(counts_matrix), as.character(isoform_annotation$isoform_harmonized_id))) {
  stop("Module 02 count rows and annotation rows are not aligned.")
}
if (!setequal(colnames(counts_matrix), sample_manifest$sample)) {
  stop("Module 02 count columns do not match the frozen sample manifest.")
}
counts_matrix <- counts_matrix[, sample_manifest$sample, drop = FALSE]
if (anyNA(counts_matrix) || any(counts_matrix < 0L) ||
    any(abs(counts_matrix - round(counts_matrix)) > 1e-8)) {
  stop("The remodeling analysis requires nonnegative integer FL counts.")
}
storage.mode(counts_matrix) <- "integer"

required_annotation_columns <- c(
  "isoform_harmonized_id", "associated_gene", "structural_category",
  "mapping_resolution", "cross_sample_comparable", "gene_mapped",
  "analyzed_for_gain_loss"
)
if (!all(required_annotation_columns %in% names(isoform_annotation))) {
  stop("Module 02 isoform annotation lacks: ",
       paste(setdiff(required_annotation_columns, names(isoform_annotation)), collapse = ", "))
}

analysis_row <- which(
  isoform_annotation$analyzed_for_gain_loss %in% TRUE &
    isoform_annotation$cross_sample_comparable %in% TRUE &
    isoform_annotation$gene_mapped %in% TRUE
)
if (!length(analysis_row)) stop("No cross-sample-comparable isoforms passed the locked universe.")

analysis_counts <- counts_matrix[analysis_row, , drop = FALSE]
analysis_annotation <- isoform_annotation[analysis_row, , drop = FALSE] %>%
  mutate(
    associated_gene = as.character(associated_gene),
    isoform_harmonized_id = as.character(isoform_harmonized_id),
    pooled_FL = as.numeric(rowSums(analysis_counts)),
    retained_for_remodeling = pooled_FL >= cfg$thresholds$min_isoform_total_fl
  )

gene_rows <- split(seq_len(nrow(analysis_annotation)), analysis_annotation$associated_gene)
patient_ids <- unique(sample_manifest$patient)
sample_by_patient <- sample_manifest %>%
  select(patient, timepoint, sample) %>%
  pivot_wider(names_from = timepoint, values_from = sample) %>%
  arrange(patient)

make_gene_object <- function(gene_id) {
  rows_all <- gene_rows[[gene_id]]
  retained_local <- which(analysis_annotation$retained_for_remodeling[rows_all])
  all_counts <- analysis_counts[rows_all, , drop = FALSE]
  retained_counts <- if (length(retained_local)) {
    all_counts[retained_local, , drop = FALSE]
  } else {
    matrix(0L, nrow = 0L, ncol = ncol(all_counts),
           dimnames = list(character(), colnames(all_counts)))
  }
  retained_ids <- analysis_annotation$isoform_harmonized_id[rows_all][retained_local]
  rownames(retained_counts) <- retained_ids
  original_gene_count <- colSums(all_counts)
  retained_gene_count <- if (nrow(retained_counts)) colSums(retained_counts) else rep(0, ncol(all_counts))
  retained_mass <- ifelse(original_gene_count > 0, retained_gene_count / original_gene_count, NA_real_)
  names(retained_mass) <- colnames(all_counts)
  sample_quantifiable <- original_gene_count >= cfg$thresholds$min_gene_sample_count &
    retained_gene_count > 0 & retained_mass >= cfg$thresholds$min_retained_mass
  names(sample_quantifiable) <- colnames(all_counts)
  pair_quantifiable <- vapply(seq_len(nrow(sample_by_patient)), function(index) {
    sample_quantifiable[[sample_by_patient$primary[[index]]]] &&
      sample_quantifiable[[sample_by_patient$recurrent[[index]]]]
  }, logical(1))
  names(pair_quantifiable) <- sample_by_patient$patient
  list(
    gene_id = gene_id,
    rows_all = rows_all,
    retained_local = retained_local,
    all_counts = all_counts,
    counts = t(retained_counts),
    transcript_ids = retained_ids,
    original_gene_count = original_gene_count,
    retained_gene_count = retained_gene_count,
    retained_mass = retained_mass,
    sample_quantifiable = sample_quantifiable,
    pair_quantifiable = pair_quantifiable,
    n_isoforms_total = length(rows_all),
    n_isoforms_retained = length(retained_ids),
    pooled_gene_count = sum(all_counts),
    pooled_retained_count = sum(retained_counts)
  )
}

gene_objects <- lapply(sort(names(gene_rows)), make_gene_object)
names(gene_objects) <- vapply(gene_objects, `[[`, character(1), "gene_id")

sample_gene_evaluability <- rbindlist(lapply(gene_objects, function(object) {
  tibble(
    gene_id = object$gene_id,
    sample = names(object$original_gene_count),
    original_gene_FL = as.numeric(object$original_gene_count),
    retained_gene_FL = as.numeric(object$retained_gene_count),
    retained_mass = as.numeric(object$retained_mass),
    sample_quantifiable = as.logical(object$sample_quantifiable)
  )
}), fill = TRUE) %>%
  left_join(sample_manifest %>% select(sample, patient, timepoint), by = "sample")

gene_quantifiability <- tibble(
  gene_id = names(gene_objects),
  pooled_gene_FL = vapply(gene_objects, `[[`, numeric(1), "pooled_gene_count"),
  pooled_retained_FL = vapply(gene_objects, `[[`, numeric(1), "pooled_retained_count"),
  n_isoforms_total = vapply(gene_objects, `[[`, integer(1), "n_isoforms_total"),
  n_isoforms_retained = vapply(gene_objects, `[[`, integer(1), "n_isoforms_retained"),
  n_evaluable_pairs = vapply(gene_objects, function(x) sum(x$pair_quantifiable), integer(1))
) %>%
  mutate(
    quantifiable = n_isoforms_retained >= 2L &
      n_evaluable_pairs >= cfg$thresholds$min_evaluable_pairs,
    quantifiability_reason = case_when(
      n_isoforms_retained < 2L ~ "fewer_than_two_retained_isoforms",
      n_evaluable_pairs < cfg$thresholds$min_evaluable_pairs ~ "insufficient_complete_quantifiable_pairs",
      TRUE ~ "quantifiable"
    ),
    gene_depth_bin = cut(
      pooled_retained_FL,
      breaks = c(-Inf, as.numeric(cfg$split_null_public$pooled_gene_depth_bin_lower_bounds %||% c(100,500,2000,10000)) - 1, Inf),
      labels = depth_bin_labels
    ),
    isoform_complexity_bin = cut(
      n_isoforms_retained,
      breaks = c(-Inf, as.numeric(cfg$split_null_public$isoform_complexity_breaks %||% c(2,4,8,20)), Inf),
      labels = complexity_bin_labels
    )
  )

universe_audit <- tibble(
  metric = c(
    "isoforms_in_Module02_matrix", "locked_comparable_isoforms",
    "retained_isoforms_after_pooled_FL_filter", "genes_in_comparable_universe",
    "quantifiable_genes", "nonquantifiable_genes", "samples", "matched_pairs",
    "total_FL_in_locked_comparable_universe"
  ),
  value = c(
    nrow(counts_matrix), nrow(analysis_annotation),
    sum(analysis_annotation$retained_for_remodeling), length(gene_objects),
    sum(gene_quantifiability$quantifiable), sum(!gene_quantifiability$quantifiable),
    nrow(sample_manifest), n_distinct(sample_manifest$patient), sum(analysis_counts)
  )
)

write_tsv(analysis_annotation, file.path(m04$tables, "remodeling_isoform_universe.tsv.gz"))
write_tsv(sample_gene_evaluability, file.path(m04$tables, "sample_gene_evaluability.tsv.gz"))
write_tsv(gene_quantifiability, file.path(m04$tables, "gene_quantifiability.tsv"))
write_tsv(universe_audit, file.path(m04$qc, "remodeling_universe_audit.tsv"))
atomic_save_rds(gene_objects, file.path(m04$objects, "gene_objects.rds"))

log_msg(
  "Remodeling universe prepared: ", comma(sum(gene_quantifiability$quantifiable)),
  " quantifiable genes and ", comma(sum(analysis_annotation$retained_for_remodeling)),
  " retained isoforms."
)

# ---- m04-observed-metrics (source-derived computation) ----
log_msg("Remodeling step: observed-metrics.")
analysis_cohorts <- cohorts

observed_results <- imap(analysis_cohorts, function(included_patients, cohort_name) {
  log_msg("Calculating observed remodeling metrics for ", cohort_name, ".")
  lapply(gene_objects, analyze_gene_observed,
         included_patients = included_patients, cohort = cohort_name)
})

observed_gene <- rbindlist(lapply(observed_results, function(cohort_result) {
  rbindlist(lapply(cohort_result, `[[`, "gene"), fill = TRUE)
}), fill = TRUE) %>% as_tibble()

observed_isoform <- rbindlist(lapply(observed_results, function(cohort_result) {
  rbindlist(lapply(cohort_result, `[[`, "isoform"), fill = TRUE)
}), fill = TRUE) %>% as_tibble()

observed_patient <- rbindlist(lapply(observed_results, function(cohort_result) {
  rbindlist(lapply(cohort_result, `[[`, "patient"), fill = TRUE)
}), fill = TRUE) %>% as_tibble()

primary_observed <- observed_gene %>% filter(cohort == primary_cohort_name)
if (!any(primary_observed$quantifiable)) {
  # A valid dataset can have no genes meeting the configured eligibility rules.
  # Preserve the full screened universe and report absent inference explicitly.
  no_gene_table <- observed_gene %>% mutate(
    architecture_class="NQ",architecture_description="non_quantifiable",
    is_remodeled_AE=FALSE,is_high_confidence_remodeled=FALSE,is_high_confidence_class=FALSE,
    downstream_primary_remodeling_set=FALSE,downstream_class_stratified_set=FALSE,
    gene_model_state="outside_registered_testing_family",gene_model_changes_class=FALSE,
    authoritative_status="no_quantifiable_genes",effect_threshold_gate_pass=NA)
  gates <- tibble(check="No eligible genes; resampling and inferential checks not applicable",
                  enabled=TRUE,blocking=FALSE,pass=TRUE)
  write_tsv(no_gene_table,file.path(m04$tables,"gene_remodeling_classes_all_cohorts.tsv.gz"))
  write_tsv(gates,file.path(m04$qc,"completion_gates.tsv"))
  write_tsv(observed_gene,file.path(m04$tables,"observed_gene_remodeling_metrics.tsv.gz"))
  empty <- tibble(cohort=character(),gene_id=character())
  result <- list(gene_table=no_gene_table,
    primary_gene_table=no_gene_table %>% filter(cohort==primary_cohort_name),
    isoform_table=observed_isoform,patient_table=observed_patient,observed_gene=observed_gene,
    patient_inference=empty,bootstrap=empty,cohort_model=empty,model_family=empty,
    config=cfg,fingerprint=input_fingerprint,paths=m04,completion_gates=gates,complete=TRUE,
    execution_profile=cfg$analysis_profile,analysis_status="no_quantifiable_genes",
    authoritative_classes=TRUE)
  atomic_save_rds(result,file.path(m04$objects,"remodeling_results.rds"))
  return(result)
}


safe_correlation <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(tibble(n = sum(keep), estimate = NA_real_, p_value = NA_real_))
  }
  result <- suppressWarnings(cor.test(x[keep], y[keep], method = method, exact = FALSE))
  tibble(n = sum(keep), estimate = unname(result$estimate), p_value = result$p.value)
}

metric_correlations <- bind_rows(
  safe_correlation(primary_observed$jsd, primary_observed$max_abs_delta_if, "pearson") %>%
    mutate(comparison = "JSD_vs_max_abs_delta_IF", method = "Pearson"),
  safe_correlation(primary_observed$jsd, primary_observed$max_abs_delta_if, "spearman") %>%
    mutate(comparison = "JSD_vs_max_abs_delta_IF", method = "Spearman"),
  safe_correlation(primary_observed$jsd, primary_observed$total_variation_distance, "pearson") %>%
    mutate(comparison = "JSD_vs_total_variation", method = "Pearson"),
  safe_correlation(primary_observed$jsd, primary_observed$total_variation_distance, "spearman") %>%
    mutate(comparison = "JSD_vs_total_variation", method = "Spearman")
) %>% select(comparison, method, everything())

effect_definition_gene <- primary_observed %>%
  filter(quantifiable) %>%
  transmute(
    gene_id,
    jsd = jsd,
    max_abs_delta_if = max_abs_delta_if,
    combined_JSD_or_deltaIF = jsd >= cfg$thresholds$jsd |
      max_abs_delta_if >= cfg$thresholds$delta_if,
    JSD_only_definition = jsd >= cfg$thresholds$jsd,
    deltaIF_only_definition = max_abs_delta_if >= cfg$thresholds$delta_if,
    effect_definition_pattern = case_when(
      JSD_only_definition & deltaIF_only_definition ~ "Both metrics pass",
      JSD_only_definition ~ "JSD threshold only",
      deltaIF_only_definition ~ "Delta-IF threshold only",
      TRUE ~ "Neither threshold"
    )
  )

effect_definition_summary <- bind_rows(
  effect_definition_gene %>% summarise(
    definition = "Combined: JSD or maximum |delta IF|",
    n_genes = sum(combined_JSD_or_deltaIF)
  ),
  effect_definition_gene %>% summarise(
    definition = "JSD threshold alone",
    n_genes = sum(JSD_only_definition)
  ),
  effect_definition_gene %>% summarise(
    definition = "Maximum |delta IF| threshold alone",
    n_genes = sum(deltaIF_only_definition)
  )
) %>%
  mutate(
    quantifiable_genes = nrow(effect_definition_gene),
    fraction = n_genes / quantifiable_genes
  )

effect_definition_overlap <- effect_definition_gene %>%
  count(effect_definition_pattern, name = "n_genes") %>%
  mutate(fraction = n_genes / sum(n_genes))

metric_identity_qc <- primary_observed %>%
  filter(quantifiable) %>%
  summarise(
    n_quantifiable = n(),
    n_jsd_outside_unit_interval = sum(jsd < -1e-12 | jsd > 1 + 1e-12),
    n_delta_if_outside_unit_interval = sum(max_abs_delta_if < -1e-12 |
                                            max_abs_delta_if > 1 + 1e-12),
    n_tvd_outside_unit_interval = sum(total_variation_distance < -1e-12 |
                                       total_variation_distance > 1 + 1e-12),
    n_max_delta_exceeds_tvd = sum(max_abs_delta_if > total_variation_distance + 1e-10),
    maximum_log_base_identity_error = max(
      abs(jsd - jsd_natural_log / log(2)), na.rm = TRUE
    ),
    n_nonfinite_metrics = sum(!is.finite(jsd) | !is.finite(max_abs_delta_if) |
                                !is.finite(total_variation_distance))
  )

isoform_contribution_qc <- observed_isoform %>%
  filter(cohort == primary_cohort_name) %>%
  group_by(gene_id) %>%
  summarise(
    sum_jsd_contribution = sum(jsd_contribution),
    sum_primary_if = sum(mean_if_primary),
    sum_recurrent_if = sum(mean_if_recurrent),
    .groups = "drop"
  ) %>%
  left_join(primary_observed %>% select(gene_id, jsd), by = "gene_id") %>%
  mutate(
    jsd_contribution_error = abs(sum_jsd_contribution - jsd),
    primary_composition_error = abs(sum_primary_if - 1),
    recurrent_composition_error = abs(sum_recurrent_if - 1)
  )

write_tsv(observed_gene, file.path(m04$tables, "observed_gene_remodeling_metrics.tsv.gz"))
write_tsv(observed_isoform, file.path(m04$tables, "observed_isoform_remodeling_metrics.tsv.gz"))
write_tsv(observed_patient, file.path(m04$tables, "observed_patient_remodeling_metrics.tsv.gz"))
write_tsv(metric_correlations, file.path(m04$tables, "metric_correlations.tsv"))
write_tsv(effect_definition_gene, file.path(m04$tables, "effect_definition_gene_comparison.tsv.gz"))
write_tsv(effect_definition_summary, file.path(m04$tables, "effect_definition_summary.tsv"))
write_tsv(effect_definition_overlap, file.path(m04$tables, "effect_definition_overlap.tsv"))
write_tsv(metric_identity_qc, file.path(m04$qc, "metric_identity_qc.tsv"))
write_tsv(isoform_contribution_qc, file.path(m04$qc, "isoform_contribution_qc.tsv.gz"))


# ---- m04-patient-level-inference (source-derived computation) ----
log_msg("Remodeling step: patient-level-inference.")
patient_pair_lrt <- function(object, patient_id) {
  pair_row <- sample_by_patient %>%
    filter(patient == patient_id)
  if (nrow(pair_row) != 1L) {
    stop("Expected one primary-recurrent sample pair for ", patient_id, ".")
  }
  primary_sample <- pair_row$primary[[1]]
  recurrent_sample <- pair_row$recurrent[[1]]
  primary_count <- as.numeric(object$counts[primary_sample, , drop = TRUE])
  recurrent_count <- as.numeric(object$counts[recurrent_sample, , drop = TRUE])
  names(primary_count) <- names(recurrent_count) <- object$transcript_ids
  primary_total <- sum(primary_count)
  recurrent_total <- sum(recurrent_count)

  output <- tibble(
    gene_id = object$gene_id,
    patient = patient_id,
    primary_sample = primary_sample,
    recurrent_sample = recurrent_sample,
    n_isoforms_retained = object$n_isoforms_retained,
    primary_retained_FL = primary_total,
    recurrent_retained_FL = recurrent_total,
    primary_retained_mass = as.numeric(object$retained_mass[[primary_sample]]),
    recurrent_retained_mass = as.numeric(object$retained_mass[[recurrent_sample]]),
    pair_quantifiable = object$pair_quantifiable[[patient_id]] %in% TRUE,
    patient_jsd = NA_real_,
    patient_max_abs_delta_if = NA_real_,
    patient_driver_isoform = NA_character_,
    patient_driver_delta_if = NA_real_,
    patient_driver_tie_n = NA_integer_,
    patient_entropy_primary = NA_real_,
    patient_entropy_recurrent = NA_real_,
    patient_delta_entropy = NA_real_,
    patient_entropy_direction = NA_character_,
    patient_effect_pass = FALSE,
    patient_test_in_family = FALSE,
    patient_test_available = FALSE,
    patient_test_reason = "not_evaluated",
    patient_test_components_before_pooling = NA_integer_,
    patient_test_components_after_pooling = NA_integer_,
    patient_test_pooled_rare_components = NA_integer_,
    patient_test_minimum_expected_count = NA_real_,
    patient_test_fraction_expected_below_5 = NA_real_,
    patient_test_statistic = NA_real_,
    patient_test_df = NA_integer_,
    patient_test_p = NA_real_
  )

  if (object$n_isoforms_retained < 2L) {
    output$patient_test_reason <- "fewer_than_two_retained_isoforms"
    return(output)
  }
  if (!isTRUE(output$pair_quantifiable[[1]])) {
    output$patient_test_reason <- "pair_not_quantifiable_under_locked_mass_and_count_rules"
    return(output)
  }
  if (primary_total <= 0 || recurrent_total <= 0) {
    output$patient_test_reason <- "zero_retained_count_in_one_condition"
    return(output)
  }

  primary_composition <- primary_count / primary_total
  recurrent_composition <- recurrent_count / recurrent_total
  delta_if <- recurrent_composition - primary_composition
  maximum_delta <- max(abs(delta_if))
  driver_candidates <- which(abs(abs(delta_if) - maximum_delta) <= 1e-12)
  driver_index <- driver_candidates[[1]]
  delta_entropy <- entropy_bits(recurrent_composition) -
    entropy_bits(primary_composition)
  output <- output %>%
    mutate(
      patient_jsd = jsd_base2(primary_composition, recurrent_composition),
      patient_max_abs_delta_if = maximum_delta,
      patient_driver_isoform = object$transcript_ids[[driver_index]],
      patient_driver_delta_if = delta_if[[driver_index]],
      patient_driver_tie_n = length(driver_candidates),
      patient_entropy_primary = entropy_bits(primary_composition),
      patient_entropy_recurrent = entropy_bits(recurrent_composition),
      patient_delta_entropy = delta_entropy,
      patient_entropy_direction = case_when(
        delta_entropy >= cfg$thresholds$entropy_delta ~
          "more_diffuse_in_recurrence",
        delta_entropy <= -cfg$thresholds$entropy_delta ~
          "more_focused_in_recurrence",
        TRUE ~ "entropy_stable"
      ),
      patient_effect_pass = patient_jsd >= cfg$thresholds$jsd |
        patient_max_abs_delta_if >= cfg$thresholds$delta_if
    )

  if (!cfg$patient_test$enabled) {
    output$patient_test_reason <- "patient_test_disabled"
    return(output)
  }
  if (object$n_isoforms_retained > cfg$patient_test$maximum_isoforms) {
    output$patient_test_reason <- "isoform_complexity_above_registered_patient_test_scope"
    return(output)
  }
  if (primary_total < cfg$patient_test$minimum_condition_FL ||
      recurrent_total < cfg$patient_test$minimum_condition_FL) {
    output$patient_test_reason <- "condition_depth_below_registered_patient_test_floor"
    return(output)
  }

  pair_total <- primary_count + recurrent_count
  present <- pair_total > 0
  primary_present <- primary_count[present]
  recurrent_present <- recurrent_count[present]
  pair_total_present <- pair_total[present]
  common <- pair_total_present >= cfg$patient_test$minimum_pooled_component_FL
  test_primary <- primary_present[common]
  test_recurrent <- recurrent_present[common]
  rare_primary <- sum(primary_present[!common])
  rare_recurrent <- sum(recurrent_present[!common])
  rare_component_n <- sum(!common)
  if ((rare_primary + rare_recurrent) > 0L) {
    test_primary <- c(test_primary, rare_primary)
    test_recurrent <- c(test_recurrent, rare_recurrent)
  }
  component_n <- length(test_primary)
  output$patient_test_components_before_pooling <- sum(present)
  output$patient_test_components_after_pooling <- component_n
  output$patient_test_pooled_rare_components <- rare_component_n
  if (component_n < 2L) {
    output$patient_test_reason <- "fewer_than_two_components_after_rare_pooling"
    return(output)
  }

  observed <- rbind(primary = test_primary, recurrent = test_recurrent)
  row_total <- rowSums(observed)
  column_total <- colSums(observed)
  grand_total <- sum(observed)
  expected <- outer(row_total, column_total) / grand_total
  minimum_expected <- min(expected)
  fraction_below_5 <- mean(expected < 5)
  output$patient_test_minimum_expected_count <- minimum_expected
  output$patient_test_fraction_expected_below_5 <- fraction_below_5
  if (minimum_expected < cfg$patient_test$minimum_expected_count) {
    output$patient_test_reason <- "minimum_expected_count_below_registered_floor"
    return(output)
  }
  if (fraction_below_5 > cfg$patient_test$maximum_fraction_expected_below_5) {
    output$patient_test_reason <- "too_many_expected_cells_below_five"
    return(output)
  }

  positive <- observed > 0
  statistic <- 2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
  statistic <- max(0, statistic)
  degrees_freedom <- component_n - 1L
  p_value <- pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  output %>%
    mutate(
      patient_test_in_family = TRUE,
      patient_test_available = is.finite(p_value),
      patient_test_reason = if_else(
        patient_test_available,
        "available_multinomial_likelihood_ratio",
        "conservative_failure_nonfinite_likelihood_ratio"
      ),
      patient_test_statistic = if_else(patient_test_available, statistic, NA_real_),
      patient_test_df = if_else(
        patient_test_available, as.integer(degrees_freedom), NA_integer_
      ),
      patient_test_p = if_else(patient_test_available, p_value, 1)
    )
}

patient_specific_inference <- if (cfg$patient_test$enabled) {
  map_dfr(patient_ids, function(patient_id) {
    map_dfr(gene_objects, patient_pair_lrt, patient_id = patient_id)
  }) %>%
    arrange(patient, gene_id) %>%
    group_by(patient) %>%
    mutate(
      patient_test_bh = {
        adjusted <- rep(NA_real_, n())
        family_index <- which(patient_test_in_family)
        adjusted[family_index] <- p.adjust(
          patient_test_p[family_index], method = "BH"
        )
        adjusted
      }
    ) %>%
    ungroup()
} else {
  patient_pair_lrt(gene_objects[[1]], patient_ids[[1]])[0, ] %>%
    mutate(
      patient_test_bh = numeric(),
      patient_test_global_bh_audit = numeric(),
      patient_nominal_p_lt_0p05 = logical(),
      patient_bh_le_0p05 = logical(),
      patient_bh_le_0p10 = logical(),
      patient_supported_bh05_effect = logical(),
      patient_supported_bh10_effect = logical(),
      patient_call_primary = logical(),
      patient_call_sensitivity = logical(),
      patient_call_rule = character(),
      raw_patient_p_downstream_use = character()
    )
}

if (nrow(patient_specific_inference)) {
  global_patient_q <- rep(NA_real_, nrow(patient_specific_inference))
  global_family_index <- which(patient_specific_inference$patient_test_in_family)
  global_patient_q[global_family_index] <- p.adjust(
    patient_specific_inference$patient_test_p[global_family_index],
    method = "BH"
  )
  patient_specific_inference <- patient_specific_inference %>%
    mutate(
      patient_test_global_bh_audit = global_patient_q,
      patient_nominal_p_lt_0p05 = patient_test_in_family &
        patient_test_p < 0.05,
      patient_bh_le_0p05 = patient_test_in_family &
        patient_test_bh <= 0.05,
      patient_bh_le_0p10 = patient_test_in_family &
        patient_test_bh <= 0.10,
      patient_supported_bh05_effect = patient_effect_pass &
        patient_bh_le_0p05,
      patient_supported_bh10_effect = patient_effect_pass &
        patient_bh_le_0p10,
      patient_call_primary = patient_effect_pass & patient_test_in_family & patient_test_bh <= cfg$patient_test$primary_fdr,
      patient_call_sensitivity = patient_effect_pass & patient_test_in_family & patient_test_bh <= cfg$patient_test$sensitivity_fdr,
      patient_call_rule = paste0(
        "within_patient_BH<=", cfg$patient_test$primary_fdr,
        " AND (JSD>=", cfg$thresholds$jsd,
        " OR max_abs_delta_IF>=", cfg$thresholds$delta_if, ")"
      ),
      raw_patient_p_downstream_use = "forbidden"
    )
}

patient_test_qc <- if (nrow(patient_specific_inference)) {
  patient_specific_inference %>%
    group_by(patient) %>%
    summarise(
      genes_screened = n(),
      frozen_testing_family = sum(patient_test_in_family),
      available_tests = sum(patient_test_available),
      nominal_p_lt_0p05 = sum(patient_nominal_p_lt_0p05),
      bh_le_0p05 = sum(patient_bh_le_0p05),
      bh_le_0p10 = sum(patient_bh_le_0p10),
      bh05_and_effect = sum(patient_supported_bh05_effect),
      bh10_and_effect = sum(patient_supported_bh10_effect),
      finite_family_p = sum(
        patient_test_in_family & is.finite(patient_test_p)
      ),
      finite_family_bh = sum(
        patient_test_in_family & is.finite(patient_test_bh)
      ),
      pass = genes_screened == length(gene_objects) &
        frozen_testing_family == available_tests &
        finite_family_p == frozen_testing_family &
        finite_family_bh == frozen_testing_family,
      .groups = "drop"
    )
} else {
  tibble(
    patient = patient_ids,
    genes_screened = 0L,
    frozen_testing_family = 0L,
    available_tests = 0L,
    nominal_p_lt_0p05 = 0L,
    bh_le_0p05 = 0L,
    bh_le_0p10 = 0L,
    bh05_and_effect = 0L,
    bh10_and_effect = 0L,
    finite_family_p = 0L,
    finite_family_bh = 0L,
    pass = !cfg$patient_test$enabled
  )
}

patient_test_exclusions <- if (nrow(patient_specific_inference)) {
  patient_specific_inference %>%
    filter(!patient_test_in_family) %>%
    count(patient, patient_test_reason, name = "n_genes")
} else {
  tibble(patient = character(), patient_test_reason = character(), n_genes = integer())
}

write_tsv(
  patient_specific_inference,
  file.path(m04$tables, "patient_specific_remodeling_inference.tsv.gz")
)
write_tsv(patient_test_qc, file.path(m04$qc, "patient_specific_test_qc.tsv"))
write_tsv(
  patient_test_exclusions,
  file.path(m04$qc, "patient_specific_test_exclusions.tsv")
)


# ---- m04-geometry-simulation (source-derived computation) ----
log_msg("Remodeling step: geometry-simulation.")
rdirichlet_one <- function(k, concentration) {
  draw <- rgamma(k, shape = concentration, rate = 1)
  draw / sum(draw)
}

simulate_geometry_stratum <- function(k, depth, replicates, seed) {
  with_seed(seed, {
    rows <- vector("list", replicates)
    for (replicate_index in seq_len(replicates)) {
      concentration <- 10^runif(1, -0.5, 1)
      primary_true <- rdirichlet_one(k, concentration)
      independent <- rdirichlet_one(k, concentration)
      mixing <- runif(1)
      recurrent_true <- probability_vector((1 - mixing) * primary_true + mixing * independent)
      primary_count <- as.numeric(rmultinom(1, depth, primary_true))
      recurrent_count <- as.numeric(rmultinom(1, depth, recurrent_true))
      primary_observed <- primary_count / sum(primary_count)
      recurrent_observed <- recurrent_count / sum(recurrent_count)
      rows[[replicate_index]] <- tibble(
        n_isoforms = k,
        depth = depth,
        replicate = replicate_index,
        mixing_fraction = mixing,
        jsd = jsd_base2(primary_observed, recurrent_observed),
        max_abs_delta_if = max(abs(recurrent_observed - primary_observed)),
        total_variation_distance = total_variation(primary_observed, recurrent_observed),
        delta_entropy = entropy_bits(recurrent_observed) - entropy_bits(primary_observed)
      )
    }
    bind_rows(rows)
  })
}

geometry_design <- tidyr::crossing(
  n_isoforms = cfg$geometry$isoform_numbers,
  depth = cfg$geometry$depths
)

geometry_cache_path <- file.path(m04$cache, "feasible_geometry_simulation.rds")
geometry_cache_valid <- FALSE
if (cfg$resume_cached_results && file.exists(geometry_cache_path)) {
  geometry_cache <- readRDS(geometry_cache_path)
  geometry_cache_valid <- geometry_cache$input_fingerprint %in%
    c(cache_input_fingerprint, legacy_input_fingerprints_v0_5) &&
    identical(geometry_cache$geometry_configuration, cfg$geometry)
}

if (geometry_cache_valid) {
  geometry_simulation <- geometry_cache$result
  log_msg("Loaded matching feasible-geometry simulation cache.")
} else {
  geometry_simulation <- pmap_dfr(
    geometry_design,
    function(n_isoforms, depth) simulate_geometry_stratum(
      n_isoforms, depth, cfg$geometry$replicates_per_stratum,
      deterministic_seed(cfg$seed, "geometry", n_isoforms, depth)
    )
  )
  atomic_save_rds(
    list(
      input_fingerprint = cache_input_fingerprint,
      geometry_configuration = cfg$geometry,
      result = geometry_simulation
    ),
    geometry_cache_path
  )
}

geometry_envelope <- geometry_simulation %>%
  mutate(
    delta_bin = cut(
      max_abs_delta_if,
      breaks = seq(0, 1, by = 0.025), include.lowest = TRUE,
      labels = FALSE
    )
  ) %>%
  filter(!is.na(delta_bin)) %>%
  group_by(n_isoforms, depth, delta_bin) %>%
  summarise(
    max_abs_delta_if_midpoint = mean(max_abs_delta_if),
    jsd_q025 = quantile(jsd, 0.025, names = FALSE),
    jsd_q50 = quantile(jsd, 0.50, names = FALSE),
    jsd_q975 = quantile(jsd, 0.975, names = FALSE),
    n = n(),
    .groups = "drop"
  )

geometry_correlations <- geometry_simulation %>%
  group_by(n_isoforms, depth) %>%
  summarise(
    n = n(),
    pearson_JSD_deltaIF = cor(jsd, max_abs_delta_if, method = "pearson"),
    spearman_JSD_deltaIF = cor(jsd, max_abs_delta_if, method = "spearman"),
    pearson_JSD_TVD = cor(jsd, total_variation_distance, method = "pearson"),
    .groups = "drop"
  )

geometry_qc <- geometry_simulation %>%
  summarise(
    n_simulations = n(),
    n_nonfinite = sum(!is.finite(jsd) | !is.finite(max_abs_delta_if) |
                        !is.finite(total_variation_distance)),
    n_jsd_outside_unit_interval = sum(jsd < -1e-12 | jsd > 1 + 1e-12),
    n_max_delta_exceeds_tvd = sum(max_abs_delta_if > total_variation_distance + 1e-10),
    strata_completed = n_distinct(interaction(n_isoforms, depth)),
    strata_expected = nrow(geometry_design),
    minimum_replicates_per_stratum = min(table(interaction(n_isoforms, depth)))
  )

write_tsv(geometry_simulation, file.path(m04$tables, "feasible_geometry_simulation.tsv.gz"))
write_tsv(geometry_envelope, file.path(m04$tables, "feasible_geometry_envelope.tsv.gz"))
write_tsv(geometry_correlations, file.path(m04$tables, "feasible_geometry_correlations.tsv"))
write_tsv(geometry_qc, file.path(m04$qc, "feasible_geometry_qc.tsv"))

# ---- m04-empirical-null (source-derived computation) ----
log_msg("Remodeling step: empirical-null.")
null_architecture_levels <- c(
  "A", "B", "C", "D", "E", "Below_effect", "Unresolved_effect", "NQ"
)
null_gene_ids <- primary_observed %>%
  filter(quantifiable) %>%
  arrange(gene_id) %>%
  pull(gene_id)
null_gene_objects <- gene_objects[null_gene_ids]
null_gene_strata <- primary_observed %>%
  filter(gene_id %in% null_gene_ids) %>%
  select(gene_id, pooled_retained_FL, n_isoforms_retained) %>%
  left_join(
    gene_quantifiability %>%
      select(gene_id, gene_depth_bin, isoform_complexity_bin),
    by = "gene_id"
  ) %>%
  arrange(match(gene_id, null_gene_ids))

analyze_pseudo_pair <- function(object, source_samples) {
  all_count <- t(object$all_counts[, source_samples, drop = FALSE])
  half_a_all <- matrix(
    rbinom(length(all_count), size = as.vector(all_count), prob = 0.5),
    nrow = nrow(all_count), ncol = ncol(all_count), dimnames = dimnames(all_count)
  )
  half_b_all <- all_count - half_a_all
  retained <- object$retained_local
  if (length(retained) < 2L) {
    return(list(
      jsd = NA_real_, max_abs_delta_if = NA_real_, tvd = NA_real_,
      architecture = "NQ", n_evaluable_pairs = 0L,
      original_total = sum(all_count), half_a_total = sum(half_a_all),
      half_b_total = sum(half_b_all), conservation_error = 0
    ))
  }
  half_a <- half_a_all[, retained, drop = FALSE]
  half_b <- half_b_all[, retained, drop = FALSE]
  original_a <- rowSums(half_a_all)
  original_b <- rowSums(half_b_all)
  retained_a <- rowSums(half_a)
  retained_b <- rowSums(half_b)
  mass_a <- ifelse(original_a > 0, retained_a / original_a, NA_real_)
  mass_b <- ifelse(original_b > 0, retained_b / original_b, NA_real_)
  evaluable <- original_a >= cfg$thresholds$min_gene_sample_count &
    original_b >= cfg$thresholds$min_gene_sample_count &
    retained_a > 0 & retained_b > 0 &
    mass_a >= cfg$thresholds$min_retained_mass &
    mass_b >= cfg$thresholds$min_retained_mass
  if (sum(evaluable) < cfg$thresholds$min_evaluable_pairs) {
    return(list(
      jsd = NA_real_, max_abs_delta_if = NA_real_, tvd = NA_real_,
      architecture = "NQ", n_evaluable_pairs = as.integer(sum(evaluable)),
      original_total = sum(all_count), half_a_total = sum(half_a_all),
      half_b_total = sum(half_b_all),
      conservation_error = sum(abs(half_a_all + half_b_all - all_count))
    ))
  }
  composition_a <- half_a[evaluable, , drop = FALSE] / retained_a[evaluable]
  composition_b <- half_b[evaluable, , drop = FALSE] / retained_b[evaluable]
  mean_a <- colMeans(composition_a)
  mean_b <- colMeans(composition_b)
  transcript_ids <- object$transcript_ids
  dominance_a <- dominance_summary(mean_a, composition_a, transcript_ids)
  dominance_b <- dominance_summary(mean_b, composition_b, transcript_ids)
  secondary <- secondary_remodeling(
    mean_a, mean_b,
    dominance_a$dominant_id, dominance_b$dominant_id,
    dominance_a$has_dominant, dominance_b$has_dominant,
    transcript_ids
  )
  jsd_value <- jsd_base2(mean_a, mean_b)
  max_delta <- max(abs(mean_b - mean_a))
  effect_pass <- jsd_value >= cfg$thresholds$jsd ||
    max_delta >= cfg$thresholds$delta_if
  architecture <- classify_architecture(
    TRUE, effect_pass, dominance_a, dominance_b, secondary
  )
  list(
    jsd = jsd_value,
    max_abs_delta_if = max_delta,
    tvd = total_variation(mean_a, mean_b),
    architecture = architecture,
    n_evaluable_pairs = as.integer(sum(evaluable)),
    original_total = sum(all_count),
    half_a_total = sum(half_a_all),
    half_b_total = sum(half_b_all),
    conservation_error = sum(abs(half_a_all + half_b_all - all_count))
  )
}

run_null_replicate <- function(replicate_index, source_timepoint) {
  source_samples <- sample_by_patient[[source_timepoint]]
  seed <- deterministic_seed(cfg$seed, "empirical_null", source_timepoint, replicate_index)
  result <- with_seed(seed, lapply(null_gene_objects, analyze_pseudo_pair,
                                   source_samples = source_samples))
  tibble(
    gene_id = null_gene_ids,
    jsd = vapply(result, `[[`, numeric(1), "jsd"),
    max_abs_delta_if = vapply(result, `[[`, numeric(1), "max_abs_delta_if"),
    total_variation_distance = vapply(result, `[[`, numeric(1), "tvd"),
    architecture = vapply(result, `[[`, character(1), "architecture"),
    n_evaluable_pairs = vapply(result, `[[`, integer(1), "n_evaluable_pairs"),
    original_total = vapply(result, `[[`, numeric(1), "original_total"),
    half_a_total = vapply(result, `[[`, numeric(1), "half_a_total"),
    half_b_total = vapply(result, `[[`, numeric(1), "half_b_total"),
    conservation_error = vapply(result, `[[`, numeric(1), "conservation_error")
  )
}

null_cache_path <- file.path(m04$cache, "empirical_split_null_cache.rds")
n_null_gene <- length(null_gene_ids)
n_null_replicate <- cfg$empirical_null$replicates
null_sources <- c("primary", "recurrent")
null_configuration_hash <- digest::digest(
  list(
    input_fingerprint = cache_input_fingerprint,
    empirical_null = cfg$empirical_null,
    thresholds = cache_thresholds,
    genes = null_gene_ids
  ), algo = "sha256"
)
legacy_null_configuration_hashes <- vapply(
  legacy_input_fingerprints_v0_5,
  function(legacy_fingerprint) {
    digest::digest(
      list(
        input_fingerprint = legacy_fingerprint,
        empirical_null = cfg$empirical_null,
        thresholds = legacy_thresholds_v0_5,
        genes = null_gene_ids
      ),
      algo = "sha256"
    )
  },
  character(1)
)

null_state <- NULL
if (cfg$resume_cached_results && file.exists(null_cache_path)) {
  candidate_state <- readRDS(null_cache_path)
  cache_hash_accepted <- candidate_state$configuration_hash %in%
    c(null_configuration_hash, legacy_null_configuration_hashes)
  if (cache_hash_accepted &&
      identical(dim(candidate_state$jsd), c(n_null_gene, n_null_replicate, 2L))) {
    null_state <- candidate_state
    null_state$configuration_hash <- null_configuration_hash
    log_msg("Loaded compatible empirical-null cache with ",
            sum(candidate_state$completed), " completed source-replicates.")
  }
}

if (is.null(null_state)) {
  null_state <- list(
    configuration_hash = null_configuration_hash,
    gene_ids = null_gene_ids,
    jsd = array(NA_real_, dim = c(n_null_gene, n_null_replicate, 2L),
                dimnames = list(null_gene_ids, NULL, null_sources)),
    max_abs_delta_if = array(NA_real_, dim = c(n_null_gene, n_null_replicate, 2L),
                             dimnames = list(null_gene_ids, NULL, null_sources)),
    tvd = array(NA_real_, dim = c(n_null_gene, n_null_replicate, 2L),
                dimnames = list(null_gene_ids, NULL, null_sources)),
    architecture = array(NA_integer_, dim = c(n_null_gene, n_null_replicate, 2L),
                         dimnames = list(null_gene_ids, NULL, null_sources)),
    evaluable_pairs = array(NA_integer_, dim = c(n_null_gene, n_null_replicate, 2L),
                            dimnames = list(null_gene_ids, NULL, null_sources)),
    completed = matrix(FALSE, nrow = n_null_replicate, ncol = 2L,
                       dimnames = list(NULL, null_sources)),
    replicate_audit = list()
  )
}

completed_since_checkpoint <- 0L
for (source_index in seq_along(null_sources)) {
  source_timepoint <- null_sources[[source_index]]
  for (replicate_index in seq_len(n_null_replicate)) {
    if (isTRUE(null_state$completed[replicate_index, source_index])) next
    result <- run_null_replicate(replicate_index, source_timepoint)
    null_state$jsd[, replicate_index, source_index] <- result$jsd
    null_state$max_abs_delta_if[, replicate_index, source_index] <- result$max_abs_delta_if
    null_state$tvd[, replicate_index, source_index] <- result$total_variation_distance
    null_state$architecture[, replicate_index, source_index] <- match(
      result$architecture, null_architecture_levels
    )
    null_state$evaluable_pairs[, replicate_index, source_index] <- result$n_evaluable_pairs
    null_state$completed[replicate_index, source_index] <- TRUE
    null_state$replicate_audit[[paste(source_timepoint, replicate_index, sep = "_")]] <- tibble(
      source_timepoint = source_timepoint,
      replicate = replicate_index,
      random_seed = deterministic_seed(cfg$seed, "empirical_null", source_timepoint, replicate_index),
      genes = nrow(result),
      quantifiable_genes = sum(result$architecture != "NQ"),
      total_original_FL_across_gene_blocks = sum(result$original_total),
      total_half_A_FL_across_gene_blocks = sum(result$half_a_total),
      total_half_B_FL_across_gene_blocks = sum(result$half_b_total),
      maximum_conservation_error = max(result$conservation_error),
      median_half_A_fraction = median(
        result$half_a_total / pmax(1, result$original_total)
      )
    )
    completed_since_checkpoint <- completed_since_checkpoint + 1L
    if (completed_since_checkpoint >= cfg$empirical_null$checkpoint_every) {
      atomic_save_rds(null_state, null_cache_path)
      completed_since_checkpoint <- 0L
      log_msg(
        "Empirical null: completed ", sum(null_state$completed), " / ",
        length(null_state$completed), " source-replicates."
      )
    }
  }
}
atomic_save_rds(null_state, null_cache_path)

log_msg("Split-null simulations complete; summarizing calibration and strata.")
null_replicate_audit <- bind_rows(null_state$replicate_audit)
write_tsv(null_replicate_audit, file.path(m04$qc, "empirical_null_split_audit.tsv.gz"))

wilson_interval <- function(success, total, confidence = 0.95) {
  if (!is.finite(total) || total <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - confidence) / 2)
  p <- success / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(max(0, center - half), min(1, center + half))
}

summarise_false_counts <- function(false, total) {
  interval <- wilson_interval(false, total)
  tibble(
    n_null_comparisons = total,
    n_false_effect_calls = false,
    false_effect_rate = if (total > 0) false / total else NA_real_,
    wilson_lower_95 = interval[[1]],
    wilson_upper_95 = interval[[2]]
  )
}

quantile_safe <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  quantile(x, probability, names = FALSE, type = cfg$bootstrap$quantile_type)
}

threshold_grid <- tidyr::crossing(
  jsd_threshold = sort(unique(c(as.numeric(cfg$split_null_public$jsd_threshold_grid %||% c(.10,.15,.20,.25,.30)), cfg$thresholds$jsd))),
  delta_if_threshold = sort(unique(c(as.numeric(cfg$split_null_public$max_delta_if_threshold_grid %||% c(.15,.20,.25,.30,.35)), cfg$thresholds$delta_if)))
)

strata_index <- null_gene_strata %>%
  mutate(
    stratum_id = interaction(
      gene_depth_bin, isoform_complexity_bin, drop = TRUE, lex.order = TRUE
    )
  )
stratum_rows <- split(seq_len(nrow(strata_index)), strata_index$stratum_id)
stratum_lookup <- strata_index %>%
  distinct(stratum_id, gene_depth_bin, isoform_complexity_bin)

lookup_null_stratum <- function(current_stratum_id) {
  stratum_metadata <- stratum_lookup %>%
    filter(
      as.character(.data$stratum_id) ==
        as.character(.env$current_stratum_id)
    )
  if (nrow(stratum_metadata) != 1L) {
    stop(
      "Expected exactly one empirical-null stratum for '",
      current_stratum_id, "'; found ", nrow(stratum_metadata), "."
    )
  }
  stratum_metadata
}

source_rate_rows <- vector("list", length(null_sources))
stratum_rate_rows <- vector("list", length(null_sources))
class_rate_rows <- vector("list", length(null_sources))
gene_stability_rows <- vector("list", length(null_sources))
plot_sample_rows <- vector("list", length(null_sources))
convergence_rows <- vector("list", length(null_sources))
threshold_false_counts <- numeric(nrow(threshold_grid))
threshold_valid_total <- 0
overall_false_total <- 0
overall_valid_total <- 0

for (source_index in seq_along(null_sources)) {
  source_timepoint <- null_sources[[source_index]]
  jsd_matrix <- matrix(null_state$jsd[, , source_index], nrow=n_null_gene, ncol=n_null_replicate)
  delta_matrix <- matrix(null_state$max_abs_delta_if[, , source_index], nrow=n_null_gene, ncol=n_null_replicate)
  tvd_matrix <- matrix(null_state$tvd[, , source_index], nrow=n_null_gene, ncol=n_null_replicate)
  architecture_matrix <- matrix(null_state$architecture[, , source_index], nrow=n_null_gene, ncol=n_null_replicate)
  valid_matrix <- !is.na(architecture_matrix) &
    architecture_matrix != match("NQ", null_architecture_levels) &
    is.finite(jsd_matrix) & is.finite(delta_matrix)
  effect_matrix <- valid_matrix &
    (jsd_matrix >= cfg$thresholds$jsd |
       delta_matrix >= cfg$thresholds$delta_if)

  valid_by_replicate <- colSums(valid_matrix)
  false_by_replicate <- colSums(effect_matrix)
  convergence_rows[[source_index]] <- tibble(
    source_timepoint = source_timepoint,
    replicate = seq_len(n_null_replicate),
    valid_comparisons_in_replicate = as.numeric(valid_by_replicate),
    false_calls_in_replicate = as.numeric(false_by_replicate),
    cumulative_valid_comparisons = cumsum(valid_comparisons_in_replicate),
    cumulative_false_calls = cumsum(false_calls_in_replicate),
    cumulative_false_effect_rate = cumulative_false_calls /
      pmax(1, cumulative_valid_comparisons)
  )

  valid_n <- sum(valid_matrix)
  false_n <- sum(effect_matrix)
  overall_valid_total <- overall_valid_total + valid_n
  overall_false_total <- overall_false_total + false_n
  threshold_valid_total <- threshold_valid_total + valid_n
  source_rate_rows[[source_index]] <- summarise_false_counts(false_n, valid_n) %>%
    mutate(source_timepoint = source_timepoint, .before = 1)

  stratum_rate_rows[[source_index]] <- imap_dfr(
    stratum_rows,
    function(rows, stratum_id) {
      stratum_valid <- sum(valid_matrix[rows, , drop = FALSE])
      stratum_false <- sum(effect_matrix[rows, , drop = FALSE])
      stratum_metadata <- lookup_null_stratum(stratum_id)
      bind_cols(
        tibble(source_timepoint = source_timepoint),
        stratum_metadata %>% select(gene_depth_bin, isoform_complexity_bin),
        summarise_false_counts(stratum_false, stratum_valid)
      )
    }
  )

  class_rate_rows[[source_index]] <- imap_dfr(
    stratum_rows,
    function(rows, stratum_id) {
      valid_values <- architecture_matrix[rows, , drop = FALSE][
        valid_matrix[rows, , drop = FALSE]
      ]
      class_counts <- tabulate(valid_values, nbins = length(null_architecture_levels))
      total <- sum(class_counts)
      stratum_metadata <- lookup_null_stratum(stratum_id)
      tidyr::crossing(
        tibble(source_timepoint = source_timepoint),
        stratum_metadata %>% select(gene_depth_bin, isoform_complexity_bin),
        tibble(
          architecture = null_architecture_levels,
          n_class_calls = as.integer(class_counts),
          n_null_comparisons = as.integer(total),
          class_fraction = if (total > 0) class_counts / total else NA_real_
        )
      )
    }
  ) %>%
    filter(architecture != "NQ")

  valid_per_gene <- rowSums(valid_matrix)
  jsd_for_quantile <- jsd_matrix
  delta_for_quantile <- delta_matrix
  jsd_for_quantile[!valid_matrix] <- NA_real_
  delta_for_quantile[!valid_matrix] <- NA_real_
  modal_architecture <- apply(architecture_matrix, 1L, function(values) {
    values <- values[values %in% seq_len(5L)]
    if (!length(values)) return(NA_character_)
    null_architecture_levels[[which.max(tabulate(values, nbins = 5L))]]
  })
  modal_architecture_probability <- apply(architecture_matrix, 1L, function(values) {
    effect_values <- values[values %in% seq_len(5L)]
    valid_total <- sum(!is.na(values) &
      values != match("NQ", null_architecture_levels))
    if (!length(effect_values) || valid_total == 0L) return(0)
    max(tabulate(effect_values, nbins = 5L)) / valid_total
  })
  gene_stability_rows[[source_index]] <- tibble(
    source_timepoint = source_timepoint,
    gene_id = null_gene_ids,
    valid_fraction = valid_per_gene / n_null_replicate,
    false_effect_rate = ifelse(
      valid_per_gene > 0, rowSums(effect_matrix) / valid_per_gene, NA_real_
    ),
    jsd_q975 = apply(jsd_for_quantile, 1L, quantile_safe, probability = 0.975),
    max_abs_delta_if_q975 = apply(
      delta_for_quantile, 1L, quantile_safe, probability = 0.975
    ),
    modal_effect_architecture = modal_architecture,
    modal_effect_architecture_probability = modal_architecture_probability
  )

  for (threshold_index in seq_len(nrow(threshold_grid))) {
    threshold_false_counts[[threshold_index]] <-
      threshold_false_counts[[threshold_index]] + sum(
        valid_matrix &
          (jsd_matrix >= threshold_grid$jsd_threshold[[threshold_index]] |
             delta_matrix >= threshold_grid$delta_if_threshold[[threshold_index]])
      )
  }

  valid_position <- which(valid_matrix)
  sample_n <- min(100000L, length(valid_position))
  selected_position <- with_seed(
    deterministic_seed(cfg$seed, "null_plot_sample", source_timepoint),
    if (sample_n > 0L) {
      valid_position[sample.int(length(valid_position), sample_n)]
    } else integer()
  )
  gene_index <- if (length(selected_position)) {
    ((selected_position - 1L) %% n_null_gene) + 1L
  } else integer()
  replicate_index <- if (length(selected_position)) {
    ((selected_position - 1L) %/% n_null_gene) + 1L
  } else integer()
  plot_sample_rows[[source_index]] <- tibble(
    gene_id = null_gene_ids[gene_index],
    replicate = replicate_index,
    source_timepoint = source_timepoint,
    jsd = jsd_matrix[selected_position],
    max_abs_delta_if = delta_matrix[selected_position],
    total_variation_distance = tvd_matrix[selected_position],
    architecture = null_architecture_levels[architecture_matrix[selected_position]]
  ) %>% left_join(null_gene_strata, by = "gene_id")

  rm(
    jsd_matrix, delta_matrix, tvd_matrix, architecture_matrix,
    valid_matrix, effect_matrix, jsd_for_quantile, delta_for_quantile
  )
  invisible(gc())
}

null_overall_rate <- summarise_false_counts(
  overall_false_total, overall_valid_total
) %>% mutate(stratum = "all_sources_all_strata")
null_source_rate <- bind_rows(source_rate_rows)
null_stratum_rate <- bind_rows(stratum_rate_rows) %>%
  mutate(
    included_in_threshold_gate = n_null_comparisons >=
      cfg$empirical_null$minimum_stratum_n
  )
null_class_rate <- bind_rows(class_rate_rows) %>%
  rowwise() %>%
  mutate(
    wilson_lower_95 = wilson_interval(n_class_calls, n_null_comparisons)[[1]],
    wilson_upper_95 = wilson_interval(n_class_calls, n_null_comparisons)[[2]]
  ) %>%
  ungroup()

expected_null_stratum_rows <- length(null_sources) * nrow(stratum_lookup)
expected_null_class_rows <- expected_null_stratum_rows *
  length(setdiff(null_architecture_levels, "NQ"))
null_stratum_mapping_qc <- tibble(
  check = c(
    "one false-rate row per source and stratum",
    "unique false-rate source-stratum keys",
    "one class-rate row per source, stratum, and non-NQ class",
    "unique class-rate source-stratum-class keys"
  ),
  observed = c(
    nrow(null_stratum_rate),
    null_stratum_rate %>%
      distinct(source_timepoint, gene_depth_bin, isoform_complexity_bin) %>%
      nrow(),
    nrow(null_class_rate),
    null_class_rate %>%
      distinct(
        source_timepoint, gene_depth_bin, isoform_complexity_bin,
        architecture
      ) %>%
      nrow()
  ),
  expected = c(
    expected_null_stratum_rows,
    expected_null_stratum_rows,
    expected_null_class_rows,
    expected_null_class_rows
  )
) %>%
  mutate(pass = observed == expected)
null_stratum_mapping_pass <- nrow(null_stratum_mapping_qc) == 4L &&
  isTRUE(all(null_stratum_mapping_qc$pass))

null_gene_stability <- bind_rows(gene_stability_rows) %>%
  mutate(
    null_stability_class = case_when(
      valid_fraction < cfg$thresholds$stability_probability ~ "NQ",
      is.finite(jsd_q975) & jsd_q975 < cfg$thresholds$jsd &
        is.finite(max_abs_delta_if_q975) &
        max_abs_delta_if_q975 < cfg$thresholds$delta_if ~ "F",
      TRUE ~ "U"
    ),
    null_full_framework_class = case_when(
      valid_fraction < cfg$thresholds$stability_probability ~ "NQ",
      null_stability_class == "F" ~ "F",
      false_effect_rate >= cfg$thresholds$stability_probability &
        modal_effect_architecture_probability >= cfg$thresholds$stability_probability &
        modal_effect_architecture %in% c("A", "B", "C", "D", "E") ~
          modal_effect_architecture,
      TRUE ~ "U"
    )
  )
null_full_framework_summary <- null_gene_stability %>%
  left_join(
    null_gene_strata %>%
      select(gene_id, gene_depth_bin, isoform_complexity_bin),
    by = "gene_id"
  ) %>%
  count(
    source_timepoint, gene_depth_bin, isoform_complexity_bin,
    null_full_framework_class, name = "n_genes"
  ) %>%
  group_by(source_timepoint, gene_depth_bin, isoform_complexity_bin) %>%
  mutate(fraction = n_genes / sum(n_genes)) %>%
  ungroup()
null_threshold_sensitivity <- threshold_grid %>%
  mutate(
    n_null_comparisons = threshold_valid_total,
    n_effect_gate_calls = threshold_false_counts,
    effect_gate_rate = n_effect_gate_calls / n_null_comparisons
  ) %>%
  rowwise() %>%
  mutate(
    wilson_lower_95 = wilson_interval(
      n_effect_gate_calls, n_null_comparisons
    )[[1]],
    wilson_upper_95 = wilson_interval(
      n_effect_gate_calls, n_null_comparisons
    )[[2]]
  ) %>%
  ungroup()
null_plot_sample <- bind_rows(plot_sample_rows)
null_convergence <- bind_rows(convergence_rows)
null_convergence_qc <- null_convergence %>%
  group_by(source_timepoint) %>%
  mutate(final_rate = last(cumulative_false_effect_rate)) %>%
  filter(replicate >= ceiling(0.75 * cfg$empirical_null$replicates)) %>%
  summarise(
    requested_replicates = cfg$empirical_null$replicates,
    completed_replicates = max(replicate),
    final_false_effect_rate = first(final_rate),
    maximum_final_quarter_deviation = max(
      abs(cumulative_false_effect_rate - first(final_rate)), na.rm = TRUE
    ),
    convergence_tolerance = cfg$empirical_null$convergence_tolerance,
    pass = completed_replicates == requested_replicates &
      is.finite(final_false_effect_rate) &
      maximum_final_quarter_deviation <= convergence_tolerance,
    .groups = "drop"
  )

overall_gate_pass <- is.finite(null_overall_rate$false_effect_rate[[1]]) &&
  null_overall_rate$false_effect_rate[[1]] <= cfg$empirical_null$overall_target
adequate_null_strata_present <- any(null_stratum_rate$included_in_threshold_gate)
stratum_gate_pass <- adequate_null_strata_present && all(
  is.finite(null_stratum_rate$false_effect_rate[
    null_stratum_rate$included_in_threshold_gate
  ]) &
    null_stratum_rate$false_effect_rate[null_stratum_rate$included_in_threshold_gate] <=
      cfg$empirical_null$stratum_target
)
null_precision_pass <- all(null_state$completed) &&
  nrow(null_replicate_audit) == 2L * cfg$empirical_null$replicates &&
  all(null_replicate_audit$maximum_conservation_error == 0)
null_convergence_pass <- nrow(null_convergence_qc) == length(null_sources) &&
  all(null_convergence_qc$pass)
threshold_core_gate_pass <- overall_gate_pass && null_precision_pass &&
  null_convergence_pass && null_stratum_mapping_pass
null_stratum_stress_warning <- !stratum_gate_pass

threshold_gate <- tibble(
  check = c(
    "overall false effect rate", "adequately populated stratum false effect rates",
    "requested split-null replicates and exact count conservation",
    "cumulative false-rate convergence",
    "one-to-one null stratum and class mapping"
  ),
  observed = c(
    null_overall_rate$false_effect_rate[[1]],
    if (any(null_stratum_rate$included_in_threshold_gate)) {
      max(null_stratum_rate$false_effect_rate[null_stratum_rate$included_in_threshold_gate],
          na.rm = TRUE)
    } else NA_real_,
    mean(null_state$completed),
    if (nrow(null_convergence_qc)) {
      max(null_convergence_qc$maximum_final_quarter_deviation, na.rm = TRUE)
    } else NA_real_,
    as.numeric(null_stratum_mapping_pass)
  ),
  target = c(
    cfg$empirical_null$overall_target,
    cfg$empirical_null$stratum_target,
    1,
    cfg$empirical_null$convergence_tolerance,
    1
  ),
  pass = c(
    overall_gate_pass, stratum_gate_pass, null_precision_pass,
    null_convergence_pass, null_stratum_mapping_pass
  ),
  blocking = c(
    TRUE, FALSE, TRUE, TRUE, TRUE
  ),
  interpretation = c(
    "Primary calibration of the combined effect gate",
    paste(
      "Conservative half-depth stress screen; an excess triggers an observed-gene",
      "edge-case audit and blocks only if it reaches the locked primary downstream set"
    ),
    "Required simulation integrity check",
    "Required simulation stability check",
    "Required reporting integrity check"
  )
)

write_tsv(null_overall_rate, file.path(m04$tables, "empirical_null_overall_false_rate.tsv"))
write_tsv(null_source_rate, file.path(m04$tables, "empirical_null_false_rate_by_source.tsv"))
write_tsv(null_stratum_rate, file.path(m04$tables, "empirical_null_false_rate_by_stratum.tsv"))
write_tsv(null_class_rate, file.path(m04$tables, "empirical_null_class_rates.tsv"))
write_tsv(null_gene_stability, file.path(m04$tables, "empirical_null_gene_stability.tsv.gz"))
write_tsv(
  null_full_framework_summary,
  file.path(m04$tables, "empirical_null_full_framework_classes.tsv")
)
write_tsv(null_threshold_sensitivity, file.path(m04$tables, "empirical_null_threshold_sensitivity.tsv"))
write_tsv(null_convergence, file.path(m04$tables, "empirical_null_convergence.tsv.gz"))
write_tsv(threshold_gate, file.path(m04$qc, "effect_threshold_gate.tsv"))
write_tsv(null_convergence_qc, file.path(m04$qc, "empirical_null_convergence_qc.tsv"))
write_tsv(
  null_stratum_mapping_qc,
  file.path(m04$qc, "empirical_null_stratum_mapping_qc.tsv")
)

write_tsv(null_plot_sample, file.path(m04$tables, "empirical_null_plot_sample.tsv.gz"))

if (!exists("figure_records")) figure_records <- list()

if (plotting_enabled) {
compact_null_stratum <- function(depth, complexity) {
  depth <- as.character(depth)
  complexity <- as.character(complexity)
  depth_label <- case_when(
    depth == ">=10000" ~ ">=10k",
    depth == "2000-9999" ~ "2k-9.9k",
    depth == "500-1999" ~ "500-1.9k",
    TRUE ~ depth
  )
  complexity_label <- if_else(
    complexity == ">20", "K>20", paste0("K=", complexity)
  )
  paste(depth_label, complexity_label, sep = " | ")
}

null_stratum_levels <- stratum_lookup %>%
  arrange(gene_depth_bin, isoform_complexity_bin) %>%
  transmute(
    stratum = compact_null_stratum(gene_depth_bin, isoform_complexity_bin)
  ) %>%
  pull(stratum) %>%
  unique()

p_null_metrics <- null_plot_sample %>%
  select(source_timepoint, jsd, max_abs_delta_if) %>%
  pivot_longer(c(jsd, max_abs_delta_if), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, jsd = "JSD (base 2)",
                    max_abs_delta_if = "Max |delta IF|"),
    threshold = if_else(metric == "JSD (base 2)",
                        cfg$thresholds$jsd, cfg$thresholds$delta_if)
  ) %>%
  ggplot(aes(x = value, fill = source_timepoint)) +
  geom_histogram(bins = 70, alpha = 0.55, position = "identity") +
  geom_vline(aes(xintercept = threshold), linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(
    values = unlist(public_cfg$reporting$palettes$timepoints %||% c(primary = "#4C78A8", recurrent = "#E45756")),
    labels = c(primary = "Primary", recurrent = "Recurrent")
  ) +
  labs(
    title = "Within-library empirical null distributions",
    subtitle = "Complementary half-library comparisons; dashed lines are prespecified effect thresholds",
    x = "Null effect metric", y = "Pseudo-comparisons", fill = "Source library"
  ) +
  paper_theme()

p_null_rates <- null_stratum_rate %>%
  filter(included_in_threshold_gate) %>%
  mutate(
    stratum = factor(
      compact_null_stratum(gene_depth_bin, isoform_complexity_bin),
      levels = null_stratum_levels
    )
  ) %>%
  ggplot(aes(x = false_effect_rate, y = stratum, color = source_timepoint)) +
  geom_vline(xintercept = cfg$empirical_null$stratum_target,
             linetype = "dashed", color = "#555555") +
  geom_errorbarh(aes(xmin = wilson_lower_95, xmax = wilson_upper_95),
                 height = 0.18, linewidth = 0.7, position = position_dodge(width = 0.45)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.45)) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(
    values = unlist(public_cfg$reporting$palettes$timepoints %||% c(primary = "#4C78A8", recurrent = "#E45756")),
    labels = c(primary = "Primary", recurrent = "Recurrent")
  ) +
  labs(
    title = "False remodeling calls under the empirical null",
    subtitle = "Points and 95% Wilson intervals; dashed line is the prespecified stratum target",
    x = "False effect-gate rate", y = "Pooled FL reads | isoforms",
    color = "Source library"
  ) +
  paper_theme()

p_threshold_grid <- null_threshold_sensitivity %>%
  mutate(
    jsd_threshold = factor(jsd_threshold),
    delta_if_threshold = factor(delta_if_threshold)
  ) %>%
  ggplot(aes(x = delta_if_threshold, y = jsd_threshold, fill = effect_gate_rate)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = percent(effect_gate_rate, accuracy = 0.1)), size = 5.0) +
  scale_fill_viridis_c(labels = percent_format(accuracy = 1), option = "C") +
  labs(
    title = "Empirical-null sensitivity of the combined effect gate",
    subtitle = "A call passes when JSD or maximum |delta IF| reaches its threshold",
    x = "Max |delta IF| threshold", y = "JSD threshold",
    fill = "Null call rate"
  ) +
  paper_theme()

null_class_plot_data <- null_class_rate %>%
  filter(architecture %in% c("A", "B", "C", "D", "E")) %>%
  left_join(
    null_stratum_rate %>%
      select(
        source_timepoint, gene_depth_bin, isoform_complexity_bin,
        included_in_threshold_gate
      ),
    by = c("source_timepoint", "gene_depth_bin", "isoform_complexity_bin")
  ) %>%
  filter(included_in_threshold_gate) %>%
  mutate(
    architecture = factor(architecture, levels = c("A", "B", "C", "D", "E")),
    stratum = factor(
      compact_null_stratum(gene_depth_bin, isoform_complexity_bin),
      levels = null_stratum_levels
    )
  )
p_null_classes <- ggplot(
  null_class_plot_data,
  aes(x = architecture, y = stratum, fill = class_fraction)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_wrap(
    ~source_timepoint, ncol = 2,
    labeller = as_labeller(c(primary = "Primary", recurrent = "Recurrent"))
  ) +
  scale_fill_viridis_c(labels = percent_format(accuracy = 0.1), option = "C") +
  labs(
    title = "Architecture-specific calls under the empirical null",
    subtitle = "Rates are stratified by pooled expression depth and retained isoform complexity",
    x = "False class", y = "Pooled FL reads | isoforms",
    fill = "Null call rate"
  ) +
  paper_theme()

p_null_convergence <- ggplot(
  null_convergence,
  aes(x = replicate, y = cumulative_false_effect_rate, color = source_timepoint)
) +
  geom_hline(
    yintercept = cfg$empirical_null$overall_target,
    linetype = "dashed", color = "#555555", linewidth = 0.8
  ) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(
    values = unlist(public_cfg$reporting$palettes$timepoints %||% c(primary = "#4C78A8", recurrent = "#E45756")),
    labels = c(primary = "Primary", recurrent = "Recurrent")
  ) +
  labs(
    title = "Convergence of the empirical-null false-call rate",
    subtitle = "Cumulative rates across complementary half-library comparisons; dashed line is the overall target",
    x = "Split-null replicate", y = "Cumulative false effect-gate rate",
    color = "Source library"
  ) +
  paper_theme()

figure_records[["empirical_null_metric_distributions"]] <- save_publication_plot(
  p_null_metrics, "empirical_null_metric_distributions", 11, 6.5
)
figure_records[["empirical_null_false_call_rates"]] <- save_publication_plot(
  p_null_rates, "empirical_null_false_call_rates", 11.5,
  max(9, 0.42 * n_distinct(null_stratum_rate$gene_depth_bin,
                            null_stratum_rate$isoform_complexity_bin))
)
figure_records[["effect_threshold_sensitivity"]] <- save_publication_plot(
  p_threshold_grid, "effect_threshold_sensitivity", 9.5, 7.5
)
figure_records[["empirical_null_class_composition"]] <- save_publication_plot(
  p_null_classes, "empirical_null_class_composition", 13,
  max(9, 0.40 * n_distinct(null_class_plot_data$stratum))
)
figure_records[["empirical_null_convergence"]] <- save_publication_plot(
  p_null_convergence, "empirical_null_convergence", 10.5, 7
)


} # optional split-null plot construction

log_msg("Split-null calibration summaries complete.")
if (!threshold_core_gate_pass && cfg$empirical_null$stop_on_failure) {
  stop(
    "The JSD/max-|delta IF| gate failed a blocking empirical-null check. ",
    "No definitive A-F/U/NQ table was written. Inspect ",
    safe_norm(file.path(m04$qc, "effect_threshold_gate.tsv")),
    " and the empirical-null figures; revise thresholds only through a documented analysis decision."
  )
}
if (!threshold_core_gate_pass && !cfg$empirical_null$stop_on_failure) {
  log_msg(
    "The effect gate did not pass every blocking empirical-null criterion. ",
    "The workflow will finish in audit-only mode so all diagnostic outputs can be reviewed; ",
    "no M04_complete checkpoint can be issued.",
    level = "WARN"
  )
}
if (null_stratum_stress_warning) {
  log_msg(
    "At least one adequately populated half-depth split-null stratum exceeded the ",
    "registered stress-screen target. Observed genes in that stratum will be ",
    "cross-audited against effect stability, patient support, and cohort P values.",
    level = "WARN"
  )
}

# ---- m04-bootstrap-run (source-derived computation) ----
log_msg("Remodeling step: bootstrap-run.")
bootstrap_qc_gene_id <- primary_observed %>%
  filter(quantifiable) %>%
  arrange(desc(n_isoforms_retained), gene_id) %>%
  slice_head(n = 1L) %>%
  pull(gene_id)
if (!length(bootstrap_qc_gene_id)) {
  stop("No quantifiable gene was available for the vectorized-bootstrap equivalence check.")
}

bootstrap_vectorization_qc <- with_seed(
  deterministic_seed(cfg$seed, "bootstrap_vectorization_qc", bootstrap_qc_gene_id),
  {
    qc_object <- gene_objects[[bootstrap_qc_gene_id]]
    qc_view <- cohort_gene_view(qc_object, patient_ids)
    qc_b <- min(250L, cfg$bootstrap$replicates)
    qc_draws <- matrix(
      sample.int(nrow(qc_view$pair_table), qc_b * nrow(qc_view$pair_table),
                 replace = TRUE),
      nrow = qc_b, ncol = nrow(qc_view$pair_table)
    )
    map_dfr(c("primary", "recurrent"), function(side) {
      sample_composition <- qc_view[[paste0(side, "_composition")]]
      group_boot <- vapply(seq_len(ncol(sample_composition)), function(index) {
        row_means_from_draws(sample_composition[, index], qc_draws)
      }, numeric(qc_b)) %>%
        matrix(nrow = qc_b, ncol = ncol(sample_composition))
      vectorized <- bootstrap_dominance_matrix(
        group_boot, sample_composition, qc_draws, qc_object$transcript_ids
      )
      reference <- lapply(seq_len(qc_b), function(index) {
        dominance_summary(
          group_boot[index, ],
          sample_composition[qc_draws[index, ], , drop = FALSE],
          qc_object$transcript_ids
        )
      })
      reference_has <- vapply(reference, `[[`, logical(1), "has_dominant")
      reference_id <- vapply(reference, function(value) {
        value$dominant_id %||% NA_character_
      }, character(1))
      reference_support <- vapply(reference, `[[`, integer(1), "support_n")
      id_match <- ifelse(
        is.na(vectorized$dominant_id) & is.na(reference_id), TRUE,
        ifelse(
          is.na(vectorized$dominant_id) | is.na(reference_id), FALSE,
          vectorized$dominant_id == reference_id
        )
      )
      tibble(
        gene_id = bootstrap_qc_gene_id,
        side = side,
        checked_replicates = qc_b,
        has_dominant_disagreements = sum(vectorized$has_dominant != reference_has),
        dominant_id_disagreements = sum(!id_match, na.rm = TRUE),
        support_count_disagreements = sum(vectorized$support_n != reference_support),
        maximum_margin_error = max(abs(
          vectorized$dominance_margin -
            vapply(reference, `[[`, numeric(1), "dominance_margin")
        )),
        pass = has_dominant_disagreements == 0L &
          dominant_id_disagreements == 0L &
          support_count_disagreements == 0L & maximum_margin_error < 1e-12
      )
    })
  }
)
write_tsv(
  bootstrap_vectorization_qc,
  file.path(m04$qc, "bootstrap_vectorization_equivalence.tsv")
)

bootstrap_results <- imap(analysis_cohorts, function(included_patients, cohort_name) {
  run_bootstrap_cohort(cohort_name, included_patients)
})

bootstrap_summary <- bind_rows(lapply(bootstrap_results, `[[`, "summary")) %>%
  mutate(
    across(
      c(
        effect_pass_bootstrap_probability,
        observed_class_bootstrap_probability,
        observed_driver_bootstrap_probability,
        primary_dominance_state_probability,
        recurrent_dominance_state_probability
      ),
      ~sqrt(.x * (1 - .x) / pmax(1, bootstrap_valid_replicates)),
      .names = "{.col}_mcse"
    )
  )
bootstrap_class_probability <- bind_rows(lapply(bootstrap_results, `[[`, "class_probability"))

remodeling_classes_interim <- observed_gene %>%
  left_join(bootstrap_summary, by = c("cohort", "gene_id")) %>%
  mutate(
    bootstrap_valid = quantifiable &
      bootstrap_valid_fraction >= cfg$bootstrap$minimum_valid_fraction,
    is_observed_remodeled_AE = observed_architecture %in% c("A", "B", "C", "D", "E"),
    effect_stable = case_when(
      is_observed_remodeled_AE & bootstrap_valid ~
        effect_pass_bootstrap_probability >= cfg$thresholds$stability_probability,
      TRUE ~ NA
    ),
    architecture_class_stable = case_when(
      is_observed_remodeled_AE & bootstrap_valid ~
        observed_class_bootstrap_probability >= cfg$thresholds$stability_probability,
      TRUE ~ NA
    ),
    driver_stable = case_when(
      is_observed_remodeled_AE & bootstrap_valid ~
        observed_driver_bootstrap_probability >= cfg$thresholds$stability_probability,
      TRUE ~ NA
    ),
    primary_dominance_stable = case_when(
      is_observed_remodeled_AE & bootstrap_valid ~
        primary_dominance_state_probability >= cfg$thresholds$stability_probability,
      TRUE ~ NA
    ),
    recurrent_dominance_stable = case_when(
      is_observed_remodeled_AE & bootstrap_valid ~
        recurrent_dominance_state_probability >= cfg$thresholds$stability_probability,
      TRUE ~ NA
    ),
    supported_stability = quantifiable & observed_architecture == "Below_effect" &
      bootstrap_valid_fraction >= cfg$bootstrap$minimum_valid_fraction &
      is.finite(jsd_bootstrap_upper) & jsd_bootstrap_upper < cfg$thresholds$jsd &
      is.finite(max_abs_delta_if_bootstrap_upper) &
      max_abs_delta_if_bootstrap_upper < cfg$thresholds$delta_if,
    architecture_class = case_when(
      !quantifiable ~ "NQ",
      observed_architecture %in% c("A", "B", "C", "D", "E") ~ observed_architecture,
      supported_stability ~ "F",
      TRUE ~ "U"
    ),
    architecture_class = factor(architecture_class, levels = cfg$class_levels),
    architecture_description = recode(
      as.character(architecture_class),
      A = "dominant_isoform_switch",
      B = "same_dominant_verified_secondary_remodeling",
      C = "clear_primary_dominance_lost",
      D = "clear_recurrent_dominance_gained",
      E = "remodeling_without_clear_dominance",
      F = "bootstrap_supported_stability",
      U = "quantifiable_inconclusive",
      NQ = "non_quantifiable"
    ),
    is_remodeled_AE = architecture_class %in% c("A", "B", "C", "D", "E"),
    is_supported_stable_F = architecture_class == "F",
    is_uncertain_U = architecture_class == "U",
    is_quantifiable = architecture_class != "NQ",
    direction_support_ge1 = driver_patient_support_n >= 1L,
    direction_support_ge2 = driver_patient_support_n >= 2L,
    direction_support_ge3 = driver_patient_support_n >= 3L,
    direction_support_ge4 = driver_patient_support_n >= 4L,
    direction_support_ge5 = driver_patient_support_n >= 5L,
    direction_support_ge6 = driver_patient_support_n >= 6L,
    direction_support_primary = driver_patient_support_n >=
      cfg$thresholds$primary_patient_support,
    is_high_confidence_remodeled = is_remodeled_AE &
      effect_stable %in% TRUE & direction_support_primary,
    is_high_confidence_class = is_high_confidence_remodeled &
      architecture_class_stable %in% TRUE & driver_stable %in% TRUE &
      primary_dominance_stable %in% TRUE & recurrent_dominance_stable %in% TRUE,
    effect_class_stability = case_when(
      is_remodeled_AE & effect_stable %in% TRUE & architecture_class_stable %in% TRUE ~
        "effect_and_class_stable",
      is_remodeled_AE & effect_stable %in% TRUE ~ "effect_stable_class_uncertain",
      is_remodeled_AE & effect_stable %in% FALSE & architecture_class_stable %in% TRUE ~
        "effect_uncertain_class_stable",
      is_remodeled_AE & effect_stable %in% FALSE ~ "effect_and_class_uncertain",
      is_remodeled_AE ~ "bootstrap_not_evaluable",
      architecture_class == "F" ~ "supported_stability",
      architecture_class == "U" ~ "inconclusive",
      TRUE ~ "non_quantifiable"
    )
  )

bootstrap_class_sum_qc <- bootstrap_class_probability %>%
  group_by(cohort, gene_id) %>%
  summarise(probability_sum = sum(probability), .groups = "drop") %>%
  mutate(error = abs(probability_sum - 1))

bootstrap_qc <- remodeling_classes_interim %>%
  group_by(cohort) %>%
  summarise(
    genes = n(),
    quantifiable_genes = sum(quantifiable),
    genes_with_bootstrap = sum(is.finite(bootstrap_valid_fraction)),
    minimum_valid_fraction = min(bootstrap_valid_fraction, na.rm = TRUE),
    maximum_reported_probability_mcse = max(
      effect_pass_bootstrap_probability_mcse,
      observed_class_bootstrap_probability_mcse,
      observed_driver_bootstrap_probability_mcse,
      primary_dominance_state_probability_mcse,
      recurrent_dominance_state_probability_mcse,
      na.rm = TRUE
    ),
    n_probability_sum_errors = sum(
      bootstrap_class_sum_qc$error[match(
        paste(cohort, gene_id),
        paste(bootstrap_class_sum_qc$cohort, bootstrap_class_sum_qc$gene_id)
      )] > 1e-10,
      na.rm = TRUE
    ),
    n_invalid_intervals = sum(
      jsd_bootstrap_lower > jsd_bootstrap_upper |
        max_abs_delta_if_bootstrap_lower > max_abs_delta_if_bootstrap_upper,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

patient_support_grid <- imap_dfr(
  analysis_cohorts,
  function(included_patients, cohort_name) {
    tibble(
      cohort = cohort_name,
      support_threshold = seq_len(length(included_patients))
    )
  }
)

patient_support_sensitivity <- patient_support_grid %>%
  pmap_dfr(function(cohort, support_threshold) {
    data <- remodeling_classes_interim %>%
      filter(.data$cohort == .env$cohort)

    combined <- data %>%
      summarise(
        cohort = .env$cohort,
        support_threshold = as.integer(.env$support_threshold),
        architecture_class = "A-E combined",
        n_observed = sum(is_remodeled_AE),
        n_with_direction_support = sum(
          is_remodeled_AE &
            driver_patient_support_n >= .env$support_threshold,
          na.rm = TRUE
        ),
        n_effect_stable_with_support = sum(
          is_remodeled_AE & effect_stable %in% TRUE &
            driver_patient_support_n >= .env$support_threshold,
          na.rm = TRUE
        )
      )

    by_class <- data %>%
      filter(is_remodeled_AE) %>%
      mutate(architecture_class = as.character(architecture_class)) %>%
      group_by(architecture_class) %>%
      summarise(
        n_observed = n(),
        n_with_direction_support = sum(
          driver_patient_support_n >= .env$support_threshold,
          na.rm = TRUE
        ),
        n_effect_stable_with_support = sum(
          effect_stable %in% TRUE &
            driver_patient_support_n >= .env$support_threshold,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      complete(
        architecture_class = c("A", "B", "C", "D", "E"),
        fill = list(
          n_observed = 0L,
          n_with_direction_support = 0L,
          n_effect_stable_with_support = 0L
        )
      ) %>%
      mutate(
        cohort = .env$cohort,
        support_threshold = as.integer(.env$support_threshold),
        .before = 1
      )

    bind_rows(combined, by_class)
  })

patient_support_exact <- imap_dfr(
  analysis_cohorts,
  function(included_patients, cohort_name) {
    cohort_data <- remodeling_classes_interim %>%
      filter(.data$cohort == .env$cohort_name, is_remodeled_AE) %>%
      mutate(architecture_class = as.character(architecture_class))
    by_class <- cohort_data %>%
      count(
        architecture_class,
        exact_direction_support_n = driver_patient_support_n,
        name = "n_observed_AE"
      ) %>%
      complete(
        architecture_class = c("A", "B", "C", "D", "E"),
        exact_direction_support_n = 0:length(included_patients),
        fill = list(n_observed_AE = 0L)
      )
    combined <- cohort_data %>%
      count(
        exact_direction_support_n = driver_patient_support_n,
        name = "n_observed_AE"
      ) %>%
      complete(
        exact_direction_support_n = 0:length(included_patients),
        fill = list(n_observed_AE = 0L)
      ) %>%
      mutate(architecture_class = "A-E combined", .before = 1L)
    bind_rows(combined, by_class) %>%
      mutate(
        cohort = cohort_name,
        available_patient_pairs = length(included_patients),
        .before = 1L
      )
  }
)

patient_support_group_qc <- patient_support_sensitivity %>%
  group_by(cohort, support_threshold) %>%
  summarise(
    n_rows = n(),
    combined_count = first(
      n_effect_stable_with_support[architecture_class == "A-E combined"],
      default = NA_integer_
    ),
    class_sum = sum(
      n_effect_stable_with_support[architecture_class %in% c("A", "B", "C", "D", "E")]
    ),
    .groups = "drop"
  )
expected_support_groups <- nrow(patient_support_grid)
primary_support_expected <- remodeling_classes_interim %>%
  filter(cohort == primary_cohort_name) %>%
  summarise(n = sum(is_high_confidence_remodeled, na.rm = TRUE)) %>%
  pull(n)
primary_support_observed <- patient_support_sensitivity %>%
  filter(
    cohort == primary_cohort_name,
    support_threshold == cfg$thresholds$primary_patient_support,
    architecture_class == "A-E combined"
  ) %>%
  pull(n_effect_stable_with_support)

patient_support_qc <- tibble(
  check = c(
    "all cohort-threshold groups are present",
    "support thresholds span one through every available pair",
    "one combined plus five A-E rows per group",
    "combined counts equal summed A-E counts",
    "exact-support distributions span zero through every available pair",
    "primary support count matches authoritative high-confidence flag"
  ),
  observed = c(
    nrow(patient_support_group_qc),
    sum(imap_lgl(analysis_cohorts, function(included_patients, cohort_name) {
      identical(
        patient_support_grid$support_threshold[
          patient_support_grid$cohort == cohort_name
        ],
        seq_len(length(included_patients))
      )
    })),
    sum(patient_support_group_qc$n_rows == 6L),
    sum(patient_support_group_qc$combined_count == patient_support_group_qc$class_sum),
    sum(imap_lgl(analysis_cohorts, function(included_patients, cohort_name) {
      identical(
        sort(unique(patient_support_exact$exact_direction_support_n[
          patient_support_exact$cohort == cohort_name
        ])),
        0:length(included_patients)
      )
    })),
    if (length(primary_support_observed) == 1L) primary_support_observed else NA_real_
  ),
  expected = c(
    expected_support_groups,
    length(analysis_cohorts),
    expected_support_groups,
    expected_support_groups,
    length(analysis_cohorts),
    primary_support_expected
  )
) %>%
  mutate(pass = observed == expected)
patient_support_qc_pass <- nrow(patient_support_qc) == 6L &&
  isTRUE(all(patient_support_qc$pass))

class_counts <- remodeling_classes_interim %>%
  count(cohort, architecture_class, name = "n_genes", .drop = FALSE) %>%
  group_by(cohort) %>%
  mutate(fraction = n_genes / sum(n_genes)) %>%
  ungroup()

write_tsv(bootstrap_summary, file.path(m04$tables, "bootstrap_gene_stability.tsv.gz"))
write_tsv(bootstrap_class_probability, file.path(m04$tables, "bootstrap_class_probabilities.tsv.gz"))
write_tsv(bootstrap_qc, file.path(m04$qc, "bootstrap_qc.tsv"))
write_tsv(bootstrap_class_sum_qc, file.path(m04$qc, "bootstrap_class_probability_sums.tsv.gz"))
write_tsv(patient_support_sensitivity, file.path(m04$tables, "patient_support_threshold_sensitivity.tsv"))
write_tsv(patient_support_exact, file.path(m04$tables, "patient_support_exact_distribution.tsv"))
write_tsv(patient_support_qc, file.path(m04$qc, "patient_support_qc.tsv"))
write_tsv(class_counts, file.path(m04$tables, "remodeling_class_counts.tsv"))
write_tsv(
  remodeling_classes_interim,
  file.path(m04$tables, "gene_remodeling_classes_v2_without_gene_model.tsv.gz")
)


# ---- m04-leave-one-out (source-derived computation) ----
log_msg("Remodeling step: leave-one-out.")
loo_options <- public_cfg$sensitivity$observed_leave_one_patient_out %||% list(enabled = TRUE, patients = "all")
loo_patients <- if (identical(loo_options$enabled, FALSE)) character() else {
  requested <- as.character(unlist(loo_options$patients %||% "all"))
  if ("all" %in% requested) patient_ids else intersect(patient_ids, requested)
}
# With no requested LOO, run no resampling or model fits; retain typed empty tables.
loo_results <- map_dfr(loo_patients, function(omitted_patient) {
  included <- setdiff(patient_ids, omitted_patient)
  cohort_name <- paste0("omit_", omitted_patient)
  log_msg("Calculating leave-one-patient-out metrics for ", omitted_patient, ".")
  rbindlist(lapply(gene_objects, function(object) {
    analyze_gene_observed(object, included, cohort_name)$gene
  }), fill = TRUE) %>%
    as_tibble() %>%
    mutate(omitted_patient = omitted_patient)
})

if (!length(loo_patients)) {
  loo_results <- observed_gene[0, , drop = FALSE] %>% mutate(omitted_patient = character())
}
primary_class_lookup <- remodeling_classes_interim %>%
  filter(cohort == primary_cohort_name) %>%
  transmute(
    gene_id,
    full_class = as.character(architecture_class),
    full_observed_architecture = observed_architecture,
    full_effect_pass = effect_pass,
    full_driver = driver_isoform,
    full_jsd = jsd,
    full_max_abs_delta_if = max_abs_delta_if
  )

loo_comparison <- loo_results %>%
  left_join(primary_class_lookup, by = "gene_id") %>%
  mutate(
    loo_effect_pass = observed_architecture %in% c("A", "B", "C", "D", "E"),
    full_is_remodeled = full_observed_architecture %in% c("A", "B", "C", "D", "E"),
    class_matches_full = quantifiable & observed_architecture == full_observed_architecture,
    effect_matches_full = loo_effect_pass == full_is_remodeled,
    driver_matches_full = quantifiable & driver_isoform == full_driver
  )

loo_gene_stability <- loo_comparison %>%
  group_by(gene_id) %>%
  summarise(
    loo_quantifiable_fraction = mean(quantifiable),
    loo_effect_agreement = mean(effect_matches_full, na.rm = TRUE),
    loo_exact_class_agreement = mean(class_matches_full, na.rm = TRUE),
    loo_driver_agreement = mean(driver_matches_full, na.rm = TRUE),
    loo_min_jsd = min(jsd, na.rm = TRUE),
    loo_max_jsd = max(jsd, na.rm = TRUE),
    loo_min_max_abs_delta_if = min(max_abs_delta_if, na.rm = TRUE),
    loo_max_max_abs_delta_if = max(max_abs_delta_if, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("loo_"), ~ifelse(is.infinite(.x), NA_real_, .x)))

loo_summary <- loo_comparison %>%
  group_by(omitted_patient) %>%
  summarise(
    quantifiable_genes = sum(quantifiable),
    effect_calls = sum(loo_effect_pass),
    full_effect_calls_retained = sum(full_is_remodeled & loo_effect_pass),
    full_effect_retention_fraction = full_effect_calls_retained /
      max(1, sum(full_is_remodeled)),
    exact_class_agreement_among_full_effects = mean(
      class_matches_full[full_is_remodeled], na.rm = TRUE
    ),
    driver_agreement_among_full_effects = mean(
      driver_matches_full[full_is_remodeled], na.rm = TRUE
    ),
    .groups = "drop"
  )


write_tsv(loo_results, file.path(m04$tables, "leave_one_patient_out_gene_metrics.tsv.gz"))
write_tsv(loo_comparison, file.path(m04$tables, "leave_one_patient_out_comparison.tsv.gz"))
write_tsv(loo_gene_stability, file.path(m04$tables, "leave_one_patient_out_gene_stability.tsv.gz"))
write_tsv(loo_summary, file.path(m04$tables, "leave_one_patient_out_summary.tsv"))

# ---- m04-reml-runner (source-derived computation) ----
log_msg("Remodeling step: reml-runner.")
prepare_gene_model <- function(object, observed_row, included_patients) {
  view <- cohort_gene_view(object, included_patients)
  unavailable <- function(reason) list(available = FALSE, reason = reason)
  if (!isTRUE(view$quantifiable)) return(unavailable("not_quantifiable"))
  k <- object$n_isoforms_retained
  if (k < 2L || k > cfg$gene_null$maximum_isoforms) {
    return(unavailable("isoform_complexity_outside_K2_to_Kmax_scope"))
  }
  counts <- rbind(view$primary_counts, view$recurrent_counts)
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) ||
      any(rowSums(counts) < 1)) {
    return(unavailable("invalid_integer_count_matrix"))
  }
  smoothed <- smoothed_composition(counts, cfg$gene_null$prior_strength)
  basis <- ilr_basis(k)
  coordinates <- ilr_transform(smoothed$composition, basis)
  n_pair <- nrow(view$pair_table)
  primary_rows <- seq_len(n_pair)
  recurrent_rows <- n_pair + primary_rows
  difference <- coordinates[recurrent_rows, , drop = FALSE] -
    coordinates[primary_rows, , drop = FALSE]
  midpoint <- 0.5 * (
    coordinates[recurrent_rows, , drop = FALSE] +
      coordinates[primary_rows, , drop = FALSE]
  )
  centered <- sweep(difference, 2L, colMeans(difference), "-")
  residual_df <- n_pair - 1L
  dimension <- ncol(difference)
  total_variance <- sum(centered^2) / (residual_df * dimension)
  list(
    available = TRUE,
    reason = "eligible_family_member",
    gene_id = object$gene_id,
    basis = basis,
    counts = counts,
    coordinates = coordinates,
    difference = difference,
    midpoint = midpoint,
    primary_totals = rowSums(view$primary_counts),
    recurrent_totals = rowSums(view$recurrent_counts),
    pair_ids = view$pair_table$patient,
    sample_ids = c(view$pair_table$primary, view$pair_table$recurrent),
    n_pairs = n_pair,
    dimension = dimension,
    residual_df = residual_df,
    total_variance = total_variance,
    median_depth = median(rowSums(counts)),
    zero_fraction = mean(counts == 0),
    rare_count_fraction = mean(counts > 0 & counts <= 2),
    observed_score = joint_score(observed_row$jsd, observed_row$max_abs_delta_if)
  )
}

fit_gene_model <- function(
    object, observed_row, included_patients, cohort_name,
    null_replicates = cfg$gene_null$replicates,
    null_seed_namespace = "reference_20000",
    compute_gpd = cfg$gpd$enabled
) {
  null_replicates <- as.integer(null_replicates)
  fit <- prepare_gene_model(object, observed_row, included_patients)
  base <- tibble(
    cohort = cohort_name,
    gene_id = object$gene_id,
    gene_tested = isTRUE(fit$available),
    gene_test_reason = fit$reason,
    gene_test_n_isoforms = object$n_isoforms_retained,
    gene_test_n_complete_pairs = if (isTRUE(fit$available)) fit$n_pairs else NA_integer_,
    gene_test_ilr_dimension = if (isTRUE(fit$available)) fit$dimension else NA_integer_,
    gene_test_residual_df = if (isTRUE(fit$available)) fit$residual_df else NA_integer_,
    gene_test_total_variance = if (isTRUE(fit$available)) fit$total_variance else NA_real_,
    gene_test_median_depth = if (isTRUE(fit$available)) fit$median_depth else NA_real_,
    gene_test_zero_fraction = if (isTRUE(fit$available)) fit$zero_fraction else NA_real_,
    gene_test_rare_count_fraction = if (isTRUE(fit$available)) fit$rare_count_fraction else NA_real_,
    gene_model_observed_score = if (isTRUE(fit$available)) fit$observed_score else NA_real_,
    gene_model_null_seed_namespace = as.character(null_seed_namespace)
  )
  failed <- function(reason, covariance_diagnostics = tibble(), additions = list()) {
    summary <- base %>% mutate(
      gene_model_p = 1,
      gene_model_exceedances = NA_integer_,
      gene_model_replicates = null_replicates,
      gene_model_minimum_p = 1 / (null_replicates + 1),
      gene_model_monte_carlo_se = NA_real_,
      gene_model_null_score_mean = NA_real_,
      gene_model_null_score_sd = NA_real_,
      gene_model_null_score_q50 = NA_real_,
      gene_model_null_score_q90 = NA_real_,
      gene_model_null_score_q95 = NA_real_,
      gene_model_null_score_q99 = NA_real_,
      gene_model_mixture_p = NA_real_,
      gene_model_boundary_safety_p = NA_real_,
      gene_model_boundary_safety_used = FALSE,
      gene_model_selected_null_branch = "conservative_p1_failure",
      gene_model_inference_reason = paste0("conservative_p1_", reason),
      gene_model_gpd_raw_p = NA_real_,
      gene_model_gpd_hybrid_p = 1,
      gene_model_gpd_hybrid_method = "conservative_p1_failure",
      gene_model_gpd_available = FALSE,
      gene_model_gpd_threshold_stable = FALSE,
      gene_model_gpd_selected_tail_fraction = NA_real_,
      gene_model_gpd_shape = NA_real_,
      gene_model_gpd_scale = NA_real_,
      gene_model_gpd_ks_distance = NA_real_,
      gene_model_gpd_reason = paste0("conservative_p1_", reason),
      gene_test_sampling_variance = NA_real_,
      gene_test_reml_estimate = NA_real_,
      gene_test_reml_boundary = NA,
      gene_test_reml_converged = FALSE,
      gene_test_reml_objective_zero = NA_real_,
      gene_test_reml_objective_hat = NA_real_,
      gene_test_reml_information = NA_real_,
      gene_test_reml_curvature = NA_real_,
      gene_test_reml_upper = NA_real_,
      gene_test_profile_tail_resolved = FALSE,
      gene_test_profile_upper = NA_real_,
      gene_test_mixture_mean = NA_real_,
      gene_test_mixture_sd = NA_real_,
      gene_test_mixture_q05 = NA_real_,
      gene_test_mixture_q50 = NA_real_,
      gene_test_mixture_q95 = NA_real_,
      gene_model_variance_uncertainty_propagated = FALSE
    )
    if (length(additions)) {
      for (name in names(additions)) summary[[name]] <- additions[[name]]
    }
    list(
      summary = summary, covariance = covariance_diagnostics,
      gpd = tibble(), null_spec = NULL
    )
  }
  if (!isTRUE(fit$available)) return(failed(fit$reason))

  n_sample <- nrow(fit$counts)
  covariance <- vector("list", n_sample)
  covariance_rows <- vector("list", n_sample)
  side <- c(rep("primary", fit$n_pairs), rep("recurrent", fit$n_pairs))
  pair_index <- c(seq_len(fit$n_pairs), seq_len(fit$n_pairs))
  for (sample_index in seq_len(n_sample)) {
    covariance_seed <- deterministic_seed(
      cfg$seed, "finite_count_covariance", cohort_name, fit$gene_id,
      fit$sample_ids[[sample_index]]
    )
    estimate <- finite_count_covariance(
      fit$counts[sample_index, ], fit$basis, covariance_seed
    )
    qc <- estimate$qc %||% list()
    covariance_rows[[sample_index]] <- tibble(
      cohort = cohort_name,
      gene_id = fit$gene_id,
      sample_role = side[[sample_index]],
      pair_id = fit$pair_ids[[pair_index[[sample_index]]]],
      sample_id = fit$sample_ids[[sample_index]],
      covariance_available = isTRUE(estimate$available),
      covariance_reason = estimate$reason,
      covariance_replicates = estimate$replicates,
      covariance_trace = qc$trace %||% NA_real_,
      covariance_trace_relative_discrepancy =
        qc$trace_relative_discrepancy %||% NA_real_,
      covariance_frobenius_relative_discrepancy =
        qc$frobenius_relative_discrepancy %||% NA_real_,
      covariance_minimum_eigenvalue = qc$minimum_eigenvalue %||% NA_real_,
      covariance_maximum_eigenvalue = qc$maximum_eigenvalue %||% NA_real_,
      covariance_condition_number = qc$condition_number %||% NA_real_,
      random_seed = covariance_seed
    )
    if (!isTRUE(estimate$available)) {
      return(failed(
        paste0("finite_count_covariance_failed:", estimate$reason),
        bind_rows(covariance_rows)
      ))
    }
    covariance[[sample_index]] <- estimate$covariance
  }
  pair_covariance <- lapply(seq_len(fit$n_pairs), function(index) {
    covariance[[index]] + covariance[[fit$n_pairs + index]]
  })
  sample_trace <- vapply(covariance, function(value) {
    sum(diag(value)) / fit$dimension
  }, numeric(1))
  pair_trace <- sample_trace[seq_len(fit$n_pairs)] +
    sample_trace[fit$n_pairs + seq_len(fit$n_pairs)]
  sampling_variance <- sum((1 - 1 / fit$n_pairs) * pair_trace) / fit$residual_df

  profile <- profile_variance_mixture(fit$difference, pair_covariance)
  if (!isTRUE(profile$available)) {
    failed_fit <- profile$fit
    return(failed(
      paste0("variance_profile_failed:", profile$reason),
      bind_rows(covariance_rows),
      list(
        gene_test_sampling_variance = sampling_variance,
        gene_test_reml_estimate = failed_fit$estimate %||% NA_real_,
        gene_test_reml_boundary = failed_fit$boundary %||% NA,
        gene_test_reml_converged = failed_fit$converged %||% FALSE,
        gene_test_reml_objective_zero = failed_fit$objective_zero %||% NA_real_,
        gene_test_reml_objective_hat = failed_fit$objective_hat %||% NA_real_,
        gene_test_reml_information = failed_fit$information %||% NA_real_,
        gene_test_reml_curvature = failed_fit$curvature %||% NA_real_,
        gene_test_reml_upper = failed_fit$upper %||% NA_real_,
        gene_test_profile_tail_resolved = profile$tail_resolved %||% FALSE,
        gene_test_profile_upper = profile$upper %||% NA_real_
      )
    ))
  }

  variance_seed <- if (identical(null_seed_namespace, "reference_20000")) {
    deterministic_seed(cfg$seed, "variance_mixture", cohort_name, fit$gene_id)
  } else {
    deterministic_seed(
      cfg$seed, "variance_mixture", null_seed_namespace, cohort_name, fit$gene_id
    )
  }
  mixture_seed <- if (identical(null_seed_namespace, "reference_20000")) {
    deterministic_seed(cfg$seed, "null_mixture", cohort_name, fit$gene_id)
  } else {
    deterministic_seed(
      cfg$seed, "null_mixture", null_seed_namespace, cohort_name, fit$gene_id
    )
  }
  boundary_seed <- if (identical(null_seed_namespace, "reference_20000")) {
    deterministic_seed(cfg$seed, "null_boundary_safety", cohort_name, fit$gene_id)
  } else {
    deterministic_seed(
      cfg$seed, "null_boundary_safety", null_seed_namespace,
      cohort_name, fit$gene_id
    )
  }
  variance_draw <- with_seed(
    variance_seed,
    sample(profile$grid, null_replicates, replace = TRUE, prob = profile$weight)
  )
  mixture_scores <- simulate_joint_scores(
    fit, variance_draw, null_replicates, mixture_seed
  )
  mixture_result <- monte_carlo_p(mixture_scores, fit$observed_score)
  if (is.null(mixture_result)) {
    return(failed("null_mixture_simulation_failed", bind_rows(covariance_rows)))
  }
  final_result <- mixture_result
  boundary_result <- NULL
  selected_branch <- "variance_profile_mixture"
  if (isTRUE(profile$fit$boundary)) {
    boundary_scores <- simulate_joint_scores(
      fit, rep(fit$total_variance, null_replicates),
      null_replicates, boundary_seed
    )
    boundary_result <- monte_carlo_p(boundary_scores, fit$observed_score)
    if (is.null(boundary_result)) {
      return(failed("boundary_safety_simulation_failed", bind_rows(covariance_rows)))
    }
    if (boundary_result$p > mixture_result$p) {
      final_result <- boundary_result
      selected_branch <- "total_variance_boundary_safety"
    }
  }

  mixture_gpd <- if (isTRUE(compute_gpd)) {
    estimate_gpd_tail(
      mixture_scores, fit$observed_score, mixture_result,
      "variance_profile_mixture"
    )
  } else {
    list(
      summary = list(
        raw_p = NA_real_, hybrid_p = mixture_result$p,
        hybrid_method = "GPD_not_requested", available = FALSE,
        threshold_stable = FALSE, selected_fraction = NA_real_,
        shape = NA_real_, scale = NA_real_, ks_distance = NA_real_,
        reason = "GPD_not_requested"
      ),
      diagnostics = tibble()
    )
  }
  boundary_gpd <- NULL
  if (!is.null(boundary_result)) {
    boundary_gpd <- if (isTRUE(compute_gpd)) {
      estimate_gpd_tail(
        boundary_scores, fit$observed_score, boundary_result,
        "total_variance_boundary_safety"
      )
    } else {
      list(
        summary = list(
          raw_p = NA_real_, hybrid_p = boundary_result$p,
          hybrid_method = "GPD_not_requested", available = FALSE,
          threshold_stable = FALSE, selected_fraction = NA_real_,
          shape = NA_real_, scale = NA_real_, ks_distance = NA_real_,
          reason = "GPD_not_requested"
        ),
        diagnostics = tibble()
      )
    }
  }
  gpd_raw_p <- if (is.null(boundary_gpd)) {
    mixture_gpd$summary$raw_p
  } else if (all(is.finite(c(
    mixture_gpd$summary$raw_p, boundary_gpd$summary$raw_p
  )))) {
    max(mixture_gpd$summary$raw_p, boundary_gpd$summary$raw_p)
  } else {
    NA_real_
  }
  gpd_hybrid_p <- if (is.null(boundary_gpd)) {
    mixture_gpd$summary$hybrid_p
  } else {
    max(mixture_gpd$summary$hybrid_p, boundary_gpd$summary$hybrid_p)
  }
  gpd_selected <- if (!is.null(boundary_gpd) &&
      boundary_gpd$summary$hybrid_p >= mixture_gpd$summary$hybrid_p) {
    boundary_gpd$summary
  } else {
    mixture_gpd$summary
  }
  gpd_diagnostics <- bind_rows(
    mixture_gpd$diagnostics,
    if (is.null(boundary_gpd)) tibble() else boundary_gpd$diagnostics
  ) %>%
    mutate(
      cohort = cohort_name,
      gene_id = fit$gene_id,
      null_seed_namespace = as.character(null_seed_namespace),
      .before = 1L
    )

  empirical_p <- if (isTRUE(profile$fit$boundary)) {
    max(mixture_result$p, boundary_result$p)
  } else mixture_result$p
  summary <- base %>% mutate(
    gene_model_p = empirical_p,
    gene_model_exceedances = final_result$exceedances,
    gene_model_replicates = final_result$replicates,
    gene_model_minimum_p = final_result$minimum_p,
    gene_model_monte_carlo_se = if (isTRUE(profile$fit$boundary)) NA_real_ else
      final_result$monte_carlo_se,
    gene_model_null_score_mean = final_result$mean,
    gene_model_null_score_sd = final_result$sd,
    gene_model_null_score_q50 = final_result$q50,
    gene_model_null_score_q90 = final_result$q90,
    gene_model_null_score_q95 = final_result$q95,
    gene_model_null_score_q99 = final_result$q99,
    gene_model_mixture_p = mixture_result$p,
    gene_model_boundary_safety_p = if (is.null(boundary_result)) NA_real_ else boundary_result$p,
    gene_model_boundary_safety_used = isTRUE(profile$fit$boundary),
    gene_model_selected_null_branch = selected_branch,
    gene_model_inference_reason = if (isTRUE(profile$fit$boundary)) {
      "tested_with_boundary_safety_rule"
    } else "tested_with_variance_uncertainty_mixture",
    gene_model_gpd_raw_p = gpd_raw_p,
    gene_model_gpd_hybrid_p = gpd_hybrid_p,
    gene_model_gpd_hybrid_method = gpd_selected$hybrid_method,
    gene_model_gpd_available = gpd_selected$available,
    gene_model_gpd_threshold_stable = gpd_selected$threshold_stable,
    gene_model_gpd_selected_tail_fraction = gpd_selected$selected_fraction,
    gene_model_gpd_shape = gpd_selected$shape,
    gene_model_gpd_scale = gpd_selected$scale,
    gene_model_gpd_ks_distance = gpd_selected$ks_distance,
    gene_model_gpd_reason = gpd_selected$reason,
    gene_test_sampling_variance = sampling_variance,
    gene_test_reml_estimate = profile$fit$estimate,
    gene_test_reml_boundary = profile$fit$boundary,
    gene_test_reml_converged = profile$fit$converged,
    gene_test_reml_objective_zero = profile$fit$objective_zero,
    gene_test_reml_objective_hat = profile$fit$objective_hat,
    gene_test_reml_information = profile$fit$information,
    gene_test_reml_curvature = profile$fit$curvature,
    gene_test_reml_upper = profile$fit$upper,
    gene_test_profile_tail_resolved = profile$tail_resolved,
    gene_test_profile_upper = profile$upper,
    gene_test_mixture_mean = profile$mean,
    gene_test_mixture_sd = profile$sd,
    gene_test_mixture_q05 = profile$q05,
    gene_test_mixture_q50 = profile$q50,
    gene_test_mixture_q95 = profile$q95,
    gene_model_variance_uncertainty_propagated = TRUE
  )
  retain_null_spec <- identical(null_seed_namespace, "reference_20000") &&
    cfg$refinement$enabled && is.finite(empirical_p) &&
    empirical_p < cfg$refinement$trigger_p
  null_spec <- if (retain_null_spec) {
    list(
      fit = fit,
      profile_grid = profile$grid,
      profile_weight = profile$weight,
      boundary = isTRUE(profile$fit$boundary),
      total_variance = fit$total_variance
    )
  } else NULL
  list(
    summary = summary,
    covariance = bind_rows(covariance_rows),
    gpd = gpd_diagnostics,
    null_spec = null_spec
  )
}

build_gene_model_schedule <- function(
    observed_lookup, family_ids, included_patients, cohort_name) {
  if (!length(family_ids)) {
    return(tibble(
      gene_id = character(), scheduling_complexity_stratum = character(),
      scheduling_depth_stratum = character(), scheduling_sparsity_stratum = character(),
      scheduling_effect_stratum = character(), scheduling_combined_stratum = character(),
      scheduling_median_sample_depth = numeric(), scheduling_zero_fraction = numeric(),
      scheduling_rare_count_fraction = numeric(), scheduling_priority = integer(),
      scheduled_for_computation = logical(), scheduling_design = character()
    ))
  }

  count_metrics <- map_dfr(family_ids, function(gene_id) {
    view <- cohort_gene_view(gene_objects[[gene_id]], included_patients)
    if (!isTRUE(view$quantifiable)) {
      stop("A frozen gene-model family member was not quantifiable: ", gene_id)
    }
    counts <- rbind(view$primary_counts, view$recurrent_counts)
    tibble(
      gene_id = gene_id,
      scheduling_median_sample_depth = median(rowSums(counts)),
      scheduling_zero_fraction = mean(counts == 0),
      scheduling_rare_count_fraction = mean(counts > 0 & counts < 5)
    )
  })

  frame <- observed_lookup %>%
    filter(gene_id %in% family_ids) %>%
    select(gene_id, n_isoforms_retained, effect_pass, joint_effect_score) %>%
    left_join(count_metrics, by = "gene_id") %>%
    mutate(
      effect_pass = effect_pass %in% TRUE,
      scheduling_complexity_stratum = cut(
        n_isoforms_retained,
        breaks = c(1, 2, 4, 9, 14, Inf),
        labels = c("K2", "K3-4", "K5-9", "K10-14", "K15-20"),
        include.lowest = TRUE, right = TRUE
      ),
      scheduling_depth_stratum = cut(
        scheduling_median_sample_depth,
        breaks = c(-Inf, 19, 49, 199, 999, Inf),
        labels = c("N<20", "N20-49", "N50-199", "N200-999", "N>=1000"),
        include.lowest = TRUE, right = TRUE
      ),
      scheduling_sparsity_stratum = cut(
        scheduling_zero_fraction,
        breaks = c(-Inf, 0.25, 0.50, 0.75, Inf),
        labels = c("zero<25%", "zero25-49%", "zero50-74%", "zero>=75%"),
        include.lowest = TRUE, right = FALSE
      ),
      scheduling_effect_stratum = if_else(
        effect_pass, "effect_gate_positive", "effect_gate_negative"
      ),
      scheduling_combined_stratum = interaction(
        scheduling_complexity_stratum,
        scheduling_depth_stratum,
        scheduling_sparsity_stratum,
        drop = TRUE, sep = "|"
      ),
      scheduling_gene_seed = vapply(gene_id, function(value) {
        deterministic_seed(cfg$seed, "gene_model_schedule", cohort_name, value)
      }, integer(1)),
      scheduling_stratum_seed = vapply(
        as.character(scheduling_combined_stratum),
        function(value) {
          deterministic_seed(cfg$seed, "gene_model_stratum", cohort_name, value)
        },
        integer(1)
      )
    ) %>%
    group_by(scheduling_effect_stratum, scheduling_combined_stratum) %>%
    arrange(scheduling_gene_seed, .by_group = TRUE) %>%
    mutate(scheduling_round = row_number()) %>%
    ungroup() %>%
    group_by(scheduling_effect_stratum) %>%
    arrange(
      scheduling_round, scheduling_stratum_seed, scheduling_gene_seed,
      .by_group = TRUE
    ) %>%
    mutate(scheduling_effect_group_rank = row_number()) %>%
    ungroup()

  maximum_genes <- min(cfg$gene_null$maximum_genes, nrow(frame))
  if (cfg$gene_null$maximum_genes == 0L) {
    selected <- frame %>%
      arrange(scheduling_gene_seed) %>%
      mutate(scheduling_priority = row_number())
    scheduling_design <- "all_frozen_family_members"
  } else {
    positive_available <- sum(frame$effect_pass)
    negative_available <- nrow(frame) - positive_available
    positive_quota <- if (positive_available > 0L) {
      min(positive_available, max(1L, ceiling(0.25 * maximum_genes)))
    } else 0L
    negative_quota <- min(negative_available, maximum_genes - positive_quota)
    if (positive_quota + negative_quota < maximum_genes) {
      positive_quota <- min(
        positive_available,
        maximum_genes - negative_quota
      )
    }

    selected <- bind_rows(
      frame %>%
        filter(effect_pass) %>%
        arrange(scheduling_effect_group_rank) %>%
        slice_head(n = positive_quota),
      frame %>%
        filter(!effect_pass) %>%
        arrange(scheduling_effect_group_rank) %>%
        slice_head(n = negative_quota)
    )
    if (nrow(selected) < maximum_genes) {
      selected <- bind_rows(
        selected,
        frame %>%
          filter(!gene_id %in% selected$gene_id) %>%
          arrange(scheduling_round, scheduling_stratum_seed, scheduling_gene_seed) %>%
          slice_head(n = maximum_genes - nrow(selected))
      )
    }
    selected <- selected %>%
      arrange(scheduling_effect_group_rank, scheduling_gene_seed) %>%
      mutate(scheduling_priority = row_number())
    scheduling_design <- paste0(
      "deterministic_stratified_cap", maximum_genes,
      "_25percent_effect_quota_round_robin_complexity_depth_sparsity"
    )
  }

  frame %>%
    mutate(
      scheduled_for_computation = gene_id %in% selected$gene_id,
      scheduling_priority = match(gene_id, selected$gene_id),
      scheduling_design = scheduling_design
    ) %>%
    select(
      gene_id, starts_with("scheduling_"), scheduled_for_computation
    )
}

run_gene_model_cohort <- function(cohort_name, included_patients) {
  observed_lookup <- observed_gene %>% filter(cohort == cohort_name) %>% arrange(gene_id)
  family <- observed_lookup %>%
    transmute(
      cohort = cohort_name,
      gene_id,
      included_in_testing_family = quantifiable &
        n_isoforms_retained >= 2L &
        n_isoforms_retained <= cfg$gene_null$maximum_isoforms &
        n_evaluable_pairs >= cfg$thresholds$min_evaluable_pairs,
      testing_family_reason = case_when(
        !quantifiable ~ quantifiability_reason,
        n_isoforms_retained > cfg$gene_null$maximum_isoforms ~
          "isoform_complexity_above_registered_scope",
        TRUE ~ "included_by_deterministic_paired_count_scope"
      )
    )
  family_ids <- family$gene_id[family$included_in_testing_family]
  schedule <- build_gene_model_schedule(
    observed_lookup, family_ids, included_patients, cohort_name
  )
  scheduled_ids <- schedule %>%
    filter(scheduled_for_computation) %>%
    arrange(scheduling_priority) %>%
    pull(gene_id)
  family <- family %>%
    left_join(schedule, by = "gene_id") %>%
    mutate(
      scheduled_for_computation = replace_na(scheduled_for_computation, FALSE),
      testing_family_rule = paste0(
      "paired_integer_count_quantifiable_K2_to_K", cfg$gene_null$maximum_isoforms,
      "_family_frozen_before_covariance_REML_or_null_inference;failures_P_equals_1"
    )
    )

  cache_path <- file.path(m04$cache, paste0("gene_null_", cohort_name, ".rds"))
  configuration_hash <- digest::digest(
    list(
      model_engine_version = "finite_count_REML_v0.5.0_unchanged",
      input_fingerprint = cache_input_fingerprint,
      cohort = cohort_name,
      included_patients = included_patients,
      gene_null = cfg$gene_null,
      gpd = cfg$gpd,
      thresholds = cache_thresholds,
      family_ids = family_ids,
      scheduled_ids = scheduled_ids
    ), algo = "sha256"
  )
  legacy_configuration_hashes <- vapply(
    legacy_input_fingerprints_v0_5,
    function(legacy_fingerprint) {
      digest::digest(
        list(
          workflow_version = "0.5.0",
          input_fingerprint = legacy_fingerprint,
          cohort = cohort_name,
          included_patients = included_patients,
          gene_null = cfg$gene_null,
          gpd = cfg$gpd,
          thresholds = legacy_thresholds_v0_5,
          family_ids = family_ids,
          scheduled_ids = scheduled_ids
        ),
        algo = "sha256"
      )
    },
    character(1)
  )
  state <- list(configuration_hash = configuration_hash, completed = list())
  if (cfg$resume_cached_results && file.exists(cache_path)) {
    candidate <- readRDS(cache_path)
    if (candidate$configuration_hash %in%
        c(configuration_hash, legacy_configuration_hashes)) {
      state <- candidate
      state$configuration_hash <- configuration_hash
      log_msg("Loaded compatible gene-null cache for ", cohort_name, ".")
    }
  }
  pending <- setdiff(scheduled_ids, names(state$completed))
  for (position in seq_along(pending)) {
    gene_id <- pending[[position]]
    observed_row <- observed_lookup[match(gene_id, observed_lookup$gene_id), , drop = FALSE]
    state$completed[[gene_id]] <- fit_gene_model(
      gene_objects[[gene_id]], observed_row, included_patients, cohort_name,
      null_replicates = cfg$gene_null$replicates,
      null_seed_namespace = "reference_20000",
      compute_gpd = cfg$gpd$enabled
    )
    if (position %% cfg$checkpoint_every_genes == 0L || position == length(pending)) {
      atomic_save_rds(state, cache_path)
      log_msg(
        "Gene null ", cohort_name, ": ", length(state$completed), " / ",
        length(scheduled_ids), " scheduled genes cached."
      )
    }
  }
  summary <- if (length(state$completed)) {
    rbindlist(lapply(state$completed, `[[`, "summary"), fill = TRUE) %>% as_tibble()
  } else fit_gene_model(gene_objects[[1L]], observed_lookup[1L, , drop=FALSE], character(), cohort_name)$summary[0L, , drop=FALSE]
  covariance <- if (length(state$completed)) {
    rbindlist(lapply(state$completed, `[[`, "covariance"), fill = TRUE) %>% as_tibble()
  } else tibble()
  gpd <- if (length(state$completed)) {
    rbindlist(lapply(state$completed, function(value) {
      value$gpd %||% tibble()
    }), fill = TRUE) %>% as_tibble()
  } else tibble()
  null_specs <- if (length(state$completed)) {
    output <- lapply(state$completed, function(value) value$null_spec %||% NULL)
    output[!vapply(output, is.null, logical(1))]
  } else list()
  unscheduled <- setdiff(family_ids, scheduled_ids)
  if (length(unscheduled)) {
    unscheduled_rows <- observed_lookup %>%
      filter(gene_id %in% unscheduled) %>%
      transmute(
        cohort = cohort_name, gene_id, gene_tested = FALSE,
        gene_test_reason = "eligible_but_not_run_under_development_cap",
        gene_test_n_isoforms = n_isoforms_retained,
        gene_test_n_complete_pairs = n_evaluable_pairs,
        gene_model_observed_score = joint_effect_score,
        gene_model_p = 1,
        gene_model_gpd_hybrid_p = 1,
        gene_model_inference_reason = "conservative_p1_development_cap"
      )
    summary <- bind_rows(summary, unscheduled_rows)
  }
  if (length(family_ids)) {
    missing_family <- setdiff(family_ids, summary$gene_id)
    if (length(missing_family)) stop("Gene-model cache did not cover the frozen testing family.")
    summary <- summary %>%
      arrange(match(gene_id, family_ids)) %>%
      mutate(
        gene_model_p = replace_na(gene_model_p, 1),
        gene_model_empirical_p = gene_model_p,
        gene_model_q = p.adjust(gene_model_p, method = "BH"),
        gene_model_empirical_bh = gene_model_q,
        gene_model_bh_le_0p05 = gene_model_q <= 0.05,
        gene_model_bh_le_0p10 = gene_model_q <= 0.10,
        gene_model_gpd_hybrid_p = replace_na(gene_model_gpd_hybrid_p, 1),
        gene_model_gpd_hybrid_q = p.adjust(
          gene_model_gpd_hybrid_p, method = "BH"
        ),
        gene_model_gpd_bh_le_0p05 = gene_model_gpd_hybrid_q <= 0.05,
        gene_model_gpd_bh_le_0p10 = gene_model_gpd_hybrid_q <= 0.10,
        gene_model_supported = gene_model_q <= cfg$gene_null$fdr,
        gene_model_effect_supported = gene_model_supported &
          gene_model_observed_score >= 1,
        gene_model_inference_reason = case_when(
          startsWith(gene_model_inference_reason, "conservative_p1_") ~
            gene_model_inference_reason,
          gene_model_q > cfg$gene_null$fdr ~ "gene_model_q_above_configured_fdr",
          gene_model_observed_score < 1 ~ "statistically_detectable_subthreshold_difference",
          TRUE ~ "gene_level_remodeling_supported"
        )
      )
  }
  if (!length(family_ids)) summary <- summary %>% mutate(
    gene_model_empirical_p = numeric(), gene_model_q = numeric(), gene_model_empirical_bh = numeric(),
    gene_model_bh_le_0p05 = logical(), gene_model_bh_le_0p10 = logical(),
    gene_model_gpd_hybrid_q = numeric(), gene_model_gpd_bh_le_0p05 = logical(),
    gene_model_gpd_bh_le_0p10 = logical(), gene_model_supported = logical(),
    gene_model_effect_supported = logical())
  list(
    family = family, summary = summary, covariance = covariance,
    gpd = gpd, null_specs = null_specs
  )
}

# ---- m04-reml-analysis (source-derived computation) ----
log_msg("Remodeling step: reml-analysis.")
if (cfg$gene_null$enabled) {
if (cfg$gene_null$enabled) {
  gene_model_results <- imap(analysis_cohorts, function(included_patients, cohort_name) {
    run_gene_model_cohort(cohort_name, included_patients)
  })
  gene_model_family <- bind_rows(lapply(gene_model_results, `[[`, "family"))
  gene_model_summary <- bind_rows(lapply(gene_model_results, `[[`, "summary"))
  gene_model_covariance <- bind_rows(lapply(gene_model_results, `[[`, "covariance"))
  gene_model_gpd_diagnostics <- bind_rows(lapply(gene_model_results, `[[`, "gpd"))
} else {
  gene_model_results <- list()
  gene_model_family <- tibble(
    cohort = character(), gene_id = character(),
    included_in_testing_family = logical(), testing_family_reason = character(),
    scheduled_for_computation = logical(), testing_family_rule = character()
  )
  gene_model_summary <- tibble(
    cohort = character(), gene_id = character(), gene_tested = logical(),
    gene_test_reason = character(), gene_model_p = numeric(), gene_model_q = numeric(),
    gene_model_supported = logical(), gene_model_effect_supported = logical(),
    gene_model_inference_reason = character(),
    gene_model_observed_score = numeric(), gene_model_exceedances = integer(),
    gene_model_gpd_hybrid_p = numeric(), gene_model_gpd_hybrid_q = numeric(),
    gene_model_gpd_hybrid_method = character(),
    gene_model_gpd_threshold_stable = logical()
  )
  gene_model_covariance <- tibble(
    cohort = character(), gene_id = character(), sample_role = character(),
    pair_id = character(), sample_id = character(), covariance_available = logical(),
    covariance_reason = character(), covariance_replicates = integer()
  )
  gene_model_gpd_diagnostics <- tibble(
    cohort = character(), gene_id = character(), null_seed_namespace = character(),
    gpd_branch = character(), gpd_tail_fraction = numeric(),
    gpd_threshold = numeric(), gpd_tail_size = integer(),
    gpd_scale = numeric(), gpd_shape = numeric(), gpd_tail_p = numeric(),
    gpd_ks_distance = numeric(), gpd_ks_limit = numeric(),
    gpd_gof_pass = logical(), gpd_available = logical(), gpd_reason = character()
  )
}
if (!all(c(
  "cohort", "gene_id", "gpd_available", "gpd_shape",
  "gpd_ks_distance", "gpd_ks_limit", "gpd_gof_pass", "gpd_tail_fraction"
) %in% names(gene_model_gpd_diagnostics))) {
  gene_model_gpd_diagnostics <- tibble(
    cohort = character(), gene_id = character(), null_seed_namespace = character(),
    gpd_branch = character(), gpd_tail_fraction = numeric(),
    gpd_threshold = numeric(), gpd_tail_size = integer(),
    gpd_scale = numeric(), gpd_shape = numeric(),
    gpd_tail_probability = numeric(), gpd_conditional_survival = numeric(),
    gpd_tail_p = numeric(), gpd_ks_distance = numeric(),
    gpd_ks_limit = numeric(), gpd_gof_pass = logical(),
    gpd_optimizer_converged = logical(), gpd_available = logical(),
    gpd_reason = character()
  )
}

gene_model_execution <- if (cfg$gene_null$enabled) {
  gene_model_family %>%
    filter(included_in_testing_family) %>%
    left_join(gene_model_summary, by = c("cohort", "gene_id")) %>%
    mutate(
      gene_tested = replace_na(gene_tested, FALSE),
      inference_reason_for_qc = replace_na(
        gene_model_inference_reason,
        "missing_gene_model_result"
      ),
      conservative_p1 = startsWith(inference_reason_for_qc, "conservative_p1_"),
      model_success = scheduled_for_computation & gene_tested & !conservative_p1,
      scheduled_failure = scheduled_for_computation & !model_success,
      scheduled_not_attempted = scheduled_for_computation & !gene_tested,
      development_cap_unrun = !scheduled_for_computation
    )
} else tibble()

gene_model_qc <- if (cfg$gene_null$enabled) {
  gene_model_execution %>%
    group_by(cohort) %>%
    summarise(
      frozen_testing_family = n(),
      scheduled_models = sum(scheduled_for_computation),
      attempted_models = sum(scheduled_for_computation & gene_tested),
      scheduled_not_attempted = sum(.data$scheduled_not_attempted),
      development_cap_unrun = sum(.data$development_cap_unrun),
      finite_p_values = sum(is.finite(gene_model_p)),
      conservative_p1_all_family = sum(conservative_p1),
      scheduled_conservative_failures = sum(scheduled_failure),
      successful_models = sum(model_success),
      model_success_fraction =
        sum(model_success) / max(1L, sum(scheduled_for_computation)),
      boundary_models = sum(
        model_success & gene_test_reml_boundary %in% TRUE,
        na.rm = TRUE
      ),
      raw_p_below_0p05 = sum(gene_model_p < 0.05, na.rm = TRUE),
      bh_below_0p05 = sum(gene_model_bh_le_0p05 %in% TRUE, na.rm = TRUE),
      bh_below_0p10 = sum(gene_model_bh_le_0p10 %in% TRUE, na.rm = TRUE),
      gpd_hybrid_bh_below_0p05 =
        sum(gene_model_gpd_bh_le_0p05 %in% TRUE, na.rm = TRUE),
      gpd_hybrid_bh_below_0p10 =
        sum(gene_model_gpd_bh_le_0p10 %in% TRUE, na.rm = TRUE),
      q_below_fdr = sum(gene_model_supported %in% TRUE, na.rm = TRUE),
      q_below_fdr_and_effect_gate =
        sum(gene_model_effect_supported %in% TRUE, na.rm = TRUE),
      minimum_attainable_p = 1 / (cfg$gene_null$replicates + 1),
      .groups = "drop"
    )
} else {
  tibble(
    cohort = names(analysis_cohorts), frozen_testing_family = NA_integer_,
    scheduled_models = NA_integer_, attempted_models = NA_integer_,
    scheduled_not_attempted = NA_integer_, development_cap_unrun = NA_integer_,
    finite_p_values = NA_integer_, conservative_p1_all_family = NA_integer_,
    scheduled_conservative_failures = NA_integer_, successful_models = NA_integer_,
    model_success_fraction = NA_real_, boundary_models = NA_integer_,
    raw_p_below_0p05 = NA_integer_, bh_below_0p05 = NA_integer_,
    bh_below_0p10 = NA_integer_, gpd_hybrid_bh_below_0p05 = NA_integer_,
    gpd_hybrid_bh_below_0p10 = NA_integer_,
    q_below_fdr = NA_integer_, q_below_fdr_and_effect_gate = NA_integer_,
    minimum_attainable_p = NA_real_
  )
}

gene_model_execution_by_stratum <- if (cfg$gene_null$enabled) {
  bind_rows(
    gene_model_execution %>% transmute(
      cohort, scheduled_for_computation, model_success,
      stratum_axis = "isoform_complexity",
      stratum = as.character(scheduling_complexity_stratum)
    ),
    gene_model_execution %>% transmute(
      cohort, scheduled_for_computation, model_success,
      stratum_axis = "sample_depth",
      stratum = as.character(scheduling_depth_stratum)
    ),
    gene_model_execution %>% transmute(
      cohort, scheduled_for_computation, model_success,
      stratum_axis = "isoform_sparsity",
      stratum = as.character(scheduling_sparsity_stratum)
    ),
    gene_model_execution %>% transmute(
      cohort, scheduled_for_computation, model_success,
      stratum_axis = "observed_effect_gate",
      stratum = as.character(scheduling_effect_stratum)
    )
  ) %>%
    filter(!is.na(stratum), nzchar(stratum)) %>%
    group_by(cohort, stratum_axis, stratum) %>%
    summarise(
      frozen_family_models = n(),
      scheduled_models = sum(scheduled_for_computation),
      scheduling_fraction =
        sum(scheduled_for_computation) / n(),
      successful_models = sum(model_success),
      model_success_fraction =
        sum(model_success) / max(1L, sum(scheduled_for_computation)),
      stratum_represented = sum(scheduled_for_computation) > 0L,
      adequately_populated =
        sum(scheduled_for_computation) >= cfg$gene_null$minimum_stratum_n,
      pass = sum(scheduled_for_computation) < cfg$gene_null$minimum_stratum_n |
        (sum(model_success) / max(1L, sum(scheduled_for_computation))) >=
          cfg$gene_null$minimum_stratum_success_fraction,
      .groups = "drop"
    )
} else tibble(adequately_populated = logical(), model_success_fraction = numeric())

minimum_overall_success <- if (cfg$gene_null$enabled && nrow(gene_model_qc)) {
  min(gene_model_qc$model_success_fraction)
} else 1
minimum_attempt_fraction <- if (cfg$gene_null$enabled && nrow(gene_model_qc)) {
  min(gene_model_qc$attempted_models / pmax(1L, gene_model_qc$scheduled_models))
} else 1
minimum_stratum_representation <- if (
    cfg$gene_null$enabled && nrow(gene_model_execution_by_stratum)) {
  min(as.numeric(gene_model_execution_by_stratum$stratum_represented))
} else 1
populated_execution_strata <- gene_model_execution_by_stratum %>%
  filter(adequately_populated)
minimum_populated_stratum_success <- if (nrow(populated_execution_strata)) {
  min(populated_execution_strata$model_success_fraction)
} else 1

gene_model_execution_gate <- tibble(
  check = c(
    "all_scheduled_models_attempted",
    "all_observed_schedule_strata_represented",
    "overall_scheduled_model_success_fraction",
    "adequately_populated_stratum_success_fraction"
  ),
  observed = c(
    minimum_attempt_fraction,
    minimum_stratum_representation,
    minimum_overall_success,
    minimum_populated_stratum_success
  ),
  target = c(
    1,
    1,
    cfg$gene_null$minimum_success_fraction,
    cfg$gene_null$minimum_stratum_success_fraction
  ),
  pass = c(
    minimum_attempt_fraction == 1,
    minimum_stratum_representation == 1,
    minimum_overall_success >= cfg$gene_null$minimum_success_fraction,
    minimum_populated_stratum_success >=
      cfg$gene_null$minimum_stratum_success_fraction
  )
)

gene_model_failure_reasons <- if (cfg$gene_null$enabled) {
  gene_model_execution %>%
    filter(scheduled_failure) %>%
    count(cohort, inference_reason_for_qc, name = "n_genes") %>%
    left_join(
      gene_model_qc %>% select(cohort, scheduled_models),
      by = "cohort"
    ) %>%
    mutate(
      fraction_of_scheduled = n_genes / pmax(1L, scheduled_models)
    ) %>%
    select(-scheduled_models)
} else tibble()

gene_model_method_status <- tibble(
  component = c(
    "method", "finite-count covariance", "between-patient variance",
    "variance uncertainty", "boundary safety", "test statistic",
    "reference Monte Carlo resolution", "resolution refinement",
    "GPD tail estimator", "cohort pairing", "directional support",
    "patient-specific test", "multiplicity", "role in class assignment",
    "locked analytical role", "validation reporting"
  ),
  specification = c(
    "DEIR-inspired standalone paired finite-count sensitivity model; manuscript implementation v0.6.0",
    "parametric multinomial covariance after 0.5-count smoothing and full ILR transformation",
    "per-gene additive REML variance in ILR difference space",
    "profile-information variance mixture",
    "maximum of profile-mixture and total-variance null P values at the zero boundary",
    "max(JSD/0.20, max|delta IF|/0.25)",
    "20,000 draws per tested gene in the production profile",
    paste(
      "fixed independent second-stage simulation for genes with reference P<0.05;",
      "the fitted covariance/REML/profile null is reused and non-selected genes receive P=1"
    ),
    "diagnostics-gated generalized-Pareto tail estimate for low-exceedance REML nulls with empirical fallback",
    "primary and recurrent samples are paired within patient before cohort estimation",
    "reported independently from the cohort P value; it is not an input to the REML null",
    paste(
      "two-condition multinomial likelihood-ratio test after rare-component pooling;",
      "technical count-composition evidence, not biological-replicate inference"
    ),
    paste(
      "cohort BH within each frozen cohort family; patient BH separately within",
      "each patient's frozen family; report 0.05 and 0.10; failures are conservative"
    ),
    "none; P/Q values are a separate inferential axis",
    paste(
      "20k paired REML is the reference cohort test; independent refinement resolves",
      "reference-P<0.05 candidates; GPD is diagnostic; patient LRT is patient-specific"
    ),
    "oracle and production truth-known comparisons are retained for the reviewer response, not required in the manuscript body"
  )
)

write_tsv(gene_model_family, file.path(m04$tables, "gene_model_testing_family.tsv.gz"))
write_tsv(gene_model_summary, file.path(m04$tables, "gene_model_results.tsv.gz"))
write_tsv(gene_model_covariance, file.path(m04$qc, "gene_model_covariance_diagnostics.tsv.gz"))
write_tsv(
  gene_model_gpd_diagnostics,
  file.path(m04$qc, "gene_model_gpd_tail_diagnostics.tsv.gz")
)
write_tsv(gene_model_qc, file.path(m04$qc, "gene_model_qc.tsv"))
write_tsv(
  gene_model_execution_by_stratum,
  file.path(m04$qc, "gene_model_execution_by_stratum.tsv")
)
write_tsv(
  gene_model_execution_gate,
  file.path(m04$qc, "gene_model_execution_gate.tsv")
)
write_tsv(
  gene_model_failure_reasons,
  file.path(m04$qc, "gene_model_failure_reasons.tsv")
)
write_tsv(gene_model_method_status, file.path(m04$tables, "gene_model_method_status.tsv"))

if (nrow(gene_model_failure_reasons)) {
}

# ---- m04-resolution-refinement (source-derived computation) ----
log_msg("Remodeling step: resolution-refinement.")
safe_binomial_interval <- function(successes, trials, confidence = 0.95) {
  if (!is.finite(successes) || !is.finite(trials) || trials < 1L ||
      successes < 0L || successes > trials) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  interval <- binom.test(
    as.integer(successes), as.integer(trials), conf.level = confidence
  )$conf.int
  c(lower = interval[[1]], upper = interval[[2]])
}

simulate_independent_refinement <- function(null_spec, cohort_name, gene_id) {
  null_seed_namespace <- "independent_resolution_refinement"
  failed <- function(reason) {
    list(
      summary = tibble(
        gene_id = gene_id,
        gene_tested = FALSE,
        gene_model_p = 1,
        gene_model_exceedances = NA_integer_,
        gene_model_replicates = cfg$refinement$replicates,
        gene_model_selected_null_branch = "conservative_p1_failure",
        gene_model_null_seed_namespace = null_seed_namespace,
        gene_model_inference_reason = paste0("conservative_p1_", reason),
        refinement_reused_frozen_null_spec = FALSE
      )
    )
  }

  required_fields <- c(
    "fit", "profile_grid", "profile_weight", "boundary", "total_variance"
  )
  if (is.null(null_spec) || !all(required_fields %in% names(null_spec))) {
    return(failed("missing_or_incomplete_reference_null_spec"))
  }
  fit <- null_spec$fit
  if (!isTRUE(fit$available) || !is.finite(fit$observed_score) ||
      length(null_spec$profile_grid) != length(null_spec$profile_weight) ||
      !length(null_spec$profile_grid) ||
      any(!is.finite(null_spec$profile_grid)) ||
      any(!is.finite(null_spec$profile_weight)) ||
      any(null_spec$profile_weight < 0) ||
      sum(null_spec$profile_weight) <= 0) {
    return(failed("invalid_reference_null_spec"))
  }

  replicates <- cfg$refinement$replicates
  variance_seed <- deterministic_seed(
    cfg$seed, "variance_mixture", null_seed_namespace, cohort_name, gene_id
  )
  mixture_seed <- deterministic_seed(
    cfg$seed, "null_mixture", null_seed_namespace, cohort_name, gene_id
  )
  boundary_seed <- deterministic_seed(
    cfg$seed, "null_boundary_safety", null_seed_namespace, cohort_name, gene_id
  )
  variance_draw <- with_seed(
    variance_seed,
    sample(
      null_spec$profile_grid, replicates, replace = TRUE,
      prob = null_spec$profile_weight
    )
  )
  mixture_scores <- simulate_joint_scores(
    fit, variance_draw, replicates, mixture_seed
  )
  mixture_result <- monte_carlo_p(mixture_scores, fit$observed_score)
  if (is.null(mixture_result)) {
    return(failed("independent_null_mixture_simulation_failed"))
  }

  final_result <- mixture_result
  boundary_result <- NULL
  selected_branch <- "variance_profile_mixture"
  if (isTRUE(null_spec$boundary)) {
    boundary_scores <- simulate_joint_scores(
      fit, rep(null_spec$total_variance, replicates),
      replicates, boundary_seed
    )
    boundary_result <- monte_carlo_p(boundary_scores, fit$observed_score)
    if (is.null(boundary_result)) {
      return(failed("independent_boundary_safety_simulation_failed"))
    }
    if (boundary_result$p > mixture_result$p) {
      final_result <- boundary_result
      selected_branch <- "total_variance_boundary_safety"
    }
  }

  list(
    summary = tibble(
      gene_id = gene_id,
      gene_tested = TRUE,
      gene_model_p = final_result$p,
      gene_model_exceedances = final_result$exceedances,
      gene_model_replicates = final_result$replicates,
      gene_model_selected_null_branch = selected_branch,
      gene_model_null_seed_namespace = null_seed_namespace,
      gene_model_inference_reason =
        "independent_fixed_second_stage_reusing_frozen_reference_null_spec",
      refinement_reused_frozen_null_spec = TRUE
    )
  )
}

run_refinement_cohort <- function(cohort_name, included_patients) {
  family_ids <- gene_model_family %>%
    filter(cohort == cohort_name, included_in_testing_family) %>%
    pull(gene_id)
  reference <- gene_model_summary %>%
    filter(cohort == cohort_name, gene_id %in% family_ids) %>%
    select(
      gene_id,
      reference_empirical_p = gene_model_p,
      reference_empirical_bh = gene_model_q,
      reference_gene_tested = gene_tested,
      reference_inference_reason = gene_model_inference_reason,
      gene_model_observed_score
    ) %>%
    arrange(gene_id)
  candidate_table <- reference %>%
    mutate(
      refinement_candidate = reference_gene_tested %in% TRUE &
        is.finite(reference_empirical_p) &
        reference_empirical_p < cfg$refinement$trigger_p
    ) %>%
    arrange(reference_empirical_p, desc(gene_model_observed_score), gene_id)
  candidate_ids <- candidate_table$gene_id[candidate_table$refinement_candidate]
  if (cfg$refinement$maximum_genes > 0L) {
    candidate_ids <- head(candidate_ids, cfg$refinement$maximum_genes)
  }

  cache_path <- file.path(
    m04$cache, paste0("gene_null_refinement_", cohort_name, ".rds")
  )
  configuration_hash <- digest::digest(
    list(
      refinement_engine_version = "independent_fixed_stage_v0.5.0_unchanged",
      input_fingerprint = cache_input_fingerprint,
      cohort = cohort_name,
      included_patients = included_patients,
      candidate_ids = candidate_ids,
      refinement = cfg$refinement,
      gene_null = cfg$gene_null,
      thresholds = cache_thresholds
    ),
    algo = "sha256"
  )
  legacy_configuration_hashes <- vapply(
    legacy_input_fingerprints_v0_5,
    function(legacy_fingerprint) {
      digest::digest(
        list(
          input_fingerprint = legacy_fingerprint,
          cohort = cohort_name,
          included_patients = included_patients,
          candidate_ids = candidate_ids,
          refinement = cfg$refinement,
          gene_null = cfg$gene_null,
          thresholds = legacy_thresholds_v0_5,
          workflow_version = "0.5.0"
        ),
        algo = "sha256"
      )
    },
    character(1)
  )
  state <- list(configuration_hash = configuration_hash, completed = list())
  if (cfg$resume_cached_results && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (cached$configuration_hash %in%
        c(configuration_hash, legacy_configuration_hashes)) {
      state <- cached
      state$configuration_hash <- configuration_hash
      log_msg("Loaded compatible resolution-refinement cache for ", cohort_name, ".")
    }
  }
  pending <- setdiff(candidate_ids, names(state$completed))
  reference_null_specs <- gene_model_results[[cohort_name]]$null_specs %||% list()
  for (position in seq_along(pending)) {
    gene_id <- pending[[position]]
    state$completed[[gene_id]] <- simulate_independent_refinement(
      reference_null_specs[[gene_id]] %||% NULL, cohort_name, gene_id
    )
    if (position %% cfg$checkpoint_every_genes == 0L ||
        position == length(pending)) {
      atomic_save_rds(state, cache_path)
      log_msg(
        "Resolution refinement ", cohort_name, ": ",
        length(state$completed), " / ", length(candidate_ids),
        " candidate genes cached."
      )
    }
  }

  refined_rows <- if (length(state$completed)) {
    rbindlist(lapply(state$completed, `[[`, "summary"), fill = TRUE) %>%
      as_tibble() %>%
      transmute(
        gene_id,
        refinement_gene_tested = gene_tested,
        refinement_inference_reason = gene_model_inference_reason,
        refinement_independent_p_observed = gene_model_p,
        refinement_exceedances = gene_model_exceedances,
        refinement_replicates_observed = gene_model_replicates,
        refinement_selected_null_branch = gene_model_selected_null_branch,
        refinement_null_seed_namespace = gene_model_null_seed_namespace,
        refinement_reused_frozen_null_spec =
          refinement_reused_frozen_null_spec
      )
  } else {
    tibble(
      gene_id = character(), refinement_gene_tested = logical(),
      refinement_inference_reason = character(),
      refinement_independent_p_observed = numeric(),
      refinement_exceedances = integer(),
      refinement_replicates_observed = integer(),
      refinement_selected_null_branch = character(),
      refinement_null_seed_namespace = character(),
      refinement_reused_frozen_null_spec = logical()
    )
  }

  result <- candidate_table %>%
    mutate(
      refinement_selected_for_run = gene_id %in% candidate_ids,
      refinement_selection_rule = paste0(
        "reference_empirical_P<", cfg$refinement$trigger_p
      )
    ) %>%
    left_join(refined_rows, by = "gene_id") %>%
    mutate(
      refinement_available = refinement_selected_for_run &
        refinement_gene_tested %in% TRUE &
        !startsWith(
          replace_na(refinement_inference_reason, ""), "conservative_p1_"
        ),
      refinement_independent_p = case_when(
        refinement_selected_for_run ~ replace_na(
          refinement_independent_p_observed, 1
        ),
        TRUE ~ 1
      ),
      refinement_bh = p.adjust(refinement_independent_p, method = "BH"),
      refinement_bh_le_0p05 = refinement_bh <= 0.05,
      refinement_bh_le_0p10 = refinement_bh <= 0.10,
      refinement_family_treatment = case_when(
        refinement_selected_for_run & refinement_available ~
          "independent_fixed_second_stage",
        refinement_selected_for_run ~ "selected_but_conservative_P1_failure",
        refinement_candidate ~ "candidate_not_run_under_development_cap_P1",
        TRUE ~ "not_selected_reference_P_at_or_above_trigger_P1"
      )
    )
  intervals <- map2(
    result$refinement_exceedances,
    result$refinement_replicates_observed,
    safe_binomial_interval
  )
  result$refinement_mc_lower_95 <- vapply(
    intervals, `[[`, numeric(1), "lower"
  )
  result$refinement_mc_upper_95 <- vapply(
    intervals, `[[`, numeric(1), "upper"
  )
  result %>% mutate(cohort = cohort_name, .before = 1L)
}

if (cfg$gene_null$enabled && cfg$refinement$enabled) {
  refinement_results <- imap_dfr(
    analysis_cohorts,
    function(included_patients, cohort_name) {
      run_refinement_cohort(cohort_name, included_patients)
    }
  )
} else {
  refinement_results <- tibble(
    cohort = character(), gene_id = character(),
    refinement_candidate = logical(), refinement_selected_for_run = logical(),
    refinement_available = logical(), refinement_independent_p = numeric(),
    refinement_bh = numeric(), refinement_bh_le_0p05 = logical(),
    refinement_bh_le_0p10 = logical(), refinement_exceedances = integer(),
    refinement_replicates_observed = integer(),
    refinement_mc_lower_95 = numeric(), refinement_mc_upper_95 = numeric(),
    refinement_family_treatment = character(),
    refinement_selection_rule = character(),
    refinement_inference_reason = character(),
    refinement_null_seed_namespace = character(),
    refinement_reused_frozen_null_spec = logical()
  )
}

refinement_join <- refinement_results %>%
  select(
    cohort, gene_id, refinement_candidate, refinement_selected_for_run,
    refinement_available, refinement_independent_p, refinement_bh,
    refinement_bh_le_0p05, refinement_bh_le_0p10,
    refinement_exceedances, refinement_replicates_observed,
    refinement_mc_lower_95, refinement_mc_upper_95,
    refinement_family_treatment, refinement_selection_rule,
    refinement_inference_reason, refinement_null_seed_namespace,
    refinement_reused_frozen_null_spec
  )
gene_model_summary <- gene_model_summary %>%
  left_join(refinement_join, by = c("cohort", "gene_id"))

method_long <- bind_rows(
  gene_model_summary %>% transmute(
    cohort, gene_id, gene_model_observed_score,
    method = "REML empirical (20k)", p_value = gene_model_p,
    bh_value = gene_model_q, method_available = gene_tested %in% TRUE
  ),
  gene_model_summary %>% transmute(
    cohort, gene_id, gene_model_observed_score,
    method = "REML independent refinement",
    p_value = refinement_independent_p,
    bh_value = refinement_bh,
    method_available = refinement_available %in% TRUE
  ),
  gene_model_summary %>% transmute(
    cohort, gene_id, gene_model_observed_score,
    method = "REML GPD hybrid", p_value = gene_model_gpd_hybrid_p,
    bh_value = gene_model_gpd_hybrid_q,
    method_available = gene_tested %in% TRUE &
      is.finite(gene_model_gpd_hybrid_p)
  )
) %>%
  mutate(
    raw_p_lt_0p05 = is.finite(p_value) & p_value < 0.05,
    bh_le_0p05 = is.finite(bh_value) & bh_value <= 0.05,
    bh_le_0p10 = is.finite(bh_value) & bh_value <= 0.10,
    passes_effect_gate = gene_model_observed_score >= 1
  )

method_summary <- method_long %>%
  group_by(cohort, method) %>%
  summarise(
    frozen_family = n(),
    available_results = sum(method_available),
    availability_fraction = mean(method_available),
    raw_p_lt_0p05 = sum(raw_p_lt_0p05),
    bh_le_0p05 = sum(bh_le_0p05),
    bh_le_0p10 = sum(bh_le_0p10),
    bh_le_0p05_and_effect_gate = sum(bh_le_0p05 & passes_effect_gate),
    bh_le_0p10_and_effect_gate = sum(bh_le_0p10 & passes_effect_gate),
    minimum_p = if (any(method_available)) {
      min(p_value[method_available], na.rm = TRUE)
    } else NA_real_,
    minimum_bh = if (any(method_available)) {
      min(bh_value[method_available], na.rm = TRUE)
    } else NA_real_,
    .groups = "drop"
  )

method_names <- unique(method_long$method)
method_pairs <- if (length(method_names) >= 2L) {
  combn(method_names, 2L, simplify = FALSE)
} else list()
method_concordance <- map_dfr(names(analysis_cohorts), function(cohort_name) {
  map_dfr(method_pairs, function(method_pair) {
    left <- method_long %>%
      filter(
        cohort == cohort_name, method == method_pair[[1]], method_available
      ) %>%
      select(gene_id, p_left = p_value, q_left = bh_value)
    right <- method_long %>%
      filter(
        cohort == cohort_name, method == method_pair[[2]], method_available
      ) %>%
      select(gene_id, p_right = p_value, q_right = bh_value)
    joined <- inner_join(left, right, by = "gene_id") %>%
      filter(is.finite(p_left), is.finite(p_right))
    tibble(
      cohort = cohort_name,
      method_1 = method_pair[[1]],
      method_2 = method_pair[[2]],
      n_compared = nrow(joined),
      spearman_raw_p = if (nrow(joined) >= 3L) {
        cor(joined$p_left, joined$p_right, method = "spearman")
      } else NA_real_,
      overlap_bh_0p05 = sum(joined$q_left <= 0.05 & joined$q_right <= 0.05),
      union_bh_0p05 = sum(joined$q_left <= 0.05 | joined$q_right <= 0.05),
      overlap_bh_0p10 = sum(joined$q_left <= 0.10 & joined$q_right <= 0.10),
      union_bh_0p10 = sum(joined$q_left <= 0.10 | joined$q_right <= 0.10)
    )
  })
})

tail_method_benchmark <- gene_model_summary %>%
  filter(refinement_selected_for_run %in% TRUE,
         refinement_available %in% TRUE,
         is.finite(refinement_independent_p), refinement_independent_p > 0) %>%
  transmute(
    cohort, gene_id,
    reference_empirical_p = gene_model_p,
    gpd_hybrid_p = gene_model_gpd_hybrid_p,
    independent_refinement_p = refinement_independent_p,
    empirical_absolute_log10_error = abs(
      log10(reference_empirical_p) - log10(independent_refinement_p)
    ),
    gpd_absolute_log10_error = abs(
      log10(gpd_hybrid_p) - log10(independent_refinement_p)
    ),
    gpd_closer_to_independent_refinement =
      gpd_absolute_log10_error < empirical_absolute_log10_error
  )

tail_method_benchmark_summary <- tail_method_benchmark %>%
  group_by(cohort) %>%
  summarise(
    refined_candidates = n(),
    median_empirical_absolute_log10_error =
      median(empirical_absolute_log10_error),
    median_gpd_absolute_log10_error = median(gpd_absolute_log10_error),
    fraction_GPD_closer_to_refinement =
      mean(gpd_closer_to_independent_refinement),
    .groups = "drop"
  )

refinement_qc <- refinement_results %>%
  group_by(cohort) %>%
  summarise(
    frozen_family = n(),
    reference_p_below_trigger = sum(refinement_candidate),
    selected_for_refinement = sum(refinement_selected_for_run),
    completed_refinements = sum(refinement_available),
    conservative_failures = sum(
      refinement_selected_for_run & !refinement_available
    ),
    independent_replicates = cfg$refinement$replicates,
    independent_minimum_p = 1 / (cfg$refinement$replicates + 1),
    selection_and_refinement_draws_separate = all(
      refinement_null_seed_namespace[refinement_selected_for_run] ==
        "independent_resolution_refinement",
      na.rm = TRUE
    ),
    frozen_reference_null_reused = all(
      refinement_reused_frozen_null_spec[refinement_selected_for_run] %in% TRUE
    ),
    pass = selected_for_refinement == completed_refinements &
      selection_and_refinement_draws_separate &
      frozen_reference_null_reused,
    .groups = "drop"
  )

gpd_qc <- gene_model_summary %>%
  group_by(cohort) %>%
  summarise(
    frozen_family = n(),
    empirical_low_exceedance_genes = sum(
      gene_model_exceedances < cfg$gpd$empirical_exceedance_trigger,
      na.rm = TRUE
    ),
    stable_GPD_tail_genes = sum(
      gene_model_gpd_threshold_stable %in% TRUE, na.rm = TRUE
    ),
    GPD_tail_used = sum(
      gene_model_gpd_hybrid_method == "GPD_tail", na.rm = TRUE
    ),
    empirical_fallbacks = sum(
      str_detect(
        replace_na(gene_model_gpd_hybrid_method, ""),
        "fallback"
      ),
      na.rm = TRUE
    ),
    finite_hybrid_p_values = sum(is.finite(gene_model_gpd_hybrid_p)),
    pass = finite_hybrid_p_values == frozen_family,
    .groups = "drop"
  )

method_lock_status <- tibble(
  method = c(
    "REML empirical (20k)", "REML independent refinement",
    "REML GPD hybrid", "Patient multinomial LRT"
  ),
  statistical_role = c(
    "locked reference cohort test for paired primary-recurrent remodeling",
    "locked resolution stage for genes with reference P < 0.05",
    "locked diagnostic tail estimator for the same REML null",
    "locked patient-specific technical count-composition test"
  ),
  enabled = c(
    cfg$gene_null$enabled,
    cfg$gene_null$enabled && cfg$refinement$enabled,
    cfg$gene_null$enabled && cfg$gpd$enabled,
    cfg$patient_test$enabled
  ),
  method_locked = TRUE,
  gene_set_selection_authorized = c(
    FALSE, FALSE, FALSE, TRUE
  ),
  authorized_selection_rule = c(
    "none; cohort P and BH annotate but do not define A-E/F/U/NQ or downstream sets",
    "none; refinement provides resolution and does not define a downstream set",
    "none; diagnostic only",
    paste0(
      "patient_test_BH<=", cfg$patient_test$primary_fdr,
      " AND patient effect gate; patient-specific analyses only"
    )
  ),
  raw_p_selection_authorized = FALSE
)

gene_inference_method_comparison <- if (nrow(gene_model_summary)) {
  gene_model_summary %>%
    select(
    cohort, gene_id, gene_model_observed_score,
    gene_model_p, gene_model_q, gene_model_bh_le_0p05,
    gene_model_bh_le_0p10, gene_model_inference_reason,
    refinement_candidate, refinement_selected_for_run, refinement_available,
    refinement_independent_p, refinement_bh, refinement_bh_le_0p05,
    refinement_bh_le_0p10, refinement_family_treatment,
    refinement_reused_frozen_null_spec,
    gene_model_gpd_raw_p, gene_model_gpd_hybrid_p,
    gene_model_gpd_hybrid_q, gene_model_gpd_bh_le_0p05,
    gene_model_gpd_bh_le_0p10, gene_model_gpd_hybrid_method,
    gene_model_gpd_threshold_stable, gene_model_gpd_reason
    ) %>%
    mutate(
      inference_method_locked = TRUE,
      downstream_selection_authorized = FALSE,
      raw_p_is_descriptive_not_a_selection_rule = TRUE
    )
} else {
  tibble(
    cohort = character(), gene_id = character(),
    inference_method_locked = logical(),
    downstream_selection_authorized = logical(),
    raw_p_is_descriptive_not_a_selection_rule = logical()
  )
}

write_tsv(refinement_results, file.path(m04$tables, "gene_model_resolution_refinement.tsv.gz"))
write_tsv(refinement_qc, file.path(m04$qc, "gene_model_resolution_refinement_qc.tsv"))
write_tsv(method_long, file.path(m04$tables, "gene_inference_methods_long.tsv.gz"))
write_tsv(
  gene_inference_method_comparison,
  file.path(m04$tables, "gene_inference_method_comparison.tsv.gz")
)
write_tsv(method_summary, file.path(m04$tables, "gene_inference_method_summary.tsv"))
write_tsv(method_concordance, file.path(m04$tables, "gene_inference_method_concordance.tsv"))
write_tsv(tail_method_benchmark, file.path(m04$tables, "tail_method_independent_benchmark.tsv.gz"))
write_tsv(tail_method_benchmark_summary, file.path(m04$tables, "tail_method_independent_benchmark_summary.tsv"))
write_tsv(gpd_qc, file.path(m04$qc, "gene_model_gpd_qc.tsv"))
write_tsv(method_lock_status, file.path(m04$tables, "inference_method_lock_status.tsv"))
write_tsv(gene_model_summary, file.path(m04$tables, "gene_model_results.tsv.gz"))


# ---- m04-truth-known-validation (source-derived computation) ----
log_msg("Remodeling step: truth-known-validation.")
draw_simplex <- function(k, concentration = 1) {
  value <- rgamma(k, shape = concentration, rate = 1)
  value / sum(value)
}

simulate_validation_object <- function(
    scenario_id, replicate_id, k, depth, biological_sd, effect_size
) {
  seed <- deterministic_seed(
    cfg$seed, "truth_known_data", scenario_id, replicate_id
  )
  with_seed(seed, {
    basis <- ilr_basis(k)
    dimension <- k - 1L
    base_probability <- draw_simplex(k, concentration = cfg$validation_public$dirichlet_concentration %||% 1.5)
    base_coordinate <- as.numeric(ilr_transform(
      matrix(base_probability, nrow = 1L), basis
    ))
    true_midpoint <- matrix(
      rep(base_coordinate, each = length(patient_ids)),
      nrow = length(patient_ids), byrow = FALSE
    ) + matrix(
      rnorm(length(patient_ids) * dimension, sd = cfg$validation_public$midpoint_ilr_sd %||% .20),
      nrow = length(patient_ids), ncol = dimension
    )
    fixed_difference <- rep(0, dimension)
    fixed_difference[[1]] <- effect_size
    true_difference <- matrix(
      rnorm(
        length(patient_ids) * dimension,
        mean = rep(fixed_difference, each = length(patient_ids)),
        sd = biological_sd
      ),
      nrow = length(patient_ids), ncol = dimension
    )
    primary_probability <- ilr_inverse(true_midpoint - 0.5 * true_difference, basis)
    recurrent_probability <- ilr_inverse(true_midpoint + 0.5 * true_difference, basis)
    primary_counts <- t(vapply(seq_along(patient_ids), function(index) {
      as.numeric(rmultinom(1, depth, primary_probability[index, ]))
    }, numeric(k)))
    recurrent_counts <- t(vapply(seq_along(patient_ids), function(index) {
      as.numeric(rmultinom(1, depth, recurrent_probability[index, ]))
    }, numeric(k)))
    counts <- matrix(
      0L, nrow = nrow(sample_manifest), ncol = k,
      dimnames = list(sample_manifest$sample, paste0("synthetic_isoform_", seq_len(k)))
    )
    counts[sample_by_patient$primary, ] <- primary_counts
    counts[sample_by_patient$recurrent, ] <- recurrent_counts
    storage.mode(counts) <- "integer"
    object <- list(
      gene_id = paste0("validation_", scenario_id, "_rep", replicate_id),
      transcript_ids = colnames(counts),
      counts = counts,
      all_counts = counts,
      pair_quantifiable = setNames(rep(TRUE, length(patient_ids)), patient_ids),
      n_isoforms_total = k,
      n_isoforms_retained = k,
      pooled_gene_count = sum(counts),
      pooled_retained_count = sum(counts)
    )
    list(
      object = object,
      true_midpoint = true_midpoint,
      true_difference = true_difference,
      true_variance = biological_sd^2,
      fixed_difference = fixed_difference,
      seed = seed
    )
  })
}

validation_design <- tidyr::crossing(
  isoform_number = as.integer(cfg$validation_public$isoform_counts %||% c(2L,4L)),
  per_sample_depth = as.integer(cfg$validation_public$sample_depths %||% c(50L,500L)),
  biological_sd = as.numeric(cfg$validation_public$biological_ilr_sd %||% c(.05,.25)),
  truth = c("null", "alternative"),
  validation_replicate = seq_len(cfg$validation$replicates_per_scenario)
) %>%
  mutate(
    effect_size = if_else(truth == "null", as.numeric(cfg$validation_public$null_mean_shift %||% 0), as.numeric(cfg$validation_public$alternative_mean_shift %||% .60)),
    scenario_id = paste0(
      "K", isoform_number, "_N", per_sample_depth,
      "_SD", str_replace(format(biological_sd, trim = TRUE), "\\.", "p"),
      "_", truth
    )
  )

run_validation_row <- function(row) {
  synthetic <- simulate_validation_object(
    row$scenario_id, row$validation_replicate, row$isoform_number,
    row$per_sample_depth, row$biological_sd, row$effect_size
  )
  observed <- analyze_gene_observed(
    synthetic$object, patient_ids, "truth_known_validation"
  )$gene
  fit <- prepare_gene_model(synthetic$object, observed, patient_ids)
  base <- row %>% transmute(
    scenario_id, validation_replicate, truth, isoform_number,
    per_sample_depth, biological_sd, true_variance = biological_sd^2,
    effect_size, simulated_gene_id = synthetic$object$gene_id,
    observed_jsd = observed$jsd,
    observed_max_abs_delta_if = observed$max_abs_delta_if,
    observed_joint_score = observed$joint_effect_score,
    effect_gate_pass = observed$effect_pass,
    data_seed = synthetic$seed
  )
  if (!isTRUE(fit$available)) {
    return(tidyr::crossing(base, method = c(
      "full_oracle", "oracle_variance", "production_REML",
      "production_REML_GPD_hybrid"
    )) %>% mutate(
      p_value = 1, called_at_0p05 = FALSE,
      inference_status = paste0("conservative_p1_", fit$reason)
    ))
  }

  full_oracle_fit <- fit
  full_oracle_fit$midpoint <- synthetic$true_midpoint
  oracle_variance <- rep(synthetic$true_variance, cfg$validation$null_replicates)
  full_oracle_scores <- simulate_joint_scores(
    full_oracle_fit, oracle_variance, cfg$validation$null_replicates,
    deterministic_seed(
      cfg$seed, "truth_full_oracle", row$scenario_id,
      row$validation_replicate
    )
  )
  variance_oracle_scores <- simulate_joint_scores(
    fit, oracle_variance, cfg$validation$null_replicates,
    deterministic_seed(
      cfg$seed, "truth_variance_oracle", row$scenario_id,
      row$validation_replicate
    )
  )
  full_oracle <- monte_carlo_p(full_oracle_scores, fit$observed_score)
  variance_oracle <- monte_carlo_p(variance_oracle_scores, fit$observed_score)
  production <- fit_gene_model(
    synthetic$object, observed, patient_ids, "truth_known_validation",
    null_replicates = cfg$validation$null_replicates
  )$summary

  bind_rows(
    base %>% mutate(
      method = "full_oracle",
      p_value = full_oracle$p %||% 1,
      inference_status = if (is.null(full_oracle)) "conservative_p1_failure" else "available"
    ),
    base %>% mutate(
      method = "oracle_variance",
      p_value = variance_oracle$p %||% 1,
      inference_status = if (is.null(variance_oracle)) "conservative_p1_failure" else "available"
    ),
    base %>% mutate(
      method = "production_REML",
      p_value = production$gene_model_p,
      inference_status = production$gene_model_inference_reason
    ),
    base %>% mutate(
      method = "production_REML_GPD_hybrid",
      p_value = production$gene_model_gpd_hybrid_p,
      inference_status = case_when(
        !is.finite(production$gene_model_gpd_hybrid_p) ~
          "conservative_p1_GPD_unavailable",
        startsWith(
          coalesce(production$gene_model_gpd_hybrid_method, ""),
          "conservative_p1"
        ) ~ coalesce(
          production$gene_model_gpd_hybrid_method,
          "conservative_p1_GPD_unavailable"
        ),
        TRUE ~ paste0(
          "available:", production$gene_model_gpd_hybrid_method
        )
      )
    )
  ) %>% mutate(called_at_0p05 = p_value <= 0.05)
}

validation_cache_path <- file.path(m04$cache, "truth_known_validation.rds")
validation_hash <- digest::digest(
  list(
    validation_status_semantics = "gpd_failure_availability_v2",
    input_fingerprint = cache_input_fingerprint,
    validation_design = validation_design,
    validation = cfg$validation,
    gene_null = cfg$gene_null,
    gpd = cfg$gpd,
    thresholds = cache_thresholds
  ), algo = "sha256"
)
legacy_validation_hashes <- unique(c(
  vapply(
    legacy_input_fingerprints_v0_5,
    function(legacy_fingerprint) {
      digest::digest(
        list(
          input_fingerprint = legacy_fingerprint,
          validation_design = validation_design,
          validation = cfg$validation,
          gene_null = cfg$gene_null,
          gpd = cfg$gpd,
          thresholds = legacy_thresholds_v0_5
        ),
        algo = "sha256"
      )
    },
    character(1)
  ),
  vapply(
    legacy_input_fingerprints_v0_5,
    function(legacy_fingerprint) {
      digest::digest(
        list(
          validation_status_semantics = "gpd_failure_availability_v2",
          input_fingerprint = legacy_fingerprint,
          validation_design = validation_design,
          validation = cfg$validation,
          gene_null = cfg$gene_null,
          gpd = cfg$gpd,
          thresholds = legacy_thresholds_v0_5
        ),
        algo = "sha256"
      )
    },
    character(1)
  )
))

if (cfg$gene_null$enabled && cfg$validation$enabled) {
  validation_state <- list(configuration_hash = validation_hash, completed = list())
  if (cfg$resume_cached_results && file.exists(validation_cache_path)) {
    candidate <- readRDS(validation_cache_path)
    if (candidate$configuration_hash %in%
        c(validation_hash, legacy_validation_hashes)) {
      validation_state <- candidate
      validation_state$configuration_hash <- validation_hash
      log_msg("Loaded compatible truth-known validation cache.")
    }
  }
  validation_keys <- paste(
    validation_design$scenario_id, validation_design$validation_replicate,
    sep = "__"
  )
  pending_index <- which(!validation_keys %in% names(validation_state$completed))
  for (position in seq_along(pending_index)) {
    index <- pending_index[[position]]
    key <- validation_keys[[index]]
    validation_state$completed[[key]] <- run_validation_row(
      validation_design[index, , drop = FALSE]
    )
    if (position %% 10L == 0L || position == length(pending_index)) {
      atomic_save_rds(validation_state, validation_cache_path)
      log_msg(
        "Truth-known validation: ", length(validation_state$completed), " / ",
        nrow(validation_design), " synthetic datasets cached."
      )
    }
  }
  truth_known_validation <- bind_rows(validation_state$completed) %>%
    mutate(
      inference_status = case_when(
        method == "production_REML_GPD_hybrid" &
          startsWith(coalesce(inference_status, ""), "available:conservative_p1") ~
            str_remove(inference_status, "^available:"),
        TRUE ~ inference_status
      )
    )
} else {
  truth_known_validation <- tibble(
    scenario_id = character(), validation_replicate = integer(), truth = character(),
    isoform_number = integer(), per_sample_depth = integer(), biological_sd = numeric(),
    true_variance = numeric(), effect_size = numeric(), simulated_gene_id = character(),
    observed_jsd = numeric(), observed_max_abs_delta_if = numeric(),
    observed_joint_score = numeric(), effect_gate_pass = logical(), data_seed = integer(),
    method = character(), p_value = numeric(), inference_status = character(),
    called_at_0p05 = logical()
  )
}

truth_known_validation <- truth_known_validation %>%
  mutate(
    inference_available = !is.na(inference_status) &
      !startsWith(inference_status, "conservative_p1_")
  )

safe_binomial_bound <- function(successes, trials, index) {
  if (!is.finite(trials) || trials < 1L) return(NA_real_)
  binom.test(as.integer(successes), as.integer(trials))$conf.int[[index]]
}

validation_summary <- truth_known_validation %>%
  group_by(method, truth, isoform_number, per_sample_depth, biological_sd) %>%
  summarise(
    simulated_datasets = n(),
    available_inferences = sum(inference_available),
    conservative_failures = n() - sum(inference_available),
    availability_fraction = mean(inference_available),
    calls_in_all_datasets = sum(called_at_0p05),
    call_rate_in_all_datasets = mean(called_at_0p05),
    calls_among_available = sum(called_at_0p05 & inference_available),
    call_rate_among_available = if_else(
      sum(inference_available) > 0L,
      sum(called_at_0p05 & inference_available) / sum(inference_available),
      NA_real_
    ),
    median_p_value_among_available = if (any(inference_available)) {
      median(p_value[inference_available])
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    estimand = if_else(truth == "null", "type_I_error", "power"),
    caution = "Simulation strata are diagnostic; aggregate calibration is the primary summary."
  )

validation_aggregate <- truth_known_validation %>%
  group_by(method, truth) %>%
  summarise(
    simulated_datasets = n(),
    available_inferences = sum(inference_available),
    conservative_failures = n() - sum(inference_available),
    availability_fraction = mean(inference_available),
    calls_in_all_datasets = sum(called_at_0p05),
    call_rate_in_all_datasets = mean(called_at_0p05),
    calls_among_available = sum(called_at_0p05 & inference_available),
    call_rate_among_available = if_else(
      sum(inference_available) > 0L,
      sum(called_at_0p05 & inference_available) / sum(inference_available),
      NA_real_
    ),
    binomial_lower_95_among_available = safe_binomial_bound(
      sum(called_at_0p05 & inference_available), sum(inference_available), 1L
    ),
    binomial_upper_95_among_available = safe_binomial_bound(
      sum(called_at_0p05 & inference_available), sum(inference_available), 2L
    ),
    .groups = "drop"
  ) %>%
  mutate(estimand = if_else(truth == "null", "type_I_error", "power"))

validation_active <- cfg$gene_null$enabled && cfg$validation$enabled
production_validation_aggregate <- validation_aggregate %>%
  filter(method == "production_REML")
production_validation_strata <- validation_summary %>%
  filter(method == "production_REML")

production_truth_availability_pass <- !validation_active ||
  (setequal(production_validation_aggregate$truth, c("null", "alternative")) &&
     all(
       production_validation_aggregate$availability_fraction >=
         cfg$validation$minimum_production_availability
     ))
production_stratum_availability_pass <- !validation_active ||
  (nrow(production_validation_strata) > 0L &&
     all(
       production_validation_strata$availability_fraction >=
         cfg$validation$minimum_stratum_availability
     ))
production_null_row <- production_validation_aggregate %>%
  filter(truth == "null")
production_null_calibration_pass <- !validation_active ||
  (nrow(production_null_row) == 1L &&
     is.finite(production_null_row$call_rate_among_available[[1]]) &&
     production_null_row$call_rate_among_available[[1]] <=
       cfg$validation$maximum_null_call_rate)

validation_qc <- tibble(
  check = c(
    "validation_requested", "all_design_rows_completed",
    "four_methods_per_dataset", "finite_unit_interval_p_values",
    "null_and_alternative_present",
    "production_REML_truth_level_availability",
    "production_REML_design_stratum_availability",
    "production_REML_null_calibration_among_available_inferences"
  ),
  pass = c(
    TRUE,
    !validation_active ||
      n_distinct(truth_known_validation$simulated_gene_id) == nrow(validation_design),
    !validation_active ||
      nrow(truth_known_validation) == 4L * nrow(validation_design),
    !validation_active || all(
      is.finite(truth_known_validation$p_value) &
        truth_known_validation$p_value >= 0 & truth_known_validation$p_value <= 1
    ),
    !validation_active || setequal(unique(truth_known_validation$truth),
                                        c("null", "alternative")),
    production_truth_availability_pass,
    production_stratum_availability_pass,
    production_null_calibration_pass
  ),
  blocking = rep(TRUE, 8L)
)

write_tsv(validation_design, file.path(m04$tables, "truth_known_validation_design.tsv"))
write_tsv(truth_known_validation, file.path(m04$tables, "truth_known_validation_results.tsv.gz"))
write_tsv(validation_summary, file.path(m04$tables, "truth_known_validation_by_stratum.tsv"))
write_tsv(validation_aggregate, file.path(m04$tables, "truth_known_validation_aggregate.tsv"))
write_tsv(validation_qc, file.path(m04$qc, "truth_known_validation_qc.tsv"))

} else {
# Disabled inferential branches have explicit empty results, never fabricated P values.
gene_model_results <- list()
gene_model_family <- tibble(cohort=character(),gene_id=character(),
  included_in_testing_family=logical(),testing_family_reason=character(),
  scheduled_for_computation=logical(),testing_family_rule=character())
gene_model_summary <- tibble(cohort=character(),gene_id=character(),gene_tested=logical(),
  gene_model_p=numeric(),gene_model_q=numeric(),gene_model_supported=logical(),
  gene_model_effect_supported=logical(),gene_model_inference_reason=character(),
  gene_model_bh_le_0p05=logical(),gene_model_bh_le_0p10=logical(),
  refinement_available=logical(),refinement_independent_p=numeric(),refinement_bh=numeric(),
  refinement_bh_le_0p05=logical(),refinement_bh_le_0p10=logical())
gene_model_covariance <- tibble(cohort=character(),gene_id=character(),covariance_available=logical())
refinement_results <- tibble(cohort=character(),gene_id=character(),refinement_available=logical())
method_long <- tibble(cohort=character(),gene_id=character(),method=character(),p_value=numeric(),bh_value=numeric())
truth_known_validation <- tibble(scenario_id=character(),method=character(),p_value=numeric())
validation_aggregate <- tibble(method=character(),truth=character(),call_rate_among_available=numeric())
gene_model_execution_gate <- tibble(check=character(),pass=logical())
refinement_qc <- tibble(cohort=character(),pass=logical())
gpd_qc <- tibble(cohort=character(),pass=logical())
validation_qc <- tibble(check=character(),pass=logical(),blocking=logical())

}

# ---- m04-authoritative-table (source-derived computation) ----
log_msg("Remodeling step: authoritative-table.")
# Add typed absent optional-model fields before source-derived overlap auditing.
for (nm in c("gene_model_bh_le_0p05", "gene_model_bh_le_0p10", "refinement_available", "refinement_bh_le_0p05", "refinement_bh_le_0p10")) {
  if (!nm %in% names(gene_model_summary)) gene_model_summary[[nm]] <- rep(NA, nrow(gene_model_summary))
}
for (nm in c("refinement_independent_p", "refinement_bh")) {
  if (!nm %in% names(gene_model_summary)) gene_model_summary[[nm]] <- rep(NA_real_, nrow(gene_model_summary))
}
if (!all(c("cohort", "gene_id") %in% names(gene_model_family))) {
  gene_model_family <- tibble(
    cohort = character(), gene_id = character(), included_in_testing_family = logical(),
    testing_family_reason = character(), scheduled_for_computation = logical(),
    testing_family_rule = character()
  )
}
if (!all(c("cohort", "gene_id") %in% names(gene_model_summary))) {
  gene_model_summary <- tibble(
    cohort = character(), gene_id = character(), gene_tested = logical(),
    gene_test_reason = character(), gene_model_p = numeric(), gene_model_q = numeric(),
    gene_model_supported = logical(), gene_model_effect_supported = logical(),
    gene_model_inference_reason = character()
  )
}

problematic_null_strata <- null_stratum_rate %>%
  filter(
    included_in_threshold_gate,
    is.finite(false_effect_rate),
    false_effect_rate > cfg$empirical_null$stratum_target
  ) %>%
  arrange(desc(false_effect_rate), source_timepoint,
          gene_depth_bin, isoform_complexity_bin)

complexity_order <- complexity_bin_labels
problematic_key_table <- problematic_null_strata %>%
  distinct(gene_depth_bin, isoform_complexity_bin)
adjacent_problematic_strata <- if (nrow(problematic_key_table)) {
  map2_dfr(
    as.character(problematic_key_table$gene_depth_bin),
    as.character(problematic_key_table$isoform_complexity_bin),
    function(depth_bin, complexity_bin) {
      complexity_index <- match(complexity_bin, complexity_order)
      neighbor_index <- c(complexity_index - 1L, complexity_index + 1L)
      neighbor_index <- neighbor_index[
        neighbor_index >= 1L & neighbor_index <= length(complexity_order)
      ]
      tibble(
        gene_depth_bin = depth_bin,
        isoform_complexity_bin = complexity_order[neighbor_index]
      )
    }
  ) %>%
    distinct()
} else {
  tibble(
    gene_depth_bin = character(),
    isoform_complexity_bin = character()
  )
}

primary_class_for_risk <- remodeling_classes_interim %>%
  filter(cohort == primary_cohort_name) %>%
  transmute(
    gene_id,
    architecture_class = as.character(architecture_class),
    is_remodeled_AE,
    effect_stable,
    driver_patient_support_n,
    is_high_confidence_remodeled,
    is_high_confidence_class,
    jsd,
    max_abs_delta_if,
    delta_entropy
  )
primary_model_for_risk <- gene_model_summary %>%
  filter(cohort == primary_cohort_name) %>%
  select(
    gene_id, gene_model_p, gene_model_q, gene_model_bh_le_0p05,
    gene_model_bh_le_0p10, refinement_available,
    refinement_independent_p, refinement_bh,
    refinement_bh_le_0p05, refinement_bh_le_0p10
  )
primary_family_for_risk <- gene_model_family %>%
  filter(cohort == primary_cohort_name) %>%
  select(gene_id, included_in_testing_family)
patient_risk_gene_summary <- if (nrow(patient_specific_inference)) {
  patient_specific_inference %>%
    group_by(gene_id) %>%
    summarise(
      patient_test_family_any = any(patient_test_in_family),
      patient_BH05_and_effect_any = any(patient_supported_bh05_effect),
      patient_BH10_and_effect_any = any(patient_supported_bh10_effect),
      .groups = "drop"
    )
} else {
  tibble(
    gene_id = character(), patient_test_family_any = logical(),
    patient_BH05_and_effect_any = logical(),
    patient_BH10_and_effect_any = logical()
  )
}

null_risk_observed_gene_audit <- gene_quantifiability %>%
  left_join(primary_class_for_risk, by = "gene_id") %>%
  left_join(primary_family_for_risk, by = "gene_id") %>%
  left_join(primary_model_for_risk, by = "gene_id") %>%
  left_join(patient_risk_gene_summary, by = "gene_id") %>%
  mutate(
    depth_complexity_key = paste(
      as.character(gene_depth_bin),
      as.character(isoform_complexity_bin),
      sep = "|"
    ),
    exact_problematic_stratum = depth_complexity_key %in% paste(
      as.character(problematic_null_strata$gene_depth_bin),
      as.character(problematic_null_strata$isoform_complexity_bin),
      sep = "|"
    ),
    adjacent_problematic_stratum = depth_complexity_key %in% paste(
      as.character(adjacent_problematic_strata$gene_depth_bin),
      as.character(adjacent_problematic_strata$isoform_complexity_bin),
      sep = "|"
    ),
    null_risk_group = case_when(
      exact_problematic_stratum ~ "exact_flagged_half_depth_stratum",
      adjacent_problematic_stratum ~ "same_depth_adjacent_complexity",
      TRUE ~ "other_strata"
    ),
    eligible_for_cohort_monte_carlo =
      replace_na(included_in_testing_family, FALSE),
    excluded_from_cohort_monte_carlo_by_complexity = quantifiable &
      n_isoforms_retained > cfg$gene_null$maximum_isoforms,
    cohort_reference_p_lt_0p05 = is.finite(gene_model_p) &
      gene_model_p < 0.05,
    cohort_refined_p_lt_0p05 = refinement_available %in% TRUE &
      is.finite(refinement_independent_p) & refinement_independent_p < 0.05,
    exact_overlap_with_reference_p_lt_0p05 = exact_problematic_stratum &
      cohort_reference_p_lt_0p05,
    exact_overlap_with_refined_p_lt_0p05 = exact_problematic_stratum &
      cohort_refined_p_lt_0p05
  )

null_risk_overlap_summary <- null_risk_observed_gene_audit %>%
  group_by(null_risk_group) %>%
  summarise(
    genes = n(),
    quantifiable_genes = sum(quantifiable),
    observed_A_to_E = sum(is_remodeled_AE %in% TRUE),
    effect_stable_A_to_E = sum(is_remodeled_AE %in% TRUE &
      effect_stable %in% TRUE),
    primary_downstream_genes = sum(is_high_confidence_remodeled %in% TRUE),
    strict_class_downstream_genes = sum(is_high_confidence_class %in% TRUE),
    cohort_monte_carlo_family = sum(eligible_for_cohort_monte_carlo),
    cohort_reference_p_lt_0p05 = sum(cohort_reference_p_lt_0p05),
    cohort_refined_p_lt_0p05 = sum(cohort_refined_p_lt_0p05),
    patient_test_family_any = sum(patient_test_family_any %in% TRUE),
    patient_BH05_and_effect_any = sum(patient_BH05_and_effect_any %in% TRUE),
    patient_BH10_and_effect_any = sum(patient_BH10_and_effect_any %in% TRUE),
    .groups = "drop"
  )

null_risk_primary_set_overlap_n <- sum(
  null_risk_observed_gene_audit$exact_problematic_stratum &
    null_risk_observed_gene_audit$is_high_confidence_remodeled %in% TRUE
)
null_risk_reference_p_overlap_n <- sum(
  null_risk_observed_gene_audit$exact_overlap_with_reference_p_lt_0p05
)
null_risk_refined_p_overlap_n <- sum(
  null_risk_observed_gene_audit$exact_overlap_with_refined_p_lt_0p05
)
null_risk_primary_set_pass <- nrow(problematic_null_strata) == 0L ||
  null_risk_primary_set_overlap_n == 0L
threshold_gate_pass <- threshold_core_gate_pass && null_risk_primary_set_pass

threshold_gate <- bind_rows(
  threshold_gate,
  tibble(
    check = "flagged half-depth stratum absent from primary downstream set",
    observed = as.numeric(null_risk_primary_set_overlap_n),
    target = 0,
    pass = null_risk_primary_set_pass,
    blocking = TRUE,
    interpretation = paste(
      "A split-null stratum excess is nonblocking only when no gene in the locked",
      "bootstrap-stable, >=4-patient downstream set belongs to that exact stratum"
    )
  )
)

null_risk_overlap_interpretation <- tibble(
  problematic_strata_detected = nrow(problematic_null_strata),
  worst_false_effect_rate = if (nrow(problematic_null_strata)) {
    max(problematic_null_strata$false_effect_rate)
  } else NA_real_,
  exact_problematic_observed_genes = sum(
    null_risk_observed_gene_audit$exact_problematic_stratum
  ),
  exact_problematic_primary_downstream_genes =
    null_risk_primary_set_overlap_n,
  exact_overlap_reference_p_lt_0p05 = null_risk_reference_p_overlap_n,
  exact_overlap_refined_p_lt_0p05 = null_risk_refined_p_overlap_n,
  exact_p_overlap_structurally_possible = any(
    null_risk_observed_gene_audit$exact_problematic_stratum &
      null_risk_observed_gene_audit$eligible_for_cohort_monte_carlo
  ),
  exact_patient_test_overlap_structurally_possible = any(
    null_risk_observed_gene_audit$exact_problematic_stratum &
      null_risk_observed_gene_audit$patient_test_family_any %in% TRUE
  ),
  conclusion = case_when(
    nrow(problematic_null_strata) == 0L ~
      "no adequately populated split-null stratum exceeded the registered stress target",
    null_risk_primary_set_overlap_n > 0L ~
      "flagged stratum reaches the primary downstream set; completion is blocked",
    TRUE ~ paste(
      "flagged half-depth simulation stratum is absent from the primary downstream set;",
      "the excess remains a documented nonblocking stress warning"
    )
  )
)

patient_gene_inference_summary <- if (nrow(patient_specific_inference)) {
  patient_specific_inference %>%
    group_by(gene_id) %>%
    summarise(
      patient_tests_available_n = sum(patient_test_available),
      patient_nominal_p_lt_0p05_n = sum(patient_nominal_p_lt_0p05),
      patient_supported_bh05_effect_n = sum(patient_supported_bh05_effect),
      patient_supported_bh10_effect_n = sum(patient_supported_bh10_effect),
      patient_supported_bh05_ids = paste(
        patient[patient_supported_bh05_effect], collapse = ","
      ),
      patient_supported_bh10_ids = paste(
        patient[patient_supported_bh10_effect], collapse = ","
      ),
      .groups = "drop"
    )
} else {
  tibble(
    gene_id = character(),
    patient_tests_available_n = integer(),
    patient_nominal_p_lt_0p05_n = integer(),
    patient_supported_bh05_effect_n = integer(),
    patient_supported_bh10_effect_n = integer(),
    patient_supported_bh05_ids = character(),
    patient_supported_bh10_ids = character()
  )
}

write_tsv(
  problematic_null_strata,
  file.path(m04$tables, "empirical_null_problematic_strata.tsv")
)
write_tsv(
  null_risk_observed_gene_audit,
  file.path(m04$tables, "null_risk_observed_gene_audit.tsv.gz")
)
write_tsv(
  null_risk_overlap_summary,
  file.path(m04$tables, "null_risk_overlap_summary.tsv")
)
write_tsv(
  null_risk_overlap_interpretation,
  file.path(m04$qc, "null_risk_overlap_interpretation.tsv")
)
write_tsv(threshold_gate, file.path(m04$qc, "effect_threshold_gate.tsv"))

remodeling_analysis_all_cohorts_v2 <- remodeling_classes_interim %>%
  mutate(architecture_class = as.character(architecture_class)) %>%
  left_join(
    gene_model_family %>% select(
      cohort, gene_id, included_in_testing_family, testing_family_reason,
      scheduled_for_computation, testing_family_rule
    ),
    by = c("cohort", "gene_id")
  ) %>%
  left_join(gene_model_summary, by = c("cohort", "gene_id")) %>%
  mutate(
    included_in_testing_family = replace_na(included_in_testing_family, FALSE),
    scheduled_for_computation = replace_na(scheduled_for_computation, FALSE),
    gene_model_state = case_when(
      !cfg$gene_null$enabled ~ "not_requested",
      !included_in_testing_family ~ "outside_registered_testing_family",
      !scheduled_for_computation ~ "not_run_under_development_cap",
      startsWith(replace_na(gene_model_inference_reason, ""), "conservative_p1_") ~
        "conservative_P1_failure",
      gene_model_effect_supported %in% TRUE ~ "FDR_and_effect_supported",
      gene_model_supported %in% TRUE ~ "FDR_supported_below_effect_gate",
      TRUE ~ "not_FDR_supported"
    ),
    gene_model_changes_class = FALSE,
    refinement_changes_class = FALSE,
    gpd_tail_changes_class = FALSE,
    patient_test_changes_class = FALSE,
    inference_method_locked = TRUE,
    inference_selection_authorized = FALSE,
    effect_threshold_gate_pass = threshold_gate_pass,
    authoritative_status = if_else(
      threshold_gate_pass,
      "eligible_subject_to_remaining_completion_checks",
      "audit_only_effect_threshold_gate_failed"
    ),
    class_source = "Module04_gene_remodeling_classes_v2",
    class_framework = "A-F-U-NQ_v2",
    effect_gate_rule = paste0(
      "JSD_base2>=", cfg$thresholds$jsd,
      " OR max_abs_delta_IF>=", cfg$thresholds$delta_if
    ),
    primary_direction_support_rule = paste0(
      "driver direction supported in >=", cfg$thresholds$primary_patient_support,
      " paired patients"
    )
  )

# Public result contract: one row per cohort/gene; inferential evidence never changes class.
remodeling_analysis_all_cohorts_v2 <- remodeling_analysis_all_cohorts_v2 %>%
  mutate(
    downstream_primary_remodeling_set = is_high_confidence_remodeled,
    downstream_class_stratified_set = is_high_confidence_class,
    primary_cohort = cohort == primary_cohort_name
  )
gene_remodeling_classes_v2 <- remodeling_analysis_all_cohorts_v2 %>%
  filter(cohort == primary_cohort_name) %>%
  left_join(loo_gene_stability, by = "gene_id") %>%
  left_join(patient_gene_inference_summary, by = "gene_id")

class_codebook <- tribble(
  ~architecture_class, ~description,
  "A", "Dominant isoform changes identity between conditions; effect gate passes",
  "B", "Same dominant isoform; verified secondary-isoform remodeling; effect gate passes",
  "C", "Clear reference-condition dominance is lost; effect gate passes",
  "D", "Clear comparison-condition dominance is gained; effect gate passes",
  "E", "Remodeling without clear dominance in either condition; effect gate passes",
  "F", "Both upper bootstrap effect limits remain below the prespecified thresholds",
  "U", "Quantifiable gene with inconclusive remodeling or stability evidence",
  "NQ", "Insufficient retained isoforms or complete quantifiable pairs"
)
class_invariant_qc <- remodeling_analysis_all_cohorts_v2 %>%
  summarise(
    unique_gene_cohort_keys = !anyDuplicated(paste(cohort, gene_id)),
    allowed_classes = all(architecture_class %in% cfg$class_levels),
    nq_not_quantifiable = all(quantifiable[architecture_class == "NQ"] %in% FALSE),
    remodeling_has_effect = all(effect_pass[architecture_class %in% LETTERS[1:5]] %in% TRUE),
    f_supported_stable = all(supported_stability[architecture_class == "F"] %in% TRUE),
    primary_set_correct = all(downstream_primary_remodeling_set ==
      (is_remodeled_AE & effect_stable %in% TRUE & direction_support_primary)),
    strict_set_subset = all(!downstream_class_stratified_set | downstream_primary_remodeling_set)
  )

# Sensitivity comparisons are metadata-driven. Each cohort is a separate omission/refit.
base_lookup <- gene_remodeling_classes_v2 %>% select(
  gene_id, architecture_class, observed_architecture, jsd, max_abs_delta_if,
  driver_isoform, is_high_confidence_remodeled
)
sensitivity_comparison <- imap_dfr(analysis_cohorts[-1L], function(ids, cohort_name) {
  followup <- remodeling_analysis_all_cohorts_v2 %>% filter(cohort == cohort_name) %>%
    select(all_of(names(base_lookup)))
  left_join(base_lookup, followup, by = "gene_id", suffix = c("_primary", "_sensitivity")) %>%
    mutate(cohort = cohort_name, excluded_patients = paste(setdiff(patient_ids, ids), collapse = ","),
      exact_class_agreement = architecture_class_primary == architecture_class_sensitivity,
      effect_agreement = (observed_architecture_primary %in% LETTERS[1:5]) ==
        (observed_architecture_sensitivity %in% LETTERS[1:5]),
      driver_agreement = driver_isoform_primary == driver_isoform_sensitivity)
})
if (!length(analysis_cohorts[-1L])) sensitivity_comparison <- tibble(
  cohort = character(), gene_id = character(), excluded_patients = character(),
  exact_class_agreement = logical(), effect_agreement = logical(), driver_agreement = logical()
)
sensitivity_summary <- sensitivity_comparison %>% group_by(cohort, excluded_patients) %>%
  summarise(genes_compared = n(), exact_class_agreement = mean(exact_class_agreement, na.rm = TRUE),
    effect_agreement = mean(effect_agreement, na.rm = TRUE),
    driver_agreement = mean(driver_agreement, na.rm = TRUE), .groups = "drop")

# Configurable alpha lists: retain source compatibility columns and export full long thresholds.
patient_inference_fdr <- map_dfr(cfg$inference$fdr_thresholds, function(alpha) {
  patient_specific_inference %>% mutate(fdr_threshold = alpha,
    fdr_supported = patient_test_available & is.finite(patient_test_bh) & patient_test_bh <= alpha,
    fdr_and_effect_supported = fdr_supported & patient_effect_pass)
})
gene_inference_fdr <- if (nrow(gene_model_summary)) map_dfr(cfg$inference$fdr_thresholds, function(alpha) {
  gene_model_summary %>% mutate(fdr_threshold = alpha,
    reference_fdr_supported = is.finite(gene_model_q) & gene_model_q <= alpha,
    refinement_fdr_supported = refinement_available %in% TRUE & is.finite(refinement_bh) & refinement_bh <= alpha,
    gpd_fdr_supported = is.finite(gene_model_gpd_hybrid_q) & gene_model_gpd_hybrid_q <= alpha)
}) else tibble(cohort = character(), gene_id = character(), fdr_threshold = numeric())

write_tsv(remodeling_analysis_all_cohorts_v2, file.path(m04$tables, "gene_remodeling_classes_all_cohorts.tsv.gz"))
write_tsv(gene_remodeling_classes_v2, file.path(m04$tables, "gene_remodeling_classes_primary.tsv.gz"))
write_tsv(class_codebook, file.path(m04$tables, "remodeling_class_codebook.tsv"))
write_tsv(class_invariant_qc, file.path(m04$qc, "class_invariant_qc.tsv"))
write_tsv(sensitivity_comparison, file.path(m04$tables, "exclusion_sensitivity_comparison.tsv.gz"))
write_tsv(sensitivity_summary, file.path(m04$tables, "exclusion_sensitivity_summary.tsv"))
write_tsv(patient_inference_fdr, file.path(m04$tables, "patient_inference_fdr_thresholds.tsv.gz"))
write_tsv(gene_inference_fdr, file.path(m04$tables, "gene_inference_fdr_thresholds.tsv.gz"))

metric_qc_pass <- metric_identity_qc$n_nonfinite_metrics[[1]] == 0L &&
  metric_identity_qc$n_jsd_outside_unit_interval[[1]] == 0L &&
  metric_identity_qc$n_delta_if_outside_unit_interval[[1]] == 0L &&
  metric_identity_qc$n_tvd_outside_unit_interval[[1]] == 0L &&
  metric_identity_qc$n_max_delta_exceeds_tvd[[1]] == 0L &&
  metric_identity_qc$maximum_log_base_identity_error[[1]] < 1e-10 &&
  all(isoform_contribution_qc$jsd_contribution_error < 1e-9) &&
  all(isoform_contribution_qc$primary_composition_error < 1e-10) &&
  all(isoform_contribution_qc$recurrent_composition_error < 1e-10)
geometry_qc_pass <- geometry_qc$n_nonfinite[[1]] == 0L &&
  geometry_qc$n_jsd_outside_unit_interval[[1]] == 0L &&
  geometry_qc$n_max_delta_exceeds_tvd[[1]] == 0L &&
  geometry_qc$strata_completed[[1]] == geometry_qc$strata_expected[[1]] &&
  geometry_qc$minimum_replicates_per_stratum[[1]] >= cfg$geometry$replicates_per_stratum
bootstrap_qc_pass <- all(bootstrap_vectorization_qc$pass) &&
  all(bootstrap_qc$n_probability_sum_errors == 0L) && all(bootstrap_qc$n_invalid_intervals == 0L) &&
  all(bootstrap_qc$maximum_reported_probability_mcse <= cfg$bootstrap$mcse_target)
completion_gates <- bind_rows(
  tibble(check = c("metric identities", "geometry simulation integrity", "paired-bootstrap integrity and Monte Carlo precision",
                    "patient support accounting", "class invariants"),
    pass = c(metric_qc_pass, geometry_qc_pass, bootstrap_qc_pass, patient_support_qc_pass,
             all(unlist(class_invariant_qc))), blocking = TRUE, enabled = TRUE),
  threshold_gate %>% transmute(check = paste0("split-null: ", check), pass, blocking, enabled = TRUE),
  tibble(check = "complete frozen cohort-testing family scheduled", blocking = TRUE,
    enabled = cfg$gene_null$enabled && cfg$analysis_profile == "production",
    pass = all(gene_model_family$scheduled_for_computation[gene_model_family$included_in_testing_family] %in% TRUE)),
  patient_test_qc %>% transmute(check = paste0("patient count-composition test: ", patient), pass,
    blocking = TRUE, enabled = cfg$patient_test$enabled),
  gene_model_execution_gate %>% transmute(check = paste0("cohort-model: ", check), pass,
    blocking = TRUE, enabled = cfg$gene_null$enabled),
  refinement_qc %>% transmute(check = paste0("independent refinement: ", cohort), pass,
    blocking = TRUE, enabled = cfg$gene_null$enabled && cfg$refinement$enabled),
  gpd_qc %>% transmute(check = paste0("tail results finite in frozen family: ", cohort), pass,
    blocking = TRUE, enabled = cfg$gene_null$enabled && cfg$gpd$enabled),
  validation_qc %>% transmute(check = paste0("truth-known validation: ", check), pass, blocking,
    enabled = cfg$gene_null$enabled && cfg$validation$enabled)
) %>% mutate(pass = as.logical(pass), blocking = as.logical(blocking))
complete <- all(completion_gates$pass[completion_gates$enabled & completion_gates$blocking] %in% TRUE)
write_tsv(completion_gates, file.path(m04$qc, "completion_gates.tsv"))

# Publication figures are reader-facing. Every facet and composite is also exported separately.
if (!identical(public_cfg$reporting$make_plots, FALSE) && !identical(public_cfg$reporting$enabled, FALSE)) {
  class_colors <- unlist(public_cfg$reporting$palettes$remodeling_classes %||% c(A="#B2182B",B="#EF8A62",C="#7B3294",D="#2166AC",E="#1B7837",F="#636363",U="#BDBDBD",NQ="#F0F0F0"))
  add_plot <- function(p, stem, w=7, h=5) {
    figure_records[[stem]] <<- save_publication_plot(p + paper_theme(), stem, w, h)
  }
  add_plot(ggplot(gene_remodeling_classes_v2 %>% filter(quantifiable),
    aes(max_abs_delta_if, jsd, color = architecture_class)) + geom_point(alpha=.65, size=1.8) +
    geom_vline(xintercept=cfg$thresholds$delta_if, linetype=2) +
    geom_hline(yintercept=cfg$thresholds$jsd, linetype=2) +
    scale_color_manual(values=class_colors, drop=FALSE) +
    labs(title="Isoform remodeling", x="Maximum absolute change in isoform fraction", y="JSD (base 2)", color="Class"),
    "observed_remodeling_geometry")
  add_plot(ggplot(class_counts, aes(architecture_class,n_genes,fill=architecture_class)) +
    geom_col() + facet_wrap(~cohort, scales="free_y") + scale_fill_manual(values=class_colors,guide="none") +
    labs(title="Gene remodeling classes",x="Class",y="Genes"), "remodeling_class_counts",8,5)
  for (cohort_name in names(analysis_cohorts)) {
    cd <- class_counts %>% filter(cohort == cohort_name)
    add_plot(ggplot(cd,aes(architecture_class,n_genes,fill=architecture_class)) + geom_col() +
      scale_fill_manual(values=class_colors,guide="none") + labs(title="Gene remodeling classes",x="Class",y="Genes"),
      paste0("remodeling_class_counts_",cohort_name),6,4.5)
    bs <- remodeling_analysis_all_cohorts_v2 %>% filter(cohort == cohort_name, is_remodeled_AE)
    for (metric in c("effect_pass_bootstrap_probability","observed_class_bootstrap_probability","observed_driver_bootstrap_probability")) {
      if (nrow(bs)) add_plot(ggplot(bs,aes(x=.data[[metric]],fill=architecture_class)) +
        geom_histogram(bins=20) + geom_vline(xintercept=cfg$thresholds$stability_probability,linetype=2) +
        scale_fill_manual(values=class_colors) + labs(title="Paired-bootstrap support",x=gsub("_"," ",metric),y="Genes",fill="Class"),
        paste0("bootstrap_",metric,"_",cohort_name),7,5)
    }
    ps <- patient_support_sensitivity %>% filter(cohort==cohort_name)
    add_plot(ggplot(ps,aes(support_threshold,n_effect_stable_with_support,color=architecture_class)) +
      geom_line(linewidth=.8) + geom_point(size=2.2) +
      scale_color_manual(values=c(class_colors, `A-E combined`="#222222")) +
      labs(title="Cross-patient direction support",x="Minimum supporting patients",y="Genes",color="Class"),
      paste0("patient_support_sensitivity_",cohort_name),7,5)
  }
  for (metric in c("jsd","max_abs_delta_if","delta_entropy")) {
    dd <- gene_remodeling_classes_v2 %>% filter(quantifiable)
    add_plot(ggplot(dd,aes(x=.data[[metric]])) + geom_histogram(bins=40,fill="#2C7FB8") +
      labs(title="Isoform remodeling",x=gsub("_"," ",metric),y="Genes"),paste0("remodeling_",metric,"_distribution"),6,4.5)
  }
  if (nrow(loo_summary)) add_plot(ggplot(loo_summary,aes(omitted_patient,full_effect_retention_fraction)) +
    geom_col(fill="#2C7FB8") + scale_y_continuous(labels=scales::percent) +
    labs(title="Leave-one-patient-out stability",x="Omitted patient",y="Observed remodeling retained"),"leave_one_patient_out_stability",7,5)
  for (k in sort(unique(geometry_simulation$n_isoforms))) {
    add_plot(ggplot(geometry_simulation %>% filter(n_isoforms==k),aes(max_abs_delta_if,jsd,color=factor(depth))) +
      geom_point(alpha=.2,size=.7) + labs(title=paste("Metric geometry:",k,"isoforms"),
        x="Maximum absolute change in isoform fraction",y="JSD (base 2)",color="Reads"),paste0("feasible_geometry_K",k),6,5)
  }
  if (nrow(validation_aggregate)) for (truth_label in unique(validation_aggregate$truth)) {
    dd <- validation_aggregate %>% filter(truth == truth_label)
    add_plot(ggplot(dd,aes(method,call_rate_among_available)) + geom_col(fill="#2C7FB8") +
      coord_flip() + scale_y_continuous(labels=scales::percent) +
      labs(title=if(truth_label=="null") "Truth-known null calibration" else "Truth-known alternative detection",x=NULL,y="Called at P < 0.05"),
      paste0("truth_known_validation_",truth_label),8,4.8)
  }
  if (nrow(gene_model_summary)) for (cohort_name in names(analysis_cohorts)) {
    dd <- gene_model_summary %>% filter(cohort==cohort_name)
    add_plot(ggplot(dd,aes(gene_model_observed_score,-log10(pmax(gene_model_p,.Machine$double.xmin)))) +
      geom_point(alpha=.65,size=1.8) + labs(title="Cohort remodeling evidence",x="Combined effect score",y=expression(-log[10](P))),
      paste0("cohort_monte_carlo_evidence_",cohort_name),7,5)
  }
  # All patient-specific FDR figures use both configured thresholds, with neutral labels.
  for (alpha in cfg$inference$fdr_thresholds) {
    dd <- patient_inference_fdr %>% filter(fdr_threshold==alpha) %>%
      count(patient, fdr_and_effect_supported, name="genes") %>% filter(fdr_and_effect_supported)
    add_plot(ggplot(dd,aes(patient,genes)) + geom_col(fill="#2C7FB8") +
      labs(title=paste0("Patient-specific remodeling: FDR < ",format(alpha)),x="Patient",y="Genes"),
      paste0("patient_specific_calls_FDR",gsub("\\.","p",format(alpha))),6,4.5)
  }
}
figure_manifest <- bind_rows(figure_records)
write_tsv(figure_manifest, file.path(m04$tables,"figure_manifest.tsv"))
result <- list(
  gene_table = remodeling_analysis_all_cohorts_v2, primary_gene_table = gene_remodeling_classes_v2,
  isoform_table = observed_isoform, patient_table = observed_patient,
  patient_inference = patient_specific_inference, patient_inference_fdr = patient_inference_fdr,
  observed_gene = observed_gene, bootstrap = bootstrap_summary,
  bootstrap_class_probability = bootstrap_class_probability,
  cohort_model = gene_model_summary, model_family = gene_model_family,
  model_covariance = gene_model_covariance, resolution_refinement = refinement_results,
  inference_fdr = gene_inference_fdr, inference_methods = method_long,
  truth_known_validation = truth_known_validation, validation_summary = validation_aggregate,
  leave_one_patient_out = loo_comparison, sensitivity_comparison = sensitivity_comparison,
  split_null = list(overall = null_overall_rate, strata = null_stratum_rate,
                   threshold_gate = threshold_gate, risk_audit = null_risk_observed_gene_audit),
  geometry = geometry_simulation, config = cfg, fingerprint = input_fingerprint,
  paths = m04, completion_gates = completion_gates, complete = complete,
  execution_profile = cfg$analysis_profile, authoritative_classes = TRUE
)
atomic_save_rds(result, file.path(m04$objects, if(complete) "remodeling_results.rds" else "remodeling_audit_only.rds"))
if (cfg$strict_completion && !complete) stop(
  "Remodeling did not pass strict completion. Inspect qc/completion_gates.tsv. Audit outputs were retained. Failed gates: ",
  paste(completion_gates$check[completion_gates$enabled & completion_gates$blocking & !completion_gates$pass],collapse="; ")
)
result

}
