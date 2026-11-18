# Translate the public YAML to the original, explicit statistical parameter names.
rna_remodeling_parameters <- function(public_cfg) {
  defaults <- list(
  project_root = "F:/PacBio revision analysis",
  analysis_profile = "production",
  jsd_threshold = 0.2,
  max_delta_if_threshold = 0.25,
  entropy_direction_threshold = 0.2,
  direction_delta_if = 0.1,
  min_gene_sample_count = 5,
  min_isoform_total_fl = 5,
  min_retained_mass = 0.8,
  min_evaluable_pairs = 4,
  dominance_mean_if = 0.35,
  dominance_margin = 0.1,
  dominance_support_n = 4,
  dominance_sample_min_if = 0.25,
  secondary_min_mass = 0.15,
  stability_probability = 0.8,
  bootstrap_confidence = 0.95,
  bootstrap_mcse_target = 0.01,
  n_bootstrap_replicates = 5000,
  n_split_null_replicates = 1000,
  checkpoint_every_null_replicates = 25,
  n_geometry_per_stratum = 2500,
  geometry_isoform_numbers_csv = "2,3,4,6,10,20",
  geometry_depths_csv = "20,50,100,500,2000",
  patient_support_thresholds_csv = "1,2,3,4,5,6",
  primary_patient_support = 4,
  empirical_null_overall_target = 0.05,
  empirical_null_stratum_target = 0.1,
  empirical_null_min_stratum_n = 200,
  empirical_null_convergence_tolerance = 0.005,
  stop_on_threshold_failure = FALSE,
  run_gene_null = TRUE,
  gene_null_replicates = 20000,
  gene_null_max_isoforms = 20,
  gene_null_max_genes = 0,
  finite_count_pseudocount = 0.5,
  covariance_initial_replicates = 4096,
  covariance_maximum_replicates = 16384,
  covariance_trace_tolerance = 0.1,
  covariance_frobenius_tolerance = 0.25,
  profile_grid_points = 257,
  profile_tail_log_drop = 30,
  gene_null_fdr = 0.05,
  inferential_fdr_thresholds_csv = "0.05,0.10",
  run_resolution_refinement = TRUE,
  refinement_trigger_p = 0.05,
  refinement_replicates = 120000,
  refinement_max_genes = 0,
  run_gpd_tail = TRUE,
  gpd_empirical_exceedance_trigger = 10,
  gpd_tail_fractions_csv = "0.025,0.05,0.10",
  gpd_min_tail_size = 80,
  gpd_log10_stability_tolerance = 0.75,
  gpd_ks_multiplier = 1.63,
  run_patient_level_test = TRUE,
  patient_test_min_condition_fl = 30,
  patient_test_min_pooled_component_fl = 5,
  patient_test_max_isoforms = 20,
  patient_test_min_expected_count = 1,
  patient_test_max_fraction_expected_below_5 = 0.2,
  patient_test_primary_fdr = 0.05,
  patient_test_sensitivity_fdr = 0.1,
  gene_model_minimum_success_fraction = 0.95,
  gene_model_minimum_stratum_success_fraction = 0.9,
  gene_model_minimum_stratum_n = 20,
  run_truth_known_validation = TRUE,
  validation_replicates_per_scenario = 12,
  validation_null_replicates = 2000,
  validation_minimum_production_availability = 0.95,
  validation_minimum_stratum_availability = 0.75,
  validation_maximum_null_call_rate = 0.1,
  checkpoint_every_genes = 100,
  resume_cached_results = TRUE,
  seed = 20260903,
  strict_completion = FALSE
)
  flatten <- function(x) {
    out <- list()
    for (key in names(x)) {
      if (is.list(x[[key]]) && !is.data.frame(x[[key]]) && !is.null(names(x[[key]]))) {
        out <- c(out, flatten(x[[key]]))
      } else out[[key]] <- if (is.list(x[[key]])) unlist(x[[key]], use.names=FALSE) else x[[key]]
    }
    out
  }
  supplied <- flatten(public_cfg$remodeling)
  for (key in names(defaults)) {
    if (!is.null(supplied[[key]])) defaults[[key]] <- supplied[[key]]
    else if (grepl("_csv$", key)) {
      values <- supplied[[sub("_csv$", "", key)]]
      if (!is.null(values)) defaults[[key]] <- paste(unlist(values), collapse = ",")
    }
  }
  defaults$project_root <- "."
  defaults$analysis_profile <- public_cfg$execution$profile %||% "production"
  defaults$seed <- public_cfg$execution$seeds$remodeling %||% defaults$seed
  defaults$resume_cached_results <- isTRUE(public_cfg$execution$resume) &&
    !isTRUE(public_cfg$execution$force_refit)
  defaults$strict_completion <- isTRUE(public_cfg$execution$strict_completion)
  defaults$checkpoint_every_genes <- public_cfg$execution$checkpoint_every_genes %||% 100L
  alpha <- sort(unique(as.numeric(unlist(public_cfg$multiplicity$fdr_thresholds %||% c(.05,.10)))))
  defaults$gene_null_fdr <- min(alpha)
  defaults$inferential_fdr_thresholds_csv <- paste(alpha, collapse = ",")
  defaults$patient_test_primary_fdr <- min(alpha)
  defaults$patient_test_sensitivity_fdr <- max(alpha)
  defaults
}

rna_remodeling_engine_config <- function(public_cfg) {
  params <- rna_remodeling_parameters(public_cfg)
  project_root <- "."
cfg <- list(
  workflow_version = "0.6.0",
  project_root = project_root,
  analysis_profile = tolower(as.character(params$analysis_profile)),
  thresholds = list(
    jsd = as.numeric(params$jsd_threshold),
    delta_if = as.numeric(params$max_delta_if_threshold),
    entropy_delta = as.numeric(params$entropy_direction_threshold),
    direction_delta_if = as.numeric(params$direction_delta_if),
    min_gene_sample_count = as.integer(params$min_gene_sample_count),
    min_isoform_total_fl = as.integer(params$min_isoform_total_fl),
    min_retained_mass = as.numeric(params$min_retained_mass),
    min_evaluable_pairs = as.integer(params$min_evaluable_pairs),
    dominance_mean_if = as.numeric(params$dominance_mean_if),
    dominance_margin = as.numeric(params$dominance_margin),
    dominance_support_n = as.integer(params$dominance_support_n),
    dominance_sample_min_if = as.numeric(params$dominance_sample_min_if),
    secondary_min_mass = as.numeric(params$secondary_min_mass),
    stability_probability = as.numeric(params$stability_probability),
    patient_support = parse_numeric_csv(params$patient_support_thresholds_csv, TRUE),
    primary_patient_support = as.integer(params$primary_patient_support)
  ),
  bootstrap = list(
    replicates = as.integer(params$n_bootstrap_replicates),
    confidence = as.numeric(params$bootstrap_confidence),
    mcse_target = as.numeric(params$bootstrap_mcse_target)
  ),
  empirical_null = list(
    replicates = as.integer(params$n_split_null_replicates),
    checkpoint_every = as.integer(params$checkpoint_every_null_replicates),
    overall_target = as.numeric(params$empirical_null_overall_target),
    stratum_target = as.numeric(params$empirical_null_stratum_target),
    minimum_stratum_n = as.integer(params$empirical_null_min_stratum_n),
    convergence_tolerance = as.numeric(params$empirical_null_convergence_tolerance),
    stop_on_failure = isTRUE(params$stop_on_threshold_failure)
  ),
  geometry = list(
    replicates_per_stratum = as.integer(params$n_geometry_per_stratum),
    isoform_numbers = parse_numeric_csv(params$geometry_isoform_numbers_csv, TRUE),
    depths = parse_numeric_csv(params$geometry_depths_csv, TRUE)
  ),
  gene_null = list(
    enabled = isTRUE(params$run_gene_null),
    replicates = as.integer(params$gene_null_replicates),
    maximum_isoforms = as.integer(params$gene_null_max_isoforms),
    maximum_genes = as.integer(params$gene_null_max_genes),
    prior_strength = as.numeric(params$finite_count_pseudocount),
    covariance_initial_replicates = as.integer(params$covariance_initial_replicates),
    covariance_maximum_replicates = as.integer(params$covariance_maximum_replicates),
    covariance_trace_tolerance = as.numeric(params$covariance_trace_tolerance),
    covariance_frobenius_tolerance = as.numeric(params$covariance_frobenius_tolerance),
    profile_grid_points = as.integer(params$profile_grid_points),
    profile_tail_log_drop = as.numeric(params$profile_tail_log_drop),
    fdr = as.numeric(params$gene_null_fdr),
    minimum_success_fraction = as.numeric(params$gene_model_minimum_success_fraction),
    minimum_stratum_success_fraction =
      as.numeric(params$gene_model_minimum_stratum_success_fraction),
    minimum_stratum_n = as.integer(params$gene_model_minimum_stratum_n)
  ),
  inference = list(
    fdr_thresholds = parse_numeric_csv(
      params$inferential_fdr_thresholds_csv, FALSE
    ),
    primary_fdr = as.numeric(params$gene_null_fdr)
  ),
  refinement = list(
    enabled = isTRUE(params$run_resolution_refinement),
    trigger_p = as.numeric(params$refinement_trigger_p),
    requested_replicates = as.integer(params$refinement_replicates),
    replicates = as.integer(params$refinement_replicates),
    maximum_genes = as.integer(params$refinement_max_genes),
    design = paste(
      "independent_fixed_second_stage_selected_by_reference_P",
      "reusing_frozen_reference_covariance_REML_and_variance_profile",
      sep = ";"
    )
  ),
  gpd = list(
    enabled = isTRUE(params$run_gpd_tail),
    empirical_exceedance_trigger =
      as.integer(params$gpd_empirical_exceedance_trigger),
    tail_fractions = parse_numeric_csv(params$gpd_tail_fractions_csv, FALSE),
    minimum_tail_size = as.integer(params$gpd_min_tail_size),
    log10_stability_tolerance =
      as.numeric(params$gpd_log10_stability_tolerance),
    ks_multiplier = as.numeric(params$gpd_ks_multiplier),
    interpretation = "tail_estimator_for_same_REML_null"
  ),
  patient_test = list(
    enabled = isTRUE(params$run_patient_level_test),
    method = "two_condition_multinomial_likelihood_ratio",
    minimum_condition_FL = as.integer(params$patient_test_min_condition_fl),
    minimum_pooled_component_FL =
      as.integer(params$patient_test_min_pooled_component_fl),
    maximum_isoforms = as.integer(params$patient_test_max_isoforms),
    minimum_expected_count = as.numeric(params$patient_test_min_expected_count),
    maximum_fraction_expected_below_5 =
      as.numeric(params$patient_test_max_fraction_expected_below_5),
    primary_fdr = as.numeric(params$patient_test_primary_fdr),
    sensitivity_fdr = as.numeric(params$patient_test_sensitivity_fdr),
    multiplicity = "BH_separately_within_each_patient_frozen_family",
    interpretation = paste(
      "technical_count_composition_evidence; patient_BH_and_effect_gate_required;",
      "not_biological_replicate_inference"
    )
  ),
  validation = list(
    enabled = isTRUE(params$run_truth_known_validation),
    replicates_per_scenario = as.integer(params$validation_replicates_per_scenario),
    null_replicates = as.integer(params$validation_null_replicates),
    minimum_production_availability =
      as.numeric(params$validation_minimum_production_availability),
    minimum_stratum_availability =
      as.numeric(params$validation_minimum_stratum_availability),
    maximum_null_call_rate = as.numeric(params$validation_maximum_null_call_rate)
  ),
  checkpoint_every_genes = as.integer(params$checkpoint_every_genes),
  resume_cached_results = isTRUE(params$resume_cached_results),
  seed = as.integer(params$seed),
  strict_completion = isTRUE(params$strict_completion),
  class_levels = c("A", "B", "C", "D", "E", "F", "U", "NQ")
)


  cfg$workflow_version <- "0.6.0-public-0.1.0"
  cfg$bootstrap$quantile_type <- as.integer(public_cfg$remodeling$bootstrap$quantile_type %||% 8L)
  cfg$bootstrap$minimum_valid_fraction <- as.numeric(public_cfg$remodeling$bootstrap$minimum_valid_fraction %||% .80)
  cfg$split_null_public <- public_cfg$remodeling$split_null %||% list()
  cfg$validation_public <- public_cfg$remodeling$validation %||% list()
  cfg$inference$fdr_thresholds <- sort(unique(cfg$inference$fdr_thresholds))
  cfg
}
