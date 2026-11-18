# Independent units, exact label permutations, fixed universes and RNG hygiene.
(function() {
  dir <- tempfile("unpaired_methods_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  cfg <- yaml::read_yaml(file.path(.repo_root, "config", "analysis.yaml"))
  cfg$comparison$paired <- FALSE
  cfg$unpaired <- list(min_samples_per_condition = 3L, minimum_evaluable_per_condition = 3L,
    n_bootstrap = 9L, n_permutations = 100L, minimum_prevalence_difference = .5)
  cfg$reporting$enabled <- FALSE; cfg$execution$resume <- TRUE
  cfg$execution$checkpoint_every_genes <- 1L
  cfg$gain_loss$n_downsample_replicates <- 3L
  cfg$gain_loss$keep_structural_categories <- "full-splice_match"
  cfg$gain_loss$gene_evaluable_min_fl <- 10L
  cfg$gain_loss$isoform_presence_min_fl <- 5L
  cfg$gain_loss$isoform_presence_min_fraction <- .1
  cfg$remodeling$definition$min_gene_sample_count <- 10L
  cfg$remodeling$definition$min_isoform_total_fl <- 1L
  cfg$remodeling$definition$min_retained_mass <- .8
  counts <- rbind(c(100, 100, 100, 0, 0, 0), c(0, 0, 0, 100, 100, 100),
                 rep(50, 6), rep(50, 6))
  rownames(counts) <- paste0("JUNC_", 1:4); colnames(counts) <- paste0("sample", 1:6)
  samples <- data.frame(sample = colnames(counts), patient = paste0("independent", 1:6),
    timepoint = rep(c("primary", "recurrent"), each = 3),
    condition = rep(c("control", "treated"), each = 3))
  annotation <- data.frame(isoform_harmonized_id = rownames(counts), associated_gene = rep(c("shift", "same"), each = 2),
    structural_category = "full-splice_match", gene_mapped = TRUE, cross_sample_comparable = TRUE)
  inputs <- list(counts_matrix = counts, sample_manifest = samples, isoform_annotation = annotation)
  cohorts <- list(all_independent = samples$sample)
  ctx <- rna_unpaired_context(inputs, cfg, cohorts)
  prevalence <- rna_unpaired_prevalence(counts, ctx, cfg, cohorts[[1]])
  stopifnot(identical(prevalence$prevalence_class,
    c("comparison_depleted", "comparison_enriched", "no_large_prevalence_difference", "no_large_prevalence_difference")),
    all(prevalence$reference_evaluable == 3), all(prevalence$comparison_evaluable == 3))
  state <- dtu_rng_snapshot()
  gain <- run_unpaired_gain_loss(inputs, cfg, cohorts, file.path(dir, "gain"))
  remodel <- run_unpaired_remodeling(inputs, cfg, cohorts, file.path(dir, "remodel"))
  qc <- run_unpaired_qc(inputs, cfg, cohorts, gain, file.path(dir, "qc"))
  stopifnot(identical(state, dtu_rng_snapshot()), gain$complete, remodel$complete, qc$complete,
    gain$rarefaction_target == 200, all(gain$tables$prevalence$depth_robust[1:2]))
  shift <- remodel$gene_table[remodel$gene_table$gene_id == "shift", ]
  same <- remodel$gene_table[remodel$gene_table$gene_id == "same", ]
  stopifnot(shift$jsd == 1, shift$total_variation == 1, shift$max_abs_delta_if == 1,
    shift$gene_model_p == .1, shift$permutation_draws == 20,
    shift$minimum_attainable_permutation_p == .1,
    shift$bootstrap_confidence == cfg$remodeling$bootstrap$bootstrap_confidence,
    same$jsd == 0, same$gene_model_p == 1,
    !"architecture_class" %in% names(remodel$gene_table), is.null(remodel$patient_inference),
    all(remodel$gene_table$design == "unpaired"), all(qc$tests$p_value == 1))
  # Reading matching cache must preserve numerical values and caller RNG state.
  again <- run_unpaired_remodeling(inputs, cfg, cohorts, file.path(dir, "remodel"))
  stopifnot(identical(remodel$gene_table, again$gene_table), identical(state, dtu_rng_snapshot()))
  # An interrupted cohort can resume previously committed gene blocks without
  # restarting all genome-wide bootstrap/permutation calculations.
  unlink(file.path(dir, "remodel", "cache", "all_independent_unpaired_remodeling.rds"))
  block_paths <- list.files(file.path(dir, "remodel", "cache", "unpaired_gene_blocks"), full.names = TRUE, recursive = TRUE)
  stopifnot(length(block_paths) == 2L)
  block_hashes <- tools::md5sum(block_paths)
  message_log <- capture.output(resumed <- run_unpaired_remodeling(inputs, cfg, cohorts, file.path(dir, "remodel")), type = "message")
  stopifnot(identical(resumed$gene_table, remodel$gene_table),
    identical(block_hashes, tools::md5sum(block_paths)),
    any(grepl("Reusing unpaired remodeling", message_log, fixed = TRUE)))
  narrower <- cfg; narrower$remodeling$bootstrap$bootstrap_confidence <- .8
  changed_ci <- run_unpaired_remodeling(inputs, narrower, cohorts, file.path(dir, "narrow_ci"))
  stopifnot(all(changed_ci$gene_table$bootstrap_confidence == .8))
  # A moderate redistribution can cross the fraction threshold without crossing
  # the JSD threshold. It still meets the shared prespecified effect rule.
  moderate <- inputs
  moderate$counts_matrix[3:4, 4:6] <- c(80, 20)
  moderate_result <- run_unpaired_remodeling(moderate, cfg, cohorts, file.path(dir, "moderate"))
  m <- moderate_result$gene_table[moderate_result$gene_table$gene_id == "same", ]
  stopifnot(m$jsd < cfg$remodeling$definition$jsd_threshold,
    m$max_abs_delta_if > cfg$remodeling$definition$max_delta_if_threshold,
    m$unpaired_effect_class == "usage_shift")
  rare <- dtu_with_rng(NULL, 127L, rna_unpaired_rarefy(counts, 43))
  stopifnot(all(colSums(rare) == 43), all(rare <= counts), all(rare == floor(rare)))
  low <- counts; low[, 1] <- 0
  lowctx <- ctx; lowctx$counts <- low
  z <- rna_unpaired_prevalence(low, lowctx, cfg, cohorts[[1]])
  stopifnot(all(z$reference_evaluable == 2), all(z$prevalence_class == "not_quantifiable"))
  repeated <- inputs; repeated$sample_manifest$patient[2] <- repeated$sample_manifest$patient[1]
  stopifnot(inherits(try(rna_unpaired_context(repeated, cfg, cohorts), silent = TRUE), "try-error"))
  cfg$unpaired$min_samples_per_condition <- 4L
  stopifnot(inherits(try(rna_unpaired_context(inputs, cfg, cohorts), silent = TRUE), "try-error"))
  cfg$unpaired$min_samples_per_condition <- 2L
  stopifnot(inherits(try(rna_unpaired_settings(cfg), silent = TRUE), "try-error"))
})()
message("Independent-sample prevalence, remodeling, exact permutations and QC checks passed.")

# Raw Pigeon import through the three independent-sample engines, including
# every leave-one-sample-out cohort. DTU fitting has its own dependency tests.
(function() {
  work <- tempfile("unpaired_pipeline_"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixture <- make_rna_fixture(file.path(work, "pigeon"), n_pairs = 4L)
  metadata <- utils::read.csv(fixture$metadata_file, stringsAsFactors = FALSE)
  metadata$subject_id <- metadata$sample_id
  utils::write.csv(metadata, fixture$metadata_file, row.names = FALSE)
  cfg <- yaml::read_yaml(file.path(.repo_root, "config", "analysis.yaml"))
  cfg$inputs$sample_metadata <- fixture$metadata_file
  cfg$inputs$root_directory <- dirname(fixture$metadata_file)
  cfg$output$directory <- file.path(work, "workflow")
  cfg$comparison$paired <- FALSE
  cfg$cohorts <- list(list(id = "all_samples", exclude_samples = character()))
  cfg$sensitivity$mode <- "all"
  cfg$sensitivity$samples <- NULL
  cfg$unpaired <- list(min_samples_per_condition = 3L, n_bootstrap = 9L,
    n_permutations = 100L, minimum_prevalence_difference = .5)
  cfg$gain_loss$n_downsample_replicates <- 2L
  cfg$gain_loss$n_depth_curve_replicates <- 2L
  cfg$execution$profile <- "smoke"
  cfg$execution$checkpoint_every_genes <- 2L
  cfg$execution$resume <- TRUE
  cfg$execution$force_refit <- FALSE
  cfg$reporting$enabled <- FALSE
  cfg$enrichment$enabled <- FALSE
  cfg$references <- list(annotation_gtf = NULL, gene_id_map = NULL,
    pathway_term2gene = NULL, pathway_term_names = NULL)
  cfg$enrichment$reference_cache_directory <- file.path(work, "references")
  cfg$stages <- list(import = TRUE, gain_loss = TRUE, qc = TRUE, remodeling = TRUE, dtu = FALSE)
  config_file <- file.path(work, "analysis.yaml")
  yaml::write_yaml(cfg, config_file)
  old_root <- getOption("IsoRemodel.code_root")
  options(IsoRemodel.code_root = .repo_root)
  on.exit(options(IsoRemodel.code_root = old_root), add = TRUE)
  resolved <- read_rna_config(config_file, validate_inputs = FALSE)
  defaults <- rna_unpaired_settings(list())
  stopifnot(resolved$unpaired$minimum_evaluable_per_condition == 3L,
    defaults$min_samples_per_condition == 3L, defaults$n_bootstrap == 1000L,
    defaults$n_permutations == 1000L)
  log <- capture.output(out <- run_rna_analysis(config_file), type = "message")
  stopifnot(identical(out$run_record$status, "diagnostic_complete"),
    identical(out$run_record$design, "unpaired"), length(out$run_record$cohorts) == 9L,
    all(lengths(out$run_record$cohorts[-1L]) == 7L),
    all(vapply(out$run_record$stage_records, function(x) isTRUE(x$complete), logical(1))),
    !isTRUE(out$run_record$completed), !file.exists(file.path(cfg$output$directory, "COMPLETE")),
    length(out$results$gain_loss$cohort_results) == 9L,
    length(out$results$remodeling$cohort_results) == 9L,
    length(unique(out$results$qc$tests$cohort)) == 9L,
    all(out$results$gain_loss$tables$prevalence$reference_samples %in% 3:4),
    all(out$results$gain_loss$tables$prevalence$comparison_samples %in% 3:4),
    out$results$gain_loss$rarefaction_target == 381L,
    !"architecture_class" %in% names(out$results$remodeling$gene_table))
})()
message("Raw-input unpaired smoke passed for the full cohort and all eight sample exclusions.")
