# Synthetic, bounded execution check; never a production result.
rna_smoke_test <- function(repo_root = getOption("IsoRemodel.code_root", "."),
                           output_root = tempfile("RNA_smoke_"),
                           include_dtu = FALSE, make_plots = FALSE,
                           paired = TRUE, leave_one_out = FALSE) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo_root, "R", "load.R"), env)
  env$load_rna_code(repo_root, envir = env)
  sys.source(file.path(repo_root, "tests", "make_fixture.R"), env)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  fixture <- env$make_rna_fixture(file.path(output_root, "synthetic_inputs"))
  cfg <- env$read_rna_config(file.path(repo_root, "config", "analysis.yaml"), validate_inputs = FALSE)
  cfg <- cfg[!startsWith(names(cfg), ".")]
  cfg$comparison$paired <- isTRUE(paired)
  if (!isTRUE(paired)) {
    samples <- utils::read.csv(fixture$metadata_file, stringsAsFactors = FALSE)
    samples$subject_id <- samples$sample_id
    utils::write.csv(samples, fixture$metadata_file, row.names = FALSE, na = "")
    cfg$cohorts <- list(list(id = "all_samples", exclude_samples = character()))
    cfg$sensitivity$mode <- if (isTRUE(leave_one_out)) "all" else "none"
    cfg$sensitivity$patients <- NULL
    cfg$sensitivity$samples <- character()
    cfg$sensitivity$observed_leave_one_patient_out$enabled <- FALSE
    cfg$unpaired$n_bootstrap <- 50L
    cfg$unpaired$n_permutations <- 99L
  } else if (isTRUE(leave_one_out)) cfg$sensitivity$mode <- "all"
  cfg$inputs$root_directory <- dirname(fixture$metadata_file)
  cfg$inputs$sample_metadata <- fixture$metadata_file
  cfg$references$annotation_gtf <- fixture$gtf_file
  cfg$references$gene_id_map <- cfg$references$pathway_term2gene <- cfg$references$pathway_term_names <- NULL
  cfg$output$directory <- file.path(output_root, "results")
  cfg$enrichment$enabled <- FALSE
  cfg$enrichment$reference_cache_directory <- file.path(output_root, "references")
  cfg$execution$profile <- "smoke"
  cfg$execution$strict_completion <- FALSE
  cfg$execution$checkpoint_every_genes <- 2L
  cfg$reporting$enabled <- isTRUE(make_plots)
  cfg$gain_loss$n_downsample_replicates <- 4L
  cfg$gain_loss$n_depth_curve_replicates <- 2L
  cfg$gain_loss$checkpoint_every <- 2L
  cfg$remodeling$bootstrap$n_bootstrap_replicates <- 50L
  cfg$remodeling$bootstrap$bootstrap_mcse_target <- .2
  cfg$remodeling$split_null$n_split_null_replicates <- 4L
  cfg$remodeling$split_null$checkpoint_every_null_replicates <- 2L
  cfg$remodeling$geometry$n_geometry_per_stratum <- 10L
  cfg$remodeling$cohort_test$gene_null_replicates <- 40L
  cfg$remodeling$cohort_test$covariance_initial_replicates <- 32L
  cfg$remodeling$cohort_test$covariance_maximum_replicates <- 64L
  cfg$remodeling$cohort_test$profile_grid_points <- 33L
  cfg$remodeling$refinement$run_resolution_refinement <- FALSE
  cfg$remodeling$tail_diagnostic$run_gpd_tail <- FALSE
  cfg$remodeling$validation$run_truth_known_validation <- FALSE
  cfg$stages$dtu <- isTRUE(include_dtu)
  path <- file.path(output_root, "smoke_analysis.yaml")
  yaml::write_yaml(cfg, path)
  result <- env$run_rna_analysis(path)
  message("Synthetic smoke output: ", output_root,
          "\nThis run uses reduced simulations and does not certify study-result reproduction.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  root <- if (length(arg)) dirname(dirname(normalizePath(sub("^--file=", "", arg[[1L]]), winslash = "/"))) else "."
  rna_smoke_test(root, include_dtu = "--include-dtu" %in% commandArgs(TRUE),
                 make_plots = "--plots" %in% commandArgs(TRUE),
                 paired = !"--unpaired" %in% commandArgs(TRUE),
                 leave_one_out = "--leave-one-out" %in% commandArgs(TRUE))
}
