# Full plotting-only integration fixture: saved parents and corrected DTU,
# strict RNG provenance gates, source immutability, review PDF/ZIP, and no fits.
(function() {
  root <- .repo_root
  temp <- tempfile("isoremodel_release_replot_")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  previous_root <- getOption("IsoRemodel.code_root")
  options(IsoRemodel.code_root = root)
  on.exit(options(IsoRemodel.code_root = previous_root), add = TRUE)
  source <- file.path(temp, "previous")
  dtu <- file.path(temp, "confirmed", "M06_dtu")
  dir.create(file.path(source, "M03_qc", "figures"), recursive = TRUE)
  dir.create(file.path(source, "M06_dtu", "figures"), recursive = TRUE)
  dir.create(file.path(dtu, "objects"), recursive = TRUE)
  p <- ggplot2::ggplot(data.frame(sample = c("ptn1_primary", "ptn1_recurrent"),
      reads = c(10, 12)), ggplot2::aes(sample, reads)) + ggplot2::geom_col()
  saveRDS(list(plot = p, width = 7, height = 5),
    file.path(source, "M03_qc", "figures", "full_length_reads_by_sample.plot.rds"))
  # An invalid old DTU plot must never be read, much less exported.
  saveRDS(list(invalid_old_DTU = TRUE),
    file.path(source, "M06_dtu", "figures", "stale_DTU.plot.rds"))
  cohorts <- c("all_6_pairs", "excluding_ptn2_5_pairs")
  results <- stats::setNames(lapply(cohorts, function(cohort) {
    patients <- if (cohort == "all_6_pairs") paste0("ptn", 1:6) else paste0("ptn", c(1, 3:6))
    genes <- data.frame(gene_id = c("GENEA", "GENEB", "GENEC"),
      dtu_fdr = c(.001, .02, .08), raw_dtu_p_value = c(.0001, .003, .02),
      fit_status = "finite_test", architecture_class = c("A", "B", "E"),
      DTU_driver_abs_mean_delta_if = c(.2, .3, .1), DTU_driver_directional_support = c(4, 4, 3),
      jsd = c(.2, .3, .15), max_abs_delta_if = c(.2, .3, .2),
      dominant_isoform_switched = c(TRUE, FALSE, TRUE),
      primary_dominant_structure = c("full-splice_match", "novel_in_catalog", "incomplete-splice_match"),
      recurrent_dominant_structure = c("incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"))
    samples <- data.frame(patient = rep(patients, each = 2),
      timepoint = rep(c("primary", "recurrent"), length(patients)))
    feature <- expand.grid(gene_id = genes$gene_id, isoform_n = 1:2, stringsAsFactors = FALSE)
    feature$feature_id <- paste(feature$gene_id, feature$isoform_n, sep = "|")
    feature$structural_category <- ifelse(feature$isoform_n == 1, "full-splice_match", "incomplete-splice_match")
    feature$mean_delta_if <- ifelse(feature$isoform_n == 1, -.2, .2)
    feature$stageR_confirmed_alpha05 <- feature$gene_id != "GENEC"
    feature$stageR_confirmed_alpha10 <- TRUE
    usage <- merge(feature, samples)
    usage$isoform_fraction <- ifelse(usage$isoform_n == 1,
      ifelse(usage$timepoint == "primary", .7, .5), ifelse(usage$timepoint == "primary", .3, .5))
    patient <- expand.grid(gene_id = genes$gene_id, patient = patients, stringsAsFactors = FALSE)
    patient$cohort_driver_patient_delta_if <- rep(c(.2, -.1, .3), length.out = nrow(patient))
    patient$patient_direction_supported <- patient$cohort_driver_patient_delta_if >= .05
    list(gene = genes, patient = patient, effects = list(feature = feature, usage = usage),
      broad_structure = data.frame(), prepared = list(samples = samples,
        flow = data.frame(stage = c("Comparable universe", ">=2 observed isoforms", "Proportion filter", "Count filter"),
          genes = c(10, 8, 6, 3), isoforms = c(40, 36, 20, 6))))
  }), cohorts)
  contract <- dtu_rng_contract(NULL, 260605L, cohorts)
  full <- list(complete = TRUE, cohort_results = results, rng_contract = contract,
    gates = data.frame(gate = "fixture", pass = TRUE))
  complete_path <- file.path(dtu, "objects", "DTU_complete.rds")
  saveRDS(full, complete_path)
  status_path <- file.path(dirname(dtu), "refit_status.rds")
  saveRDS(list(completed = TRUE, original_results_reproduced = TRUE,
    status = "original_results_reproduced"), status_path)
  mr <- do.call(rbind, lapply(cohorts, function(cohort) data.frame(cohort = cohort,
    gene_id = c("GENEA", "GENEB", "GENEC", "OUTSIDE_DTU"),
    downstream_primary_remodeling_set = c(FALSE, TRUE, FALSE, TRUE))))
  mp <- file.path(source, "M04_remodeling", "cache")
  dir.create(mp, recursive = TRUE)
  saveRDS(list(complete = TRUE, result = list(gene_table = mr)), file.path(mp, "stage_result.rds"))
  ora <- data.frame(pathway = c("KEGG_LONG_REFERENCE_PATHWAY_LABEL_FOR_VARIABLE_WRAPPING_ONE",
      "KEGG_SECOND_REFERENCE_PATHWAY_WITH_A_LONG_LABEL"), collection = "KEGG",
    selected_fraction = c(.1, .15), observed_genes = c(4, 5),
    matched_p_value = c(.001, .02), matched_bh = c(.2, .3))
  for (cohort in cohorts) for (alpha in c(.05, .10)) {
    dd <- file.path(dtu, cohort, "enrichment", "tables")
    dir.create(dd, recursive = TRUE, showWarnings = FALSE)
    rna_write_tsv(ora, file.path(dd, paste0("matched_ORA_", dtu_fdr_tag(alpha), "_KEGG.tsv.gz")))
  }
  cfg <- yaml::read_yaml(file.path(root, "config", "analysis.yaml"), eval.expr = FALSE)
  cfg$reporting$separate_panels_for_composites <- TRUE
  cfg$reporting$formats <- "png"; cfg$reporting$dpi <- 72
  cfg$dtu$run_upset <- TRUE
  cfg$dtu$maximum_gene_plots <- 2
  cfg$enrichment$enabled <- TRUE
  file <- file.path(temp, "analysis.yaml")
  yaml::write_yaml(cfg, file)
  original <- rna_replot_validate_dtu(dtu)
  stopifnot(identical(original$contract, contract))
  bad <- full; bad$rng_contract$kind <- "Mersenne-Twister"
  saveRDS(bad, complete_path)
  stopifnot(inherits(try(rna_replot_validate_dtu(dtu), silent = TRUE), "try-error"))
  saveRDS(full, complete_path)
  file.rename(status_path, paste0(status_path, ".hold"))
  stopifnot(inherits(try(rna_replot_validate_dtu(dtu), silent = TRUE), "try-error"))
  file.rename(paste0(status_path, ".hold"), status_path)
  inputs <- list.files(temp, recursive = TRUE, full.names = TRUE)
  before <- tools::md5sum(inputs)
  state <- dtu_rng_snapshot()
  opt_names <- c("IsoRemodel.log_file", "IsoRemodel.progress_file",
    "IsoRemodel.progress_latest_file", "IsoRemodel.run_id", "IsoRemodel.run_start")
  options_before <- stats::setNames(lapply(opt_names, getOption), opt_names)
  res <- replot_isoremodel(file, source, dtu, file.path(temp, "replot"), formats = "png", dpi = 72)
  stopifnot(identical(before, tools::md5sum(inputs)),
    identical(state, dtu_rng_snapshot()), nrow(res$figure_index) >= 20,
    identical(options_before, stats::setNames(lapply(opt_names, getOption), opt_names)),
    all(res$figure_index$status == "complete"), file.exists(res$review_zip),
    length(res$review_pdfs) >= 1, all(file.exists(res$review_pdfs)),
    !any(grepl("stale_DTU", res$figure_index$relative_stem)),
    any(grepl("matched_ORA_FDR05_KEGG", res$figure_index$relative_stem)),
    any(grepl("DTU_effect_evidence", res$figure_index$relative_stem)),
    any(grepl("DTU_remodeling_intersections", res$figure_index$relative_stem)),
    !length(list.files(res$output_directory, pattern = "[.]rds$", recursive = TRUE)))
  zip <- zip::zip_list(res$review_zip)
  stopifnot("figure_index.tsv" %in% zip$filename,
    !any(grepl("^/|^[A-Za-z]:", zip$filename)),
    any(grepl("^review_volumes/", zip$filename)),
    isTRUE(res$summary$source_files_unchanged), !isTRUE(res$summary$analytical_modules_run))

  # Native independent-group results carry their own design and completion
  # manifest. They contain no paired architecture or patient-support columns.
  native <- file.path(temp, "independent")
  native_dtu <- file.path(native, "M06_dtu")
  dir.create(file.path(native_dtu, "objects"), recursive = TRUE)
  independent_samples <- data.frame(sample_id = paste0("sample", 1:8),
    timepoint = rep(c("primary", "recurrent"), each = 4), stringsAsFactors = FALSE)
  independent_cohorts <- list(all_samples = independent_samples$sample_id,
    excluding_sample8_7_samples = independent_samples$sample_id[-8])
  independent_results <- lapply(independent_cohorts, function(ids) {
    g <- data.frame(gene_id = c("GENEA", "GENEB", "GENEC"),
      dtu_fdr = c(.001, .02, .08), raw_dtu_p_value = c(.0001, .003, .02),
      fit_status = "finite_test", DTU_driver_abs_mean_delta_if = c(.2, .3, .1))
    features <- expand.grid(gene_id = g$gene_id, isoform_n = 1:2,
      stringsAsFactors = FALSE)
    features$feature_id <- paste(features$gene_id, features$isoform_n, sep = "|")
    features$structural_category <- ifelse(features$isoform_n == 1,
      "full-splice_match", "incomplete-splice_match")
    features$mean_delta_if <- ifelse(features$isoform_n == 1, -.2, .2)
    samples <- independent_samples[independent_samples$sample_id %in% ids, , drop = FALSE]
    usage <- merge(features, samples)
    usage$isoform_fraction <- ifelse(usage$isoform_n == 1,
      ifelse(usage$timepoint == "primary", .7, .5),
      ifelse(usage$timepoint == "primary", .3, .5))
    list(gene = g, effects = list(feature = features, usage = usage),
      prepared = list(samples = samples, flow = data.frame(stage = c("Input", "Tested"),
        genes = c(8, 3), isoforms = c(16, 6))))
  })
  independent <- list(complete = TRUE, analysis_design = "unpaired",
    cohort_results = independent_results,
    rng_contract = dtu_rng_contract(NULL, 260605L, names(independent_results)),
    gates = data.frame(gate = "fixture", pass = TRUE))
  independent_path <- file.path(native_dtu, "objects", "DTU_complete.rds")
  saveRDS(independent, independent_path)
  native_manifest <- list(completed = TRUE, status = "complete", package_version = "0.2.0",
    design = "unpaired", stage_records = list(dtu = list(complete = TRUE)))
  native_manifest_path <- file.path(native, "latest_run_manifest.rds")
  saveRDS(native_manifest, native_manifest_path)
  stopifnot(identical(rna_replot_validate_dtu(native_dtu)$analysis_design, "unpaired"),
    identical(original$analysis_design, "paired"))
  wrong <- native_manifest; wrong$design <- "paired"
  saveRDS(wrong, native_manifest_path)
  stopifnot(inherits(try(rna_replot_validate_dtu(native_dtu), silent = TRUE), "try-error"))
  saveRDS(native_manifest, native_manifest_path)
  native_status <- file.path(native, "refit_status.rds")
  file.copy(status_path, native_status)
  stopifnot(inherits(try(rna_replot_validate_dtu(native_dtu), silent = TRUE), "try-error"))
  unlink(native_status)
  stopifnot(inherits(try(replot_isoremodel(file, native, native_dtu,
    file.path(temp, "wrong_design")), silent = TRUE), "try-error"),
    !dir.exists(file.path(temp, "wrong_design")))
  native_mr <- file.path(native, "M04_remodeling", "cache")
  dir.create(native_mr, recursive = TRUE)
  saveRDS(list(complete = TRUE, result = list(gene_table = data.frame(
    cohort = "all_samples", gene_id = c("GENEA", "GENEB", "GENEC"),
    downstream_primary_remodeling_set = c(TRUE, FALSE, FALSE)))),
    file.path(native_mr, "stage_result.rds"))
  independent_cfg <- yaml::read_yaml(file.path(root, "config", "analysis_unpaired.yaml"),
    eval.expr = FALSE)
  independent_cfg$reporting$formats <- "png"; independent_cfg$reporting$dpi <- 72
  independent_cfg$dtu$maximum_gene_plots <- 2
  independent_cfg_path <- file.path(temp, "independent.yaml")
  yaml::write_yaml(independent_cfg, independent_cfg_path)
  native_files <- list.files(native, recursive = TRUE, full.names = TRUE)
  native_before <- tools::md5sum(native_files)
  native_result <- replot_isoremodel(independent_cfg_path, native, native_dtu,
    file.path(temp, "independent_replot"), formats = "png", dpi = 72)
  stopifnot(identical(native_before, tools::md5sum(native_files)),
    identical(native_result$summary$analysis_design, "unpaired"),
    isTRUE(native_result$summary$source_files_unchanged),
    !isTRUE(native_result$summary$analytical_modules_run),
    file.exists(native_result$review_zip),
    nrow(native_result$figure_index) >= 8L,
    all(native_result$figure_index$status == "complete"),
    any(grepl("DTU_independent_effect_evidence", native_result$figure_index$relative_stem)),
    any(grepl("DTU_independent_isoform_usage", native_result$figure_index$relative_stem)),
    any(grepl("DTU_sensitivity", native_result$figure_index$relative_stem)),
    !any(grepl("DTU_effect_evidence|DTU_remodeling_intersections|patient_support",
      native_result$figure_index$relative_stem)))
})()
message("Paired/native unpaired replot fixtures, design/RNG gates and unchanged sources passed.")
