# Executed by tests/run_all.R. Pure identity/reference tests always run;
# package-dependent numerical and stageR tests run when their dependencies exist.
stopifnot(identical(dtu_alpha_tag(.05), "alpha05"), identical(dtu_alpha_tag(.10), "alpha10"),
          identical(dtu_fdr_tag(.025), "FDR025"))

.dtu_test_env <- new.env(parent = globalenv())
sys.source(file.path(.repo_root, "inst", "engines", "dtu_core.R"), envir = .dtu_test_env)
.dtu_samples <- data.frame(sample_id = c("A_p", "A_r", "B_p", "B_r", "C_p", "C_r"),
    patient = factor(rep(c("A", "B", "C"), each = 2)), timepoint = factor(rep(c("primary", "recurrent"), 3), levels = c("primary", "recurrent")))
.dtu_counts <- data.frame(gene_id = c("g1", "g1", "g2", "g2"), feature_id = c("g1|supported", "g1|unsupported", "g2|a", "g2|b"),
    A_p = c(10, 0, 10, 20), A_r = c(10, 0, 20, 10), B_p = c(10, 2, 10, 20), B_r = c(10, 2, 20, 10), C_p = c(10, 2, 10, 20), C_r = c(10, 2, 20, 10), check.names = FALSE)
.dtu_order <- .dtu_test_env$m06_reference_order(.dtu_counts, .dtu_samples)
stopifnot(identical(.dtu_order$audit$selected_reference, c("g1|supported", "g2|b")),
  setequal(.dtu_order$counts$feature_id, .dtu_counts$feature_id),
  all(colSums(.dtu_order$counts[, .dtu_samples$sample_id]) == colSums(.dtu_counts[, .dtu_samples$sample_id])),
  all(.dtu_order$audit$reference_is_best_available))

.dtu_dependencies <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "scales", "forcats", "ggplot2")
if (all(vapply(.dtu_dependencies, requireNamespace, logical(1), quietly = TRUE))) {
  .dtu_test_env <- dtu_engine_environment()
  sys.source(file.path(.repo_root, "inst", "engines", "dtu_core.R"), envir = .dtu_test_env)
  sys.source(file.path(.repo_root, "inst", "engines", "dtu_enrichment_core.R"), envir = .dtu_test_env)
  .dtu_test_env$empty_dtu_prepared <- dtu_empty_prepared
  .dtu_test_env$analysis_cohorts <- list(fixture = c("A", "B", "C"))
  .dtu_test_env$sample_manifest <- dplyr::mutate(tibble::as_tibble(.dtu_samples), sample = sample_id)
  .dtu_test_env$counts_matrix <- matrix(10L, 2L, 6L, dimnames = list(c("single1", "single2"), .dtu_samples$sample_id))
  .dtu_test_env$isoform_annotation <- tibble::tibble(isoform_harmonized_id = c("single1", "single2"), associated_gene = c("g1", "g2"),
      structural_category = "full-splice_match", gene_mapped = TRUE, cross_sample_comparable = TRUE, analyzed_for_gain_loss = TRUE)
  .dtu_test_env$cfg <- list(filters = list(minimum_isoforms_per_gene = 2L))
  .dtu_zero_family <- .dtu_test_env$prepare_dtu_family("fixture")
  stopifnot(isTRUE(.dtu_zero_family$empty), nrow(.dtu_zero_family$testing_family) == 0L)
  .dtu_test_env$cfg <- list(ora = list(minimum_genes = 2L, minimum_term_genes = 2L, maximum_term_genes = 1000L))
  .dtu_test_env$map_gene_symbols <- function(genes) tibble::tibble(gene_id = genes, SYMBOL = genes, mapping_type = "fixture_identity")
  .dtu_family <- tibble::tibble(gene_id = paste0("gene", 1:16), pooled_mean_gene_CPM = 1:16,
      detected_isoforms = rep(3:6, 4), tested_isoforms = rep(2:3, 8), gencode_annotated_transcripts = rep(6:9, 4))
  .dtu_targets <- c("gene1", "gene5", "gene9", "gene13")
  .dtu_match <- .dtu_test_env$make_matched_design(.dtu_family, .dtu_targets, 3L)
  .dtu_match_reordered <- .dtu_test_env$make_matched_design(.dtu_family[16:1, ], rev(.dtu_targets), 3L)
  stopifnot(isTRUE(.dtu_match$available), nrow(.dtu_match$pairs) == 12L,
      !anyDuplicated(.dtu_match$pairs$control_SYMBOL), !any(.dtu_match$pairs$control_SYMBOL %in% .dtu_targets),
      identical(.dtu_match$pairs, .dtu_match_reordered$pairs))
  .dtu_missing <- .dtu_family; .dtu_missing$gencode_annotated_transcripts[1] <- NA_integer_
  .dtu_unavailable <- .dtu_test_env$make_matched_design(.dtu_missing, .dtu_targets, 3L)
  stopifnot(!.dtu_unavailable$available, identical(.dtu_unavailable$reason, "incomplete_target_matching_covariates"))
  .dtu_terms <- tibble::tibble(collection = "fixture", pathway = rep(c("one", "two"), each = 4), gene = paste0("gene", 1:8))
  .dtu_ora <- .dtu_test_env$standard_ora(c("gene1", "gene2"), paste0("gene", 1:8), .dtu_terms, "fixture")
  stopifnot(nrow(.dtu_ora) == 2L, isTRUE(all.equal(.dtu_ora$ora_p_value[.dtu_ora$pathway == "one"], stats::phyper(1, 4, 4, 2, lower.tail = FALSE))),
    all(.dtu_ora$ora_bh == stats::p.adjust(.dtu_ora$ora_p_value, "BH")))

  if (requireNamespace("stageR", quietly = TRUE)) {
    .dtu_test_env$dtu_fits <- list(fixture = list(gene = tibble::tibble(gene_id = c("g1", "g2", "g3"), dtu_bh = stats::p.adjust(c(.001, .02, 1), "BH")),
      feature = tibble::tibble(gene_id = rep(c("g1", "g2", "g3"), each = 2), feature_id = paste0("tx", 1:6), isoform_dtu_p_value = c(.001, .004, .01, .04, NA, NA))))
    .dtu_test_env$dtu_prepared <- list(fixture = list(feature_map = .dtu_test_env$dtu_fits$fixture$feature[, c("gene_id", "feature_id")]))
    .dtu_stage05 <- .dtu_test_env$run_stageR_alpha("fixture", .05)
    .dtu_stage10 <- .dtu_test_env$run_stageR_alpha("fixture", .10)
    stopifnot(.dtu_stage05$screening_genes == 3L, .dtu_stage05$confirmation_family_genes == 3L,
      setequal(.dtu_stage05$feature$feature_id, paste0("tx", 1:6)), .dtu_stage05$conservative_feature_P1 == 2L,
      isTRUE(all.equal(.dtu_stage05$gene$stageR_gene_adjusted_alpha05, .dtu_stage10$gene$stageR_gene_adjusted_alpha10)))
  } else message("SKIP DTU stageR execution test: stageR unavailable")
} else message("SKIP DTU matching/stageR tests: plotting/data packages unavailable")

if (requireNamespace("DRIMSeq", quietly = TRUE) && requireNamespace("MASS", quietly = TRUE) &&
    as.character(utils::packageVersion("DRIMSeq")) == "1.38.0") {
  .dtu_hessian_before <- body(get("dm_Hessian_regG_prop", asNamespace("DRIMSeq")))
  .dtu_adapter <- .dtu_test_env$make_m06_drimseq_engine()
  .dtu_private <- environment(.dtu_adapter$dmPrecision)
  .dtu_y <- matrix(c(12, 8, 11, 9, 4, 16, 5, 15), nrow = 2L)
  .dtu_x <- cbind(1, c(0, 0, 1, 1)); .dtu_b <- c(.2, -.4); .dtu_eps <- 1e-6
  .dtu_gradient <- vapply(seq_along(.dtu_b), function(i) {
    a <- b <- .dtu_b; a[i] <- a[i] + .dtu_eps; b[i] <- b[i] - .dtu_eps
    (.dtu_private$stable_dm_lik(a, .dtu_x, 30, .dtu_y) - .dtu_private$stable_dm_lik(b, .dtu_x, 30, .dtu_y)) / (2 * .dtu_eps)
  }, numeric(1))
  stopifnot(max(abs(.dtu_gradient - .dtu_private$stable_dm_score(.dtu_b, .dtu_x, 30, .dtu_y))) < 1e-5,
    identical(.dtu_hessian_before, body(get("dm_Hessian_regG_prop", asNamespace("DRIMSeq")))))
  .dtu_boundary <- .dtu_private$dm_fitOneGroup(matrix(c(20, 0, 10, 0), nrow = 2L), 30)
  stopifnot(identical(.dtu_boundary$prop, c(1, 0)), .dtu_boundary$lik == 0,
    all(abs(.dtu_private$stable_bb_likelihood(matrix(c(20, 0, 10, 0), nrow = 2L), 30, matrix(c(1, 0, 1, 0), nrow = 2L))) < 1e-10))
} else message("SKIP DTU adapter numerical tests: exact DRIMSeq 1.38.0 unavailable")
message("DTU tests passed (available dependencies).")
