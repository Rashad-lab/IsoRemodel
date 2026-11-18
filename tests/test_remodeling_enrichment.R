.re_genes <- data.frame(gene_id = paste0("g", 1:6), architecture_class = c("A", "B", "C", "F", "U", "NQ"), quantifiable = c(rep(TRUE, 5), FALSE),
  included_in_testing_family = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), downstream_primary_remodeling_set = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  is_high_confidence_class = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), gene_model_p = c(.001, .04, .02, 1, NA, NA), gene_model_q = c(.004, .08, .04, 1, NA, NA))
.re_sets <- remodeling_enrichment_sets(.re_genes, c(.05, .10))
stopifnot(setequal(.re_sets$authoritative_biological$target, c("g1", "g2")),
  setequal(.re_sets$authoritative_biological$background, paste0("g", 1:5)),
  setequal(.re_sets$nominal_monte_carlo_all$background, paste0("g", 1:4)),
  "g4" %in% .re_sets$nominal_monte_carlo_all$background,
  !"g4" %in% .re_sets$nominal_monte_carlo_all$target,
  setequal(.re_sets$cohort_FDR05$target, c("g1", "g3")), setequal(.re_sets$cohort_FDR10$target, c("g1", "g2", "g3")),
  identical(.re_sets$strict_architecture$target, "g1"))
.re_empty <- remodeling_enrichment_sets(.re_genes[0, ], c(.05, .10))
stopifnot(all(vapply(.re_empty, function(x) !length(x$target) && !length(x$background), logical(1))))
if (requireNamespace("tibble", quietly = TRUE)) {
  .re_mapping <- data.frame(input_gene = c("g1", "g2"), SYMBOL = c("A", "B"), ENTREZID = c("1", "2"))
  .re_small <- remodeling_ora_set("g1", c("g1", "g2"), list(mapping = .re_mapping), list(enrichment = list(minimum_foreground_genes = 8L)))
  stopifnot(identical(.re_small$status$status, "skipped_too_few_mapped_target_genes"), nrow(.re_small$result) == 0L, .re_small$status$input_background_genes == 2L)
}
message("Remodeling enrichment set/background tests passed.")
