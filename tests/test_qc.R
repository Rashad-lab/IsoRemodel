# ORF-free quality checks: explicit denominators, annotation missingness,
# canonical-count conservation, paired tests and direction-aware event burden.
test_qc <- function() {
  required <- c("data.table", "ggplot2", "scales", "readr", "digest")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("QC tests require: ", paste(missing, collapse = ", "))
  stopifnot(identical(.qc_flag(c("TRUE", "false", "", "unavailable", NA)), c(TRUE, FALSE, NA, NA, NA)))
  stopifnot(.qc_wmedian(c(100, 1000, 2000), c(1, 1, 8)) == 2000,
    .qc_wmedian(c(100, 1000), c(1, 1)) == 100,
    .qc_wfraction(c(TRUE, FALSE, NA), c(2, 3, 100)) == .4,
    .qc_coverage(c(1, NA), c(1, 3)) == .25,
    is.na(.qc_wfraction(c(NA, NA), c(2, 3))))
  stopifnot(.qc_paired_test(1:4, 1:4)$p_value == 1,
    is.na(.qc_paired_test(1:3, 1:3)$p_value),
    .qc_correlation(rep(1, 4), 1:4)$status == "constant_variable",
    .qc_correlation(c(1, 2, NA), 1:3)$n == 2L)
  colors <- .qc_structural_palette(list(reporting = list(palettes = list(structural_categories = c(genic = "#123456")))), c("genic", "Unassigned"))
  stopifnot(colors[["genic"]] == "#123456", colors[["Unassigned"]] == "#B3B3B3")
  subjects <- paste0("subject", 1:4)
  manifest <- do.call(rbind, lapply(subjects, function(p) data.frame(
    sample = paste(p, c("baseline", "followup"), sep = "_"), patient = p,
    timepoint = c("primary", "recurrent"), stringsAsFactors = FALSE)))
  ids <- c("JUNC_A", "JUNC_B", "JUNC_C")
  counts <- matrix(c(20L, 20L, 5L), nrow = 3L, ncol = nrow(manifest),
    dimnames = list(ids, manifest$sample))
  annotation <- data.frame(isoform_harmonized_id = ids,
    associated_gene = c("gene1", "gene1", "gene2"),
    structural_category = c("full-splice_match", "incomplete-splice_match", NA),
    gene_mapped = TRUE, analyzed_for_gain_loss = c(TRUE, TRUE, FALSE), stringsAsFactors = FALSE)
  local <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(j) data.frame(
    manifest[rep(j, 3L), ], isoform_local = c("PB.1.1", "PB.1.2", "PB.2.1"),
    isoform_harmonized_id = ids, FL = counts[, j], associated_gene = annotation$associated_gene,
    structural_category = annotation$structural_category, length = c(2000, 5000, NA),
    subcategory = c(NA, "5prime_fragment", NA), chrom = c("chr1", "chr1", "chrM"),
    within_CAGE_peak = c("true", "false", NA), within_polyA_site = c("false", "true", NA),
    polyA_motif_found = c(NA, "false", NA), diff_to_gene_TSS = c(40, 100, NA),
    diff_to_gene_TTS = c(10, 200, NA), stringsAsFactors = FALSE)))
  inp <- list(local_transcripts = local, counts_matrix = counts,
    isoform_annotation = annotation, sample_manifest = manifest, fingerprint = "qc_synthetic_no_ORF")
  cfg <- list(technical_qc = list(terminal_proximity_nt = 50, terminus_bin_width_nt = 20,
    minimum_correlation_n = 4L, long_transcript_threshold_nt = 4000,
    minimum_chromosome_annotation_fl_coverage = .95,
    mitochondrial_contigs = c("M", "MT", "chrM", "chrMT")),
    global_qc = list(detected_min_fl = c(1, 5), complexity_support_min_fl = 3,
      complexity_support_min_subjects_in_either_condition = 4L),
    remodeling = list(definition = list(min_isoform_total_fl = 5L)),
    gain_loss = list(gene_evaluable_min_fl = 20),
    reporting = list(enabled = FALSE), comparison = list(reference = "Baseline", comparison = "Follow-up"))
  local_with_old_bins <- inp; local_with_old_bins$local_transcripts$tss_bin <- 999
  stopifnot(all(.qc_prepare_local(local_with_old_bins, cfg)$data[isoform_local == "PB.1.1"]$tss_bin == 40))
  prepared <- .qc_prepare_local(inp, cfg)
  first <- prepared$data[sample == manifest$sample[1]]
  stopifnot(first[isoform_local == "PB.1.1"]$tss_near,
    first[isoform_local == "PB.1.1"]$tss_bin == 40,
    is.na(first[isoform_local == "PB.1.1"]$both_end_flag),
    !first[isoform_local == "PB.1.2"]$both_end_flag)
  sq <- .qc_sample_metrics(prepared$data, inp, cfg)
  stopifnot(all(sq$count_total_difference == 0), all(sq$detected_isoforms_FLge1 == 3L),
    all(sq$detected_genes_FLge1 == 2L), all(sq$mitochondrial_chromosome_FL_fraction == 5 / 45),
    all(sq$transcript_length_FL_coverage == 40 / 45))
  # Absent chromosome annotations cannot be read as zero mitochondrial content.
  incomplete <- inp; incomplete$local_transcripts$chrom <- NULL
  missing_chr <- .qc_sample_metrics(.qc_prepare_local(incomplete, cfg)$data, incomplete, cfg)
  stopifnot(all(is.na(missing_chr$mitochondrial_chromosome_FL_fraction)),
    all(missing_chr$chromosome_annotation_FL_coverage == 0))
  term <- .qc_terminus(prepared$data, manifest, 20)
  stopifnot(nrow(term$events) == 4L, all(term$events$isoform_local_primary == "PB.1.1"),
    all(term$summary$tss_bin_changed_fraction == 0))
  # Unknown category remains an explicit denominator category; missing categories
  # in a library are completed to zero before condition means are calculated.
  complete <- .qc_complete_composition(data.table::data.table(sample = manifest$sample[1],
    structural_category = "A", n_isoforms = 2L, FL = 20), manifest[1:2, ], c("A", "B"))
  stopifnot(nrow(complete) == 4L, complete[sample == manifest$sample[1] & structural_category == "B"]$FL == 0,
    all(is.na(complete[sample == manifest$sample[2]]$FL_fraction)))
  # Mismatched harmonized counts fail rather than silently changing denominators.
  broken <- inp; broken$counts_matrix[1, 1] <- broken$counts_matrix[1, 1] + 1L
  stopifnot(inherits(try(.qc_sample_metrics(prepared$data, broken, cfg), silent = TRUE), "try-error"))
  # Exercise full ORF-free stage and direction-aware associations without graphics.
  inp$gain_loss <- list(full_state = list(patient_summary = data.frame(patient = subjects,
    n_gained = 1:4, n_lost = 4:1)), patient_summary = do.call(rbind, lapply(subjects, function(p)
      data.frame(patient = p, direction = c("recurrent_gained", "primary_lost"), median_events = c(1, 2)))))
  folder <- tempfile("qc_test_"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  result <- run_qc(inp, cfg, list(all_pairs = subjects, omit_last = subjects[1:3]), folder)
  stopifnot(result$complete, length(result$cohorts) == 2L,
    nrow(result$cohorts$all_pairs$paired_tests) == 16L,
    result$cohorts$all_pairs$gain_loss$available,
    nrow(result$cohorts$all_pairs$gain_loss$associations) == 84L,
    all(result$cohorts$omit_last$paired_tests$status == "insufficient_complete_pairs"),
    all(result$cohorts$all_pairs$global$support_filtered_ids %in% ids[1:2]),
    file.exists(file.path(folder, "all_pairs", "tables", "structure_by_complexity_condition.tsv")))
  message("RNA QC tests passed.")
}
test_qc()
