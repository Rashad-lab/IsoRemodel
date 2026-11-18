# RNA-only global, structure and technical quality checks.
# Derived from the original Pigeon QC and revised RNA Module 03/supplement.
# Protein/CDS/NMD fields are deliberately not prerequisites or outcomes.

.qc_numeric <- function(x) suppressWarnings(as.numeric(x))

.qc_character <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  x
}

.qc_flag <- function(x) {
  x <- tolower(.qc_character(x))
  out <- rep(NA, length(x))
  out[!is.na(x) & x %in% c("true", "t", "1", "yes", "y", "inside", "within")] <- TRUE
  out[!is.na(x) & x %in% c("false", "f", "0", "no", "n", "outside")] <- FALSE
  out
}

.qc_token <- function(x) {
  x <- tolower(.qc_character(x))
  gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", x))
}

.qc_wmean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

.qc_wmedian <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  x <- x[keep]; w <- w[keep]
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1L]]
}

.qc_wfraction <- function(flag, w) {
  keep <- !is.na(flag) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(w[keep] * as.numeric(flag[keep])) / sum(w[keep])
}

.qc_coverage <- function(value, w) {
  keep <- is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(w[keep & !is.na(value)]) / sum(w[keep])
}

.qc_safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

.qc_correlation <- function(x, y, minimum_n = 4L, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  n <- sum(keep)
  status <- if (n < minimum_n) "insufficient_complete_pairs" else if (
    length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L
  ) "constant_variable" else "estimated"
  if (status != "estimated") return(data.frame(n = n, estimate = NA_real_, p_value = NA_real_, status = status))
  test <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = method, exact = FALSE))
  data.frame(n = n, estimate = unname(test$estimate), p_value = test$p.value, status = status)
}

.qc_paired_test <- function(primary, recurrent, minimum_n = 4L) {
  keep <- is.finite(primary) & is.finite(recurrent)
  delta <- recurrent[keep] - primary[keep]
  n <- length(delta)
  p <- if (n < minimum_n) NA_real_ else if (all(delta == 0)) 1 else
    suppressWarnings(stats::wilcox.test(recurrent[keep], primary[keep], paired = TRUE, exact = FALSE)$p.value)
  data.frame(n_pairs = n,
    median_primary = .qc_safe_median(primary[keep]),
    median_recurrent = .qc_safe_median(recurrent[keep]),
    median_delta = .qc_safe_median(delta),
    mean_delta = if (n) mean(delta) else NA_real_,
    n_positive = sum(delta > 0), n_negative = sum(delta < 0), n_zero = sum(delta == 0),
    p_value = p, status = if (n < minimum_n) "insufficient_complete_pairs" else "estimated")
}

.qc_dictionary <- function(proximity = 50, long_nt = 4000) {
  metric <- c("total_full_length_reads", "log10_total_full_length_reads", "detected_genes_FLge1",
    "detected_isoforms_FLge1", "median_transcript_length_FL_weighted", "long_transcript_FL_fraction",
    "incomplete_splice_match_FL_fraction", "ism_5prime_fragment_FL_fraction",
    "ism_3prime_fragment_FL_fraction", "ism_internal_fragment_FL_fraction",
    "within_CAGE_peak_FL_fraction", "tss_within_threshold_FL_fraction",
    "within_polyA_site_FL_fraction", "tts_within_threshold_FL_fraction",
    "both_end_supported_FL_fraction", "mitochondrial_chromosome_FL_fraction")
  labels <- c("Full-length library depth", "Log10 full-length library depth", "Detected genes",
    "Detected splice chains", "FL-weighted median transcript length",
    paste0("FL fraction from transcripts >=", long_nt, " nt"), "Incomplete-splice-match fraction",
    "5-prime-fragment ISM fraction", "3-prime-fragment ISM fraction", "Internal-fragment ISM fraction",
    "Within-CAGE-peak fraction", paste0("Gene-TSS proximity <=", proximity, " nt"),
    "Within-poly(A)-site fraction", paste0("Gene-TTS proximity <=", proximity, " nt"),
    "Both-end-support fraction", "Mitochondrial chromosome read fraction")
  role <- c("technical_predictor", "technical_predictor", "descriptive_qc", "descriptive_qc",
    "biological_outcome", "biological_outcome", "technical_predictor", "technical_predictor",
    "secondary_technical", "secondary_technical", "technical_predictor", "technical_predictor",
    "secondary_technical", "secondary_technical", "secondary_technical", "technical_predictor")
  denominator <- c("All positive local FL counts; conserved by harmonized matrix",
    "All positive local FL counts", "Mapped genes with positive harmonized FL", "Harmonized rows with positive FL",
    "Positive FL with a valid transcript length", "Positive FL with a valid transcript length",
    "Positive FL with a known structural category", rep("Positive FL with a known fragment status", 3),
    "Positive FL with a known CAGE flag", "Positive FL with a known gene-TSS distance",
    "Positive FL with a known poly(A)-site flag", "Positive FL with a known gene-TTS distance",
    "Positive FL with a determinate combined end-support status",
    "All positive FL; unavailable below minimum chromosome-annotation coverage")
  data.frame(metric = metric, display_label = labels, role = role,
    native_unit = c("reads", "log10_reads", "genes", "splice_chains", "nt", rep("fraction", 11)),
    denominator = denominator,
    plot_multiplier = c(1e-6, 1, 1, 1, .001, rep(100, 11)), stringsAsFactors = FALSE)
}

.qc_prepare_local <- function(inputs, cfg) {
  x <- data.table::as.data.table(data.table::copy(inputs$local_transcripts))
  required <- c("sample", "patient", "timepoint", "isoform_local", "isoform_harmonized_id", "associated_gene", "FL")
  if (length(setdiff(required, names(x)))) rna_stop("QC local table lacks: ", paste(setdiff(required, names(x)), collapse = ", "))
  optional_numeric <- c("length", "diff_to_gene_TSS", "diff_to_gene_TTS", "dist_to_CAGE_peak", "dist_to_polyA_site", "tss_bin", "tts_bin")
  optional_character <- c("chrom", "structural_category", "subcategory", "within_CAGE_peak", "within_polyA_site", "polyA_motif_found", "all_canonical", "RTS_stage")
  availability <- data.frame(column = c(optional_numeric, optional_character),
    present_in_input = c(optional_numeric, optional_character) %in% names(x), stringsAsFactors = FALSE)
  for (name in setdiff(optional_numeric, names(x))) x[, (name) := NA_real_]
  for (name in setdiff(optional_character, names(x))) x[, (name) := NA_character_]
  for (name in c(optional_numeric, "FL")) x[, (name) := .qc_numeric(get(name))]
  for (name in c(optional_character, "sample", "patient", "timepoint", "isoform_local", "isoform_harmonized_id", "associated_gene")) x[, (name) := .qc_character(get(name))]
  if (any(!is.finite(x$FL) | x$FL < 0)) rna_stop("QC requires complete nonnegative FL counts.")
  availability$n_positive_FL_rows_with_value <- vapply(availability$column,
    function(name) sum(x$FL > 0 & !is.na(x[[name]])), integer(1))
  availability$positive_FL_coverage <- vapply(availability$column,
    function(name) .qc_coverage(x[[name]], x$FL), numeric(1))
  if (anyNA(x$sample) || anyNA(x$isoform_local) || anyDuplicated(x[, c("sample", "isoform_local"), with = FALSE])) rna_stop("QC local transcript keys are missing or duplicated.")
  q <- cfg$technical_qc
  proximity <- q$terminal_proximity_nt %||% 50
  bin_width <- q$terminus_bin_width_nt %||% 50
  x[!is.finite(length) | length <= 0, length := NA_real_]
  x[, structural_category_normalized := .qc_token(structural_category)]
  x[, subcategory_normalized := .qc_token(subcategory)]
  x[, is_ism := ifelse(is.na(structural_category_normalized), NA,
    structural_category_normalized == "incomplete_splice_match")]
  for (term in c("5prime", "3prime", "internal")) {
    pattern <- switch(term, `5prime` = "(^|_)5_?prime(_|$)|(^|_)5p(_|$)|5prime_fragment|five_prime_fragment",
      `3prime` = "(^|_)3_?prime(_|$)|(^|_)3p(_|$)|3prime_fragment|three_prime_fragment", internal = "internal_fragment")
    value <- ifelse(is.na(x$is_ism), NA, ifelse(!x$is_ism, FALSE,
      ifelse(is.na(x$subcategory_normalized), NA, grepl(pattern, x$subcategory_normalized))))
    x[, (paste0("is_ism_", term, "_fragment")) := value]
  }
  x[, cage_flag := .qc_flag(within_CAGE_peak)]
  x[, polya_flag := .qc_flag(within_polyA_site)]
  x[, motif_flag := .qc_flag(polyA_motif_found)]
  x[, tss_near := ifelse(is.finite(diff_to_gene_TSS), abs(diff_to_gene_TSS) <= proximity, NA)]
  x[, tts_near := ifelse(is.finite(diff_to_gene_TTS), abs(diff_to_gene_TTS) <= proximity, NA)]
  # Preserve three-valued logic: TRUE OR NA is TRUE, FALSE OR NA is NA.
  x[, both_end_flag := cage_flag & (polya_flag | motif_flag)]
  mito_contigs <- tolower(as.character(q$mitochondrial_contigs %||% c("M", "MT", "chrM", "chrMT")))
  x[, mito_flag := ifelse(is.na(chrom), NA, tolower(chrom) %in% mito_contigs)]
  # QC terminus resolution is editable independently of import end-ID bins.
  x[is.finite(diff_to_gene_TSS), tss_bin := round(diff_to_gene_TSS / bin_width) * bin_width]
  x[is.finite(diff_to_gene_TTS), tts_bin := round(diff_to_gene_TTS / bin_width) * bin_width]
  list(data = x[FL > 0], availability = availability)
}

.qc_sample_metrics <- function(x, inputs, cfg) {
  long_nt <- cfg$technical_qc$long_transcript_threshold_nt %||% 4000
  minimum_chrom_coverage <- cfg$technical_qc$minimum_chromosome_annotation_fl_coverage %||% .95
  y <- x[, {
    valid_gene <- !is.na(associated_gene) & associated_gene != "unknown"
    ism_weight <- sum(FL[is_ism %in% TRUE])
    total <- sum(FL)
    chrom_coverage <- .qc_coverage(chrom, FL)
    list(total_full_length_reads = total,
      log10_total_full_length_reads = if (total > 0) log10(total) else NA_real_,
      detected_genes_local_annotations = data.table::uniqueN(associated_gene[valid_gene]),
      detected_isoforms_local = data.table::uniqueN(isoform_local),
      median_transcript_length_FL_weighted = .qc_wmedian(length, FL),
      mean_transcript_length_FL_weighted = .qc_wmean(length, FL),
      median_detected_isoform_length_unweighted = .qc_safe_median(length[length > 0]),
      transcript_length_FL_coverage = .qc_coverage(ifelse(is.finite(length) & length > 0, length, NA_real_), FL),
      long_transcript_FL_fraction = .qc_wfraction(ifelse(is.finite(length) & length > 0, length >= long_nt, NA), FL),
      incomplete_splice_match_FL_fraction = .qc_wfraction(is_ism, FL),
      structural_annotation_FL_coverage = .qc_coverage(structural_category, FL),
      ISM_subcategory_FL_coverage = if (ism_weight > 0) sum(FL[is_ism %in% TRUE & !is.na(subcategory)]) / ism_weight else NA_real_,
      ism_5prime_fragment_FL_fraction = .qc_wfraction(is_ism_5prime_fragment, FL),
      ism_3prime_fragment_FL_fraction = .qc_wfraction(is_ism_3prime_fragment, FL),
      ism_internal_fragment_FL_fraction = .qc_wfraction(is_ism_internal_fragment, FL),
      within_CAGE_peak_FL_fraction = .qc_wfraction(cage_flag, FL),
      CAGE_annotation_FL_coverage = .qc_coverage(cage_flag, FL),
      median_abs_CAGE_distance_FL_weighted = .qc_wmedian(abs(dist_to_CAGE_peak), FL),
      tss_within_threshold_FL_fraction = .qc_wfraction(tss_near, FL),
      TSS_annotation_FL_coverage = .qc_coverage(diff_to_gene_TSS, FL),
      within_polyA_site_FL_fraction = .qc_wfraction(polya_flag, FL),
      polyA_site_annotation_FL_coverage = .qc_coverage(polya_flag, FL),
      polyA_motif_FL_fraction = .qc_wfraction(motif_flag, FL),
      polyA_motif_annotation_FL_coverage = .qc_coverage(motif_flag, FL),
      median_abs_polyA_distance_FL_weighted = .qc_wmedian(abs(dist_to_polyA_site), FL),
      tts_within_threshold_FL_fraction = .qc_wfraction(tts_near, FL),
      TTS_annotation_FL_coverage = .qc_coverage(diff_to_gene_TTS, FL),
      both_end_supported_FL_fraction = .qc_wfraction(both_end_flag, FL),
      both_end_annotation_FL_coverage = .qc_coverage(both_end_flag, FL),
      mitochondrial_full_length_reads = if (any(!is.na(mito_flag))) sum(FL[mito_flag %in% TRUE]) else NA_real_,
      mitochondrial_chromosome_FL_fraction = if (is.finite(chrom_coverage) && chrom_coverage >= minimum_chrom_coverage)
        sum(FL[mito_flag %in% TRUE]) / total else NA_real_,
      chromosome_annotation_FL_coverage = chrom_coverage)
  }, by = sample]
  manifest <- data.table::as.data.table(inputs$sample_manifest)
  y <- merge(manifest[, .(sample, patient, timepoint)], y, by = "sample", all.x = TRUE, sort = FALSE)
  counts <- inputs$counts_matrix
  annotation <- inputs$isoform_annotation
  gene_ok <- annotation$gene_mapped %in% TRUE
  genes <- unique(annotation$associated_gene[gene_ok])
  if (length(genes)) {
    gene_counts <- rowsum(counts[gene_ok, , drop = FALSE], group = annotation$associated_gene[gene_ok], reorder = FALSE)
    detected_genes <- colSums(gene_counts > 0)
  } else detected_genes <- setNames(rep(0, ncol(counts)), colnames(counts))
  y[, detected_genes_FLge1 := as.integer(detected_genes[match(sample, names(detected_genes))])]
  y[, detected_isoforms_FLge1 := as.integer(colSums(counts > 0)[match(sample, colnames(counts))])]
  y[, count_matrix_total_FL := as.numeric(colSums(counts)[match(sample, colnames(counts))])]
  y[, count_total_difference := total_full_length_reads - count_matrix_total_FL]
  if (any(!is.finite(y$count_total_difference)) || any(y$count_total_difference != 0)) rna_stop("QC local/master count totals disagree.")
  data.table::setorderv(y, c("patient", "timepoint"))
  y
}

.qc_terminus <- function(x, samples, gene_floor) {
  valid <- x[!is.na(associated_gene) & associated_gene != "unknown"]
  if (!nrow(valid)) return(list(events = data.table::data.table(), summary = data.frame(patient = unique(samples$patient))))
  local <- valid[, {
    w <- FL
    list(FL = sum(w), transcript_length_nt = .qc_wmedian(length, w),
      tss_position_nt = .qc_wmedian(diff_to_gene_TSS, w), tts_position_nt = .qc_wmedian(diff_to_gene_TTS, w),
      tss_bin_nt = .qc_wmedian(tss_bin, w), tts_bin_nt = .qc_wmedian(tts_bin, w))
  }, by = .(sample, patient, timepoint, associated_gene, isoform_local)]
  gene <- local[, .(gene_FL = sum(FL)), by = .(patient, timepoint, associated_gene)]
  totals <- data.table::dcast(gene, patient + associated_gene ~ timepoint, value.var = "gene_FL", fill = 0)
  for (name in setdiff(c("primary", "recurrent"), names(totals))) totals[, (name) := 0]
  totals <- totals[primary >= gene_floor & recurrent >= gene_floor]
  data.table::setorderv(local, c("patient", "associated_gene", "timepoint", "FL", "isoform_local"), c(1L, 1L, 1L, -1L, 1L))
  dominant <- local[, .SD[1L], by = .(patient, associated_gene, timepoint)]
  p <- data.table::copy(dominant[timepoint == "primary"]); p[, timepoint := NULL]
  r <- data.table::copy(dominant[timepoint == "recurrent"]); r[, timepoint := NULL]
  other <- setdiff(names(p), c("patient", "associated_gene"))
  data.table::setnames(p, other, paste0(other, "_primary")); data.table::setnames(r, other, paste0(other, "_recurrent"))
  events <- merge(p, r, by = c("patient", "associated_gene"))
  events <- merge(events, totals, by = c("patient", "associated_gene"))
  events[, `:=`(tss_comparable = is.finite(tss_bin_nt_primary) & is.finite(tss_bin_nt_recurrent),
    tts_comparable = is.finite(tts_bin_nt_primary) & is.finite(tts_bin_nt_recurrent))]
  events[, `:=`(tss_bin_changed = ifelse(tss_comparable, tss_bin_nt_primary != tss_bin_nt_recurrent, NA),
    tts_bin_changed = ifelse(tts_comparable, tts_bin_nt_primary != tts_bin_nt_recurrent, NA),
    signed_tss_shift_nt = tss_position_nt_recurrent - tss_position_nt_primary,
    signed_tts_shift_nt = tts_position_nt_recurrent - tts_position_nt_primary,
    signed_transcript_length_change_nt = transcript_length_nt_recurrent - transcript_length_nt_primary)]
  events[, any_terminus_comparable := tss_comparable | tts_comparable]
  events[, any_terminus_changed := ifelse(any_terminus_comparable,
    tss_bin_changed %in% TRUE | tts_bin_changed %in% TRUE, NA)]
  summary <- data.table::rbindlist(lapply(unique(samples$patient), function(id) {
    z <- events[patient == id]
    n_tss <- sum(z$tss_comparable); n_tts <- sum(z$tts_comparable); n_any <- sum(z$any_terminus_comparable)
    data.frame(patient = id, n_evaluable_genes = nrow(z), n_tss_comparable = n_tss,
      n_tss_bin_changed = sum(z$tss_bin_changed %in% TRUE),
      tss_bin_changed_fraction = if (n_tss) sum(z$tss_bin_changed %in% TRUE) / n_tss else NA_real_,
      median_abs_tss_shift_nt = .qc_safe_median(abs(z$signed_tss_shift_nt)),
      n_tts_comparable = n_tts, n_tts_bin_changed = sum(z$tts_bin_changed %in% TRUE),
      tts_bin_changed_fraction = if (n_tts) sum(z$tts_bin_changed %in% TRUE) / n_tts else NA_real_,
      median_abs_tts_shift_nt = .qc_safe_median(abs(z$signed_tts_shift_nt)),
      n_any_terminus_comparable = n_any, n_any_terminus_changed = sum(z$any_terminus_changed %in% TRUE),
      any_terminus_changed_fraction = if (n_any) sum(z$any_terminus_changed %in% TRUE) / n_any else NA_real_,
      median_signed_dominant_transcript_length_change_nt = .qc_safe_median(z$signed_transcript_length_change_nt))
  }), fill = TRUE)
  list(events = events, summary = summary)
}

.qc_complete_composition <- function(observed, samples, categories) {
  grid <- data.table::CJ(sample = as.character(samples$sample), structural_category = categories, unique = TRUE)
  x <- merge(grid, observed, by = c("sample", "structural_category"), all.x = TRUE, sort = FALSE)
  x[is.na(n_isoforms), n_isoforms := 0L]; x[is.na(FL), FL := 0]
  x[, `:=`(isoform_fraction = if (sum(n_isoforms) > 0) n_isoforms / sum(n_isoforms) else NA_real_,
    FL_fraction = if (sum(FL) > 0) FL / sum(FL) else NA_real_), by = sample]
  merge(x, data.table::as.data.table(samples)[, .(sample, patient, timepoint)], by = "sample", all.x = TRUE)
}

.qc_chain_universe <- function(inputs, keep_rows = rep(TRUE, nrow(inputs$counts_matrix))) {
  counts <- inputs$counts_matrix
  anno <- data.table::as.data.table(inputs$isoform_annotation)
  if (!identical(rownames(counts), as.character(anno$isoform_harmonized_id))) rna_stop("QC count rows and annotations are not aligned.")
  if (!any(keep_rows)) return(data.table::data.table(sample = character(), isoform_harmonized_id = character(), gene_id = character(), structural_category = character(), FL = numeric()))
  out <- data.table::rbindlist(lapply(seq_len(ncol(counts)), function(j) {
    idx <- which(keep_rows & counts[, j] > 0)
    data.table::data.table(sample = colnames(counts)[j],
      isoform_harmonized_id = rownames(counts)[idx], gene_id = anno$associated_gene[idx],
      structural_category = as.character(anno$structural_category[idx]), FL = as.numeric(counts[idx, j]))
  }), use.names = TRUE)
  out[is.na(structural_category) | !nzchar(structural_category), structural_category := "Unassigned"]
  out
}

.qc_global <- function(inputs, local, cfg, samples, dirs) {
  annotation <- inputs$isoform_annotation
  counts <- inputs$counts_matrix
  threshold <- cfg$remodeling$definition$min_isoform_total_fl %||% 5
  keep <- annotation$analyzed_for_gain_loss %in% TRUE & rowSums(counts) >= threshold
  universes <- list(all_harmonized = .qc_chain_universe(inputs),
    support_filtered = .qc_chain_universe(inputs, keep))
  universes$raw_local <- local[, .(sample, isoform_harmonized_id = isoform_local,
    gene_id = associated_gene, structural_category = ifelse(is.na(structural_category), "Unassigned", structural_category), FL)]
  composition <- complexity <- distribution <- category_complexity <- list()
  for (name in names(universes)) {
    z <- universes[[name]][sample %in% samples$sample]
    category <- sort(unique(c("Unassigned", z$structural_category)))
    raw_comp <- z[, .(n_isoforms = .N, FL = sum(FL)), by = .(sample, structural_category)]
    comp <- .qc_complete_composition(raw_comp, samples, category); comp[, universe := name]
    composition[[name]] <- comp
    g <- z[!is.na(gene_id) & nzchar(gene_id) & gene_id != "unknown", .(n_isoforms = .N, gene_FL = sum(FL)), by = .(sample, gene_id)]
    g[, complexity := cut(n_isoforms, breaks = c(0, 1, 2, 4, 8, 20, Inf), labels = c("1", "2", "3-4", "5-8", "9-20", ">20"))]
    g <- merge(g, data.table::as.data.table(samples)[, .(sample, patient, timepoint)], by = "sample", all.x = TRUE)
    g[, universe := name]; complexity[[name]] <- g
    dist <- g[, .(n_genes = .N), by = .(sample, complexity)]
    grid <- data.table::CJ(sample = as.character(samples$sample), complexity = c("1", "2", "3-4", "5-8", "9-20", ">20"), unique = TRUE)
    dist[, complexity := as.character(complexity)]
    dist <- merge(grid, dist, by = c("sample", "complexity"), all.x = TRUE)
    dist[is.na(n_genes), n_genes := 0L]
    dist[, fraction := if (sum(n_genes)) n_genes / sum(n_genes) else NA_real_, by = sample]
    dist <- merge(dist, data.table::as.data.table(samples)[, .(sample, patient, timepoint)], by = "sample")
    dist[, universe := name]; distribution[[name]] <- dist
    if (nrow(g)) {
      with_comp <- merge(z, g[, .(sample, gene_id, complexity)], by = c("sample", "gene_id"))
      cc <- with_comp[, .(n_isoforms = .N, FL = sum(FL)), by = .(sample, complexity, structural_category)]
      existing <- unique(g[, .(sample, complexity)]); existing[, qc_cross_key := 1L]
      grid_cc <- merge(existing, data.table::data.table(structural_category = category, qc_cross_key = 1L), by = "qc_cross_key", allow.cartesian = TRUE)
      grid_cc[, qc_cross_key := NULL]
      cc <- merge(grid_cc, cc, by = c("sample", "complexity", "structural_category"), all.x = TRUE)
      cc[is.na(n_isoforms), n_isoforms := 0L]; cc[is.na(FL), FL := 0]
      cc[, `:=`(isoform_fraction = n_isoforms / sum(n_isoforms), FL_fraction = FL / sum(FL)), by = .(sample, complexity)]
      cc <- merge(cc, data.table::as.data.table(samples)[, .(sample, patient, timepoint)], by = "sample")
      cc[, universe := name]; category_complexity[[name]] <- cc
    }
  }
  comp <- data.table::rbindlist(composition, fill = TRUE)
  group <- comp[, .(mean_isoform_fraction = mean(isoform_fraction, na.rm = TRUE),
    sd_isoform_fraction = stats::sd(isoform_fraction, na.rm = TRUE),
    mean_FL_fraction = mean(FL_fraction, na.rm = TRUE), sd_FL_fraction = stats::sd(FL_fraction, na.rm = TRUE),
    pooled_isoform_observations = sum(n_isoforms), pooled_FL = sum(FL), n_samples = .N),
    by = .(universe, timepoint, structural_category)]
  group[, `:=`(pooled_isoform_fraction = if (sum(pooled_isoform_observations)) pooled_isoform_observations / sum(pooled_isoform_observations) else NA_real_,
    pooled_FL_fraction = if (sum(pooled_FL)) pooled_FL / sum(pooled_FL) else NA_real_), by = .(universe, timepoint)]
  # Original 06B support/complexity QC intentionally uses only local rows with
  # the canonical winning original gene label. Core counts remain unchanged.
  map <- data.table::as.data.table(inputs$isoform_annotation)[gene_mapped %in% TRUE, .(isoform_harmonized_id, associated_gene, structural_category)]
  matched <- merge(local, map[, .(isoform_harmonized_id, associated_gene)], by = c("isoform_harmonized_id", "associated_gene"))
  iso_sample <- matched[, .(FL = sum(FL)), by = .(isoform_harmonized_id, associated_gene, sample, patient, timepoint)]
  min_fl <- cfg$global_qc$complexity_support_min_fl %||% 3
  min_patients <- cfg$global_qc$complexity_support_min_subjects_in_either_condition %||% 4
  support <- iso_sample[sample %in% samples$sample, .(n_patients_supported = data.table::uniqueN(patient[FL >= min_fl]), total_FL = sum(FL)),
    by = .(isoform_harmonized_id, associated_gene, timepoint)]
  support[, passes_condition_support := n_patients_supported >= min_patients]
  supported <- support[, .(passes_support = any(passes_condition_support)), by = .(isoform_harmonized_id, associated_gene)]
  legacy_all <- unique(iso_sample[sample %in% samples$sample & FL > 0, .(isoform_harmonized_id, associated_gene)])
  legacy_complexity <- legacy_all[, .(n_isoforms_all_detected = .N), by = associated_gene]
  supported_gene <- supported[passes_support == TRUE, .(n_isoforms_supported = .N), by = associated_gene]
  legacy_complexity <- merge(legacy_complexity, supported_gene, by = "associated_gene", all.x = TRUE)
  legacy_complexity[is.na(n_isoforms_supported), n_isoforms_supported := 0L]
  legacy_gene_fl <- iso_sample[sample %in% samples$sample, .(gene_FL = sum(FL)), by = associated_gene]
  legacy_complexity <- merge(legacy_complexity, legacy_gene_fl, by = "associated_gene", all.x = TRUE)
  correlations <- data.table::rbindlist(lapply(c("n_isoforms_all_detected", "n_isoforms_supported"), function(col) {
    positive <- legacy_complexity[gene_FL > 0 & get(col) > 0]
    data.table::rbindlist(lapply(c("spearman", "pearson_log10_plus1"), function(method) {
      res <- .qc_correlation(if (method == "spearman") positive$gene_FL else log10(positive$gene_FL + 1),
        if (method == "spearman") positive[[col]] else log10(positive[[col]] + 1), 3L,
        if (method == "spearman") "spearman" else "pearson")
      data.frame(universe = col, method = method, res)
    }))
  }))
  legacy_structure <- data.table::rbindlist(lapply(c("all_detected", "supported"), function(universe) {
    iso <- merge(legacy_all, map, by = c("isoform_harmonized_id", "associated_gene"))
    if (universe == "supported") iso <- merge(iso, supported[passes_support == TRUE], by = c("isoform_harmonized_id", "associated_gene"))
    number <- iso[, .(n_isoforms = .N), by = associated_gene]
    number[, complexity := cut(n_isoforms, breaks = c(0, 1, 2, 5, 10, 20, 50, 100, Inf),
      labels = c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"))]
    iso <- merge(iso, number, by = "associated_gene")
    iso[is.na(structural_category) | !nzchar(structural_category), structural_category := "Unassigned"]
    out <- iso[, .(n_isoforms = .N, n_genes_with_category = data.table::uniqueN(associated_gene)), by = .(complexity, structural_category)]
    out[, isoform_fraction := n_isoforms / sum(n_isoforms), by = complexity]
    out[, universe := universe]
    out
  }), fill = TRUE)
  rna_write_tsv(legacy_structure, file.path(dirs$tables, "legacy_structure_by_gene_complexity.tsv"))
  rna_write_tsv(comp, file.path(dirs$tables, "structural_composition_by_sample.tsv.gz"))
  rna_write_tsv(group, file.path(dirs$tables, "structural_composition_by_condition.tsv"))
  rna_write_tsv(data.table::rbindlist(complexity, fill = TRUE), file.path(dirs$tables, "gene_sample_complexity_and_support.tsv.gz"))
  rna_write_tsv(data.table::rbindlist(distribution, fill = TRUE), file.path(dirs$tables, "complexity_distribution.tsv"))
  cc_all <- data.table::rbindlist(category_complexity, fill = TRUE)
  if (!length(names(cc_all))) cc_all <- data.table::data.table(sample = character(), complexity = character(), structural_category = character(), n_isoforms = integer(), FL = numeric(), isoform_fraction = numeric(), FL_fraction = numeric(), patient = character(), timepoint = character(), universe = character())
  rna_write_tsv(cc_all, file.path(dirs$tables, "structure_by_complexity_sample.tsv.gz"))
  cc_group <- cc_all[, .(mean_isoform_fraction = mean(isoform_fraction), mean_FL_fraction = mean(FL_fraction), n_samples_with_complexity = .N), by = .(universe, timepoint, complexity, structural_category)]
  rna_write_tsv(cc_group, file.path(dirs$tables, "structure_by_complexity_condition.tsv"))
  rna_write_tsv(support, file.path(dirs$tables, "legacy_support_by_condition.tsv.gz"))
  rna_write_tsv(legacy_complexity, file.path(dirs$tables, "legacy_gene_complexity.tsv"))
  rna_write_tsv(correlations, file.path(dirs$tables, "legacy_complexity_count_correlations.tsv"))
  rna_write_tsv(data.frame(metric = c("local_FL_total", "winning_original_gene_row_FL", "excluded_from_legacy_complexity_FL"),
    value = c(sum(local[sample %in% samples$sample]$FL), sum(matched[sample %in% samples$sample]$FL),
      sum(local[sample %in% samples$sample]$FL) - sum(matched[sample %in% samples$sample]$FL))),
    file.path(dirs$qc, "legacy_complexity_annotation_selection_audit.tsv"))
  list(composition = comp, condition_composition = group, gene_sample_complexity = data.table::rbindlist(complexity, fill = TRUE),
    complexity_distribution = data.table::rbindlist(distribution, fill = TRUE), legacy_complexity = legacy_complexity, legacy_structure = legacy_structure,
    support = support, structural_by_complexity = cc_all, structural_by_complexity_condition = cc_group, support_filtered_ids = rownames(counts)[keep])
}

.qc_association_grid <- function(data, predictors, outcomes, cohort, scope, minimum_n) {
  if (!length(predictors) || !length(outcomes)) return(data.table::data.table())
  result <- data.table::rbindlist(lapply(predictors, function(pred) {
    data.table::rbindlist(lapply(outcomes, function(outcome) {
      data.frame(cohort = cohort, scope = scope, predictor = pred, outcome = outcome,
        .qc_correlation(data[[pred]], data[[outcome]], minimum_n))
    }))
  }))
  result[, FDR_within_outcome := stats::p.adjust(p_value, method = "BH"), by = outcome]
  result[, FDR_global := stats::p.adjust(p_value, method = "BH")]
  result
}

.qc_gain_loss_associations <- function(gain_loss, sample_qc, cohort, subjects, cfg, dirs) {
  # A missing optional upstream result is reported, never interpreted as no events.
  if (is.null(gain_loss) || is.null(gain_loss$full_state$patient_summary)) {
    availability <- data.frame(cohort = cohort, available = FALSE,
      reason = "Gain/loss result not supplied to RNA QC; event associations unavailable.")
    rna_write_tsv(availability, file.path(dirs$qc, "gain_loss_association_availability.tsv"))
    return(list(available = FALSE, burden = data.table::data.table(),
      design = data.table::data.table(), associations = data.table::data.table()))
  }
  raw <- data.table::as.data.table(gain_loss$full_state$patient_summary)
  required <- c("patient", "n_gained", "n_lost")
  if (length(setdiff(required, names(raw)))) rna_stop("Gain/loss patient summary lacks RNA QC fields: ", paste(setdiff(required, names(raw)), collapse = ", "))
  raw <- raw[patient %in% subjects]
  burden <- data.table::rbindlist(list(
    raw[, .(patient, direction = "recurrent_gained", n_events = as.numeric(n_gained), burden_basis = "full_depth")],
    raw[, .(patient, direction = "primary_lost", n_events = as.numeric(n_lost), burden_basis = "full_depth")]))
  common <- data.table::as.data.table(gain_loss$patient_summary)
  if (all(c("patient", "direction", "median_events") %in% names(common))) {
    common <- common[patient %in% subjects]
    burden <- data.table::rbindlist(list(burden,
      common[, .(patient, direction, n_events = as.numeric(median_events), burden_basis = "common_depth_median")]), use.names = TRUE)
  }
  burden[, `:=`(cohort = cohort,
    detection_timepoint = ifelse(direction == "recurrent_gained", "recurrent", "primary"),
    absence_timepoint = ifelse(direction == "recurrent_gained", "primary", "recurrent"))]
  predictors <- c("log10_total_full_length_reads", "median_transcript_length_FL_weighted",
    "incomplete_splice_match_FL_fraction", "ism_5prime_fragment_FL_fraction",
    "within_CAGE_peak_FL_fraction", "tss_within_threshold_FL_fraction", "mitochondrial_chromosome_FL_fraction")
  predictors <- intersect(predictors, names(sample_qc))
  metrics <- data.table::melt(data.table::copy(sample_qc),
    id.vars = c("sample", "patient", "timepoint"), measure.vars = predictors,
    variable.name = "metric", value.name = "metric_value", variable.factor = FALSE)
  burden[, qc_cross_key := 1L]
  design <- merge(burden, data.table::data.table(metric = predictors, qc_cross_key = 1L), by = "qc_cross_key", allow.cartesian = TRUE)
  burden[, qc_cross_key := NULL]; design[, qc_cross_key := NULL]
  detection <- metrics[, .(patient, detection_timepoint = timepoint, metric, detection_value = metric_value)]
  absence <- metrics[, .(patient, absence_timepoint = timepoint, metric, absence_value = metric_value)]
  design <- merge(design, detection, by = c("patient", "detection_timepoint", "metric"), all.x = TRUE)
  design <- merge(design, absence, by = c("patient", "absence_timepoint", "metric"), all.x = TRUE)
  design[, detection_minus_absence := detection_value - absence_value]
  long <- data.table::melt(design, measure.vars = c("detection_value", "absence_value", "detection_minus_absence"),
    variable.name = "technical_role", value.name = "technical_value", variable.factor = FALSE)
  tests <- long[, .qc_correlation(technical_value, n_events, cfg$technical_qc$minimum_correlation_n %||% 4L),
    by = .(cohort, burden_basis, direction, metric, technical_role)]
  tests[, FDR_within_burden_direction_role := stats::p.adjust(p_value, "BH"), by = .(burden_basis, direction, technical_role)]
  tests[, FDR_global := stats::p.adjust(p_value, "BH")]
  tests[, adjustment_family := "RNA_only_direction_aware_gain_loss_associations_within_cohort"]
  rna_write_tsv(data.frame(cohort = cohort, available = TRUE, reason = "Upstream full-depth and available common-depth event counts."),
    file.path(dirs$qc, "gain_loss_association_availability.tsv"))
  rna_write_tsv(burden, file.path(dirs$tables, "gain_loss_burden_inputs.tsv"))
  rna_write_tsv(design, file.path(dirs$tables, "gain_loss_technical_design.tsv"))
  rna_write_tsv(tests, file.path(dirs$tables, "gain_loss_technical_associations.tsv"))
  if (!isFALSE(cfg$reporting$enabled)) {
    labels <- setNames(.qc_dictionary(cfg$technical_qc$terminal_proximity_nt %||% 50)$display_label,
      .qc_dictionary(cfg$technical_qc$terminal_proximity_nt %||% 50)$metric)
    for (basis in unique(tests$burden_basis)) for (role in unique(tests$technical_role)) {
      z <- data.table::copy(tests[burden_basis == basis & technical_role == role & is.finite(estimate)])
      if (!nrow(z)) next
      z[, metric_label := vapply(labels[metric], function(v) paste(strwrap(v, 29), collapse = "\n"), character(1))]
      z[, direction_label := ifelse(direction == "recurrent_gained", "Gained", "Lost")]
      role_label <- switch(role, detection_value = "Detection-side library", absence_value = "Absence-side library",
        detection_minus_absence = "Detection minus absence")
      p <- ggplot2::ggplot(z, ggplot2::aes(x = direction_label, y = metric_label, fill = estimate)) +
        ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", estimate)), size = 4) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Technical properties and isoform events",
          subtitle = paste(if (basis == "full_depth") "Full depth" else "Common depth", role_label, sep = "; ")) + rna_theme(cfg)
      rna_save_plot(p, paste0("gain_loss_technical_", basis, "_", role), dirs, cfg, width = 7, height = 6.5)
    }
  }
  list(available = TRUE, burden = burden, design = design, associations = tests)
}

.qc_structural_palette <- function(cfg, categories) {
  fallback <- c("full-splice_match" = "#66C2A5", "incomplete-splice_match" = "#FC8D62",
    "novel_in_catalog" = "#8DA0CB", "novel_not_in_catalog" = "#E78AC3",
    genic = "#A6D854", antisense = "#FFD92F", fusion = "#E5C494", intergenic = "#B3B3B3")
  palette <- unlist(cfg$reporting$palettes$structural_categories %||% fallback)
  palette[setdiff(categories, names(palette))] <- "#B3B3B3"
  palette
}

.qc_plot_outputs <- function(sample_qc, dictionary, global, paired, associations, dirs, cfg) {
  if (isFALSE(cfg$reporting$enabled)) return(invisible(NULL))
  labels <- setNames(dictionary$display_label, dictionary$metric)
  display_conditions <- c(primary = cfg$comparison$reference, recurrent = cfg$comparison$comparison)
  temporal <- unlist(cfg$reporting$palettes$timepoints %||% cfg$reporting$palettes$timepoint %||%
    c(primary = "#4C78A8", recurrent = "#E45756"))
  timepoint_palette <- stats::setNames(as.character(temporal[c("primary", "recurrent")]), unname(display_conditions))
  structural_palette <- .qc_structural_palette(cfg, unique(c(global$composition$structural_category,
    global$legacy_structure$structural_category)))
  long <- data.table::melt(data.table::as.data.table(sample_qc), id.vars = c("sample", "patient", "timepoint"),
    measure.vars = dictionary$metric, variable.name = "metric", value.name = "value", variable.factor = FALSE)
  long[, label := labels[metric]]
  long[, condition_display := factor(display_conditions[timepoint], levels = unname(display_conditions))]
  library <- data.table::copy(sample_qc)
  library[, sample_display := factor(sample, levels = sample_qc$sample)]
  library[, condition_display := factor(display_conditions[timepoint], levels = unname(display_conditions))]
  p <- ggplot2::ggplot(library, ggplot2::aes(x = sample_display, y = total_full_length_reads / 1e6, fill = condition_display)) +
    ggplot2::geom_col(width = .8) + ggplot2::scale_fill_manual(values = timepoint_palette) + ggplot2::labs(x = NULL, y = "Full-length reads (millions)", fill = NULL, title = "Full-length read support by sample") +
    rna_theme(cfg) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  rna_save_plot(p, "full_length_reads_by_sample", dirs, cfg, width = 9, height = 5)
  library_mean <- library[, .(mean_FL = mean(total_full_length_reads) / 1e6, sd_FL = stats::sd(total_full_length_reads) / 1e6), by = condition_display]
  p <- ggplot2::ggplot(library_mean, ggplot2::aes(x = condition_display, y = mean_FL, fill = condition_display)) +
    ggplot2::geom_col(width = .65) + ggplot2::scale_fill_manual(values = timepoint_palette) + ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, mean_FL - sd_FL), ymax = mean_FL + sd_FL), width = .15) +
    ggplot2::labs(x = NULL, y = "Full-length reads (millions)", fill = NULL, title = "Full-length read support", subtitle = "Mean and standard deviation") + rna_theme(cfg)
  rna_save_plot(p, "full_length_reads_condition_mean_sd", dirs, cfg, width = 5, height = 4.5)
  coverage_names <- grep("_coverage$", names(library), value = TRUE)
  if (length(coverage_names)) {
    coverage <- data.table::melt(library, id.vars = "sample", measure.vars = coverage_names,
      variable.name = "annotation", value.name = "fraction", variable.factor = FALSE)
    coverage[, annotation := gsub("_", " ", gsub("_FL_coverage$|_annotation_fl_coverage$", "", annotation))]
    p <- ggplot2::ggplot(coverage, ggplot2::aes(x = sample, y = annotation, fill = fraction)) + ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), labels = scales::label_percent(), na.value = "#BDBDBD") +
      ggplot2::labs(x = NULL, y = NULL, fill = "Annotated reads", title = "Availability of transcript annotations") + rna_theme(cfg) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    rna_save_plot(p, "technical_annotation_coverage", dirs, cfg, width = 11, height = 7)
  }
  for (metric_name in dictionary$metric) {
    z <- long[metric == metric_name & is.finite(value)]
    if (!nrow(z)) next
    scale <- dictionary$plot_multiplier[match(metric_name, dictionary$metric)]
    z[, plotted_value := value * scale]
    axis_label <- if (scale == 100) paste0(labels[[metric_name]], " (%)") else if (scale == .001)
      "Median transcript length (kb)" else if (scale == 1e-6) "Full-length reads (millions)" else labels[[metric_name]]
    p <- ggplot2::ggplot(z, ggplot2::aes(x = condition_display, y = plotted_value, group = patient)) +
      ggplot2::geom_line(color = "#9A9A9A", linewidth = .5) +
      ggplot2::geom_point(ggplot2::aes(color = condition_display), size = 2.5) +
      ggplot2::scale_color_manual(values = timepoint_palette) +
      ggplot2::labs(x = NULL, y = paste(strwrap(axis_label, 32), collapse = "\n"),
        title = paste(strwrap(labels[[metric_name]], 40), collapse = "\n"), color = NULL) + rna_theme(cfg)
    rna_save_plot(p, paste0("paired_", metric_name), dirs, cfg, width = 6, height = 4.5)
  }
  for (universe_name in c("raw_local", "all_harmonized", "support_filtered")) {
  comp <- data.table::copy(global$composition[universe == universe_name])
  for (weight in c("isoform_fraction", "FL_fraction")) if (nrow(comp) && any(is.finite(comp[[weight]]))) {
    comp[, sample_display := factor(sample, levels = sample_qc$sample)]
    p <- ggplot2::ggplot(comp, ggplot2::aes(x = sample_display, y = .data[[weight]], fill = structural_category)) +
      ggplot2::geom_col(width = .8) + ggplot2::scale_fill_manual(values = structural_palette, na.value = "#B3B3B3") +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = NULL, y = if (weight == "FL_fraction") "Full-length reads" else "Detected splice chains",
        title = "Transcript structural categories", fill = "Structural category") +
      rna_theme(cfg) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    rna_save_plot(p, paste0("RNA_structural_composition_", universe_name, "_", weight), dirs, cfg, width = 10, height = 6)
  }
  group_comp <- data.table::copy(global$condition_composition[universe == universe_name])
  group_comp[, condition_display := factor(display_conditions[timepoint], levels = unname(display_conditions))]
  for (weight in c("isoform", "FL")) for (basis in c("pooled", "mean")) {
    column <- if (basis == "pooled") paste0("pooled_", weight, "_fraction") else paste0("mean_", weight, "_fraction")
    if (!nrow(group_comp) || !any(is.finite(group_comp[[column]]))) next
    p <- ggplot2::ggplot(group_comp, ggplot2::aes(x = condition_display, y = .data[[column]], fill = structural_category)) +
      ggplot2::geom_col(width = .65) + ggplot2::scale_fill_manual(values = structural_palette, na.value = "#B3B3B3") +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = NULL, y = if (weight == "FL") "Full-length reads" else "Detected splice chains", fill = "Structural category",
        title = "Transcript structural categories", subtitle = if (basis == "pooled") "Pooled across samples" else "Mean sample fractions") + rna_theme(cfg)
    rna_save_plot(p, paste0("RNA_structural_composition_condition_", universe_name, "_", weight, "_", basis), dirs, cfg, width = 7, height = 6)
  }
  }
  cc_group <- data.table::copy(global$structural_by_complexity_condition[universe == "support_filtered"])
  for (tp in c("primary", "recurrent")) for (weight in c("mean_isoform_fraction", "mean_FL_fraction")) {
    z <- data.table::copy(cc_group[timepoint == tp & is.finite(get(weight))])
    if (!nrow(z)) next
    z[, complexity := factor(complexity, levels = c("1", "2", "3-4", "5-8", "9-20", ">20"))]
    p <- ggplot2::ggplot(z, ggplot2::aes(x = complexity, y = .data[[weight]], fill = structural_category)) +
      ggplot2::geom_col(width = .8) + ggplot2::scale_fill_manual(values = structural_palette, na.value = "#B3B3B3") +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = "Detected splice chains per gene", y = if (weight == "mean_FL_fraction") "Mean read fraction" else "Mean splice-chain fraction",
        title = "Structural categories across isoform complexity", subtitle = display_conditions[[tp]], fill = "Structural category") + rna_theme(cfg)
    rna_save_plot(p, paste0("RNA_structure_by_complexity_", tp, "_", weight), dirs, cfg, width = 8, height = 6)
  }
  dist <- global$complexity_distribution[universe == "support_filtered"]
  if (nrow(dist) && any(is.finite(dist$fraction))) {
    dist[, complexity := factor(complexity, levels = c("1", "2", "3-4", "5-8", "9-20", ">20"))]
    dist[, condition_display := factor(display_conditions[timepoint], levels = unname(display_conditions))]
    p <- ggplot2::ggplot(dist, ggplot2::aes(x = complexity, y = fraction, color = condition_display, group = sample)) +
      ggplot2::geom_line(alpha = .6) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = timepoint_palette) +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = "Detected splice chains per gene", y = "Genes within sample", color = NULL,
        title = "Isoform complexity across samples") + rna_theme(cfg)
    rna_save_plot(p, "RNA_detected_isoforms_per_gene", dirs, cfg, width = 7, height = 5)
  }
  gene_sample <- global$gene_sample_complexity[universe == "support_filtered"]
  for (tp in c("primary", "recurrent")) {
    z <- gene_sample[timepoint == tp & gene_FL > 0 & n_isoforms > 0]
    if (!nrow(z)) next
    p <- ggplot2::ggplot(z, ggplot2::aes(x = gene_FL, y = n_isoforms)) +
      ggplot2::geom_point(alpha = .12, size = .7, color = unname(temporal[["primary"]])) +
      ggplot2::scale_x_log10(labels = scales::label_comma()) + ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Gene full-length reads", y = "Detected splice chains per gene",
        title = "Read support and isoform detection", subtitle = display_conditions[[tp]]) + rna_theme(cfg)
    rna_save_plot(p, paste0("RNA_read_support_vs_isoform_number_", tp), dirs, cfg, width = 7, height = 5)
  }
  for (column in c("n_isoforms_all_detected", "n_isoforms_supported")) {
    z <- data.table::copy(global$legacy_complexity[get(column) > 0])
    if (!nrow(z)) next
    upper <- as.numeric(stats::quantile(z[[column]], .99, names = FALSE, na.rm = TRUE))
    p <- ggplot2::ggplot(z, ggplot2::aes(x = .data[[column]])) +
      ggplot2::geom_histogram(bins = 60, fill = unname(temporal[["primary"]]), color = "white", linewidth = .2) +
      ggplot2::coord_cartesian(xlim = c(0, max(upper, 1))) +
      ggplot2::labs(x = if (column == "n_isoforms_supported") "Support-filtered splice chains per gene" else "Detected splice chains per gene",
        y = "Genes", title = "Gene isoform complexity", subtitle = "Range through the 99th percentile") + rna_theme(cfg)
    rna_save_plot(p, paste0("gene_complexity_", column), dirs, cfg, width = 7, height = 5)
    p <- ggplot2::ggplot(z[gene_FL > 0], ggplot2::aes(x = gene_FL, y = .data[[column]])) +
      ggplot2::geom_point(alpha = .25, size = 1.1, color = unname(temporal[["primary"]])) +
      ggplot2::scale_x_log10(labels = scales::label_comma()) + ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Total gene full-length reads", y = "Splice chains per gene", title = "Gene read support and isoform complexity") + rna_theme(cfg)
    rna_save_plot(p, paste0("gene_read_support_vs_complexity_", column), dirs, cfg, width = 7, height = 5)
  }
  for (universe_name in c("all_detected", "supported")) {
    z <- data.table::copy(global$legacy_structure[universe == universe_name])
    if (!nrow(z)) next
    p <- ggplot2::ggplot(z, ggplot2::aes(x = complexity, y = isoform_fraction, fill = structural_category)) +
      ggplot2::geom_col(width = .8) + ggplot2::scale_fill_manual(values = structural_palette, na.value = "#B3B3B3") +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = "Splice chains per gene", y = "Splice-chain fraction", fill = "Structural category",
        title = "Structural categories across gene complexity", subtitle = if (universe_name == "supported") "Supported splice chains" else "All detected splice chains") + rna_theme(cfg)
    rna_save_plot(p, paste0("gene_complexity_structural_categories_", universe_name), dirs, cfg, width = 9, height = 6)
  }
  if (nrow(associations)) {
    z <- associations[is.finite(estimate)]
    if (nrow(z)) {
      z[, predictor := gsub("delta__", "", predictor, fixed = TRUE)]
      z[, outcome := gsub("delta__", "", outcome, fixed = TRUE)]
      z[, predictor_label := ifelse(predictor %in% names(labels), labels[predictor], gsub("_", " ", predictor))]
      z[, outcome_label := ifelse(outcome %in% names(labels), labels[outcome], gsub("_", " ", outcome))]
      p <- ggplot2::ggplot(z, ggplot2::aes(x = outcome_label, y = predictor_label, fill = estimate)) +
        ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", estimate)), size = 3.5) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Technical properties and paired RNA changes") +
        rna_theme(cfg) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
      rna_save_plot(p, "paired_technical_association_heatmap", dirs, cfg, width = 10, height = 7)
    }
  }
  invisible(NULL)
}

#' Run ORF-independent RNA quality and technical checks.
#'
#' The RNA-only association families intentionally exclude protein outcomes.
#' Consequently their FDR family is explicitly recorded and is not claimed
#' to reproduce an earlier combined RNA/protein association screen.
#' @export
run_qc <- function(inputs, cfg, cohorts, stage_dir) {
  rna_require(c("data.table", "ggplot2", "scales"))
  dirs <- rna_dirs(stage_dir)
  prepared <- .qc_prepare_local(inputs, cfg)
  local <- prepared$data
  sample_qc <- .qc_sample_metrics(local, inputs, cfg)
  dictionary <- .qc_dictionary(cfg$technical_qc$terminal_proximity_nt %||% 50,
    cfg$technical_qc$long_transcript_threshold_nt %||% 4000)
  terminus <- .qc_terminus(local, inputs$sample_manifest, cfg$gain_loss$gene_evaluable_min_fl %||% 20)
  rna_write_tsv(prepared$availability, file.path(dirs$qc, "technical_annotation_schema.tsv"))
  rna_write_tsv(sample_qc, file.path(dirs$tables, "sample_technical_qc_metrics.tsv"))
  rna_write_tsv(dictionary, file.path(dirs$tables, "technical_metric_dictionary.tsv"))
  rna_write_tsv(terminus$events, file.path(dirs$tables, "dominant_transcript_terminus_events.tsv.gz"))
  rna_write_tsv(terminus$summary, file.path(dirs$tables, "patient_terminus_change_summary.tsv"))
  subcategory <- local[, .(full_length_reads = sum(FL), n_local_transcripts = .N),
    by = .(structural_category_normalized, subcategory_normalized)]
  rna_write_tsv(subcategory, file.path(dirs$tables, "structural_subcategory_audit.tsv"))
  count_thresholds <- as.numeric(cfg$global_qc$detected_min_fl %||% c(1, 2, 3, 5))
  detection <- data.table::rbindlist(lapply(count_thresholds, function(threshold) {
    x <- local[, .(detected_local_transcripts = sum(FL >= threshold)), by = sample]
    data.frame(sample = colnames(inputs$counts_matrix), minimum_FL = threshold,
      detected_isoforms = as.integer(colSums(inputs$counts_matrix >= threshold)),
      detected_local_transcripts = as.integer(x$detected_local_transcripts[match(colnames(inputs$counts_matrix), x$sample)]))
  }))
  detection <- merge(detection, data.table::as.data.table(inputs$sample_manifest)[, .(sample, patient, timepoint)], by = "sample")
  rna_write_tsv(detection, file.path(dirs$tables, "isoform_detection_by_FL_threshold.tsv"))
  minimum_n <- cfg$technical_qc$minimum_correlation_n %||% 4L
  predictors <- c("log10_total_full_length_reads", "incomplete_splice_match_FL_fraction",
    "ism_5prime_fragment_FL_fraction", "within_CAGE_peak_FL_fraction",
    "tss_within_threshold_FL_fraction", "mitochondrial_chromosome_FL_fraction")
  outcomes <- c("median_transcript_length_FL_weighted", "long_transcript_FL_fraction")
  results <- vector("list", length(cohorts)); names(results) <- names(cohorts)
  for (cohort in names(cohorts)) {
    rna_log("RNA QC: ", cohort)
    cd <- rna_dirs(file.path(stage_dir, cohort))
    manifest <- as.data.frame(inputs$sample_manifest)
    samples <- manifest[manifest$patient %in% cohorts[[cohort]], , drop = FALSE]
    sq <- sample_qc[sample %in% samples$sample]
    long <- data.table::melt(sq, id.vars = c("sample", "patient", "timepoint"),
      measure.vars = dictionary$metric, variable.name = "metric", value.name = "value", variable.factor = FALSE)
    paired <- data.table::dcast(long, patient + metric ~ timepoint, value.var = "value")
    for (name in setdiff(c("primary", "recurrent"), names(paired))) paired[, (name) := NA_real_]
    paired[, delta_recurrent_minus_primary := recurrent - primary]
    tests <- data.table::rbindlist(lapply(dictionary$metric, function(metric_name) {
      z <- paired[metric == metric_name]
      data.frame(cohort = cohort, metric = metric_name, .qc_paired_test(z$primary, z$recurrent, minimum_n))
    }))
    tests[, FDR := stats::p.adjust(p_value, method = "BH")]
    tests <- merge(tests, dictionary, by = "metric", all.x = TRUE)
    tests[, adjustment_family := "RNA_only_sample_metric_tests_within_cohort"]
    pair_wide <- data.table::dcast(paired, patient ~ metric, value.var = "delta_recurrent_minus_primary")
    rename <- setdiff(names(pair_wide), "patient")
    data.table::setnames(pair_wide, rename, paste0("delta__", rename))
    pair_analysis <- merge(data.table::data.table(patient = as.character(cohorts[[cohort]])), pair_wide, by = "patient", all.x = TRUE)
    pair_analysis <- merge(pair_analysis, data.table::as.data.table(terminus$summary), by = "patient", all.x = TRUE)
    paired_outcomes <- c(paste0("delta__", outcomes), "tss_bin_changed_fraction", "tts_bin_changed_fraction",
      "any_terminus_changed_fraction", "median_signed_dominant_transcript_length_change_nt")
    paired_outcomes <- intersect(paired_outcomes, names(pair_analysis))
    pa <- .qc_association_grid(pair_analysis, paste0("delta__", predictors), paired_outcomes,
      cohort, "within_patient_change", minimum_n)
    sa <- .qc_association_grid(sq, predictors, outcomes, cohort,
      "sample_level_descriptive_paired_libraries_not_independent", minimum_n)
    global <- .qc_global(inputs, local, cfg, samples, cd)
    gain_loss_qc <- .qc_gain_loss_associations(inputs$gain_loss, sq, cohort, cohorts[[cohort]], cfg, cd)
    rna_write_tsv(sq, file.path(cd$tables, "sample_technical_qc_metrics.tsv"))
    rna_write_tsv(long, file.path(cd$tables, "sample_technical_metrics_long.tsv"))
    group_summary <- long[, .(n_finite = sum(is.finite(value)),
      mean = if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_,
      sd = if (sum(is.finite(value)) >= 2) stats::sd(value[is.finite(value)]) else NA_real_), by = .(timepoint, metric)]
    rna_write_tsv(group_summary, file.path(cd$tables, "condition_technical_metric_mean_sd.tsv"))
    rna_write_tsv(paired, file.path(cd$tables, "paired_technical_metric_changes.tsv"))
    rna_write_tsv(tests, file.path(cd$tables, "paired_technical_metric_tests.tsv"))
    rna_write_tsv(pair_analysis, file.path(cd$tables, "within_patient_analysis_matrix.tsv"))
    rna_write_tsv(pa, file.path(cd$tables, "within_patient_technical_associations.tsv"))
    rna_write_tsv(sa, file.path(cd$tables, "sample_level_technical_associations_descriptive.tsv"))
    availability <- data.frame(cohort = cohort, metric = dictionary$metric,
      n_finite_samples = vapply(dictionary$metric, function(name) sum(is.finite(sq[[name]])), integer(1)),
      n_requested_samples = nrow(sq))
    availability$status <- ifelse(availability$n_finite_samples == 0, "unavailable",
      ifelse(availability$n_finite_samples == availability$n_requested_samples, "available", "partly_available"))
    rna_write_tsv(availability, file.path(cd$qc, "metric_availability.tsv"))
    .qc_plot_outputs(sq, dictionary, global, paired, pa, cd, cfg)
    term_summary <- data.table::as.data.table(terminus$summary)
    term_summary <- term_summary[patient %in% cohorts[[cohort]]]
    results[[cohort]] <- list(sample_metrics = sq, paired_metrics = paired, paired_tests = tests,
      paired_associations = pa, sample_associations = sa, terminus_summary = term_summary,
      global = global, gain_loss = gain_loss_qc, availability = availability)
  }
  result <- list(complete = TRUE, version = "0.1.0", sample_metrics = sample_qc,
    metric_dictionary = dictionary, terminus_events = terminus$events, terminus_summary = terminus$summary,
    cohorts = results, fingerprint = rna_hash(list(inputs$fingerprint, inputs$gain_loss$fingerprint,
      cfg$global_qc, cfg$technical_qc, cfg$gain_loss$gene_evaluable_min_fl,
      cfg$remodeling$definition$min_isoform_total_fl, cohorts)))
  rna_atomic_rds(result, file.path(dirs$objects, "rna_qc_results.rds"))
  rna_write_tsv(data.frame(check = c("local_count_conservation", "RNA_only_optional_annotations"),
    pass = c(all(sample_qc$count_total_difference == 0), TRUE),
    detail = c("All harmonized counts reconcile with local transcript FL totals.",
      "Unavailable optional measurements are exported explicitly; no protein fields required.")),
    file.path(dirs$qc, "completion_checks.tsv"))
  result
}
