# Small synthetic raw Pigeon fixture. Genomic locations are artificial and
# must not be used as biological results or a real annotation reference.
make_rna_fixture <- function(path, n_pairs = 6L) {
  n_pairs <- as.integer(n_pairs)
  stopifnot(length(n_pairs) == 1L, is.finite(n_pairs), n_pairs >= 2L)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  metadata <- vector("list", 2L * n_pairs)
  index <- 0L
  for (p in seq_len(n_pairs)) {
    for (condition in c("primary", "recurrent")) {
      index <- index + 1L
      sample <- paste0("ptn", p, "_", condition)
      ids <- paste0("PB.", p, ".", if (condition == "primary") "P" else "R", ".", seq_len(10L))
      counts <- if (condition == "primary") c(90L, 10L, 60L, 40L, 95L, 5L, 20L, 40L, 20L, 1L) else c(30L, 70L, 60L, 40L, 0L, 100L, 20L, 40L, 20L, 1L)
      classification <- data.frame(
        isoform = ids,
        chrom = c(rep("chr1", 2L), rep("chr2", 2L), rep("chr3", 2L), rep("chr4", 3L), "chrM"),
        strand = "+", length = c(1000L, 1050L, 1100L, 900L, 1250L, 980L, 1000L, 1150L, 1200L, 800L),
        exons = c(rep(2L, 9L), 1L),
        structural_category = c("full-splice_match", "novel_in_catalog", "full-splice_match", "incomplete-splice_match",
                                  "full-splice_match", "novel_not_in_catalog", "full-splice_match", "full-splice_match", "novel_in_catalog", "genic"),
        associated_gene = c("TP53", "TP53", "EGFR", "EGFR", "PTEN", "PTEN", "ATRX", if (p == 1L) "OTHER_LABEL" else "ATRX", "ATRX", "MT_SYNTHETIC"),
        associated_transcript = paste0("ENST_SYNTHETIC_", seq_len(10L)),
        subcategory = c(rep("multi-exon", 9L), "mono-exon"),
        all_canonical = c(rep("canonical", 9L), NA_character_),
        RTS_stage = "FALSE", diff_to_TSS = 0, diff_to_TTS = 0,
        diff_to_gene_TSS = 0, diff_to_gene_TTS = 0,
        dist_to_CAGE_peak = 0, within_CAGE_peak = "TRUE",
        dist_to_polyA_site = 0, within_polyA_site = "TRUE",
        polyA_motif = "AATAAA", polyA_dist = 15, polyA_motif_found = "TRUE",
        check.names = FALSE, stringsAsFactors = FALSE)
      classification[["FL.sample"]] <- counts
      classification[["FL_TPM.sample"]] <- counts / sum(counts) * 1e6
      # Local transcripts 7 and 8 share the same full splice chain; their ends
      # and lengths differ. This exercises within-sample count aggregation.
      junctions <- data.frame(
        isoform = ids[1:9], chrom = classification$chrom[1:9], strand = "+",
        genomic_start_coord = c(101L, 101L, 301L, 301L, 501L, 501L, 701L, 701L, 751L),
        genomic_end_coord = c(200L, 250L, 400L, 450L, 600L, 650L, 800L, 800L, 850L),
        splice_site = "GT-AG", junction_number = "junction_1", canonical = "canonical",
        junction_category = "known", RTS_junction = "FALSE", stringsAsFactors = FALSE)
      class_file <- file.path(path, paste0(sample, "_classification.filtered_lite_classification.txt"))
      junction_file <- file.path(path, paste0(sample, "_classification.filtered_lite_junctions.txt"))
      utils::write.table(classification, class_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
      utils::write.table(junctions, junction_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
      metadata[[index]] <- data.frame(
        sample_id = sample, subject_id = paste0("ptn", p), condition = condition,
        classification_file = class_file, junctions_file = junction_file,
        fl_count_column = "FL.sample", pigeon_summary_file = "", filter_reasons_file = "",
        barcode = sprintf("bc%02d", index), include = TRUE,
        sample = sample, patient = paste0("ptn", p), timepoint = condition,
        stringsAsFactors = FALSE)
    }
  }
  samples <- do.call(rbind, metadata)
  metadata_file <- file.path(path, "samples.csv")
  utils::write.csv(samples[, setdiff(names(samples), c("sample", "patient", "timepoint"))], metadata_file, row.names = FALSE, na = "")
  gtf_file <- file.path(path, "synthetic_annotation.gtf")
  gtf <- unlist(lapply(seq_along(c("TP53", "EGFR", "PTEN", "ATRX")), function(i) {
    gene <- c("TP53", "EGFR", "PTEN", "ATRX")[i]
    vapply(1:2, function(t) {
      paste(paste0("chr", i), "synthetic", "transcript", 1L, 1000L, ".", "+", ".",
              paste0('gene_id "', gene, '"; gene_name "', gene, '"; transcript_id "SYNTHETIC_', gene, '_', t, '";'), sep = "\t")
    }, character(1))
  }))
  writeLines(gtf, gtf_file)
  cfg <- list(
    comparison = list(reference = "primary", comparison = "recurrent"),
    import = list(count_column_default = "FL.sample"),
    harmonization = list(identity_version = "junction_chain_paper_v1",
                          canonical_gene_policy = "paper_support_rank_preserve_all_counts", optional_end_bin_width_nt = 50),
    gain_loss = list(keep_structural_categories = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"),
                     require_cross_sample_comparable_id = TRUE))
  list(samples = samples, metadata_file = metadata_file, gtf_file = gtf_file, cfg = cfg,
         raw_total_FL_per_sample = setNames(rep(381, nrow(samples)), samples$sample))
}
