# Raw Pigeon import and splice-chain harmonization.
# Extracted from the study's Pigeon importer and revised M01/M02 rules.
# A harmonized ID denotes a splice chain, not exact transcript-end identity.

rna_import_clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  x
}

rna_import_read_table <- function(path, required, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    rna_stop(label, " does not exist: ", path)
  }
  out <- data.table::fread(path, data.table = TRUE, colClasses = "character",
                          check.names = FALSE, na.strings = c("NA", "NaN", ""),
                          showProgress = FALSE)
  absent <- setdiff(required, names(out))
  if (length(absent)) rna_stop(label, " lacks columns: ", paste(absent, collapse = ", "))
  if (anyDuplicated(names(out))) rna_stop(label, " has duplicate column names.")
  out
}

rna_import_integer <- function(x, label, positive = FALSE) {
  y <- suppressWarnings(as.numeric(x))
  lower <- if (positive) 1 else 0
  bad <- !is.finite(y) | y < lower |
    abs(y - round(y)) > 1e-8 | y > .Machine$integer.max
  if (any(bad)) rna_stop(label, " must contain finite ",
                         if (positive) "positive" else "nonnegative",
                         " integer values; invalid rows: ",
                         paste(utils::head(which(bad), 10L), collapse = ", "))
  as.integer(round(y))
}

rna_import_canonical_annotation <- function(local, isoform_ids) {
  # IMPORTANT: annotation ranking never removes count rows. Only the ranking
  # of structural categories is conditioned on the chosen original gene label.
  dt <- data.table::as.data.table(local)
  support <- dt[!is.na(associated_gene) & nzchar(associated_gene) &
                  associated_gene != "unknown",
                .(total_FL_support = sum(FL),
                  n_samples_detected = data.table::uniqueN(sample[FL >= 1L]),
                  n_patients_detected = data.table::uniqueN(patient[FL >= 1L]),
                  n_source_rows = .N),
                by = .(isoform_harmonized_id, associated_gene)]
  support <- support[total_FL_support >= 1]
  if (nrow(support)) {
    data.table::setorderv(support,
                         c("isoform_harmonized_id", "total_FL_support",
                           "n_samples_detected", "associated_gene"),
                         c(1L, -1L, -1L, 1L))
    support[, `:=`(
      gene_support_rank = seq_len(.N),
      n_gene_candidates = .N,
      runner_up_total_FL = if (.N > 1L) total_FL_support[2L] else NA_real_,
      winning_support_fraction = total_FL_support / sum(total_FL_support),
      n_top_tied_gene_labels = sum(total_FL_support == total_FL_support[1L] &
                                    n_samples_detected == n_samples_detected[1L])
    ), by = isoform_harmonized_id]
    chosen <- support[gene_support_rank == 1L]
    chosen[, mapping_resolution := ifelse(
      n_gene_candidates == 1L, "unique_gene_annotation",
      ifelse(n_top_tied_gene_labels == 1L, "resolved_by_FL_and_sample_support",
             "resolved_by_original_lexical_tiebreak"))]
    category_rows <- merge(
      dt[, .(isoform_harmonized_id, associated_gene, structural_category, FL)],
      chosen[, .(isoform_harmonized_id, associated_gene)],
      by = c("isoform_harmonized_id", "associated_gene"),
      all = FALSE, allow.cartesian = FALSE, sort = FALSE)
    category_rows[is.na(structural_category) | !nzchar(structural_category),
                  structural_category := "unknown"]
    category <- category_rows[FL >= 1L,
                              .(category_FL_support = sum(FL),
                                category_source_rows = .N),
                              by = .(isoform_harmonized_id, associated_gene,
                                     structural_category)]
    data.table::setorderv(category,
                         c("isoform_harmonized_id", "associated_gene",
                           "category_FL_support", "category_source_rows",
                           "structural_category"), c(1L, 1L, -1L, -1L, 1L))
    category <- category[, .SD[1L], by = c("isoform_harmonized_id", "associated_gene")]
    chosen <- merge(chosen, category,
                    by = c("isoform_harmonized_id", "associated_gene"), all.x = TRUE,
                    sort = FALSE)
  } else {
    chosen <- data.table::data.table(
      isoform_harmonized_id = character(), associated_gene = character(),
      total_FL_support = numeric(), n_samples_detected = integer(),
      n_patients_detected = integer(), n_source_rows = integer(),
      gene_support_rank = integer(), n_gene_candidates = integer(),
      runner_up_total_FL = numeric(), winning_support_fraction = numeric(),
      n_top_tied_gene_labels = integer(), mapping_resolution = character(),
      structural_category = character(), category_FL_support = numeric(),
      category_source_rows = integer())
  }
  annotation <- merge(data.table::data.table(isoform_harmonized_id = isoform_ids),
                      chosen, by = "isoform_harmonized_id", all.x = TRUE, sort = FALSE)
  annotation <- annotation[match(isoform_ids, isoform_harmonized_id)]
  annotation[is.na(mapping_resolution), mapping_resolution := "unmapped_gene_annotation"]
  list(annotation = annotation, gene_candidates = support)
}

import_pigeon <- function(cfg, samples, stage_dir) {
  rna_require(c("data.table", "digest"))
  dirs <- rna_dirs(stage_dir)
  # setDT/setorder mutate by reference; never alter the caller's metadata.
  samples <- data.table::copy(as.data.frame(samples, stringsAsFactors = FALSE))
  required_metadata <- c("sample_id", "subject_id", "condition", "classification_file", "junctions_file")
  missing_metadata <- setdiff(required_metadata, names(samples))
  if (length(missing_metadata)) rna_stop("Sample metadata lacks: ", paste(missing_metadata, collapse = ", "))
  if (!nrow(samples)) rna_stop("No included samples were supplied.")
  if (anyNA(samples$sample_id) || any(!nzchar(samples$sample_id)) || anyDuplicated(samples$sample_id)) {
    rna_stop("Sample IDs must be nonempty and unique.")
  }
  if (!"sample" %in% names(samples)) samples$sample <- as.character(samples$sample_id)
  if (!"patient" %in% names(samples)) samples$patient <- as.character(samples$subject_id)
  if (!"timepoint" %in% names(samples)) {
    samples$timepoint <- ifelse(samples$condition == cfg$comparison$reference, "primary",
                                ifelse(samples$condition == cfg$comparison$comparison, "recurrent", NA_character_))
  }
  if (anyNA(samples$timepoint) || any(!samples$timepoint %in% c("primary", "recurrent"))) {
    rna_stop("Sample conditions must match the two configured comparison labels.")
  }
  if (!"barcode" %in% names(samples)) samples$barcode <- NA_character_
  if (!"fl_count_column" %in% names(samples)) {
    samples$fl_count_column <- cfg$import$count_column_default %||% "FL.sample"
  }
  data.table::setorderv(data.table::setDT(samples), c("patient", "timepoint", "sample"))
  # Internal temporal roles deliberately order primary before recurrent.
  samples <- samples[order(patient, match(timepoint, c("primary", "recurrent")), sample)]
  class_required <- c("isoform", "chrom", "strand", "length", "exons",
                       "structural_category", "associated_gene")
  junction_required <- c("isoform", "chrom", "strand", "genomic_start_coord",
                          "genomic_end_coord", "splice_site")
  numeric_annotations <- c("diff_to_TSS", "diff_to_TTS", "diff_to_gene_TSS", "diff_to_gene_TTS",
                            "dist_to_CAGE_peak", "dist_to_polyA_site", "polyA_dist",
                            "perc_A_downstream_TTS", "ratio_TSS")
  character_annotations <- c("associated_transcript", "subcategory", "all_canonical", "RTS_stage",
                              "within_CAGE_peak", "within_polyA_site", "polyA_motif", "polyA_motif_found")
  retained_class <- c("isoform", "chrom", "strand", "length", "exons", "structural_category",
                       "associated_gene", character_annotations, numeric_annotations)
  junction_character <- c("junction_number", "junction_category", "canonical", "RTS_junction",
                           "start_site_category", "end_site_category", "bite_junction", "indel_near_junct")
  local_list <- vector("list", nrow(samples))
  junction_list <- vector("list", nrow(samples))
  unmatched_list <- vector("list", nrow(samples))
  source_rows <- list()
  for (i in seq_len(nrow(samples))) {
    s <- samples[i]
    sample_id <- as.character(s$sample)
    rna_log("Importing Pigeon: ", sample_id)
    count_column <- as.character(s$fl_count_column)
    if (is.na(count_column) || !nzchar(count_column)) {
      count_column <- cfg$import$count_column_default %||% "FL.sample"
    }
    samples[i, fl_count_column := count_column]
    raw <- rna_import_read_table(s$classification_file, c(class_required, count_column),
                                 paste0("Classification for ", sample_id))
    if (!nrow(raw)) rna_stop("Classification table is empty for ", sample_id)
    raw[, isoform := rna_import_clean_character(isoform)]
    if (anyNA(raw$isoform) || anyDuplicated(raw$isoform)) {
      duplicate <- raw[is.na(isoform) | duplicated(isoform) | duplicated(isoform, fromLast = TRUE)]
      rna_write_tsv(as.data.frame(duplicate), file.path(dirs$qc, paste0(sample_id, "_invalid_local_keys.tsv")))
      rna_stop("Missing or duplicated local transcript IDs in classification for ", sample_id)
    }
    FL <- rna_import_integer(raw[[count_column]], paste0(sample_id, " ", count_column))
    if (sum(FL) <= 0) rna_stop("Sample has no positive FL counts: ", sample_id)
    for (nm in setdiff(retained_class, names(raw))) raw[, (nm) := NA_character_]
    local <- data.table::copy(raw[, ..retained_class])
    data.table::setnames(local, "isoform", "isoform_local")
    local[, `:=`(sample = sample_id, patient = as.character(s$patient),
                  timepoint = as.character(s$timepoint), group = as.character(s$timepoint), FL = FL)]
    local[, length := rna_import_integer(length, paste0(sample_id, " transcript length"), positive = TRUE)]
    local[, exons := rna_import_integer(exons, paste0(sample_id, " exon count"), positive = TRUE)]
    for (nm in c("chrom", "strand", "structural_category", "associated_gene", character_annotations)) {
      local[, (nm) := rna_import_clean_character(get(nm))]
    }
    for (nm in numeric_annotations) local[, (nm) := suppressWarnings(as.numeric(get(nm)))]
    tpm_column <- sub("^FL\\.", "FL_TPM.", count_column)
    local[, FL_TPM := if (tpm_column %in% names(raw)) suppressWarnings(as.numeric(raw[[tpm_column]])) else NA_real_]
    # Coding/ORF/CDS/NMD fields are intentionally not copied from raw Pigeon files.
    junc <- rna_import_read_table(s$junctions_file, junction_required,
                                  paste0("Junctions for ", sample_id))
    for (nm in setdiff(junction_character, names(junc))) junc[, (nm) := NA_character_]
    keep_junction <- c(junction_required, junction_character)
    junc <- data.table::copy(junc[, ..keep_junction])
    junc[, isoform := rna_import_clean_character(isoform)]
    junc[, sample := sample_id]
    unmatched_list[[i]] <- junc[is.na(isoform) | !isoform %in% local$isoform_local]
    # An unfiltered junction file may contain local IDs not retained in the
    # classification file. Audit them and retain only classification members.
    junc <- junc[!is.na(isoform) & isoform %in% local$isoform_local]
    if (nrow(junc)) {
      for (nm in c("chrom", "strand", "splice_site")) junc[, (nm) := rna_import_clean_character(get(nm))]
      if (anyNA(junc$chrom) || anyNA(junc$strand) || anyNA(junc$splice_site) ||
          any(!junc$strand %in% c("+", "-"))) {
        rna_stop("Junction chromosome, strand or motif is unavailable/invalid for ", sample_id)
      }
      for (nm in c("genomic_start_coord", "genomic_end_coord")) {
        junc[, (nm) := rna_import_integer(get(nm), paste0(sample_id, " ", nm), positive = TRUE)]
      }
      if (any(junc$genomic_end_coord < junc$genomic_start_coord)) rna_stop("Reversed junction coordinates: ", sample_id)
      duplicates <- duplicated(junc[, .(isoform, genomic_start_coord, genomic_end_coord)])
      if (any(duplicates)) {
        rna_write_tsv(as.data.frame(junc[duplicates]), file.path(dirs$qc, paste0(sample_id, "_duplicate_junctions.tsv")))
        rna_stop("Duplicated junction coordinates within a local transcript: ", sample_id)
      }
      consistency <- junc[, .(n_chrom = data.table::uniqueN(chrom), n_strand = data.table::uniqueN(strand)), by = isoform]
      if (any(consistency$n_chrom != 1L | consistency$n_strand != 1L)) {
        rna_stop("Inconsistent chromosome/strand within a local transcript: ", sample_id)
      }
      sig <- junc[, {
        ord <- order(genomic_start_coord, genomic_end_coord)
        .(chrom_signature = chrom[1L], strand_signature = strand[1L],
          n_junctions = .N,
          junction_signature = paste0(genomic_start_coord[ord], ":", genomic_end_coord[ord], ":",
                                        splice_site[ord], collapse = "|"))
      }, by = isoform]
      sig[, full_signature := paste(chrom_signature, strand_signature, junction_signature, sep = "|")]
      # Keep digest's default serialize=TRUE exactly as in the original code.
      sig[, isoform_harmonized_id := paste0("JUNC_", substr(vapply(full_signature, digest::digest,
                                                                       character(1), algo = "xxhash64"), 1L, 16L))]
      k <- match(local$isoform_local, sig$isoform)
      has <- !is.na(k)
      if (any(has & (is.na(local$chrom) | is.na(local$strand)))) rna_stop("Missing classification chromosome/strand for spliced transcript: ", sample_id)
      if (any(local$chrom[has] != sig$chrom_signature[k[has]] |
              local$strand[has] != sig$strand_signature[k[has]])) {
        rna_stop("Classification/junction chromosome or strand disagreement: ", sample_id)
      }
      local[, `:=`(isoform_harmonized_id = sig$isoform_harmonized_id[k],
                    junction_signature = sig$junction_signature[k],
                    full_signature = sig$full_signature[k], n_junctions = sig$n_junctions[k])]
    } else {
      local[, `:=`(isoform_harmonized_id = NA_character_, junction_signature = NA_character_,
                    full_signature = NA_character_, n_junctions = NA_integer_)]
    }
    mismatch <- (local$exons > 1L & is.na(local$n_junctions)) |
      (!is.na(local$n_junctions) & local$n_junctions != local$exons - 1L)
    if (any(mismatch)) {
      rna_write_tsv(as.data.frame(local[mismatch]), file.path(dirs$qc, paste0(sample_id, "_exon_junction_mismatch.tsv")))
      rna_stop("Exon/junction inconsistency or missing junctions for spliced transcript: ", sample_id)
    }
    local[is.na(isoform_harmonized_id), isoform_harmonized_id := paste0("PBONLY_", sample, "|", isoform_local)]
    local[is.na(n_junctions), n_junctions := 0L]
    bin_width <- as.numeric(cfg$harmonization$optional_end_bin_width_nt %||% 50)
    if (!is.finite(bin_width) || bin_width <= 0) rna_stop("Invalid optional_end_bin_width_nt")
    local[, `:=`(tss_bin = round(diff_to_gene_TSS / bin_width) * bin_width,
                  tts_bin = round(diff_to_gene_TTS / bin_width) * bin_width)]
    local[, isoform_junc_end_id := ifelse(startsWith(isoform_harmonized_id, "JUNC_"),
                                         paste0(isoform_harmonized_id, "_TSS", tss_bin, "_TTS", tts_bin),
                                         isoform_harmonized_id)]
    local_list[[i]] <- local
    junction_list[[i]] <- junc
    for (role in intersect(c("classification_file", "junctions_file", "pigeon_summary_file", "filter_reasons_file"), names(samples))) {
      path <- as.character(s[[role]])
      exists <- length(path) == 1L && !is.na(path) && nzchar(path) && file.exists(path)
      source_rows[[length(source_rows) + 1L]] <- data.frame(
        sample = sample_id, role = role, path = if (length(path)) path else NA_character_,
        exists = exists, required = role %in% c("classification_file", "junctions_file"),
        sha256 = if (exists) digest::digest(file = path, algo = "sha256") else NA_character_,
        stringsAsFactors = FALSE)
    }
  }
  local <- data.table::rbindlist(local_list, use.names = TRUE, fill = TRUE)
  data.table::setorderv(local, c("sample", "isoform_local"))
  signatures <- unique(local[startsWith(isoform_harmonized_id, "JUNC_"),
                             .(isoform_harmonized_id, full_signature)])
  if (anyDuplicated(signatures$isoform_harmonized_id)) rna_stop("Hash collision: different splice chains share an ID.")
  aggregate_counts <- local[, .(FL = sum(as.numeric(FL))), by = .(isoform_harmonized_id, sample)]
  if (any(aggregate_counts$FL > .Machine$integer.max)) rna_stop("Aggregated FL cell exceeds integer limit.")
  isoform_ids <- sort(unique(local$isoform_harmonized_id))
  matrix <- matrix(0L, nrow = length(isoform_ids), ncol = nrow(samples),
                     dimnames = list(isoform_ids, as.character(samples$sample)))
  matrix[cbind(match(aggregate_counts$isoform_harmonized_id, isoform_ids),
                match(aggregate_counts$sample, colnames(matrix)))] <- as.integer(aggregate_counts$FL)
  raw_totals <- local[, .(raw_FL = sum(as.numeric(FL))), by = sample]
  raw_totals[, harmonized_FL := as.numeric(colSums(matrix)[sample])]
  raw_totals[, difference := harmonized_FL - raw_FL]
  if (any(raw_totals$difference != 0)) rna_stop("FL count conservation failed during harmonization.")
  resolved <- rna_import_canonical_annotation(local, isoform_ids)
  annotation <- resolved$annotation
  keep_categories <- cfg$gain_loss$keep_structural_categories %||%
    c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")
  annotation[, `:=`(
    isoform_id_type = ifelse(startsWith(isoform_harmonized_id, "JUNC_"), "junction_chain", "sample_specific_monoexonic_fallback"),
    cross_sample_comparable = startsWith(isoform_harmonized_id, "JUNC_"),
    gene_mapped = !is.na(associated_gene) & nzchar(associated_gene) & associated_gene != "unknown",
    structural_category_retained = !is.na(structural_category) & structural_category %in% keep_categories)]
  annotation[, legacy_inclusive_gain_loss_eligible := gene_mapped & structural_category_retained]
  comparable_required <- isTRUE(cfg$gain_loss$require_cross_sample_comparable_id %||% TRUE)
  annotation[, analyzed_for_gain_loss := legacy_inclusive_gain_loss_eligible & (!comparable_required | cross_sample_comparable)]
  lost_label_count_audit <- merge(local[, .(isoform_harmonized_id, associated_gene, FL)],
                                  annotation[, .(isoform_harmonized_id, canonical_gene = associated_gene)],
                                  by = "isoform_harmonized_id", all.x = TRUE, sort = FALSE)
  # Competing-label FL is a descriptive audit, and remains in the count matrix.
  competing_FL <- sum(lost_label_count_audit$FL[
    !is.na(lost_label_count_audit$canonical_gene) &
      (is.na(lost_label_count_audit$associated_gene) |
         lost_label_count_audit$associated_gene != lost_label_count_audit$canonical_gene)])
  audit <- data.frame(metric = c("local_transcript_rows", "harmonized_isoforms", "total_FL",
                                  "sample_specific_local_transcripts", "competing_gene_label_FL_retained",
                                  "ORF_input_fields_required"),
                       value = c(nrow(local), nrow(matrix), sum(matrix),
                                   sum(startsWith(local$isoform_harmonized_id, "PBONLY_")), competing_FL, 0))
  file_manifest <- data.table::rbindlist(source_rows, fill = TRUE)
  fingerprint <- rna_hash(list(identity_version = "junction_chain_paper_v1",
                                source_files = file_manifest[, .(sample, role, sha256)],
                                sample_roles = samples[, .(sample, patient, timepoint)],
                                canonical_policy = "paper_support_rank_preserve_all_counts",
                                eligibility = list(keep_categories, comparable_required)))
  inputs <- list(
    local_transcripts = as.data.frame(local), counts_matrix = matrix,
    isoform_annotation = as.data.frame(annotation), sample_manifest = as.data.frame(samples),
    canonical_annotations = as.data.frame(annotation), gene_annotation_candidates = as.data.frame(resolved$gene_candidates),
    isoform_dictionary = as.data.frame(signatures),
    local_to_harmonized_map = as.data.frame(local[, .(sample, isoform_local, isoform_harmonized_id,
                                                    junction_signature, n_junctions)]),
    junctions = as.data.frame(data.table::rbindlist(junction_list, fill = TRUE)),
    source_manifest = as.data.frame(file_manifest), count_conservation = as.data.frame(raw_totals),
    fingerprint = fingerprint, schema_version = "1.0.0")
  rna_write_tsv(inputs$local_transcripts, file.path(dirs$tables, "pigeon_master_long.tsv.gz"))
  rna_write_tsv(as.data.frame(aggregate_counts), file.path(dirs$tables, "harmonized_counts_long.tsv.gz"))
  rna_write_tsv(inputs$isoform_annotation, file.path(dirs$tables, "canonical_isoform_annotation.tsv.gz"))
  rna_write_tsv(inputs$gene_annotation_candidates, file.path(dirs$qc, "isoform_gene_support_candidates.tsv.gz"))
  rna_write_tsv(inputs$local_to_harmonized_map, file.path(dirs$tables, "local_to_harmonized_map.tsv.gz"))
  rna_write_tsv(inputs$isoform_dictionary, file.path(dirs$tables, "splice_chain_dictionary.tsv.gz"))
  rna_write_tsv(inputs$source_manifest, file.path(dirs$qc, "input_file_manifest.tsv"))
  rna_write_tsv(inputs$count_conservation, file.path(dirs$qc, "count_conservation.tsv"))
  rna_write_tsv(as.data.frame(data.table::rbindlist(unmatched_list, fill = TRUE)),
                file.path(dirs$qc, "junction_rows_not_in_classification.tsv.gz"))
  rna_write_tsv(audit, file.path(dirs$qc, "import_audit.tsv"))
  rna_atomic_rds(inputs, file.path(dirs$objects, "imported_inputs.rds"))
  rna_log("Imported ", nrow(local), " local transcript rows; ", nrow(matrix),
            " harmonized IDs; FL count conservation verified.")
  inputs
}
