# Paired and independent-group differential transcript usage. The numerical engine
# lives in inst/engines/dtu_core.R and is evaluated in a private environment.

dtu_alpha_tag <- function(alpha) {
  paste0("alpha", sub("^FDR", "", dtu_fdr_tag(alpha)))
}

dtu_fdr_tag <- function(alpha) {
  z <- sub("^0\\.", "", format(alpha, scientific = FALSE, trim = TRUE))
  if (nchar(z) < 2L) z <- paste0(z, "0")
  paste0("FDR", gsub("[^0-9]", "", z))
}

dtu_engine_environment <- function(parent = parent.frame()) {
  e <- new.env(parent = parent)
  # Source functions use the same non-standard evaluation semantics as the
  # original scripts. Import exports into this environment, without attaching
  # packages or writing into .GlobalEnv or a package namespace.
  for (pkg in c("dplyr", "forcats", "ggplot2", "purrr", "readr", "scales", "stringr", "tibble", "tidyr")) {
    ns <- asNamespace(pkg)
    for (name in getNamespaceExports(pkg)) assign(name, getExportedValue(pkg, name), envir = e)
  }
  e
}

dtu_engine_path <- function(cfg, name) {
  candidate <- file.path(cfg$.code_root %||% "", "inst", "engines", name)
  if (!file.exists(candidate)) candidate <- system.file("engines", name, package = "IsoRemodel")
  if (!nzchar(candidate) || !file.exists(candidate)) stop("Missing packaged DTU engine: ", name, call. = FALSE)
  candidate
}

dtu_cohort_seed <- function(master, cohort) {
  # Retain the two publication namespaces, and use stable IDs for new cohorts.
  if (identical(cohort, "all_6_pairs")) return(as.integer(master + 1000L))
  if (identical(cohort, "excluding_ptn2_5_pairs")) return(as.integer(master + 2000L))
  rna_stage_seed(master, "dtu", cohort)
}

# Random-number settings are part of the DTU model contract. A numerical seed
# alone does not select the same 10% precision subset across RNG algorithms.
dtu_rng_snapshot <- function() {
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(kind = RNGkind(), seed_present = present,
       seed = if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL)
}

dtu_rng_restore <- function(state) {
  do.call(RNGkind, as.list(state$kind))
  if (isTRUE(state$seed_present)) assign(".Random.seed", state$seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv)
  invisible(NULL)
}

dtu_rng_spec <- function(rng = NULL) {
  defaults <- list(kind = "L'Ecuyer-CMRG", normal_kind = "Inversion", sample_kind = "Rejection")
  if (is.null(rng)) rng <- list()
  if (!is.list(rng) || is.null(names(rng)) && length(rng))
    stop("dtu.rng must be a named YAML mapping.", call. = FALSE)
  unknown <- setdiff(names(rng), names(defaults))
  if (length(unknown)) stop("Unknown dtu.rng setting(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  out <- utils::modifyList(defaults, rng, keep.null = TRUE)
  for (key in names(defaults)) {
    value <- out[[key]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
      stop("dtu.rng.", key, " must be one nonempty character value.", call. = FALSE)
  }
  # Validate through R itself, while preserving both the algorithm and seed,
  # including callers that have not yet created .Random.seed.
  previous <- dtu_rng_snapshot()
  on.exit(dtu_rng_restore(previous), add = TRUE)
  tryCatch(do.call(RNGkind, list(kind = out$kind, normal.kind = out$normal_kind,
                               sample.kind = out$sample_kind)),
           error = function(e) stop("Invalid dtu.rng settings: ", conditionMessage(e), call. = FALSE))
  out
}

dtu_with_rng <- function(rng, seed, expr) {
  previous <- dtu_rng_snapshot()
  on.exit(dtu_rng_restore(previous), add = TRUE)
  rng <- dtu_rng_spec(rng)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 1 || seed > .Machine$integer.max)
    stop("DTU cohort seed must be one positive integer.", call. = FALSE)
  do.call(RNGkind, list(kind = rng$kind, normal.kind = rng$normal_kind,
                       sample.kind = rng$sample_kind))
  set.seed(as.integer(seed))
  force(expr)
}

dtu_rng_contract <- function(rng, master_seed, cohorts) {
  rng <- dtu_rng_spec(rng)
  data.frame(cohort = as.character(cohorts),
             seed = vapply(as.character(cohorts), function(cohort) dtu_cohort_seed(master_seed, cohort), integer(1)),
             kind = rng$kind, normal_kind = rng$normal_kind, sample_kind = rng$sample_kind,
             contract_version = "dtu-rng-v1", stringsAsFactors = FALSE, row.names = NULL)
}

dtu_empty_prepared <- function(cohort, samples, flow) {
  list(cohort = cohort, dm = NULL, empty = TRUE,
       counts = tibble::tibble(gene_id = character(), feature_id = character()),
       samples = tibble::as_tibble(samples), design = NULL, coefficient = NA_character_,
       feature_map = tibble::tibble(gene_id = character(), feature_id = character(), isoform_harmonized_id = character(), structural_category = character()),
       testing_family = tibble::tibble(gene_id = character(), tested_isoforms = integer(), detected_isoforms = integer(), pooled_mean_gene_CPM = numeric(), gencode_annotated_transcripts = integer()),
       flow = flow)
}

dtu_empty_fit <- function(cohort) {
  list(cohort = cohort, fitted_dm = NULL, failures = tibble::tibble(),
       gene = tibble::tibble(gene_id = character(), raw_dtu_p_value = numeric(), dtu_p_value = numeric(), dtu_bh = numeric(), dtu_fdr = numeric(), fit_status = character(), numerical_failure_stage = character(), CR_boundary_associated_P1 = logical(), dtu_primary_bh05 = logical(), dtu_exploratory_bh10 = logical(), dtu_exploratory_only_bh10 = logical()),
       feature = tibble::tibble(gene_id = character(), feature_id = character(), isoform_dtu_p_value = numeric()))
}

dtu_filter_flow_plots <- function(flows, cfg, dirs) {
  if (!nrow(flows) || !isTRUE(cfg$reporting$enabled %||% TRUE)) return(invisible(NULL))
  for (cohort in unique(flows$cohort)) {
    x <- as.data.frame(flows[flows$cohort == cohort, , drop = FALSE])
    cd <- rna_dirs(file.path(dirs$root, cohort))
    long <- tidyr::pivot_longer(x, c("genes", "isoforms"), names_to = "unit", values_to = "retained")
    long$unit <- ifelse(long$unit == "genes", "Genes", "Isoforms")
    view <- rna_flow_plot(as.data.frame(long), cfg)
    subcfg <- cfg; subcfg$reporting$separate_panels_for_composites <- FALSE
    rna_save_plot(view$plot, "DTU_filtering_testing_family_flow", cd, subcfg,
                  view$width, view$height)
    # Meaningful standalone panels have explicit names, not spacer components.
    if (isTRUE(cfg$reporting$separate_panels_for_composites))
      for (name in names(view$components)) rna_save_plot(view$components[[name]],
        paste0("DTU_filtering_testing_family_flow__", name), cd, subcfg, 7, 5.2)
  }
  invisible(NULL)
}

dtu_reference_complexity <- function(cfg, dirs) {
  path <- cfg$references$annotation_gtf
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    if (isTRUE(cfg$enrichment$enabled) && isTRUE(cfg$enrichment$matched_background)) {
      stop("Matched DTU enrichment needs references.annotation_gtf for annotated transcript complexity.", call. = FALSE)
    }
    return(tibble::tibble(gene_id = character(), gencode_annotated_transcripts = integer()))
  }
  rna_require("rtracklayer")
  hash <- rna_hash(list(content = digest::digest(file = path, algo = "sha256"), method = "gene_transcript_alias_counts_v1"))
  cache_path <- file.path(dirs$cache, "reference_transcript_complexity.rds")
  cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
  if (!is.null(cached) && identical(cached$hash, hash)) {
    rna_write_tsv(cached$table, file.path(dirs$tables, "reference_transcript_complexity.tsv.gz"))
    return(cached$table)
  }
  gr <- rtracklayer::import(path, feature.type = "transcript")
  x <- as.data.frame(gr)
  if ("type" %in% names(x)) x <- x[x$type == "transcript", , drop = FALSE]
  if (!all(c("gene_id", "transcript_id") %in% names(x))) stop("Reference GTF lacks gene_id/transcript_id.")
  if (!"gene_name" %in% names(x)) x$gene_name <- NA_character_
  clean <- function(z) { z <- trimws(as.character(z)); z[z %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_; z }
  strip <- function(z) sub("\\.[0-9]+$", "", clean(z))
  aliases <- unique(rbind(data.frame(gene_id = strip(x$gene_id), transcript_id = strip(x$transcript_id)),
                          data.frame(gene_id = clean(x$gene_name), transcript_id = strip(x$transcript_id))))
  aliases <- aliases[!is.na(aliases$gene_id) & nzchar(aliases$gene_id) & !is.na(aliases$transcript_id), , drop = FALSE]
  out <- dplyr::count(tibble::as_tibble(aliases), gene_id, name = "gencode_annotated_transcripts")
  rna_atomic_rds(list(hash = hash, table = out), cache_path)
  rna_write_tsv(out, file.path(dirs$tables, "reference_transcript_complexity.tsv.gz"))
  out
}

dtu_effect_tables <- function(prepared, fit, stages, alphas, support_cutoff) {
  if (isTRUE(prepared$empty)) return(list(feature = prepared$feature_map, paired = tibble::tibble(), usage = tibble::tibble(), driver = tibble::tibble(gene_id = character())))
  env <- dtu_engine_environment()
  env$prepared <- prepared; env$fit <- fit; env$stages <- stages
  env$alphas <- alphas; env$support_cutoff <- support_cutoff
  eval(quote({
    sample_columns <- as.character(prepared$samples$sample_id)
    usage <- prepared$counts %>%
      left_join(prepared$feature_map, by = c("gene_id", "feature_id")) %>%
      pivot_longer(all_of(sample_columns), names_to = "sample_id", values_to = "FL") %>%
      left_join(prepared$samples %>% mutate(patient = as.character(patient), timepoint = as.character(timepoint)), by = "sample_id") %>%
      group_by(gene_id, sample_id) %>%
      mutate(gene_FL = sum(FL), isoform_fraction = ifelse(gene_FL > 0, FL / gene_FL, 0)) %>% ungroup()
    paired <- usage %>% select(gene_id, feature_id, isoform_harmonized_id, structural_category, patient, timepoint, FL, gene_FL, isoform_fraction) %>%
      pivot_wider(names_from = timepoint, values_from = c(FL, gene_FL, isoform_fraction), names_sep = "_") %>%
      mutate(patient_pair_evaluable = gene_FL_primary > 0 & gene_FL_recurrent > 0,
             delta_if = ifelse(patient_pair_evaluable, isoform_fraction_recurrent - isoform_fraction_primary, NA_real_), abs_delta_if = abs(delta_if))
    finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
    feature <- paired %>% group_by(gene_id, feature_id, isoform_harmonized_id, structural_category) %>%
      summarise(mean_primary_fraction = finite_mean(isoform_fraction_primary[gene_FL_primary > 0]),
                mean_recurrent_fraction = finite_mean(isoform_fraction_recurrent[gene_FL_recurrent > 0]),
                mean_delta_if = finite_mean(delta_if), median_delta_if = if (any(is.finite(delta_if))) median(delta_if[is.finite(delta_if)]) else NA_real_,
                max_patient_abs_delta_if = if (any(is.finite(abs_delta_if))) max(abs_delta_if[is.finite(abs_delta_if)]) else NA_real_,
                patients_increased = sum(delta_if >= support_cutoff, na.rm = TRUE), patients_decreased = sum(delta_if <= -support_cutoff, na.rm = TRUE),
                evaluable_patient_pairs = sum(patient_pair_evaluable), .groups = "drop") %>%
      mutate(directional_support = ifelse(mean_delta_if > 0, patients_increased, ifelse(mean_delta_if < 0, patients_decreased, 0L))) %>%
      left_join(fit$feature, by = c("gene_id", "feature_id"))
    for (alpha in alphas) {
      tag <- dtu_alpha_tag(alpha)
      current <- stages[[tag]]$feature
      if (!is.null(current) && nrow(current)) feature <- left_join(feature, current, by = c("gene_id", "feature_id"))
      value <- paste0("stageR_isoform_adjusted_", tag)
      if (!value %in% names(feature)) feature[[value]] <- NA_real_
      feature[[paste0("stageR_confirmed_", tag)]] <- is.finite(feature[[value]]) & feature[[value]] <= alpha
    }
    driver <- feature %>% group_by(gene_id) %>% arrange(desc(abs(mean_delta_if)), feature_id, .by_group = TRUE) %>% slice_head(n = 1L) %>% ungroup() %>%
      transmute(gene_id, DTU_driver_isoform = isoform_harmonized_id, DTU_driver_structure = structural_category,
                DTU_driver_mean_delta_if = mean_delta_if, DTU_driver_abs_mean_delta_if = abs(mean_delta_if), DTU_driver_directional_support = directional_support)
    list(feature = feature, paired = paired, usage = usage, driver = driver)
  }), envir = env)
}

dtu_patient_summary <- function(effects, genes, support_cutoff) {
  if (!nrow(effects$paired)) return(tibble::tibble())
  env <- dtu_engine_environment(); env$effects <- effects; env$genes <- genes; env$cutoff <- support_cutoff
  eval(quote({
    largest <- effects$paired %>% group_by(gene_id, patient) %>%
      arrange(desc(abs_delta_if), feature_id, .by_group = TRUE) %>%
      summarise(patient_max_abs_delta_if_from_DTU_family = if (any(is.finite(abs_delta_if))) max(abs_delta_if[is.finite(abs_delta_if)]) else NA_real_,
        largest_shift_isoform_in_DTU_family = if (any(is.finite(abs_delta_if))) isoform_harmonized_id[which.max(replace(abs_delta_if, !is.finite(abs_delta_if), -Inf))] else NA_character_,
        largest_shift_delta_if = if (any(is.finite(abs_delta_if))) delta_if[which.max(replace(abs_delta_if, !is.finite(abs_delta_if), -Inf))] else NA_real_, .groups = "drop")
    shifts <- effects$paired %>% inner_join(effects$driver, by = c("gene_id", "isoform_harmonized_id" = "DTU_driver_isoform")) %>%
      transmute(gene_id, patient, cohort_driver_isoform = isoform_harmonized_id, cohort_driver_mean_delta_if = DTU_driver_mean_delta_if,
                cohort_driver_patient_delta_if = delta_if, patient_pair_evaluable,
                patient_any_shift_above_threshold = is.finite(delta_if) & abs(delta_if) >= cutoff,
                patient_direction_supported = patient_pair_evaluable & is.finite(delta_if) & is.finite(DTU_driver_mean_delta_if) & abs(delta_if) >= cutoff & sign(delta_if) == sign(DTU_driver_mean_delta_if),
                patient_direction_opposed = patient_pair_evaluable & is.finite(delta_if) & is.finite(DTU_driver_mean_delta_if) & abs(delta_if) >= cutoff & sign(delta_if) == -sign(DTU_driver_mean_delta_if))
    left_join(largest, shifts, by = c("gene_id", "patient")) %>%
      mutate(patient_any_shift_above_threshold = is.finite(patient_max_abs_delta_if_from_DTU_family) & patient_max_abs_delta_if_from_DTU_family >= cutoff) %>%
      left_join(genes %>% select(gene_id, dtu_p_value, dtu_fdr, fit_status, any_of("architecture_class")), by = "gene_id")
  }), env)
}

dtu_dominance_tables <- function(effects, genes, cfg) {
  empty <- list(gene = tibble::tibble(gene_id = character()), patient = tibble::tibble(gene_id = character(), patient = character()))
  if (!nrow(effects$usage)) return(empty)
  env <- dtu_engine_environment(); env$effects <- effects; env$genes <- genes; env$cfg <- cfg
  eval(quote({
    annotation <- feature <- effects$feature %>% select(isoform_harmonized_id, structural_category) %>% distinct(isoform_harmonized_id, .keep_all = TRUE)
    dominance <- tibble(gene_id = genes$gene_id)
    if (all(c("has_dominant_primary", "has_dominant_recurrent", "dominant_primary", "dominant_recurrent") %in% names(genes))) {
      dominance <- genes %>% select(gene_id, has_dominant_primary, has_dominant_recurrent, dominant_primary, dominant_recurrent) %>%
        left_join(annotation %>% rename(dominant_primary = isoform_harmonized_id, primary_dominant_structure = structural_category), by = "dominant_primary") %>%
        left_join(annotation %>% rename(dominant_recurrent = isoform_harmonized_id, recurrent_dominant_structure = structural_category), by = "dominant_recurrent") %>%
        mutate(dominant_isoform_switched = ifelse(has_dominant_primary & has_dominant_recurrent, dominant_primary != dominant_recurrent, NA),
          dominant_structural_category_switched = ifelse(has_dominant_primary & has_dominant_recurrent & !is.na(primary_dominant_structure) & !is.na(recurrent_dominant_structure), primary_dominant_structure != recurrent_dominant_structure, NA)) %>%
        select(gene_id, primary_dominant_structure, recurrent_dominant_structure, dominant_isoform_switched, dominant_structural_category_switched)
    }
    patient <- effects$usage %>% filter(gene_FL > 0) %>% group_by(gene_id, patient, timepoint) %>%
      group_modify(~ {
        ordering <- order(-.x$isoform_fraction, .x$isoform_harmonized_id)
        top <- ordering[1L]; runner <- if (length(ordering) > 1L) .x$isoform_fraction[ordering[2L]] else 0
        fraction <- .x$isoform_fraction[top]; margin <- fraction - runner
        clear <- is.finite(fraction) && fraction >= cfg$remodeling$definition$dominance_mean_if && margin >= cfg$remodeling$definition$dominance_margin
        tibble(has_dominant = clear, dominant_isoform = if (clear) .x$isoform_harmonized_id[top] else NA_character_, top_isoform = .x$isoform_harmonized_id[top], top_fraction = fraction, dominance_margin = margin)
      }) %>% ungroup() %>%
      pivot_wider(names_from = timepoint, values_from = c(has_dominant, dominant_isoform, top_isoform, top_fraction, dominance_margin))
    for (role in c("primary", "recurrent")) {
      for (name in c("has_dominant", "dominant_isoform", "top_isoform", "top_fraction", "dominance_margin")) {
        key <- paste(name, role, sep = "_")
        if (!key %in% names(patient)) patient[[key]] <- if (name == "has_dominant") NA else if (name %in% c("dominant_isoform", "top_isoform")) NA_character_ else NA_real_
      }
    }
    patient <- patient %>% mutate(patient_dominant_switch = ifelse(has_dominant_primary & has_dominant_recurrent, dominant_isoform_primary != dominant_isoform_recurrent, NA))
    support <- patient %>% group_by(gene_id) %>% summarise(patients_with_dominant_switch = sum(patient_dominant_switch, na.rm = TRUE), patients_evaluable_for_dominance = sum(!is.na(patient_dominant_switch)), .groups = "drop")
    list(gene = left_join(dominance, support, by = "gene_id"), patient = patient)
  }), env)
}

dtu_gene_examples <- function(result, cfg, dirs) {
  if (!isTRUE(cfg$reporting$enabled %||% TRUE) || !nrow(result$gene)) return(invisible(NULL))
  alpha <- max(as.numeric(cfg$multiplicity$fdr_thresholds)); genes <- result$gene
  genes <- genes[is.finite(genes$dtu_fdr) & genes$dtu_fdr < alpha, , drop = FALSE]
  if (!nrow(genes)) return(invisible(NULL))
  priority <- genes$gene_id %in% unlist(cfg$dtu$priority_example_genes)
  genes <- genes[order(!priority, genes$dtu_fdr, -genes$DTU_driver_abs_mean_delta_if, genes$gene_id), , drop = FALSE]
  if (!isTRUE(cfg$dtu$plot_all_dtu_genes)) {
    keep <- unlist(lapply(split(seq_len(nrow(genes)), as.character(genes$architecture_class)), head, cfg$dtu$examples_per_class %||% 3L), use.names = FALSE)
    genes <- genes[sort(unique(c(keep, which(genes$gene_id %in% unlist(cfg$dtu$priority_example_genes))))), , drop = FALSE]
  }
  if ((cfg$dtu$maximum_gene_plots %||% 0L) > 0L) genes <- head(genes, cfg$dtu$maximum_gene_plots)
  for (gene in genes$gene_id) {
    row <- genes[genes$gene_id == gene, , drop = FALSE]
    features <- result$effects$feature[result$effects$feature$gene_id == gene, , drop = FALSE]
    # Rank confirmation first, then effect, using both alpha runs independently.
    confirm_columns <- intersect(paste0("stageR_confirmed_", vapply(sort(as.numeric(cfg$multiplicity$fdr_thresholds)), dtu_alpha_tag, character(1))), names(features))
    keys <- lapply(confirm_columns, function(nm) -as.integer(features[[nm]]))
    ordering <- do.call(order, c(keys, list(-abs(features$mean_delta_if), features$feature_id)))
    features <- head(features[ordering, , drop = FALSE], cfg$dtu$isoforms_per_gene_example %||% 2L)
    x <- result$effects$usage[result$effects$usage$gene_id == gene & result$effects$usage$feature_id %in% features$feature_id, , drop = FALSE]
    x$isoform <- factor(x$feature_id, levels = features$feature_id, labels = paste0(seq_len(nrow(features)), ": ", features$structural_category))
    x$timepoint <- factor(x$timepoint, levels = c("primary", "recurrent"), labels = c(cfg$comparison$reference, cfg$comparison$comparison))
    p <- ggplot2::ggplot(x, ggplot2::aes(timepoint, isoform_fraction, group = patient, color = patient)) + ggplot2::geom_line(linewidth = .8) + ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~isoform, nrow = 1L) + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + rna_theme(cfg) +
      ggplot2::labs(title = paste0(gene, ": paired isoform usage"), subtitle = paste0("Class ", row$architecture_class, "; FDR = ", format(row$dtu_fdr, digits = 3)), x = NULL, y = "Isoform proportion within gene", color = "Patient")
    pp <- cfg$reporting$palettes$patients
    if (!is.null(pp) && length(pp) >= length(unique(x$patient))) p <- p + ggplot2::scale_color_manual(values = unlist(pp))
    stem <- paste0("DTU_paired_isoform_trajectories_", gsub("[^A-Za-z0-9._-]", "_", gene))
    rna_write_tsv(x, file.path(dirs$tables, paste0(stem, "_plot_data.tsv.gz")))
    rna_write_tsv(features, file.path(dirs$tables, paste0(stem, "_isoform_key.tsv")))
    rna_save_plot(p, stem, dirs, cfg, width = max(6, 3.5 * nrow(features)), height = 5)
  }
  invisible(NULL)
}

dtu_plot_outputs <- function(result, cfg, dirs, alpha) {
  if (!isTRUE(cfg$reporting$enabled %||% TRUE) || !nrow(result$gene)) return(invisible(NULL))
  genes <- result$gene; alpha_label <- paste0("FDR < ", alpha); flag <- is.finite(genes$dtu_fdr) & genes$dtu_fdr < alpha
  if ("architecture_class" %in% names(genes)) {
    x <- dplyr::count(genes[flag, , drop = FALSE], architecture_class, name = "genes")
    if (nrow(x)) {
      p <- ggplot2::ggplot(x, ggplot2::aes(x = architecture_class, y = genes, fill = architecture_class)) + ggplot2::geom_col() +
        ggplot2::labs(title = "DTU genes by remodeling class", subtitle = alpha_label, x = "Remodeling class", y = "Genes") + rna_theme(cfg)
      cp <- cfg$reporting$palettes$remodeling_classes
      if (!is.null(cp)) p <- p + ggplot2::scale_fill_manual(values = unlist(cp), guide = "none")
      rna_save_plot(p, paste0("DTU_remodeling_classes_", dtu_fdr_tag(alpha)), dirs, cfg, width = 6, height = 4.5)
    }
  }
  if ("DTU_driver_directional_support" %in% names(genes)) {
    x <- dplyr::count(genes[flag, , drop = FALSE], DTU_driver_directional_support, name = "genes")
    if (nrow(x)) {
      p <- ggplot2::ggplot(x, ggplot2::aes(x = DTU_driver_directional_support, y = genes)) + ggplot2::geom_col(fill = "#4E79A7") +
        ggplot2::labs(title = "Patient support for DTU drivers", subtitle = alpha_label, x = "Patients supporting the driver direction", y = "Genes") + rna_theme(cfg)
      rna_save_plot(p, paste0("DTU_patient_support_", dtu_fdr_tag(alpha)), dirs, cfg, width = 6, height = 4.5)
    }
  }
  feature <- result$effects$feature
  confirmed <- paste0("stageR_confirmed_", dtu_alpha_tag(alpha))
  if (nrow(feature) && confirmed %in% names(feature)) {
    x <- dplyr::count(feature[feature[[confirmed]], , drop = FALSE], structural_category, name = "isoforms")
    if (nrow(x)) {
      p <- ggplot2::ggplot(x, ggplot2::aes(x = reorder(structural_category, isoforms), y = isoforms, fill = structural_category)) + ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::labs(title = "Structure of confirmed DTU isoforms", subtitle = alpha_label, x = NULL, y = "Isoforms") + rna_theme(cfg)
      sp <- cfg$reporting$palettes$structural_categories
      if (!is.null(sp)) p <- p + ggplot2::scale_fill_manual(values = unlist(sp), guide = "none")
      rna_save_plot(p, paste0("DTU_confirmed_structure_", dtu_fdr_tag(alpha)), dirs, cfg, width = 6, height = 4.5)
    }
  }
  if (nrow(result$patient)) {
    chosen <- genes$gene_id[flag]
    x <- result$patient[result$patient$gene_id %in% chosen, , drop = FALSE]
    if (nrow(x)) {
      p <- ggplot2::ggplot(x, ggplot2::aes(x = patient, y = cohort_driver_patient_delta_if)) + ggplot2::geom_hline(yintercept = 0, color = "grey60") +
        ggplot2::geom_boxplot(outlier.shape = NA, fill = "#D8E5F0") +
        ggplot2::labs(title = "Patient-level DTU driver changes", subtitle = alpha_label, x = "Patient", y = "Change in isoform fraction") + rna_theme(cfg)
      rna_save_plot(p, paste0("DTU_patient_driver_distributions_", dtu_fdr_tag(alpha)), dirs, cfg, width = 7, height = 4.5)
    }
  }
  for (metric in intersect(c("jsd", "max_abs_delta_if"), names(genes))) {
    x <- genes[is.finite(genes[[metric]]) & is.finite(genes$raw_dtu_p_value), , drop = FALSE]
    if (nrow(x)) {
      x$effect <- x[[metric]]; x$score <- -log10(pmax(x$raw_dtu_p_value, .Machine$double.xmin)); x$significant <- x$dtu_fdr < alpha
      label <- if (metric == "jsd") "Jensen-Shannon divergence" else "Maximum change in isoform fraction"
      p <- ggplot2::ggplot(x, ggplot2::aes(effect, score, color = significant)) + ggplot2::geom_point(alpha = .7, size = 1.7) +
        ggplot2::scale_color_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = "#B2182B"), labels = c("Other genes", alpha_label)) + rna_theme(cfg) +
        ggplot2::labs(title = "DTU evidence and isoform redistribution", x = label, y = expression(-log[10](P)), color = NULL)
      rna_save_plot(p, paste0("DTU_effect_evidence_", metric, "_", dtu_fdr_tag(alpha)), dirs, cfg, 6, 5)
    }
  }
  if (nrow(result$patient)) {
    ids <- head(genes$gene_id[order(genes$dtu_fdr, genes$gene_id)][genes$dtu_fdr[order(genes$dtu_fdr, genes$gene_id)] < alpha], cfg$dtu$maximum_heatmap_genes %||% 60L)
    x <- result$patient[result$patient$gene_id %in% ids, , drop = FALSE]
    if (nrow(x)) {
      x$gene_id <- factor(x$gene_id, levels = rev(ids))
      p <- ggplot2::ggplot(x, ggplot2::aes(patient, gene_id, fill = cohort_driver_patient_delta_if)) + ggplot2::geom_tile(color = "white", linewidth = .25) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey80") + rna_theme(cfg) +
        ggplot2::labs(title = "Paired changes in DTU driver usage", subtitle = alpha_label, x = "Patient", y = NULL, fill = "Change in\nisoform fraction")
      rna_save_plot(p, paste0("DTU_patient_driver_heatmap_", dtu_fdr_tag(alpha)), dirs, cfg, 7, max(4.5, 2 + .24 * length(ids)))
    }
  }
  if (all(c("primary_dominant_structure", "recurrent_dominant_structure") %in% names(genes))) {
    x <- genes[flag & !is.na(genes$dominant_isoform_switched) & !is.na(genes$primary_dominant_structure) & !is.na(genes$recurrent_dominant_structure), , drop = FALSE]
    if (nrow(x)) {
      tab <- dplyr::count(x, primary_dominant_structure, recurrent_dominant_structure, name = "genes")
      rna_write_tsv(tab, file.path(dirs$tables, paste0("DTU_dominant_structural_transitions_", dtu_fdr_tag(alpha), ".tsv")))
      p <- ggplot2::ggplot(tab, ggplot2::aes(primary_dominant_structure, recurrent_dominant_structure, fill = genes)) + ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = genes), size = 4.5) +
        ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") + rna_theme(cfg) + ggplot2::labs(title = "Dominant-isoform structural transitions", subtitle = alpha_label, x = "Reference condition structure", y = "Comparison condition structure", fill = "Genes")
      rna_save_plot(p, paste0("DTU_dominant_structural_transitions_", dtu_fdr_tag(alpha)), dirs, cfg, 7, 6)
    }
  }
  invisible(NULL)
}

dtu_broad_structure_outputs <- function(broad, genes, cfg, dirs, alphas) {
  if (!nrow(broad)) return(tibble::tibble(gene_id = character(), timepoint = character(), mass = numeric(), evaluable = logical(), mass_conserved = logical()))
  e <- dtu_engine_environment(); e$broad <- broad; e$genes <- genes; e$cfg <- cfg; e$dirs <- dirs; e$alphas <- alphas
  eval(quote({
    long <- broad %>% select(gene_id, structural_category_display, mean_primary_category_fraction, mean_recurrent_category_fraction) %>%
      pivot_longer(c(mean_primary_category_fraction, mean_recurrent_category_fraction), names_to = "timepoint", values_to = "category_fraction") %>%
      mutate(timepoint = ifelse(timepoint == "mean_primary_category_fraction", "primary", "recurrent"))
    audit <- long %>% group_by(gene_id, timepoint) %>% summarise(mass = sum(category_fraction, na.rm = TRUE), evaluable = any(is.finite(category_fraction)), .groups = "drop") %>%
      mutate(mass_conserved = ifelse(evaluable, abs(mass - 1) <= 1e-8, mass == 0))
    rna_write_tsv(audit, file.path(dirs$qc, "broad_structure_fraction_mass_audit.tsv.gz"))
    categories <- sort(unique(long$structural_category_display))
    for (alpha in alphas) {
      target <- genes$gene_id[is.finite(genes$dtu_fdr) & genes$dtu_fdr < alpha]
      if (!length(target)) next
      # A category absent from one gene contributes zero. The denominator is
      # all evaluable target genes, not only genes containing that category.
      included <- audit %>% filter(gene_id %in% target, evaluable) %>% select(gene_id, timepoint)
      x <- tidyr::crossing(included, structural_category_display = categories) %>%
        left_join(long, by = c("gene_id", "timepoint", "structural_category_display")) %>%
        mutate(category_fraction = coalesce(category_fraction, 0)) %>%
        group_by(timepoint, structural_category_display) %>%
        summarise(mean_gene_fraction = mean(category_fraction), evaluable_genes = n_distinct(gene_id), .groups = "drop")
      rna_write_tsv(x, file.path(dirs$tables, paste0("DTU_gene_weighted_structural_usage_", dtu_fdr_tag(alpha), ".tsv")))
      if (nrow(x) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
        p <- ggplot(x, aes(structural_category_display, mean_gene_fraction, fill = timepoint)) + geom_col(position = position_dodge(width = .8), width = .7) +
          scale_y_continuous(labels = percent) + rna_theme(cfg) +
          labs(title = "Structural usage among DTU genes", subtitle = paste0("FDR < ", alpha), x = NULL, y = "Mean proportion per gene", fill = NULL)
        palette <- cfg$reporting$palettes$timepoints %||% cfg$reporting$palettes$timepoint %||% c(primary = "#4C78A8", recurrent = "#E45756")
        p <- p + scale_fill_manual(values = unlist(palette), labels = c(primary = cfg$comparison$reference, recurrent = cfg$comparison$comparison))
        rna_save_plot(p, paste0("DTU_gene_weighted_structural_usage_", dtu_fdr_tag(alpha)), dirs, cfg, 7.5, 5)
      }
    }
    audit
  }), e)
}

dtu_upset_components <- function(tab, sets, stem, dirs, cfg, title, subtitle = NULL) {
  if (!length(sets) || !nrow(tab)) return(invisible(NULL))
  view <- rna_upset_plot(tab, sets, cfg, title, subtitle,
    max_intersections = cfg$reporting$maximum_upset_intersections %||% 40L)
  if (is.null(view$plot)) return(invisible(NULL))
  rna_write_tsv(view$counts, file.path(dirs$tables, paste0(stem, "__exact_intersections.tsv")))
  rna_write_tsv(view$sizes, file.path(dirs$tables, paste0(stem, "__set_sizes.tsv")))
  subcfg <- cfg; subcfg$reporting$separate_panels_for_composites <- FALSE
  rna_save_plot(view$plot, stem, dirs, subcfg, view$width, view$height)
  if (isTRUE(cfg$reporting$separate_panels_for_composites))
    for (name in names(view$components)) rna_save_plot(view$components[[name]],
      paste0(stem, "__", name), dirs, subcfg, max(7.5, view$width),
      if (grepl("membership", name)) max(3.6, .55 * length(sets) + 1.7) else 4.5)
  invisible(view)
}

dtu_collective_reports <- function(results, cfg, dirs, alphas, completion = NULL, remodeling = NULL) {
  if (isTRUE(cfg$dtu$write_supplement_workbook)) {
    rna_require("openxlsx")
    workbook <- openxlsx::createWorkbook()
    index <- 0L
    add <- function(x, label) {
      if (!ncol(x)) return(invisible(NULL))
      groups <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / 1000000L))
      if (!length(groups)) groups <- list(integer())
      for (rows in groups) {
        index <<- index + 1L
        name <- substr(paste0(sprintf("%03d_", index), gsub("[^A-Za-z0-9_]", "_", label)), 1, 31)
        openxlsx::addWorksheet(workbook, name)
        openxlsx::writeData(workbook, name, as.data.frame(x[rows, , drop = FALSE]))
        openxlsx::freezePane(workbook, name, firstRow = TRUE)
      }
    }
    for (cohort in names(results)) {
      add(results[[cohort]]$gene, paste0(cohort, "_genes"))
      add(results[[cohort]]$effects$feature, paste0(cohort, "_isoforms"))
      add(results[[cohort]]$patient, paste0(cohort, "_patient_support"))
    }
    if (!is.null(completion)) add(completion, "completion")
    openxlsx::saveWorkbook(workbook, file.path(dirs$tables, "DTU_supplementary_results.xlsx"), overwrite = TRUE)
  }
  if (isTRUE(cfg$dtu$run_upset) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
    for (alpha in alphas) {
      remodel <- remodeling$gene_table
      ids <- sort(unique(c(unlist(lapply(results, function(z) z$gene$gene_id)), remodel$gene_id)))
      tab <- data.frame(gene_id = ids, check.names = FALSE)
      main_patients <- as.character(results[[1L]]$prepared$samples$patient)
      for (cohort in names(results)) {
        x <- results[[cohort]]$gene
        omitted <- setdiff(main_patients, as.character(results[[cohort]]$prepared$samples$patient))
        context <- if (!length(omitted)) "all patients" else paste0("excluding ", paste(omitted, collapse = ", "))
        label <- paste0("DTU: ", context)
        tab[[label]] <- ids %in% x$gene_id[is.finite(x$dtu_fdr) & x$dtu_fdr < alpha]
        if (!is.null(remodel) && all(c("cohort", "gene_id", "downstream_primary_remodeling_set") %in% names(remodel))) {
          selected <- remodel$cohort == cohort & remodel$downstream_primary_remodeling_set %in% TRUE
          tab[[paste0("Remodeling: ", context)]] <- ids %in% remodel$gene_id[selected]
        }
      }
      sets <- setdiff(names(tab), "gene_id")
      sets <- sets[vapply(tab[sets], any, logical(1))]
      rna_write_tsv(tab, file.path(dirs$tables, paste0("DTU_intersection_membership_", dtu_fdr_tag(alpha), ".tsv.gz")))
      if (length(sets)) dtu_upset_components(tab, sets,
        paste0("DTU_remodeling_intersections_", dtu_fdr_tag(alpha)), dirs, cfg,
        "DTU and remodeling overlap", paste0("DTU FDR < ", alpha))
    }
    for (cohort in names(results)) for (alpha in alphas) {
      result <- results[[cohort]]
      targets <- result$gene$gene_id[is.finite(result$gene$dtu_fdr) & result$gene$dtu_fdr < alpha]
      patient <- result$patient
      tab <- data.frame(gene_id = targets, check.names = FALSE)
      for (id in sort(unique(patient$patient))) tab[[id]] <- targets %in% patient$gene_id[patient$patient == id & patient$patient_direction_supported %in% TRUE]
      sets <- setdiff(names(tab), "gene_id"); sets <- sets[vapply(tab[sets], any, logical(1))]
      cd <- rna_dirs(file.path(dirs$root, cohort)); stem <- paste0("DTU_patient_directional_support_", dtu_fdr_tag(alpha))
      rna_write_tsv(tab, file.path(cd$tables, paste0(stem, "_membership.tsv.gz")))
      if (!length(sets)) next
      dtu_upset_components(tab, sets, stem, cd, cfg,
        "Shared DTU driver support across patients", paste0("DTU FDR < ", alpha))
    }
  }
  invisible(NULL)
}

run_dtu <- function(inputs, cfg, cohorts, remodeling = NULL, stage_dir) {
  previous_rng <- dtu_rng_snapshot()
  on.exit(dtu_rng_restore(previous_rng), add = TRUE)
  cfg$dtu$rng <- dtu_rng_spec(cfg$dtu$rng)
  rng_contract <- dtu_rng_contract(cfg$dtu$rng, as.integer(cfg$execution$seeds$dtu), names(cohorts))
  if (!is.null(cfg$dtu$count_source) && cfg$dtu$count_source != "full_depth_integer_counts") stop("DTU count_source must be full_depth_integer_counts.", call. = FALSE)
  if (!is.null(cfg$dtu$engine) && cfg$dtu$engine != "drimseq_1_38_0_adapter_1_3_0") stop("Unsupported DTU engine; this release requires the audited DRIMSeq 1.38.0 adapter.", call. = FALSE)
  rna_require(c("DRIMSeq", "BiocParallel", "MASS", "dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "scales", "forcats", "ggplot2", "digest"))
  if (isTRUE(cfg$dtu$run_stageR)) rna_require("stageR")
  if (isTRUE(cfg$dtu$write_supplement_workbook)) rna_require("openxlsx")
  if (is.null(cfg$reporting$palettes$patients)) {
    patients <- unique(as.character(inputs$sample_manifest$patient))
    cfg$reporting$palettes$patients <- stats::setNames(scales::hue_pal()(length(patients)), patients)
  }
  dirs <- rna_dirs(stage_dir)
  e <- dtu_engine_environment(parent = environment())
  sys.source(dtu_engine_path(cfg, "dtu_core.R"), envir = e)
  e$empty_dtu_prepared <- dtu_empty_prepared
  e$make_bpparam <- function(workers) BiocParallel::SerialParam(progressbar = FALSE)
  e$write_table <- rna_write_tsv; e$atomic_save_rds <- rna_atomic_rds
  e$log_msg <- function(..., level = "INFO") rna_log(..., level = level)
  e$cohort_names <- names(cohorts); e$analysis_cohorts <- cohorts
  e$m06 <- dirs; e$m06_root <- dirs$root
  e$cohort_folder <- stats::setNames(names(cohorts), names(cohorts))
  e$sample_manifest <- tibble::as_tibble(inputs$sample_manifest)
  e$counts_matrix <- as.matrix(inputs$counts_matrix)
  e$isoform_annotation <- tibble::as_tibble(inputs$isoform_annotation)
  if (!identical(rownames(e$counts_matrix), as.character(e$isoform_annotation$isoform_harmonized_id))) stop("DTU count/annotation row mismatch.")
  if (any(!is.finite(e$counts_matrix)) || any(e$counts_matrix < 0) || any(abs(e$counts_matrix - round(e$counts_matrix)) > 1e-8)) stop("DTU requires nonnegative integer FL counts.")
  if (any(grepl("|", e$isoform_annotation$associated_gene, fixed = TRUE), na.rm = TRUE)) stop("Gene identifiers cannot contain the DTU feature separator '|'.")
  paired_design <- isTRUE(cfg$comparison$paired)
  e$cfg <- list(model_contract_version = if (paired_design) "1.3.0-public-0.1.2-rng-v1" else "independent-DTU-v1",
                paired = paired_design, workflow_version = rna_package_version(),
                analysis_profile = if ((cfg$execution$profile %||% "production") %in% c("smoke", "pilot")) "pilot" else "production",
                force_refit = isTRUE(cfg$execution$force_refit) || isTRUE(cfg$dtu$force_refit), seed = as.integer(cfg$execution$seeds$dtu), n_workers = 1L, rng = cfg$dtu$rng,
                pilot_max_genes = as.integer(cfg$dtu$pilot_max_genes %||% 100L),
                filters = cfg$dtu[c("minimum_isoform_fraction", "minimum_samples_isoform_fraction", "minimum_gene_count", "minimum_samples_gene_count", "minimum_isoform_count", "minimum_samples_isoform_count", "minimum_isoforms_per_gene")],
                fdr = list(primary = .05, exploratory = .10), run_stageR = isTRUE(cfg$dtu$run_stageR))
  alphas <- sort(unique(as.numeric(cfg$multiplicity$fdr_thresholds)))
  e$input_fingerprint <- rna_hash(list(inputs = inputs$fingerprint, counts = e$counts_matrix, samples = e$sample_manifest,
                                     parameters = rna_scientific_config(cfg)$dtu, analysis_design = if (paired_design) "paired" else "unpaired", seed = e$cfg$seed, cohorts = cohorts, code = digest::digest(file = dtu_engine_path(cfg, "dtu_core.R"), algo = "sha256"),
                                     DRIMSeq = as.character(utils::packageVersion("DRIMSeq")), rng_contract = rng_contract))
  rng_audit <- rng_contract
  rng_audit$R_version <- R.version.string
  rng_audit$DRIMSeq_version <- as.character(utils::packageVersion("DRIMSeq"))
  rng_audit$model_contract_version <- e$cfg$model_contract_version
  rng_audit$input_fingerprint <- e$input_fingerprint
  rna_write_tsv(rng_audit, file.path(dirs$qc, "dtu_rng_settings.tsv"))
  e$gencode_complexity <- dtu_reference_complexity(cfg, dirs)
  for (cohort in names(cohorts)) rna_dirs(file.path(stage_dir, cohort))
  e$dtu_prepared <- stats::setNames(lapply(names(cohorts), e$prepare_dtu_family), names(cohorts))
  # The targeted refit runner verifies every prepared cohort against saved
  # original inputs before any numerical fit starts. This private callback is
  # absent during ordinary workflow runs and is excluded from scientific keys.
  if (is.function(cfg$.dtu_preflight)) cfg$.dtu_preflight(e$dtu_prepared)
  design_audit <- dplyr::bind_rows(lapply(names(e$dtu_prepared), function(cohort) {
    prepared <- e$dtu_prepared[[cohort]]
    design <- prepared$design
    if (is.null(design)) design <- e$m06_design_matrix(as.data.frame(prepared$samples), paired_design)
    dplyr::bind_cols(tibble::tibble(cohort = cohort, sample_id = prepared$samples$sample_id,
      analysis_design = if (paired_design) "paired" else "unpaired",
      formula = if (paired_design) "~ patient + timepoint" else "~ timepoint",
      null_formula = if (paired_design) "~ patient" else "~ 1"),
      tibble::as_tibble(design, .name_repair = "minimal"))
  }))
  rna_write_tsv(design_audit, file.path(dirs$qc, "dtu_model_design.tsv"))
  flows <- dplyr::bind_rows(lapply(e$dtu_prepared, `[[`, "flow"))
  families <- dplyr::bind_rows(lapply(names(cohorts), function(x) dplyr::mutate(e$dtu_prepared[[x]]$testing_family, cohort = x)))
  rna_write_tsv(flows, file.path(dirs$qc, "dtu_filter_flow.tsv"))
  dtu_filter_flow_plots(flows, cfg, dirs)
  rna_write_tsv(families, file.path(dirs$tables, "frozen_dtu_testing_families.tsv.gz"))
  rna_write_tsv(tibble::tibble(requested_workers = cfg$execution$workers %||% 1L, effective_workers = 1L, backend = "SerialParam", reason = "Preserve complete numerical profile audits"), file.path(dirs$qc, "execution_backend.tsv"))
  e$dtu_fits <- stats::setNames(vector("list", length(cohorts)), names(cohorts))
  for (cohort in names(cohorts)) {
    e$dtu_fits[[cohort]] <- if (isTRUE(e$dtu_prepared[[cohort]]$empty)) dtu_empty_fit(cohort) else e$fit_drimseq_model(e$dtu_prepared[[cohort]])
  }
  results <- list(); fit_qc <- list(); stage_qc <- list(); structure_qc <- list()
  for (cohort in names(cohorts)) {
    cd <- rna_dirs(file.path(stage_dir, cohort)); fit <- e$dtu_fits[[cohort]]; prepared <- e$dtu_prepared[[cohort]]
    fit$gene$dtu_fdr <- fit$gene$dtu_bh
    for (alpha in alphas) fit$gene[[paste0("dtu_", dtu_fdr_tag(alpha))]] <- is.finite(fit$gene$dtu_fdr) & fit$gene$dtu_fdr < alpha
    e$dtu_fits[[cohort]] <- fit
    stages <- list(); gene <- dplyr::left_join(fit$gene, prepared$testing_family, by = "gene_id")
    for (alpha in alphas) {
      tag <- dtu_alpha_tag(alpha)
      stage <- if (isTRUE(cfg$dtu$run_stageR) && !isTRUE(prepared$empty)) e$run_stageR_alpha(cohort, alpha) else list(available = FALSE, reason = if (isTRUE(prepared$empty)) "empty_testing_family" else "disabled", gene = tibble::tibble(gene_id = character()), feature = tibble::tibble(gene_id = character(), feature_id = character()))
      stages[[tag]] <- stage
      if (isTRUE(stage$available)) {
        gene <- dplyr::left_join(gene, stage$gene, by = "gene_id")
        rna_atomic_rds(stage$object, file.path(cd$objects, paste0("stageR_", tag, ".rds")))
      }
      rna_write_tsv(stage$gene, file.path(cd$tables, paste0("stageR_gene_screening_", tag, ".tsv.gz")))
      rna_write_tsv(stage$feature, file.path(cd$tables, paste0("stageR_isoform_confirmation_", tag, ".tsv.gz")))
      value <- paste0("stageR_gene_adjusted_", tag)
      diff <- if (isTRUE(stage$available)) {
        z <- merge(stage$gene, fit$gene[, c("gene_id", "dtu_fdr")], by = "gene_id")
        if (nrow(z) && all(is.finite(z[[value]]))) max(abs(z[[value]] - z$dtu_fdr)) else Inf
      } else NA_real_
      stage_qc[[paste(cohort, tag)]] <- tibble::tibble(cohort = cohort, alpha = alpha, available = isTRUE(stage$available), reason = stage$reason,
          scheduled_genes = nrow(fit$gene), screening_genes = stage$screening_genes %||% 0L, confirmation_family_genes = stage$confirmation_family_genes %||% 0L,
          gene_screen_FDR_max_difference = diff,
          family_conserved = isTRUE(prepared$empty) || !isTRUE(cfg$dtu$run_stageR) || (isTRUE(stage$available) && stage$screening_genes == nrow(fit$gene) && stage$confirmation_family_genes == nrow(fit$gene) && setequal(stage$gene$gene_id, fit$gene$gene_id) && setequal(stage$feature$feature_id, prepared$feature_map$feature_id)))
    }
    effects <- if (paired_design) dtu_effect_tables(prepared, fit, stages, alphas, cfg$dtu$patient_support_min_abs_delta_if) else
      dtu_unpaired_effect_tables(prepared, fit, stages, alphas)
    gene <- dplyr::left_join(gene, effects$driver, by = "gene_id")
    if (!is.null(remodeling$gene_table)) {
      annotation <- tibble::as_tibble(remodeling$gene_table)
      if (all(c("cohort", "gene_id") %in% names(annotation))) {
        annotation <- annotation[as.character(annotation$cohort) == cohort, , drop = FALSE]
        wanted <- intersect(c("gene_id", "architecture_class", "jsd", "max_abs_delta_if", "delta_entropy", "has_dominant_primary", "dominant_primary", "dominant_primary_mean_if", "dominant_primary_margin", "dominant_primary_support_n", "has_dominant_recurrent", "dominant_recurrent", "dominant_recurrent_mean_if", "dominant_recurrent_margin", "dominant_recurrent_support_n", "downstream_primary_remodeling_set", "is_high_confidence_remodeled"), names(annotation))
        annotation <- annotation[, wanted, drop = FALSE]
        if (anyDuplicated(annotation$gene_id)) stop("Duplicate remodeling gene/cohort keys in DTU join.")
        gene <- dplyr::left_join(gene, annotation, by = "gene_id")
      }
    }
    if (!paired_design) gene$analysis_design <- rep("unpaired", nrow(gene))
    if (paired_design && !"architecture_class" %in% names(gene)) gene$architecture_class <- rep("NQ", nrow(gene))
    if (paired_design) {
      gene$architecture_class <- as.character(gene$architecture_class)
      gene$architecture_class[is.na(gene$architecture_class)] <- "NQ"
    }
    dominance <- if (paired_design) dtu_dominance_tables(effects, gene, cfg) else
      list(gene = tibble::tibble(gene_id = character()), patient = tibble::tibble(gene_id = character(), patient = character()))
    gene <- dplyr::left_join(gene, dominance$gene, by = "gene_id")
    patient <- if (paired_design) dtu_patient_summary(effects, gene, cfg$dtu$patient_support_min_abs_delta_if) else
      tibble::tibble(gene_id = character(), patient = character(), patient_direction_supported = logical())
    if (nrow(patient)) patient <- dplyr::left_join(patient, dominance$patient, by = c("gene_id", "patient"))
    if (!is.null(remodeling$patient_inference) && nrow(patient)) {
      pi <- tibble::as_tibble(remodeling$patient_inference)
      if ("cohort" %in% names(pi)) pi <- pi[as.character(pi$cohort) == cohort, , drop = FALSE]
      keep <- intersect(c("gene_id", "patient", "patient_test_p", "patient_test_bh", "patient_effect_pass", "patient_supported_bh05_effect", "patient_supported_bh10_effect"), names(pi))
      if (all(c("gene_id", "patient") %in% keep)) patient <- dplyr::left_join(patient, dplyr::distinct(pi[, keep, drop = FALSE]), by = c("gene_id", "patient"))
    }
    if (all(c("patient_test_bh", "patient_effect_pass") %in% names(patient))) {
      patient$patient_test_fdr <- patient$patient_test_bh
      for (alpha in alphas) patient[[paste0("patient_supported_", dtu_fdr_tag(alpha), "_effect")]] <-
        is.finite(patient$patient_test_fdr) & patient$patient_test_fdr < alpha & patient$patient_effect_pass %in% TRUE
    }
    e$dtu_prepared[[cohort]] <- prepared
    broad <- if (!isTRUE(prepared$empty)) e$build_broad_structural_category_means(cohort) else tibble::tibble()
    structure_qc[[cohort]] <- dplyr::mutate(dtu_broad_structure_outputs(broad, gene, cfg, cd, alphas), cohort = cohort)
    result <- list(cohort = cohort, gene = gene, effects = effects, patient = patient, dominance = dominance, prepared = prepared, fit = fit, stageR = stages, broad_structure = broad)
    results[[cohort]] <- result
    rna_write_tsv(dplyr::mutate(gene, cohort = cohort), file.path(cd$tables, "DTU_gene_results_integrated.tsv.gz"))
    rna_write_tsv(fit$failures, file.path(cd$qc, "drimseq_nonfinite_P1_audit.tsv"))
    rna_write_tsv(dplyr::mutate(effects$feature, cohort = cohort), file.path(cd$tables, "DTU_isoform_effects_and_confirmation.tsv.gz"))
    if (paired_design) {
      rna_write_tsv(effects$paired, file.path(cd$tables, "patient_paired_isoform_fraction_changes.tsv.gz"))
      rna_write_tsv(patient, file.path(cd$tables, "DTU_patient_gene_support.tsv.gz"))
    } else rna_write_tsv(effects$usage, file.path(cd$tables, "independent_sample_isoform_fractions.tsv.gz"))
    rna_write_tsv(broad, file.path(cd$tables, "broad_comparable_isoform_structural_category_mean_usage.tsv.gz"))
    rna_atomic_rds(result, file.path(cd$objects, "DTU_results.rds"))
    if (!is.null(fit$fitted_dm)) rna_atomic_rds(fit$fitted_dm, file.path(cd$objects, "drimseq_fitted_model.rds"))
    for (alpha in alphas) {
      rna_write_tsv(gene[is.finite(gene$dtu_fdr) & gene$dtu_fdr < alpha, , drop = FALSE], file.path(cd$tables, paste0("DTU_genes_", dtu_fdr_tag(alpha), ".tsv.gz")))
      if (paired_design) dtu_plot_outputs(result, cfg, cd, alpha) else dtu_unpaired_plot_outputs(result, cfg, cd, alpha)
    }
    if (paired_design) dtu_gene_examples(result, cfg, cd) else dtu_unpaired_gene_examples(result, cfg, cd)
    unresolved <- if (nrow(fit$gene)) sum(fit$gene$numerical_failure_stage != "none" & !fit$gene$CR_boundary_associated_P1, na.rm = TRUE) else 0L
    fit_qc[[cohort]] <- tibble::tibble(cohort = cohort, scheduled_genes = nrow(prepared$testing_family), result_genes = nrow(fit$gene),
       finite_test_genes = sum(fit$gene$fit_status == "finite_test"), unavailable_P1_genes = sum(fit$gene$fit_status != "finite_test"), unexplained_numerical_genes = unresolved,
       family_conserved = setequal(prepared$testing_family$gene_id, fit$gene$gene_id), status = if (isTRUE(prepared$empty)) "empty_testing_family" else "completed")
  }
  fit_qc <- dplyr::bind_rows(fit_qc); stage_qc <- dplyr::bind_rows(stage_qc)
  structure_qc <- dplyr::bind_rows(structure_qc)
  rna_write_tsv(fit_qc, file.path(dirs$qc, "drimseq_fit_and_discovery_summary.tsv"))
  rna_write_tsv(stage_qc, file.path(dirs$qc, "stageR_completion_summary.tsv"))
  sensitivity <- dtu_sensitivity_outputs(results, cfg, dirs, alphas)
  enrichment <- if (isTRUE(cfg$enrichment$enabled)) run_dtu_enrichment(results, cfg, stage_dir) else list(status = "disabled")
  pass <- all(fit_qc$family_conserved) && all(fit_qc$unexplained_numerical_genes == 0L) && all(stage_qc$family_conserved) &&
    all(is.na(stage_qc$gene_screen_FDR_max_difference) | stage_qc$gene_screen_FDR_max_difference <= 1e-8) && all(structure_qc$mass_conserved)
  completion <- tibble::tibble(gate = c("Frozen gene families retained", "No unexplained numerical fits", "stageR complete gene families", "stageR screening agrees with full-family FDR", "Broad structural fractions conserve mass"),
      pass = c(all(fit_qc$family_conserved), all(fit_qc$unexplained_numerical_genes == 0L), all(stage_qc$family_conserved), all(is.na(stage_qc$gene_screen_FDR_max_difference) | stage_qc$gene_screen_FDR_max_difference <= 1e-8), all(structure_qc$mass_conserved)))
  rna_write_tsv(completion, file.path(dirs$qc, "completion_gates.tsv"))
  if (paired_design) dtu_collective_reports(results, cfg, dirs, alphas, completion, remodeling)
  out <- list(analysis_design = if (paired_design) "paired" else "unpaired", completed = pass, complete = pass, gates = completion, cohort_results = results,
       gene_table = dplyr::bind_rows(lapply(names(results), function(x) dplyr::mutate(results[[x]]$gene, cohort = x))),
       feature_table = dplyr::bind_rows(lapply(names(results), function(x) dplyr::mutate(results[[x]]$effects$feature, cohort = x))),
       patient_table = dplyr::bind_rows(lapply(names(results), function(x) dplyr::mutate(results[[x]]$patient, cohort = x))),
       fit_qc = fit_qc, stageR_qc = stage_qc, structure_qc = structure_qc, sensitivity = sensitivity, enrichment = enrichment, fingerprint = e$input_fingerprint, rng_contract = rng_contract)
  rna_atomic_rds(out, file.path(dirs$objects, if (pass) "DTU_complete.rds" else "DTU_audit_only.rds"))
  if (!pass && isTRUE(cfg$execution$strict_completion)) stop("DTU failed completion; numerical/family audit outputs were retained.", call. = FALSE)
  out
}

dtu_sensitivity_outputs <- function(results, cfg, dirs, alphas) {
  if (length(results) < 2L) return(list(status = "not_requested"))
  main <- results[[1L]]$gene; summaries <- list(); comparisons <- list()
  for (cohort in names(results)[-1L]) {
    other <- results[[cohort]]$gene
    a <- main[, c("gene_id", "raw_dtu_p_value", "dtu_fdr", "fit_status"), drop = FALSE]
    b <- other[, c("gene_id", "raw_dtu_p_value", "dtu_fdr", "fit_status"), drop = FALSE]
    pair <- merge(a, b, by = "gene_id", all = TRUE, suffixes = c("_main", "_sensitivity"))
    finite <- is.finite(pair$raw_dtu_p_value_main) & is.finite(pair$raw_dtu_p_value_sensitivity)
    rho <- if (sum(finite) >= 3L && stats::sd(pair$raw_dtu_p_value_main[finite]) > 0 && stats::sd(pair$raw_dtu_p_value_sensitivity[finite]) > 0) stats::cor(pair$raw_dtu_p_value_main[finite], pair$raw_dtu_p_value_sensitivity[finite], method = "spearman") else NA_real_
    pair$cohort <- cohort; comparisons[[cohort]] <- pair
    for (alpha in alphas) {
      main_sig <- main$gene_id[main$dtu_fdr < alpha]; other_sig <- other$gene_id[other$dtu_fdr < alpha]
      summaries[[paste(cohort, alpha)]] <- tibble::tibble(cohort = cohort, alpha = alpha, shared_finite_tests = sum(finite), raw_P_spearman = rho,
             main_discoveries = length(main_sig), sensitivity_discoveries = length(other_sig), shared_discoveries = length(intersect(main_sig, other_sig)),
             main_only = length(setdiff(main_sig, other_sig)), sensitivity_only = length(setdiff(other_sig, main_sig)))
    }
    rna_write_tsv(pair, file.path(dirs$tables, paste0("DTU_sensitivity_", cohort, ".tsv.gz")))
    if (any(finite) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
      x <- pair[finite, , drop = FALSE]
      x$main_score <- -log10(pmax(x$raw_dtu_p_value_main, .Machine$double.xmin)); x$sensitivity_score <- -log10(pmax(x$raw_dtu_p_value_sensitivity, .Machine$double.xmin))
      p <- ggplot2::ggplot(x, ggplot2::aes(x = main_score, y = sensitivity_score)) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + ggplot2::geom_point(alpha = .5, size = 1.5, color = "#4E79A7") +
        ggplot2::labs(title = if (isTRUE(cfg$comparison$paired)) "DTU sensitivity to patient exclusion" else "DTU sensitivity to sample exclusion", subtitle = paste0("Spearman correlation = ", if (is.finite(rho)) round(rho, 3) else "unavailable"), x = "Main cohort -log10(P)", y = "Exclusion cohort -log10(P)") + rna_theme(cfg)
      rna_save_plot(p, paste0("DTU_sensitivity_", cohort), dirs, cfg, width = 6, height = 5)
    }
  }
  tab <- dplyr::bind_rows(summaries)
  rna_write_tsv(tab, file.path(dirs$tables, "DTU_sensitivity_summary.tsv"))
  list(status = "completed", summary = tab, comparisons = comparisons)
}

# Independent biological samples: fractions are summarized within condition.
# There is intentionally no sample matching, paired delta, or patient-support call.
dtu_unpaired_effect_tables <- function(prepared, fit, stages, alphas) {
  if (isTRUE(prepared$empty)) return(list(feature = prepared$feature_map,
    paired = tibble::tibble(), usage = tibble::tibble(),
    driver = tibble::tibble(gene_id = character())))
  sample_ids <- as.character(prepared$samples$sample_id)
  usage <- dplyr::left_join(prepared$counts, prepared$feature_map, by = c("gene_id", "feature_id"))
  usage <- tidyr::pivot_longer(usage, dplyr::all_of(sample_ids), names_to = "sample_id", values_to = "FL")
  samples <- as.data.frame(prepared$samples)
  samples$patient <- as.character(samples$patient); samples$timepoint <- as.character(samples$timepoint)
  usage <- dplyr::left_join(usage, samples, by = "sample_id")
  usage <- dplyr::group_by(usage, gene_id, sample_id)
  usage <- dplyr::mutate(usage, gene_FL = sum(FL),
    isoform_fraction = ifelse(gene_FL > 0, FL / gene_FL, NA_real_))
  usage <- dplyr::ungroup(usage)
  usage$analysis_design <- "unpaired"
  average <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  feature <- dplyr::summarise(dplyr::group_by(usage, gene_id, feature_id,
      isoform_harmonized_id, structural_category),
    mean_primary_fraction = average(isoform_fraction[timepoint == "primary"]),
    mean_recurrent_fraction = average(isoform_fraction[timepoint == "recurrent"]),
    evaluable_reference_samples = sum(timepoint == "primary" & is.finite(isoform_fraction)),
    evaluable_comparison_samples = sum(timepoint == "recurrent" & is.finite(isoform_fraction)),
    .groups = "drop")
  feature$mean_delta_if <- feature$mean_recurrent_fraction - feature$mean_primary_fraction
  feature$analysis_design <- "unpaired"
  feature <- dplyr::left_join(feature, fit$feature, by = c("gene_id", "feature_id"))
  for (alpha in alphas) {
    tag <- dtu_alpha_tag(alpha)
    current <- stages[[tag]]$feature
    if (!is.null(current) && nrow(current)) feature <- dplyr::left_join(feature, current, by = c("gene_id", "feature_id"))
    column <- paste0("stageR_isoform_adjusted_", tag)
    if (!column %in% names(feature)) feature[[column]] <- NA_real_
    feature[[paste0("stageR_confirmed_", tag)]] <- is.finite(feature[[column]]) & feature[[column]] <= alpha
  }
  ordered <- dplyr::arrange(dplyr::group_by(feature, gene_id), dplyr::desc(abs(mean_delta_if)), feature_id, .by_group = TRUE)
  driver <- dplyr::ungroup(dplyr::slice_head(ordered, n = 1L))
  driver <- dplyr::transmute(driver, gene_id, DTU_driver_isoform = isoform_harmonized_id,
    DTU_driver_structure = structural_category, DTU_driver_mean_delta_if = mean_delta_if,
    DTU_driver_abs_mean_delta_if = abs(mean_delta_if))
  list(feature = feature, paired = tibble::tibble(), usage = usage, driver = driver)
}

dtu_unpaired_plot_outputs <- function(result, cfg, dirs, alpha) {
  if (!isTRUE(cfg$reporting$enabled %||% TRUE) || !nrow(result$gene)) return(invisible(NULL))
  d <- as.data.frame(result$gene)
  d$DTU_significant <- is.finite(d$dtu_fdr) & d$dtu_fdr < alpha
  raw_p <- if ("raw_dtu_p_value" %in% names(d)) d$raw_dtu_p_value else d$dtu_p_value
  evaluable <- is.finite(raw_p) & is.finite(d$DTU_driver_abs_mean_delta_if)
  unavailable <- sum(!evaluable)
  d$log10_P <- -log10(pmax(raw_p, .Machine$double.xmin))
  d <- d[evaluable, , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  p <- ggplot2::ggplot(d, ggplot2::aes(DTU_driver_abs_mean_delta_if, log10_P, colour = DTU_significant)) +
    ggplot2::geom_point(size = 2, alpha = .7) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = "#B2182B"),
      labels = c("Other genes", paste0("FDR < ", alpha))) +
    ggplot2::labs(title = "Independent-group differential isoform usage",
      subtitle = paste0(cfg$comparison$comparison, " versus ", cfg$comparison$reference,
        "; ", nrow(d), " finite tests, ", unavailable, " unavailable"),
      x = "Absolute difference in mean isoform fraction", y = expression(-log[10](P)), colour = NULL) + rna_theme(cfg)
  rna_save_plot(p, paste0("DTU_independent_effect_evidence_", dtu_fdr_tag(alpha)), dirs, cfg, 7, 5)
  invisible(NULL)
}

dtu_unpaired_gene_examples <- function(result, cfg, dirs) {
  if (!isTRUE(cfg$reporting$enabled %||% TRUE) || !nrow(result$gene)) return(invisible(NULL))
  genes <- as.data.frame(result$gene)
  genes <- genes[is.finite(genes$dtu_fdr) & genes$dtu_fdr < max(as.numeric(cfg$multiplicity$fdr_thresholds)), , drop = FALSE]
  genes <- genes[order(genes$dtu_fdr, genes$gene_id), , drop = FALSE]
  maximum <- cfg$dtu$maximum_gene_plots %||% 12L
  if (!isTRUE(cfg$dtu$plot_all_dtu_genes)) genes <- utils::head(genes, if (maximum > 0) maximum else 12L)
  else if (maximum > 0) genes <- utils::head(genes, maximum)
  for (id in genes$gene_id) {
    feature <- as.data.frame(result$effects$feature)
    feature <- feature[feature$gene_id == id, , drop = FALSE]
    feature <- utils::head(feature[order(-abs(feature$mean_delta_if), feature$feature_id), , drop = FALSE],
      cfg$dtu$isoforms_per_gene_example %||% 2L)
    d <- as.data.frame(result$effects$usage)
    d <- d[d$gene_id == id & d$feature_id %in% feature$feature_id & is.finite(d$isoform_fraction), , drop = FALSE]
    if (!nrow(d)) next
    d$condition <- factor(d$timepoint, levels = c("primary", "recurrent"),
      labels = c(cfg$comparison$reference, cfg$comparison$comparison))
    d$isoform <- factor(d$feature_id, levels = feature$feature_id,
      labels = paste(seq_len(nrow(feature)), feature$structural_category, sep = ": "))
    p <- ggplot2::ggplot(d, ggplot2::aes(condition, isoform_fraction, colour = condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = .55) +
      ggplot2::geom_point(ggplot2::aes(group = sample_id),
        position = ggplot2::position_dodge(width = .35), size = 2.5) +
      ggplot2::facet_wrap(~isoform, nrow = 1L) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::labs(title = paste0(id, ": independent sample isoform usage"),
        subtitle = paste0("FDR = ", format(genes$dtu_fdr[genes$gene_id == id], digits = 3)),
        x = NULL, y = "Isoform fraction within gene", colour = "Condition") + rna_theme(cfg)
    rna_save_plot(p, paste0("DTU_independent_isoform_usage_", id), dirs, cfg, max(7, 3.5 * nrow(feature)), 5)
  }
  invisible(NULL)
}
