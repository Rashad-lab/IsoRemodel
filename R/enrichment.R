# DTU enrichment uses the exact publication hypergeometric and deterministic
# matching functions. Membership and identifier mapping are captured once and
# saved alongside every result, so reruns cannot silently use a new database.

dtu_enrichment_resources <- function(results, cfg, dirs, env) {
  if (!is.null(cfg$references$species) && cfg$references$species != "Homo sapiens" &&
      (is.null(cfg$references$gene_id_map) || is.null(cfg$references$pathway_term2gene)))
    stop("Nonhuman enrichment requires explicit species-matched gene_id_map and pathway_term2gene files.", call. = FALSE)
  reference_dir <- cfg$enrichment$reference_cache_directory %||% file.path(dirs$root, "references")
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  read_table <- function(path) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  term_input <- cfg$references$pathway_term2gene
  snapshot <- file.path(reference_dir, "dtu_term2gene_frozen.tsv.gz")
  metadata_path <- file.path(reference_dir, "dtu_reference_sources_frozen.rds")
  reference_metadata <- if (file.exists(metadata_path)) readRDS(metadata_path) else list()
  if (!is.null(term_input) && nzchar(term_input)) {
    terms <- read_table(term_input); term_source <- normalizePath(term_input, winslash = "/")
    reference_metadata$pathway_source <- term_source
    reference_metadata$pathway_database_version <- NA_character_
  } else if (file.exists(snapshot)) {
    terms <- read_table(snapshot); term_source <- "previously_frozen_snapshot"
  } else {
    rna_require("msigdbr")
    captured <- env$load_module06_pathways()
    terms <- captured$term2gene
    term_source <- paste0("msigdbr ", utils::packageVersion("msigdbr"), "; first capture")
    reference_metadata$pathway_source <- term_source
    reference_metadata$pathway_database_version <- captured$database_version
  }
  if (!all(c("collection", "pathway", "gene") %in% names(terms))) stop("Pathway membership requires collection,pathway,gene columns.")
  aliases <- c(GO_BP = "GO biological process", HALLMARK = "MSigDB Hallmark", KEGG = "KEGG")
  terms$collection <- ifelse(terms$collection %in% names(aliases), unname(aliases[terms$collection]), terms$collection)
  terms <- unique(terms[, c("collection", "pathway", "gene")])
  terms <- terms[!is.na(terms$gene) & !is.na(terms$pathway) & nzchar(terms$gene) & nzchar(terms$pathway), , drop = FALSE]
  missing <- setdiff(env$cfg$ora$collections, unique(terms$collection))
  if (length(missing)) stop("Frozen term2gene lacks requested collections: ", paste(missing, collapse = ", "))
  rna_write_tsv(terms, snapshot)

  genes <- sort(unique(unlist(lapply(results, function(x) x$gene$gene_id))))
  map_input <- cfg$references$gene_id_map
  map_snapshot <- file.path(reference_dir, "dtu_gene_symbol_mapping_frozen.tsv.gz")
  if (!is.null(map_input) && nzchar(map_input)) {
    mapping <- read_table(map_input); map_source <- normalizePath(map_input, winslash = "/")
    reference_metadata$mapping_package <- map_source
    if (!"SYMBOL" %in% names(mapping) && "symbol" %in% names(mapping)) names(mapping)[names(mapping) == "symbol"] <- "SYMBOL"
    if (!"mapping_type" %in% names(mapping)) mapping$mapping_type <- "user_supplied"
  } else {
    mapping <- if (file.exists(map_snapshot)) read_table(map_snapshot) else tibble::tibble(gene_id = character(), SYMBOL = character(), mapping_type = character())
    missing_genes <- setdiff(genes, mapping$gene_id)
    # Persist unavailable mappings too. They are unavailable, not negative.
    if (length(missing_genes)) {
      rna_require(c("AnnotationDbi", "org.Hs.eg.db"))
      added <- env$map_gene_symbols(missing_genes)
      added <- dplyr::left_join(tibble::tibble(gene_id = missing_genes), added, by = "gene_id")
      mapping <- dplyr::bind_rows(mapping, added)
      reference_metadata$mapping_package <- paste0("org.Hs.eg.db ", utils::packageVersion("org.Hs.eg.db"))
    }
    map_source <- "frozen symbol map; additional IDs mapped with org.Hs.eg.db"
  }
  if (!all(c("gene_id", "SYMBOL") %in% names(mapping))) stop("Gene mapping requires gene_id,SYMBOL columns.")
  if (anyDuplicated(mapping$gene_id)) stop("Gene mapping contains duplicate gene_id keys.")
  rna_write_tsv(mapping, map_snapshot)
  rna_atomic_rds(reference_metadata, metadata_path)
  rna_write_tsv(terms, file.path(dirs$tables, "pathway_memberships_used.tsv.gz"))
  rna_write_tsv(mapping, file.path(dirs$tables, "gene_symbol_mapping_used.tsv.gz"))
  provenance <- tibble::tibble(resource = c("term2gene", "gene_symbol_mapping"), source = c(term_source, map_source),
      snapshot = c(snapshot, map_snapshot), content_sha256 = c(rna_hash(terms), rna_hash(mapping)),
      original_source = c(reference_metadata$pathway_source %||% term_source, reference_metadata$mapping_package %||% map_source),
      database_version = c(reference_metadata$pathway_database_version %||% NA_character_, NA_character_), captured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  rna_write_tsv(provenance, file.path(dirs$tables, "reference_provenance.tsv"))
  # Freeze lookup inside the matching functions, too, rather than remapping on
  # each foreground and potentially changing identifier collisions.
  env$map_gene_symbols <- local({ saved <- mapping; function(genes) saved[saved$gene_id %in% genes & !is.na(saved$SYMBOL) & nzchar(saved$SYMBOL), , drop = FALSE] })
  list(term2gene = tibble::as_tibble(terms), mapping = mapping, provenance = provenance)
}

run_dtu_enrichment <- function(results, cfg, stage_dir) {
  dirs <- rna_dirs(file.path(stage_dir, "enrichment"))
  e <- dtu_engine_environment(parent = environment())
  sys.source(dtu_engine_path(cfg, "dtu_core.R"), envir = e)
  sys.source(dtu_engine_path(cfg, "dtu_enrichment_core.R"), envir = e)
  alias <- c(GO_BP = "GO biological process", HALLMARK = "MSigDB Hallmark", KEGG = "KEGG")
  collections <- as.character(unlist(cfg$enrichment$collections))
  collections <- ifelse(collections %in% names(alias), unname(alias[collections]), collections)
  e$cfg <- list(ora = list(collections = collections, minimum_genes = cfg$enrichment$minimum_foreground_genes,
      minimum_term_genes = cfg$enrichment$minimum_term_genes, maximum_term_genes = cfg$enrichment$maximum_term_genes,
      reporting_fdr = cfg$enrichment$reporting_fdr, nominal_fallback = cfg$enrichment$nominal_display_fallback_p,
      top_terms = cfg$enrichment$maximum_display_terms))
  e$publication_layout <- function(...) rna_theme(cfg)
  resources <- dtu_enrichment_resources(results, cfg, dirs, e)
  alphas <- sort(unique(as.numeric(cfg$multiplicity$fdr_thresholds)))
  standard <- matched <- patient_results <- balance <- integrity <- status <- list()
  for (cohort in names(results)) {
    result <- results[[cohort]]
    # Failed/boundary tests are retained in correction families but cannot enter
    # the eligible enrichment background: that background is finite-tested.
    family <- result$gene[result$gene$fit_status == "finite_test", , drop = FALSE]
    mapping <- e$map_gene_symbols(family$gene_id)
    symbols <- unique(mapping$SYMBOL)
    cd <- rna_dirs(file.path(stage_dir, cohort, "enrichment"))
    rna_write_tsv(mapping, file.path(cd$tables, "testing_family_gene_symbol_mapping.tsv.gz"))
    for (alpha in alphas) {
      tag <- dtu_fdr_tag(alpha)
      target <- result$gene$gene_id[is.finite(result$gene$dtu_fdr) & result$gene$dtu_fdr < alpha]
      target_symbols <- unique(mapping$SYMBOL[mapping$gene_id %in% target])
      design <- if (isTRUE(cfg$enrichment$matched_background)) e$make_matched_design(family, target, cfg$enrichment$matched_controls_per_target) else
        list(available = FALSE, reason = "disabled", pairs = tibble::tibble(), controls = tibble::tibble(), universe = tibble::tibble(), integrity = tibble::tibble(), balance = tibble::tibble())
      key <- paste(cohort, tag, sep = "__")
      integrity[[key]] <- dplyr::mutate(design$integrity, cohort = cohort, threshold = tag, available = isTRUE(design$available), reason = design$reason)
      balance[[key]] <- dplyr::mutate(design$balance, cohort = cohort, threshold = tag)
      rna_write_tsv(design$pairs, file.path(cd$tables, paste0("matched_background_pairs_", tag, ".tsv.gz")))
      rna_write_tsv(design$controls, file.path(cd$tables, paste0("matched_background_controls_", tag, ".tsv.gz")))
      rna_write_tsv(design$universe, file.path(cd$tables, paste0("matched_background_covariate_universe_", tag, ".tsv.gz")))
      for (collection in collections) {
        k <- paste(key, e$safe_name(collection), sep = "__")
        standard[[k]] <- dplyr::mutate(e$standard_ora(target_symbols, symbols, resources$term2gene, collection), cohort = cohort, threshold = tag)
        matched[[k]] <- if (isTRUE(design$available)) dplyr::mutate(e$matched_pathway_test(design, resources$term2gene, collection), cohort = cohort, threshold = tag) else tibble::tibble()
        for (method in c("standard", "matched")) {
          tab <- if (method == "standard") standard[[k]] else matched[[k]]
          stem <- paste0(method, "_ORA_", tag, "_", e$safe_name(collection))
          rna_write_tsv(tab, file.path(cd$tables, paste0(stem, ".tsv.gz")))
          if (nrow(tab) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
            p <- e$plot_pathway_dotplot(tab, paste0("DTU pathway enrichment: ", collection), paste0("FDR < ", alpha, "; ", method, " background"), matched = method == "matched")
            if (!is.null(p)) rna_save_plot(p, stem, cd, cfg, width = 8.5, height = max(5, 2.8 + .48 * nrow(e$select_pathway_plot_data(tab, method == "matched"))))
          }
        }
        status[[k]] <- tibble::tibble(cohort = cohort, alpha = alpha, collection = collection, scheduled_genes = nrow(result$gene), finite_test_genes = nrow(family),
          mapped_background_genes = length(unique(mapping$gene_id)), mapped_background_symbols = length(symbols), mapping_rate = if (nrow(family)) length(unique(mapping$gene_id)) / nrow(family) else NA_real_,
          symbol_collisions = length(unique(mapping$gene_id)) - length(symbols), target_genes = length(target), target_symbols = length(target_symbols),
          standard_terms_tested = nrow(standard[[k]]), matched_available = isTRUE(design$available), matched_reason = design$reason, matched_terms_tested = nrow(matched[[k]]))
        for (patient in sort(unique(result$patient$patient))) {
          support <- result$patient
          chosen <- support$gene_id[support$patient == patient & support$gene_id %in% target & support$patient_direction_supported %in% TRUE]
          chosen_symbols <- unique(mapping$SYMBOL[mapping$gene_id %in% chosen])
          pk <- paste(k, patient, sep = "__")
          patient_results[[pk]] <- dplyr::mutate(e$standard_ora(chosen_symbols, symbols, resources$term2gene, collection), cohort = cohort, threshold = tag, patient = patient)
          tab <- patient_results[[pk]]
          pd <- rna_dirs(file.path(stage_dir, cohort, "enrichment", paste0("patient_", e$safe_name(patient))))
          stem <- paste0("patient_supported_ORA_", tag, "_", e$safe_name(collection))
          rna_write_tsv(tab, file.path(pd$tables, paste0(stem, ".tsv.gz")))
          if (nrow(tab) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
            p <- e$plot_pathway_dotplot(tab, "Patient-supported DTU pathways", paste0(patient, "; FDR < ", alpha), matched = FALSE)
            if (!is.null(p)) rna_save_plot(p, stem, pd, cfg, width = 8.5, height = max(5, 2.8 + .48 * nrow(e$select_pathway_plot_data(tab, FALSE))))
          }
        }
      }
    }
  }
  out <- list(status = "completed", standard = dplyr::bind_rows(standard), matched = dplyr::bind_rows(matched),
       patient = dplyr::bind_rows(patient_results), balance = dplyr::bind_rows(balance), integrity = dplyr::bind_rows(integrity), runs = dplyr::bind_rows(status), resources = resources)
  for (name in c("standard", "matched", "patient", "balance", "integrity", "runs")) rna_write_tsv(out[[name]], file.path(dirs$tables, paste0("DTU_enrichment_", name, ".tsv.gz")))
  if (nrow(out$balance) && isTRUE(cfg$reporting$enabled %||% TRUE)) {
    for (cohort in unique(out$balance$cohort)) for (alpha in alphas) {
      x <- out$balance[out$balance$cohort == cohort & out$balance$threshold == dtu_fdr_tag(alpha), , drop = FALSE]
      if (!nrow(x)) next
      x$covariate <- gsub("_", " ", sub("^log_", "", x$covariate))
      p <- ggplot2::ggplot(x, ggplot2::aes(x = standardized_mean_difference, y = covariate)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") + ggplot2::geom_point(size = 3, color = "#4E79A7") + rna_theme(cfg) +
        ggplot2::labs(title = "Matched-background balance", subtitle = paste0("FDR < ", alpha), x = "Standardized mean difference", y = NULL)
      rna_save_plot(p, paste0("matched_background_balance_", cohort, "_", dtu_fdr_tag(alpha)), dirs, cfg, width = 7, height = 4.8)
    }
  }
  rna_atomic_rds(out, file.path(dirs$objects, "DTU_enrichment.rds"))
  out
}
