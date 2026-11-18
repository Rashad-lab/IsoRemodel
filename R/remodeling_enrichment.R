# ORF-free functional follow-up of the final remodeling classes. These sets,
# backgrounds, GO enrichGO and Hallmark enricher calls come from Module 05,
# which received the locked Module 04 bootstrap/direction/class definitions.

remodeling_enrichment_sets <- function(genes, alphas = c(.05, .10), nominal_p = .05) {
  n <- nrow(genes)
  logical_column <- function(name) if (name %in% names(genes)) genes[[name]] %in% TRUE else rep(FALSE, n)
  numeric_column <- function(name) if (name %in% names(genes)) as.numeric(genes[[name]]) else rep(NA_real_, n)
  ids <- as.character(genes$gene_id)
  quantifiable <- logical_column("quantifiable")
  family <- logical_column("included_in_testing_family")
  authoritative <- if ("downstream_primary_remodeling_set" %in% names(genes)) logical_column("downstream_primary_remodeling_set") else logical_column("is_high_confidence_remodeled")
  ae <- as.character(genes$architecture_class) %in% LETTERS[1:5]
  p <- numeric_column("gene_model_p"); q <- numeric_column("gene_model_q")
  nominal <- family & is.finite(p) & p < nominal_p
  make <- function(flag, background, label) list(target = unique(ids[flag %in% TRUE]), background = unique(ids[background %in% TRUE]), label = label)
  out <- list(
    authoritative_biological = make(authoritative, quantifiable, "Stable, patient-supported remodeling"),
    strict_architecture = make(logical_column("is_high_confidence_class"), quantifiable, "Stable remodeling classes"),
    nominal_monte_carlo_all = make(nominal, family, paste0("Cohort test P < ", nominal_p)),
    nominal_monte_carlo_AE = make(nominal & ae, family, paste0("Classes A-E; cohort test P < ", nominal_p)),
    concordant_authoritative_and_nominal = make(authoritative & nominal & ae, family, "Stable remodeling with cohort-test support"))
  for (alpha in alphas) out[[paste0("cohort_", dtu_fdr_tag(alpha))]] <- make(family & is.finite(q) & q <= alpha, family, paste0("Cohort remodeling: FDR < ", alpha))
  for (class in LETTERS[1:5]) out[[paste0("stable_class_", class)]] <- make(authoritative & as.character(genes$architecture_class) == class, quantifiable, paste0("Stable remodeling: class ", class))
  out
}

remodeling_enrichment_resources <- function(genes, cfg, dirs) {
  collections <- as.character(unlist(cfg$enrichment$remodeling_collections %||% cfg$enrichment$collections))
  alias <- c(GO_BP = "GO biological process", HALLMARK = "MSigDB Hallmark", KEGG = "KEGG")
  collections <- ifelse(collections %in% names(alias), unname(alias[collections]), collections)
  go_enabled <- "GO biological process" %in% collections
  backend <- cfg$enrichment$remodeling_go_backend %||% "publication_AnnotationDbi_GO"
  if (go_enabled && backend != "publication_AnnotationDbi_GO") stop("Unsupported remodeling GO backend: ", backend, ". This release retains publication_AnnotationDbi_GO.", call. = FALSE)
  if (!is.null(cfg$references$species) && cfg$references$species != "Homo sapiens") stop("The publication remodeling GO/mapping backend targets Homo sapiens.", call. = FALSE)
  rna_require(c("AnnotationDbi", "org.Hs.eg.db", "clusterProfiler"))
  cache <- cfg$enrichment$reference_cache_directory %||% file.path(dirs$root, "references")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  sqlite <- file.path(cache, "remodeling_org.Hs.eg.db_frozen.sqlite")
  metadata_path <- file.path(cache, "remodeling_reference_metadata_frozen.rds")
  metadata <- if (file.exists(metadata_path)) readRDS(metadata_path) else list()
  if (!file.exists(sqlite)) {
    if (!file.copy(AnnotationDbi::dbfile(org.Hs.eg.db::org.Hs.eg.db), sqlite, overwrite = FALSE)) stop("Could not capture the OrgDb annotation snapshot.")
    metadata$orgdb_package <- paste0("org.Hs.eg.db ", utils::packageVersion("org.Hs.eg.db"))
    metadata$orgdb_captured_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  }
  orgdb <- AnnotationDbi::loadDb(sqlite)
  if (go_enabled) {
    rna_require("GO.db")
    go_hash <- digest::digest(file = AnnotationDbi::dbfile(GO.db::GO.db), algo = "sha256")
    if (!is.null(metadata$GO_db_sha256) && !identical(metadata$GO_db_sha256, go_hash)) stop("GO.db changed after the remodeling reference capture. Restore the recorded GO.db version or deliberately choose a new reference_cache_directory.", call. = FALSE)
    metadata$GO_db_sha256 <- go_hash
    metadata$GO_db_package <- paste0("GO.db ", utils::packageVersion("GO.db"))
    # Save the actual propagated BP membership used by enrichGO for inspection.
    go_memberships <- file.path(cache, "remodeling_GO_BP_membership_frozen.tsv.gz")
    if (!file.exists(go_memberships)) {
      x <- suppressMessages(AnnotationDbi::select(orgdb, keys = AnnotationDbi::keys(orgdb, keytype = "ENTREZID"), keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
      x <- unique(x[!is.na(x$GOALL) & x$ONTOLOGYALL == "BP", c("GOALL", "ENTREZID")])
      names(x) <- c("pathway", "gene")
      rna_write_tsv(x, go_memberships)
    }
  }
  map_path <- file.path(cache, "remodeling_identifier_map_frozen.tsv.gz")
  supplied <- cfg$references$gene_id_map
  if (!is.null(supplied) && nzchar(supplied)) {
    mapping <- readr::read_tsv(supplied, show_col_types = FALSE, progress = FALSE)
    if ("gene_id" %in% names(mapping)) names(mapping)[names(mapping) == "gene_id"] <- "input_gene"
    if ("symbol" %in% names(mapping) && !"SYMBOL" %in% names(mapping)) names(mapping)[names(mapping) == "symbol"] <- "SYMBOL"
    if (!all(c("input_gene", "SYMBOL") %in% names(mapping))) stop("Remodeling gene_id_map requires gene_id/input_gene and SYMBOL/symbol.")
    if (!"ENTREZID" %in% names(mapping)) {
      added <- suppressMessages(AnnotationDbi::select(orgdb, keys = unique(stats::na.omit(mapping$SYMBOL)), keytype = "SYMBOL", columns = "ENTREZID"))
      mapping <- dplyr::left_join(mapping, added, by = "SYMBOL")
    }
    metadata$mapping_source <- normalizePath(supplied, winslash = "/")
  } else {
    mapping <- if (file.exists(map_path)) readr::read_tsv(map_path, show_col_types = FALSE, progress = FALSE) else tibble::tibble(input_gene = character(), SYMBOL = character(), ENTREZID = character())
    missing <- setdiff(unique(as.character(genes)), mapping$input_gene)
    if (length(missing)) {
      # Preserve original M05 mapping: exact symbols and leading ENSG IDs;
      # unlike the DTU mapper, compound-symbol prefixes are not guessed.
      symbols <- intersect(missing, AnnotationDbi::keys(orgdb, keytype = "SYMBOL"))
      direct <- if (length(symbols)) suppressMessages(AnnotationDbi::select(orgdb, keys = symbols, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))) else data.frame(SYMBOL = character(), ENTREZID = character())
      direct$input_gene <- direct$SYMBOL
      ens <- missing[grepl("^ENSG", missing)]
      ens_input <- data.frame(input_gene = ens, ENSEMBL = sub("\\.[0-9]+$", "", ens))
      ensembl <- if (length(ens)) suppressMessages(AnnotationDbi::select(orgdb, keys = unique(ens_input$ENSEMBL), keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID"))) else data.frame(ENSEMBL = character(), SYMBOL = character(), ENTREZID = character())
      ensembl <- merge(ens_input, ensembl, by = "ENSEMBL", all = FALSE)
      added <- dplyr::bind_rows(direct[, c("input_gene", "SYMBOL", "ENTREZID")], ensembl[, c("input_gene", "SYMBOL", "ENTREZID")])
      added <- dplyr::left_join(tibble::tibble(input_gene = missing), tibble::as_tibble(added), by = "input_gene")
      mapping <- dplyr::bind_rows(mapping, added)
    }
    metadata$mapping_source <- "Publication M05 symbol/ENSEMBL mapping against frozen OrgDb"
  }
  mapping$ENTREZID <- as.character(mapping$ENTREZID)
  mapping <- unique(mapping[, c("input_gene", "SYMBOL", "ENTREZID")])
  rna_write_tsv(mapping, map_path)
  terms <- tibble::tibble(collection = character(), pathway = character(), gene = character())
  if (any(collections %in% c("MSigDB Hallmark", "KEGG"))) {
    term_path <- file.path(cache, "remodeling_msigdb_term2gene_frozen.tsv.gz")
    supplied_terms <- cfg$references$pathway_term2gene
    if (!is.null(supplied_terms) && nzchar(supplied_terms)) {
      terms <- readr::read_tsv(supplied_terms, show_col_types = FALSE, progress = FALSE)
      metadata$msigdb_source <- normalizePath(supplied_terms, winslash = "/")
    } else if (file.exists(term_path)) terms <- readr::read_tsv(term_path, show_col_types = FALSE, progress = FALSE) else {
      rna_require("msigdbr")
      e <- dtu_engine_environment(); sys.source(dtu_engine_path(cfg, "dtu_core.R"), e); sys.source(dtu_engine_path(cfg, "dtu_enrichment_core.R"), e)
      e$cfg <- list(ora = list(collections = setdiff(collections, "GO biological process")))
      captured <- e$load_module06_pathways(); terms <- captured$term2gene
      metadata$msigdb_source <- paste0("msigdbr ", utils::packageVersion("msigdbr"))
      metadata$msigdb_database_version <- captured$database_version
    }
    if (!all(c("collection", "pathway", "gene") %in% names(terms))) stop("Remodeling term2gene requires collection,pathway,gene.")
    terms$collection <- ifelse(terms$collection %in% names(alias), unname(alias[terms$collection]), terms$collection)
    terms <- unique(terms[terms$collection %in% c("MSigDB Hallmark", "KEGG"), c("collection", "pathway", "gene")])
    if (length(setdiff(setdiff(collections, "GO biological process"), unique(terms$collection)))) stop("Missing requested Hallmark/KEGG collections in remodeling term2gene.")
    rna_write_tsv(terms, term_path)
  }
  metadata$orgdb_sha256 <- digest::digest(file = sqlite, algo = "sha256")
  rna_atomic_rds(metadata, metadata_path)
  rna_write_tsv(mapping, file.path(dirs$tables, "gene_identifier_mapping_used.tsv.gz"))
  rna_write_tsv(terms, file.path(dirs$tables, "MSigDB_membership_used.tsv.gz"))
  rna_write_tsv(tibble::tibble(resource = names(metadata), value = vapply(metadata, as.character, character(1))), file.path(dirs$tables, "reference_provenance.tsv"))
  list(orgdb = orgdb, mapping = mapping, terms = terms, collections = collections, metadata = metadata)
}

remodeling_ora_set <- function(target, background, resources, cfg) {
  background <- unique(as.character(background)); background <- background[!is.na(background) & nzchar(background)]
  target <- intersect(unique(as.character(target)), background)
  mapping <- resources$mapping
  mapping <- mapping[mapping$input_gene %in% background & !is.na(mapping$SYMBOL) & !is.na(mapping$ENTREZID), , drop = FALSE]
  foreground <- mapping[mapping$input_gene %in% target, , drop = FALSE]
  status <- list(input_target_genes = length(target), input_background_genes = length(background), mapped_target_symbols = length(unique(foreground$SYMBOL)),
    mapped_target_entrez = length(unique(foreground$ENTREZID)), mapped_background_symbols = length(unique(mapping$SYMBOL)), mapped_background_entrez = length(unique(mapping$ENTREZID)))
  if (status$mapped_target_entrez < cfg$enrichment$minimum_foreground_genes) return(list(result = tibble::tibble(), mapping = mapping, status = c(status, list(status = "skipped_too_few_mapped_target_genes"))))
  results <- list(); errors <- character()
  for (collection in resources$collections) {
    fit <- tryCatch({
      if (collection == "GO biological process") suppressMessages(clusterProfiler::enrichGO(gene = unique(foreground$ENTREZID), universe = unique(mapping$ENTREZID),
        OrgDb = resources$orgdb, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE,
        minGSSize = cfg$enrichment$remodeling_go_minimum_term_genes %||% 10L, maxGSSize = cfg$enrichment$remodeling_go_maximum_term_genes %||% 500L)) else {
        term2gene <- resources$terms[resources$terms$collection == collection, c("pathway", "gene")]
        suppressMessages(clusterProfiler::enricher(gene = unique(foreground$SYMBOL), universe = unique(mapping$SYMBOL), TERM2GENE = term2gene,
          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = cfg$enrichment$minimum_term_genes, maxGSSize = cfg$enrichment$maximum_term_genes))
      }
    }, error = function(e) e)
    if (inherits(fit, "error")) { errors <- c(errors, paste(collection, conditionMessage(fit), sep = ": ")); next }
    x <- if (is.null(fit)) data.frame() else as.data.frame(fit)
    if (!nrow(x)) next
    x$collection <- collection; x$FDR <- as.numeric(x$p.adjust)
    ratio <- strsplit(as.character(x$GeneRatio), "/", fixed = TRUE)
    x$enrichment_fraction <- vapply(ratio, function(z) as.numeric(z[1L]) / as.numeric(z[2L]), numeric(1))
    results[[collection]] <- tibble::as_tibble(x)
  }
  list(result = dplyr::bind_rows(results), mapping = mapping,
    status = c(status, list(status = if (length(errors)) "failed_collection" else "completed", errors = paste(errors, collapse = "; "))))
}

remodeling_ora_plot <- function(tab, label, cfg) {
  if (!nrow(tab)) return(NULL)
  chosen <- tab[is.finite(tab$FDR) & tab$FDR < cfg$enrichment$reporting_fdr, , drop = FALSE]
  caption <- paste0("Pathway FDR < ", cfg$enrichment$reporting_fdr)
  if (!nrow(chosen)) {
    chosen <- tab[is.finite(tab$pvalue) & tab$pvalue < cfg$enrichment$nominal_display_fallback_p, , drop = FALSE]
    caption <- paste0("Nominal P < ", cfg$enrichment$nominal_display_fallback_p, "; no pathways pass FDR < ", cfg$enrichment$reporting_fdr)
  }
  if (!nrow(chosen)) return(NULL)
  chosen <- head(chosen[order(chosen$FDR, chosen$pvalue, -chosen$Count), , drop = FALSE], cfg$enrichment$maximum_display_terms)
  chosen$label <- stringr::str_wrap(gsub("_", " ", sub("^(HALLMARK_|KEGG_)", "", chosen$Description)), 38)
  chosen$score <- -log10(pmax(chosen$pvalue, .Machine$double.xmin))
  p <- ggplot2::ggplot(chosen, ggplot2::aes(enrichment_fraction, reorder(label, enrichment_fraction), size = Count, color = score)) + ggplot2::geom_point(alpha = .95) +
    ggplot2::scale_x_continuous(labels = scales::percent) + ggplot2::scale_color_gradient(low = "#F4A582", high = "#B2182B") + rna_theme(cfg) +
    ggplot2::labs(title = label, subtitle = unique(chosen$collection), caption = caption, x = "Selected genes in pathway", y = NULL, size = "Genes", color = expression(-log[10](P)))
  attr(p, "selection_table") <- chosen
  p
}

run_remodeling_enrichment <- function(inputs, remodeling, cfg, stage_dir) {
  if (!isTRUE(cfg$enrichment$enabled)) return(list(status = "disabled", complete = TRUE))
  rna_require(c("dplyr", "readr", "tibble", "ggplot2", "stringr", "scales", "digest"))
  dirs <- rna_dirs(file.path(stage_dir, "enrichment"))
  genes <- remodeling$gene_table
  if (is.null(genes) || !nrow(genes)) return(list(status = "empty_gene_universe", complete = TRUE))
  if (!all(c("cohort", "gene_id", "architecture_class") %in% names(genes))) stop("Remodeling enrichment requires cohort,gene_id,architecture_class.")
  alphas <- sort(unique(as.numeric(cfg$multiplicity$fdr_thresholds)))
  resources <- remodeling_enrichment_resources(genes$gene_id, cfg, dirs)
  results <- statuses <- memberships <- selections <- list()
  execute <- function(spec, scope, cohort, patient, name) {
    key <- paste(scope, cohort, patient, name, sep = "__")
    folder <- file.path(dirs$root, gsub("[^A-Za-z0-9._-]", "_", key)); dd <- rna_dirs(folder)
    fit <- remodeling_ora_set(spec$target, spec$background, resources, cfg)
    result <- fit$result
    if (nrow(result)) { result$analysis_scope <- scope; result$cohort <- cohort; result$patient <- patient; result$gene_set <- name }
    results[[key]] <<- result
    statuses[[key]] <<- tibble::as_tibble(c(list(analysis_scope = scope, cohort = cohort, patient = patient, gene_set = name), fit$status))
    membership <- data.frame(gene_id = unique(spec$background), in_foreground = unique(spec$background) %in% spec$target, stringsAsFactors = FALSE)
    membership$analysis_scope <- rep(scope, nrow(membership)); membership$cohort <- rep(cohort, nrow(membership)); membership$patient <- rep(patient, nrow(membership)); membership$gene_set <- rep(name, nrow(membership))
    memberships[[key]] <<- membership
    rna_write_tsv(result, file.path(dd$tables, "ORA_all_tested_terms.tsv.gz"))
    rna_write_tsv(membership, file.path(dd$tables, "foreground_and_eligible_background.tsv.gz"))
    rna_write_tsv(fit$mapping, file.path(dd$tables, "identifier_mapping.tsv.gz"))
    if (nrow(result) && isTRUE(cfg$reporting$enabled %||% TRUE)) for (collection in unique(result$collection)) {
      p <- remodeling_ora_plot(result[result$collection == collection, , drop = FALSE], spec$label, cfg)
      if (is.null(p)) next
      selected <- attr(p, "selection_table"); selected$analysis_scope <- scope; selected$cohort <- cohort; selected$patient <- patient; selected$gene_set <- name
      selections[[paste(key, collection)]] <<- selected
      rna_save_plot(p, paste0("ORA_", gsub("[^A-Za-z0-9._-]", "_", collection)), dd, cfg, 8.5, max(5, 2.8 + .5 * nrow(selected)))
    }
  }
  for (cohort in unique(as.character(genes$cohort))) {
    data <- genes[as.character(genes$cohort) == cohort, , drop = FALSE]
    specs <- remodeling_enrichment_sets(data, alphas, cfg$enrichment$remodeling_nominal_p_cutoff %||% .05)
    for (name in names(specs)) execute(specs[[name]], "cohort", cohort, NA_character_, name)
  }
  patient <- remodeling$patient_inference
  if (!is.null(patient) && nrow(patient) && all(c("patient", "gene_id", "patient_test_in_family", "patient_test_bh", "patient_effect_pass") %in% names(patient))) {
    for (id in unique(as.character(patient$patient))) for (alpha in alphas) {
      x <- patient[as.character(patient$patient) == id, , drop = FALSE]
      background <- x$gene_id[x$patient_test_in_family %in% TRUE]
      targets <- x$gene_id[x$patient_test_in_family %in% TRUE & is.finite(x$patient_test_bh) & x$patient_test_bh <= alpha & x$patient_effect_pass %in% TRUE]
      execute(list(target = targets, background = background, label = paste0(id, ": remodeling FDR < ", alpha)), "patient_remodeling", NA_character_, id, paste0("patient_", dtu_fdr_tag(alpha), "_plus_effect"))
    }
  }
  out <- list(status = "completed", result = dplyr::bind_rows(results), runs = dplyr::bind_rows(statuses), membership = dplyr::bind_rows(memberships), plot_selection = dplyr::bind_rows(selections), references = resources$metadata)
  out$complete <- !any(out$runs$status == "failed_collection")
  for (name in c("result", "runs", "membership", "plot_selection")) rna_write_tsv(out[[name]], file.path(dirs$tables, paste0("remodeling_enrichment_", name, ".tsv.gz")))
  rna_atomic_rds(out, file.path(dirs$objects, "remodeling_enrichment.rds"))
  if (!out$complete && isTRUE(cfg$execution$strict_completion)) stop("Remodeling enrichment failed a collection; inspect enrichment run-status tables.", call. = FALSE)
  out
}
