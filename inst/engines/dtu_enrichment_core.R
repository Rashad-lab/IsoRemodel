# Source-derived deterministic matching and pathway tests, Module06 v1.3.0.

map_gene_symbols <- function(genes) {
  genes <- unique(clean_character(genes))
  genes <- genes[!is.na(genes)]
  symbol_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
  direct <- tibble(
    gene_id = intersect(genes, symbol_keys),
    SYMBOL = intersect(genes, symbol_keys),
    mapping_type = "direct_symbol"
  )

  unmatched <- setdiff(genes, direct$gene_id)
  symbol_prefix_candidates <- tibble(
    gene_id = unmatched,
    SYMBOL = str_remove(unmatched, "_.*$")
  ) %>%
    filter(SYMBOL %in% symbol_keys) %>%
    mutate(mapping_type = "Pigeon_compound_symbol_prefix") %>%
    distinct(gene_id, .keep_all = TRUE)

  unmatched <- setdiff(unmatched, symbol_prefix_candidates$gene_id)
  ensembl_candidates <- tibble(
    gene_id = unmatched,
    stripped_gene_id = strip_ensembl_version(
      str_extract(unmatched, "ENSG[0-9]+(?:\\.[0-9]+)?")
    )
  ) %>%
    filter(!is.na(stripped_gene_id))
  ensembl_map <- if (nrow(ensembl_candidates)) {
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(ensembl_candidates$stripped_gene_id),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )) %>%
      as_tibble() %>%
      transmute(
        stripped_gene_id = as.character(ENSEMBL),
        SYMBOL = clean_character(SYMBOL)
      ) %>%
      filter(!is.na(SYMBOL)) %>%
      inner_join(ensembl_candidates, by = "stripped_gene_id") %>%
      transmute(gene_id, SYMBOL, mapping_type = "ENSEMBL_to_symbol")
  } else {
    tibble(gene_id = character(), SYMBOL = character(), mapping_type = character())
  }

  bind_rows(direct, symbol_prefix_candidates, ensembl_map) %>%
    filter(!is.na(gene_id), !is.na(SYMBOL)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

load_module06_pathways <- function() {
  msig_formals <- names(formals(msigdbr::msigdbr))
  msig <- if ("db_species" %in% msig_formals) {
    msigdbr::msigdbr(db_species = "HS", species = "Homo sapiens")
  } else {
    msigdbr::msigdbr(species = "Homo sapiens")
  }
  database_version <- if ("db_version" %in% names(msig)) paste(sort(unique(as.character(msig$db_version))), collapse = ";") else NA_character_
  collection_col <- pick_column(
    msig, c("gs_collection", "gs_cat"), required = TRUE,
    label = "MSigDB collection"
  )
  subcollection_col <- pick_column(
    msig, c("gs_subcollection", "gs_subcat"), required = FALSE,
    label = "MSigDB subcollection"
  )
  term_col <- pick_column(msig, c("gs_name", "gs_id"), required = TRUE, label = "MSigDB term")
  symbol_col <- pick_column(
    msig, c("gene_symbol", "human_gene_symbol"), required = TRUE,
    label = "MSigDB gene symbol"
  )

  msig <- as_tibble(msig) %>%
    transmute(
      collection_code = as.character(.data[[collection_col]]),
      subcollection_code = if (is.null(subcollection_col)) NA_character_ else
        as.character(.data[[subcollection_col]]),
      pathway = as.character(.data[[term_col]]),
      gene = as.character(.data[[symbol_col]])
    ) %>%
    filter(!is.na(pathway), !is.na(gene)) %>%
    distinct()

  hallmark <- msig %>%
    filter(collection_code == "H") %>%
    transmute(collection = "MSigDB Hallmark", pathway, gene)

  gobp <- msig %>%
    filter(
      collection_code == "C5",
      str_detect(coalesce(subcollection_code, ""), regex("GO.?BP|GOBP", ignore_case = TRUE)) |
        str_detect(pathway, "^GOBP_")
    ) %>%
    transmute(collection = "GO biological process", pathway, gene)

  kegg <- msig %>%
    filter(
      collection_code == "C2",
      str_detect(coalesce(subcollection_code, ""), regex("KEGG", ignore_case = TRUE)) |
        str_detect(pathway, "^KEGG_")
    ) %>%
    transmute(collection = "KEGG", pathway, gene)

  pathways <- bind_rows(gobp, hallmark, kegg) %>% distinct()
  collection_qc <- pathways %>%
    count(collection, name = "gene_pathway_memberships") %>%
    left_join(
      pathways %>%
        group_by(collection) %>%
        summarise(pathways = n_distinct(pathway), .groups = "drop"),
      by = "collection"
    )
  missing_collections <- setdiff(cfg$ora$collections, collection_qc$collection)
  if (length(missing_collections)) {
    stop("No MSigDB terms recovered for: ", paste(missing_collections, collapse = ", "))
  }
  list(term2gene = pathways, qc = collection_qc, database_version = database_version)
}

standard_ora <- function(target_symbols, universe_symbols, term2gene, collection_name) {
  target_symbols <- intersect(unique(target_symbols), unique(universe_symbols))
  universe_symbols <- unique(universe_symbols)
  membership <- term2gene %>%
    filter(collection == collection_name, gene %in% universe_symbols) %>%
    distinct(pathway, gene)
  if (length(target_symbols) < cfg$ora$minimum_genes || !nrow(membership)) return(tibble())

  result <- membership %>%
    group_by(pathway) %>%
    summarise(
      background_genes = n_distinct(gene),
      overlap_genes = n_distinct(gene[gene %in% target_symbols]),
      genes = paste(sort(unique(gene[gene %in% target_symbols])), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(background_genes >= cfg$ora$minimum_term_genes, background_genes <= cfg$ora$maximum_term_genes) %>%
    mutate(
      foreground_genes = length(target_symbols),
      universe_genes = length(universe_symbols),
      selected_fraction = overlap_genes / foreground_genes,
      background_fraction = background_genes / universe_genes,
      fold_enrichment = selected_fraction / background_fraction,
      ora_p_value = phyper(
        overlap_genes - 1L,
        background_genes,
        universe_genes - background_genes,
        foreground_genes,
        lower.tail = FALSE
      ),
      ora_bh = p.adjust(ora_p_value, method = "BH"),
      collection = collection_name,
      .before = 1
    ) %>%
    arrange(ora_bh, ora_p_value, desc(overlap_genes), pathway)
  result
}

make_matched_design <- function(family_table, target_genes, controls_per_target) {
  mapping <- map_gene_symbols(family_table$gene_id)
  mapped_target_symbols <- mapping$SYMBOL[mapping$gene_id %in% target_genes] %>%
    unique()
  match_universe <- family_table %>%
    inner_join(mapping, by = "gene_id") %>%
    transmute(
      gene_id,
      SYMBOL,
      is_target_gene = gene_id %in% target_genes,
      log_pooled_gene_expression_CPM = log1p(pooled_mean_gene_CPM),
      log_detected_isoforms = log1p(detected_isoforms),
      log_tested_isoforms = log1p(tested_isoforms),
      log_annotated_transcripts = log1p(gencode_annotated_transcripts)
    ) %>%
    filter(if_all(starts_with("log_"), is.finite)) %>%
    arrange(SYMBOL, desc(is_target_gene), gene_id) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  target <- match_universe %>% filter(gene_id %in% target_genes)
  controls <- match_universe %>% filter(!(gene_id %in% target_genes))
  covariates <- c(
    "log_pooled_gene_expression_CPM", "log_detected_isoforms", "log_tested_isoforms",
    "log_annotated_transcripts"
  )
  target_covariate_coverage <- if (length(mapped_target_symbols)) {
    nrow(target) / length(mapped_target_symbols)
  } else {
    0
  }
  empty_integrity <- tibble(
    target_genes_requested = length(unique(target_genes)),
    mapped_target_symbols = length(mapped_target_symbols),
    matched_target_symbols = nrow(target),
    target_covariate_coverage = target_covariate_coverage,
    controls_per_target_requested = controls_per_target,
    controls_per_target_used = 0L,
    expected_pairs = 0L,
    observed_pairs = 0L,
    unique_control_symbols = 0L,
    duplicate_control_assignments = 0L,
    target_control_symbol_overlap = 0L
  )
  if (nrow(target) < cfg$ora$minimum_genes) {
    return(list(
      available = FALSE,
      reason = "too_few_mappable_targets_with_complete_matching_covariates",
      universe = match_universe,
      target = target,
      controls = tibble(),
      pairs = tibble(),
      balance = tibble(),
      integrity = empty_integrity
    ))
  }
  if (target_covariate_coverage < 1) {
    return(list(
      available = FALSE,
      reason = "incomplete_target_matching_covariates",
      universe = match_universe,
      target = target,
      controls = tibble(),
      pairs = tibble(),
      balance = tibble(),
      integrity = empty_integrity
    ))
  }

  controls_per_target_used <- min(
    controls_per_target,
    floor(nrow(controls) / nrow(target))
  )
  if (controls_per_target_used < 1L) {
    return(list(
      available = FALSE,
      reason = "fewer_controls_than_targets",
      universe = match_universe,
      target = target,
      controls = tibble(),
      pairs = tibble(),
      balance = tibble(),
      integrity = empty_integrity
    ))
  }

  full_covariates <- match_universe %>% select(all_of(covariates)) %>% as.matrix()
  centers <- colMeans(full_covariates)
  scales_sd <- apply(full_covariates, 2, sd)
  scales_sd[!is.finite(scales_sd) | scales_sd == 0] <- 1
  target_z <- sweep(sweep(as.matrix(target[, covariates]), 2, centers, "-"), 2, scales_sd, "/")
  control_z <- sweep(sweep(as.matrix(controls[, covariates]), 2, centers, "-"), 2, scales_sd, "/")
  distance_sq <- outer(rowSums(target_z^2), rowSums(control_z^2), "+") -
    2 * tcrossprod(target_z, control_z)
  distance_sq[distance_sq < 0] <- 0
  nearest_distance <- apply(distance_sq, 1, min)
  target_order <- order(-nearest_distance, target$SYMBOL, target$gene_id)
  available <- rep(TRUE, nrow(controls))
  pair_rows <- vector("list", nrow(target))
  selected_control_indices <- integer()
  for (order_index in seq_along(target_order)) {
    target_index <- target_order[[order_index]]
    available_indices <- which(available)
    candidate_order <- order(
      distance_sq[target_index, available_indices],
      controls$SYMBOL[available_indices],
      controls$gene_id[available_indices]
    )
    chosen <- head(available_indices[candidate_order], controls_per_target_used)
    if (length(chosen) != controls_per_target_used) {
      stop("Deterministic matching exhausted eligible controls.")
    }
    available[chosen] <- FALSE
    selected_control_indices <- c(selected_control_indices, chosen)
    pair_rows[[order_index]] <- tibble(
      target_order = order_index,
      target_gene_id = target$gene_id[[target_index]],
      target_SYMBOL = target$SYMBOL[[target_index]],
      control_rank = seq_along(chosen),
      control_gene_id = controls$gene_id[chosen],
      control_SYMBOL = controls$SYMBOL[chosen],
      standardized_distance = sqrt(distance_sq[target_index, chosen])
    )
  }
  pairs <- bind_rows(pair_rows)
  matched_controls <- controls[selected_control_indices, , drop = FALSE]

  target_means <- colMeans(as.matrix(target[, covariates]))
  control_means <- colMeans(as.matrix(matched_controls[, covariates]))
  balance <- tibble(
    covariate = covariates,
    target_mean = as.numeric(target_means),
    matched_control_mean = as.numeric(control_means),
    full_universe_sd = as.numeric(scales_sd),
    standardized_mean_difference =
      (as.numeric(target_means) - as.numeric(control_means)) / as.numeric(scales_sd),
    absolute_standardized_mean_difference = abs(standardized_mean_difference)
  )

  integrity <- tibble(
    target_genes_requested = length(unique(target_genes)),
    mapped_target_symbols = length(mapped_target_symbols),
    matched_target_symbols = nrow(target),
    target_covariate_coverage = target_covariate_coverage,
    controls_per_target_requested = controls_per_target,
    controls_per_target_used = controls_per_target_used,
    expected_pairs = nrow(target) * controls_per_target_used,
    observed_pairs = nrow(pairs),
    unique_control_symbols = n_distinct(pairs$control_SYMBOL),
    duplicate_control_assignments = nrow(pairs) - n_distinct(pairs$control_SYMBOL),
    target_control_symbol_overlap = sum(
      pairs$control_SYMBOL %in% target$SYMBOL
    )
  )

  list(
    available = TRUE,
    reason = "success",
    universe = match_universe,
    target = target,
    controls = matched_controls,
    pairs = pairs,
    balance = balance,
    integrity = integrity
  )
}

matched_pathway_test <- function(matched_design, term2gene, collection_name) {
  if (!isTRUE(matched_design$available)) return(tibble())
  target_symbols <- matched_design$target$SYMBOL
  control_symbols <- matched_design$controls$SYMBOL
  universe_symbols <- unique(c(target_symbols, control_symbols))
  membership <- term2gene %>%
    filter(collection == collection_name, gene %in% universe_symbols) %>%
    distinct(pathway, gene)
  if (!nrow(membership)) return(tibble())

  pathway_size <- membership %>% count(pathway, name = "background_genes") %>%
    filter(background_genes >= cfg$ora$minimum_term_genes, background_genes <= cfg$ora$maximum_term_genes)
  membership <- membership %>% semi_join(pathway_size, by = "pathway")
  membership %>%
    group_by(pathway) %>%
    summarise(
      background_genes = n_distinct(gene),
      observed_genes = n_distinct(gene[gene %in% target_symbols]),
      matched_control_overlap_genes = n_distinct(gene[gene %in% control_symbols]),
      genes = paste(sort(unique(gene[gene %in% target_symbols])), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(background_genes >= cfg$ora$minimum_term_genes, background_genes <= cfg$ora$maximum_term_genes) %>%
    mutate(
      collection = collection_name,
      foreground_genes = length(target_symbols),
      matched_background_genes = length(control_symbols),
      matched_universe_genes = length(universe_symbols),
      selected_fraction = observed_genes / foreground_genes,
      matched_fraction = matched_control_overlap_genes / matched_background_genes,
      fold_enrichment = if_else(
        matched_control_overlap_genes > 0,
        selected_fraction / matched_fraction,
        NA_real_
      ),
      matched_p_value = phyper(
        observed_genes - 1L,
        observed_genes + matched_control_overlap_genes,
        matched_universe_genes - observed_genes - matched_control_overlap_genes,
        foreground_genes,
        lower.tail = FALSE
      ),
      matched_bh = p.adjust(matched_p_value, method = "BH"),
      .before = 1
    ) %>%
    arrange(matched_bh, matched_p_value, desc(observed_genes), pathway)
}

format_pathway_label <- function(pathway, collection) {
  pathway <- str_remove(pathway, "^(HALLMARK_|GOBP_|KEGG_[A-Z]+_|KEGG_)")
  pathway <- str_replace_all(pathway, "_", " ")
  if (collection != "GO biological process") pathway <- str_to_title(str_to_lower(pathway))
  str_wrap(pathway, width = 34)
}

select_pathway_plot_data <- function(data, matched = TRUE) {
  if (!nrow(data)) return(data)
  p_col <- if (matched) "matched_p_value" else "ora_p_value"
  bh_col <- if (matched) "matched_bh" else "ora_bh"
  count_col <- if (matched) "observed_genes" else "overlap_genes"
  significant <- data %>% filter(.data[[bh_col]] < cfg$ora$reporting_fdr)
  if (nrow(significant)) {
    selected <- significant %>%
      arrange(.data[[bh_col]], .data[[p_col]], desc(.data[[count_col]]))
    tier <- "pathway_BH_below_0p10"
  } else {
    selected <- data %>%
      filter(.data[[p_col]] < cfg$ora$nominal_fallback) %>%
      arrange(.data[[p_col]], .data[[bh_col]])
    tier <- "raw_P_below_0p05_fallback"
  }
  if (!nrow(selected)) return(selected)
  selected %>%
    slice_head(n = cfg$ora$top_terms) %>%
    mutate(
      display_tier = tier,
      plotted_p = .data[[p_col]],
      plotted_bh = .data[[bh_col]],
      plotted_genes = .data[[count_col]],
      plot_score = -log10(pmax(plotted_p, .Machine$double.xmin)),
      pathway_label = map2_chr(pathway, collection, format_pathway_label),
      pathway_label = fct_reorder(pathway_label, selected_fraction)
    )
}

plot_pathway_dotplot <- function(data, title, subtitle, matched = TRUE) {
  plot_data <- select_pathway_plot_data(data, matched = matched)
  if (!nrow(plot_data)) return(NULL)
  nominal_only <- all(plot_data$display_tier == "raw_P_below_0p05_fallback")
  evidence_caption <- if (nominal_only) {
    paste0("Nominal P < ", cfg$ora$nominal_fallback,
           "; no pathways pass FDR < ", cfg$ora$reporting_fdr)
  } else {
    paste0("Pathway FDR < ", cfg$ora$reporting_fdr)
  }
  ggplot(
    plot_data,
    aes(x = selected_fraction, y = pathway_label, size = plotted_genes, color = plot_score)
  ) +
    geom_point(alpha = 0.94) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.10, 0.12))) +
    scale_color_gradient(low = "#F4A582", high = "#B2182B",
                         breaks = scales::breaks_pretty(n = 3),
                         labels = scales::label_number(accuracy = 0.1)) +
    scale_size_continuous(range = c(2, 7),
                          breaks = function(x) unique(round(scales::breaks_pretty(n = 3)(x)))) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = evidence_caption,
      x = "Selected genes in pathway",
      y = NULL,
      color = expression(-log[10](P)),
      size = "Genes"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18, lineheight = 1.05),
      plot.subtitle = element_text(size = 13, lineheight = 1.10),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 12, lineheight = 0.95),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(14, 24, 14, 20)
    ) + publication_layout() +
    guides(
      color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                             barwidth = grid::unit(2.6, "in"),
                             barheight = grid::unit(0.16, "in"), order = 1),
      size = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5,
                          order = 2)
    )
}
