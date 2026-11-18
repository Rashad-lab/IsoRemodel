# Plotting only: read the corrected DTU handoff and its existing result tables.
# No model fitting, hypothesis tests, stageR adjustment or enrichment tests run.

rna_saved_dtu_or <- function(x, y) if (is.null(x) || !length(x)) y else x
rna_saved_dtu_tag <- function(alpha) paste0("FDR", if (isTRUE(all.equal(alpha, .05))) "05" else if (isTRUE(all.equal(alpha, .10))) "10" else gsub("[^0-9]", "", format(alpha)))
rna_saved_dtu_wrap <- function(x, width = 38) vapply(as.character(x), function(z) paste(strwrap(z, width), collapse = "\n"), character(1))
rna_saved_dtu_structure <- function(x) {
  labels <- c(`full-splice_match` = "FSM", `incomplete-splice_match` = "ISM",
              novel_in_catalog = "NIC", novel_not_in_catalog = "NNC")
  x <- as.character(x)
  hit <- x %in% names(labels)
  x[hit] <- unname(labels[x[hit]])
  x[!hit] <- gsub("_", " ", x[!hit], fixed = TRUE)
  x
}
rna_saved_dtu_structure_key <- function(x) {
  keys <- c(FSM = "full splice match", ISM = "incomplete splice match",
            NIC = "novel in catalog", NNC = "novel not in catalog")
  selected <- intersect(unique(rna_saved_dtu_structure(x)), names(keys))
  paste(paste0(selected, ": ", keys[selected]), collapse = "; ")
}
rna_saved_dtu_read_tsv <- function(path) {
  if (!file.exists(path)) stop("Saved DTU table is missing: ", path, call. = FALSE)
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  tryCatch(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
                            quote = "\"", comment.char = "", na.strings = "NA"),
           error = function(e) {
             if (grepl("no lines available|first five rows are empty", conditionMessage(e), ignore.case = TRUE)) data.frame()
             else stop(e)
           })
}
rna_saved_dtu_require_fields <- function(x, fields, label) {
  if (!is.data.frame(x) || !all(fields %in% names(x)))
    stop(label, " lacks saved fields: ", paste(setdiff(fields, names(x)), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
rna_saved_dtu_palette <- function(cfg, patients) {
  palette <- unlist(cfg$reporting$palettes$patients)
  if (is.null(palette) || !all(patients %in% names(palette)))
    palette <- stats::setNames(scales::hue_pal()(length(patients)), patients)
  palette[patients]
}

# Selection copies the original display rule, using saved P/FDR values only.
# No pathway membership, mapping, matching or statistical value is recomputed.
rna_saved_dtu_select_ora <- function(data, matched, cfg) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  p_col <- if (matched) "matched_p_value" else "ora_p_value"
  fdr_col <- if (matched) "matched_bh" else "ora_bh"
  count_col <- if (matched) "observed_genes" else "overlap_genes"
  rna_saved_dtu_require_fields(data, c("pathway", "collection", "selected_fraction", p_col, fdr_col, count_col), "Saved pathway results")
  pathway_fdr <- rna_saved_dtu_or(cfg$enrichment$reporting_fdr, .10)
  nominal <- rna_saved_dtu_or(cfg$enrichment$nominal_display_fallback_p, .05)
  maximum <- as.integer(rna_saved_dtu_or(cfg$enrichment$maximum_display_terms, 10L))
  ok <- is.finite(data[[fdr_col]]) & data[[fdr_col]] < pathway_fdr
  if (any(ok)) {
    x <- data[ok, , drop = FALSE]
    x <- x[order(x[[fdr_col]], x[[p_col]], -x[[count_col]], na.last = TRUE), , drop = FALSE]
    tier <- "pathway_FDR"
  } else {
    x <- data[is.finite(data[[p_col]]) & data[[p_col]] < nominal, , drop = FALSE]
    if (!nrow(x)) return(x)
    x <- x[order(x[[p_col]], x[[fdr_col]], na.last = TRUE), , drop = FALSE]
    tier <- "nominal_P_only"
  }
  x <- utils::head(x, maximum)
  x$plotted_p <- x[[p_col]]; x$plotted_fdr <- x[[fdr_col]]
  x$plotted_genes <- x[[count_col]]
  x$plot_score <- -log10(pmax(x$plotted_p, .Machine$double.xmin))
  x$display_tier <- tier
  labels <- sub("^(HALLMARK_|GOBP_|KEGG_[A-Z]+_|KEGG_)", "", x$pathway)
  labels <- gsub("_", " ", labels, fixed = TRUE)
  x$pathway_label <- rna_saved_dtu_wrap(labels, 40)
  # Keep ties deterministic and preserve the selected table's relative order.
  levels <- unique(x$pathway_label[order(x$selected_fraction, seq_len(nrow(x)))])
  x$pathway_label <- factor(x$pathway_label, levels = levels)
  x
}

rna_saved_dtu_ora_plot <- function(x, alpha, method, patient = NULL, cfg = list()) {
  collection <- unique(as.character(x$collection))[[1L]]
  collection <- if (collection == "GO biological process") "GO biological process" else if (collection == "MSigDB Hallmark") "Hallmark pathways" else collection
  nominal <- identical(unique(x$display_tier), "nominal_P_only")
  cap <- if (nominal) paste0("Nominal P < ", rna_saved_dtu_or(cfg$enrichment$nominal_display_fallback_p, .05),
                             "; no pathways at FDR < ", rna_saved_dtu_or(cfg$enrichment$reporting_fdr, .10)) else
    paste0("Pathway FDR < ", rna_saved_dtu_or(cfg$enrichment$reporting_fdr, .10))
  subtitle <- paste0("DTU FDR < ", alpha, if (!is.null(patient)) paste0("; ", patient) else "",
                     if (method == "matched") "; matched background" else "")
  xmax <- min(1, max(.05, ceiling(max(x$selected_fraction, na.rm = TRUE) * 20) / 20))
  ticks <- unique(round(seq(0, xmax, length.out = 4L), 3L))
  counts <- unique(as.integer(round(scales::breaks_pretty(n = 3)(range(x$plotted_genes, na.rm = TRUE)))))
  counts <- counts[counts >= min(x$plotted_genes) & counts <= max(x$plotted_genes)]
  if (!length(counts)) counts <- unique(x$plotted_genes)
  ggplot2::ggplot(x, ggplot2::aes(selected_fraction, pathway_label, size = plotted_genes, color = plot_score)) +
    ggplot2::geom_point(alpha = .94) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, xmax),
                              breaks = ticks, expand = ggplot2::expansion(mult = c(0, .06))) +
    ggplot2::scale_color_gradient(low = "#F4A582", high = "#B2182B", breaks = scales::breaks_pretty(n = 3),
                                 labels = scales::label_number(accuracy = .1)) +
    ggplot2::scale_size_continuous(range = c(2, 7), breaks = counts) +
    ggplot2::labs(title = collection, subtitle = subtitle, caption = cap,
                  x = "Selected genes in pathway", y = NULL, color = expression(-log[10](P)), size = "Genes") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical", axis.text.y = ggplot2::element_text(lineheight = 1.05),
                   plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", barwidth = grid::unit(3.2, "in"),
                                                  barheight = grid::unit(.16, "in"), order = 1),
                    size = ggplot2::guide_legend(nrow = 1, title.position = "top", order = 2))
}

