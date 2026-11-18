# IsoRemodel 0.2 display-only layout helpers. Analytical selections and values
# remain unchanged. New .layout_* columns contain presentation coordinates only.

# Plot-only layout repairs. These helpers alter labels, guides and canvas size;
# the saved observations, confidence intervals, tests and selections stay fixed.

rna_layout_family <- function(stem, module = "") {
  stem <- as.character(stem)[1L]
  if (grepl("ORA_|enrichment", stem, ignore.case = TRUE) &&
      !grepl("structural_category_enrichment|enrichment_background", stem)) return("enrichment")
  if (grepl("paired_gain_loss_full_vs_common_depth", stem)) return("gain_loss_paired")
  if (grepl("gain_loss_gene_floor", stem)) return("gain_loss_floor")
  if (grepl("full_length_reads_by_sample|sample_depth_and_common_target", stem)) return("sample_depth")
  if (grepl("technical_annotation_coverage", stem)) return("annotation_heatmap")
  if (grepl("gain_loss_technical_|paired_technical_association", stem)) return("technical_heatmap")
  if (grepl("empirical_null_class_composition|effect_threshold_sensitivity", stem)) return("null_heatmap")
  if (grepl("DTU_patient_driver_heatmap", stem)) return("dtu_driver_heatmap")
  if (grepl("DTU_dominant_structural_transitions", stem)) return("dtu_structure")
  if (grepl("DTU_paired_isoform_trajectories", stem)) return("gene_example")
  if (grepl("filtering_testing_family_flow|filter_flow", stem)) return("filter_flow")
  if (grepl("upset|intersections|directional_support", stem, ignore.case = TRUE)) return("overlap")
  if (grepl("RNA_structural_composition|RNA_structure_by_complexity|gene_complexity_structural_categories", stem))
    return("structural_composition")
  if (grepl("^paired_", stem)) return("paired_qc")
  "other"
}

# Valid integer breaks for gene counts. Never round values outside the range.
rna_integer_legend_breaks <- function(counts) {
  r <- range(counts, finite = TRUE)
  if (any(!is.finite(r))) return(numeric())
  if (r[1L] == r[2L]) return(if (r[1L] == round(r[1L])) r[1L] else numeric())
  value <- unique(round(pretty(r, n = 3L)))
  value <- value[value >= ceiling(r[1L]) & value <= floor(r[2L])]
  if (!length(value)) value <- unique(counts[is.finite(counts) & counts == round(counts)])
  value
}

# Blank patchwork helpers/spacers are layout machinery, not separate figures.
rna_plot_meaningful <- function(plot) {
  if (inherits(plot, c("plot_spacer", "spacer", "zeroGrob", "null"))) return(FALSE)
  if (inherits(plot, "patchwork"))
    return(any(vapply(seq_len(length(plot)), function(i) rna_plot_meaningful(plot[[i]]), logical(1))))
  if (!inherits(plot, "ggplot") || !length(plot$layers)) return(FALSE)
  any(vapply(plot$layers, function(layer) !inherits(layer$geom, "GeomBlank"), logical(1)))
}

rna_layout_variant_roomy <- function(variant) {
  value <- if (is.list(variant)) paste(unlist(variant[intersect(c("id", "name", "layout"), names(variant))]), collapse = " ") else paste(variant, collapse = " ")
  grepl("roomy|wide", value, ignore.case = TRUE)
}

rna_layout_layout_family <- function(family) {
  aliases <- c(driver_heatmap = "dtu_driver_heatmap", structural_heatmap = "dtu_structure",
    intersection = "overlap", intersections = "overlap", paired_gain_loss = "gain_loss_paired")
  if (family %in% names(aliases)) unname(aliases[[family]]) else family
}

rna_layout_layout_wrap <- function(x, width = 30L) {
  if (!is.character(x)) return(x)
  vapply(x, function(z) {
    if (is.na(z)) return(NA_character_)
    paste(vapply(strsplit(z, "\n", fixed = TRUE)[[1L]], function(line)
      paste(strwrap(line, width = width), collapse = "\n"), character(1)), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

rna_layout_layout_labels <- function(x) {
  x <- as.character(x)
  dictionary <- c(
    "full-splice_match" = "FSM", "full_splice_match" = "FSM",
    "incomplete-splice_match" = "ISM", "incomplete_splice_match" = "ISM",
    "novel_in_catalog" = "NIC", "novel_not_in_catalog" = "NNC",
    "all_6_pairs" = "All six pairs", "excluding_ptn2_5_pairs" = "Patient 2 excluded")
  for (key in names(dictionary)) x <- gsub(key, unname(dictionary[[key]]), x, fixed = TRUE)
  x <- gsub("ptn([0-9]+)[_ ]primary", "P\\1 primary", x, ignore.case = TRUE)
  x <- gsub("ptn([0-9]+)[_ ]recurrent", "P\\1 recurrent", x, ignore.case = TRUE)
  x <- gsub("^ptn([0-9]+)$", "P\\1", x, ignore.case = TRUE)
  x <- gsub("[Mm]odule[[:space:]]+0*[0-9]+[[:space:]]*[-:–—]?[[:space:]]*", "", x)
  x <- gsub("\\bBH\\b", "FDR", x)
  rna_layout_layout_wrap(x, 24L)
}

rna_layout_layout_caption <- function(x) {
  if (!is.character(x)) return(x)
  x <- gsub("Nominal P < ([0-9.]+); no pathways pass FDR < ([0-9.]+)",
            "Nominal P < \\1; none meet pathway FDR < \\2", x)
  rna_layout_layout_wrap(x, 64L)
}

rna_layout_sample_labels <- function(x, compact = TRUE) {
  x <- as.character(x)
  if (compact) {
    x <- gsub("ptn([0-9]+)[_ ]primary", "P\\1-P", x, ignore.case = TRUE)
    x <- gsub("ptn([0-9]+)[_ ]recurrent", "P\\1-R", x, ignore.case = TRUE)
  } else {
    x <- rna_layout_layout_labels(x)
    x <- gsub(" (primary|recurrent)$", "\n\\1", x)
  }
  x
}

rna_layout_scale_change <- function(plot, aesthetic, change) {
  # Clone existing scales so limits, level ordering and palettes are preserved.
  # An absent position scale is introduced only for a known discrete axis.
  plot$scales <- plot$scales$clone()
  indexes <- which(vapply(plot$scales$scales, function(s) aesthetic %in% s$aesthetics, logical(1)))
  if (length(indexes)) for (i in indexes) {
    s <- plot$scales$scales[[i]]$clone()
    s <- change(s)
    plot$scales$scales[[i]] <- s
  }
  plot
}

rna_layout_discrete_labels <- function(plot, aesthetic, labels = rna_layout_layout_labels, add = FALSE) {
  found <- any(vapply(plot$scales$scales, function(s) aesthetic %in% s$aesthetics, logical(1)))
  plot <- rna_layout_scale_change(plot, aesthetic, function(s) {
    if (inherits(s, "ScaleDiscrete")) s$labels <- labels
    s
  })
  if (!found && add && aesthetic %in% c("x", "y")) {
    scale <- if (aesthetic == "x") ggplot2::scale_x_discrete else ggplot2::scale_y_discrete
    plot <- plot + scale(labels = labels)
  }
  plot
}

rna_layout_three_breaks <- function(x) {
  x <- range(x, finite = TRUE)
  if (!length(x) || any(!is.finite(x))) return(numeric())
  if (diff(x) < .Machine$double.eps * max(1, abs(x))) return(x[1L])
  z <- pretty(x, n = 3L)
  z <- z[z >= x[1L] & z <= x[2L]]
  if (!length(z)) z <- mean(x)
  if (length(z) > 3L) z <- z[unique(round(seq(1, length(z), length.out = 3L)))]
  z
}

rna_layout_numeric_labels <- function(x) {
  ifelse(is.na(x), NA_character_, format(signif(x, 3L), trim = TRUE, scientific = FALSE))
}

rna_layout_guide_repairs <- function(plot, family) {
  for (aesthetic in c("fill", "colour")) {
    plot <- rna_layout_scale_change(plot, aesthetic, function(s) {
      if (inherits(s, "ScaleDiscrete")) {
        if (family %in% c("structural_composition", "gene_example", "dtu_structure")) s$labels <- rna_layout_layout_labels
        s$guide <- ggplot2::guide_legend(nrow = if (family == "structural_composition") 2L else 1L,
          byrow = TRUE, title.position = "top", override.aes = list(alpha = 1))
      } else {
        s$breaks <- rna_layout_three_breaks
        if (!family %in% c("null_heatmap", "annotation_heatmap")) s$labels <- rna_layout_numeric_labels
        if (family == "null_heatmap") s$labels <- function(x)
          paste0(format(signif(x * 100, 3L), trim = TRUE, scientific = FALSE), "%")
        s$guide <- ggplot2::guide_colourbar(title.position = "top", title.hjust = .5,
          barwidth = grid::unit(50, "mm"), barheight = grid::unit(4, "mm"))
      }
      s
    })
  }
  plot <- rna_layout_scale_change(plot, "size", function(s) {
    if (!inherits(s, "ScaleDiscrete")) {
      s$breaks <- function(x) unique(pmax(0, round(rna_layout_three_breaks(x))))
      s$labels <- function(x) format(x, trim = TRUE, scientific = FALSE)
      s$guide <- ggplot2::guide_legend(title.position = "top", nrow = 1L)
    }
    s
  })
  plot
}

rna_layout_heatmap_text_contrast <- function(plot) {
  # Choose text ink from the existing rendered fill palette. A display mapping
  # looks up the unchanged fill values, including after individual facet export.
  built <- tryCatch(ggplot2::ggplot_build(plot), error = function(e) NULL)
  if (is.null(built)) return(plot)
  scale <- built$plot$scales$get_scales("fill")
  if (is.null(scale)) return(plot)
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (!inherits(layer$geom, "GeomText")) next
    data <- if (is.data.frame(layer$data)) layer$data else plot$data
    if (!is.data.frame(data) || !nrow(data)) next
    mapping <- layer$mapping$fill
    if (is.null(mapping)) mapping <- plot$mapping$fill
    if (is.null(mapping)) next
    values <- tryCatch(rlang::eval_tidy(mapping, data = data), error = function(e) NULL)
    if (is.null(values) || length(values) != nrow(data)) next
    colours <- tryCatch(scale$map(values), error = function(e) NULL)
    if (is.null(colours) || length(colours) != nrow(data)) next
    rgb <- tryCatch(grDevices::col2rgb(colours), error = function(e) NULL)
    if (is.null(rgb)) next
    luminance <- colSums(rgb * c(.2126, .7152, .0722)) / 255
    ink <- ifelse(luminance < .55, "white", "black")
    colour_mapping <- local({
      keys <- as.character(values); keep <- !duplicated(keys)
      saved_keys <- keys[keep]; saved_ink <- ink[keep]
      ink_for_value <- function(x) saved_ink[match(as.character(x), saved_keys)]
      ggplot2::aes(colour = I(ink_for_value(!!mapping)))$colour
    })
    revised <- local({original_layer <- layer; ggplot2::ggproto(NULL, original_layer)})
    # Keep ggplot2's mapping object intact (it is S7 in ggplot2 4.x). Assigning
    # that class to an ordinary list can create invalid recursive dispatch.
    revised$mapping <- layer$mapping
    revised$mapping$colour <- colour_mapping
    revised$aes_params <- as.list(layer$aes_params)
    revised$aes_params$colour <- NULL
    revised$show.legend <- FALSE
    plot$layers[[i]] <- revised
  }
  plot
}

rna_layout_dimensions <- function(record, variant) {
  f <- record$family
  if (is.null(f) || !length(f)) f <- rna_layout_family(record$stem, record$module)
  f <- rna_layout_layout_family(f)
  width <- if (!is.null(record$width) && is.finite(record$width)) as.numeric(record$width) else 8
  height <- if (!is.null(record$height) && is.finite(record$height)) as.numeric(record$height) else 5
  bases <- list(gain_loss_paired = c(6.4, 4.2), gain_loss_floor = c(7.2, 4.8),
    sample_depth = c(7.4, 5.1), paired_qc = c(5.4, 4.6),
    structural_composition = c(8, 5.7), technical_heatmap = c(10, 7),
    annotation_heatmap = c(9, 7), null_heatmap = c(10.5, 8),
    gene_example = c(7, 4.8), dtu_structure = c(6.8, 5.8),
    dtu_driver_heatmap = c(7.3, 7), filter_flow = c(8.2, 4.8),
    overlap = c(8.5, 5), enrichment = c(9.2, 6))
  if (f %in% names(bases)) {
    base <- bases[[f]]
    # Keep the vertical space already calculated from actual pathway label
    # lines or heatmap rows. Compact means less unused space, not fewer rows.
    width <- base[1L]
    height <- if (f %in% c("enrichment", "dtu_driver_heatmap", "technical_heatmap", "annotation_heatmap", "null_heatmap")) max(base[2L], height) else base[2L]
  }
  if (f == "enrichment" && !is.null(record$plot) && is.data.frame(record$plot$data)) {
    data <- record$plot$data
    field <- intersect(c("pathway_label", "label", "Description", "pathway"), names(data))
    if (length(field)) {
      labels <- unique(rna_layout_layout_wrap(gsub("_", " ", as.character(data[[field[1L]]])),
        if (rna_layout_variant_roomy(variant)) 48L else 40L))
      lines <- sum(vapply(strsplit(labels[!is.na(labels)], "\n", fixed = TRUE), length, integer(1)))
      height <- max(height, 2.5 + .24 * lines)
    }
  }
  if (f == "gene_example" && grepl("__panel_|__isoform_", record$stem)) width <- 5.2
  if (rna_layout_variant_roomy(variant)) c(width = width * 1.18, height = height * 1.12)
  else c(width = width, height = height)
}

rna_layout_repair_plot <- function(plot, record, variant) {
  if (!inherits(plot, "ggplot") || inherits(plot, "patchwork")) return(plot)
  family <- record$family
  if (is.null(family) || !length(family)) family <- rna_layout_family(record$stem, record$module)
  family <- rna_layout_layout_family(family)
  roomy <- rna_layout_variant_roomy(variant)
  p <- plot
  p <- rna_layout_guide_repairs(p, family)
  if (is.character(p$labels$title)) p$labels$title <- rna_layout_layout_wrap(gsub("\n", " ", rna_layout_layout_labels(p$labels$title)), if (roomy) 64L else 52L)
  if (is.character(p$labels$subtitle)) p$labels$subtitle <- rna_layout_layout_wrap(gsub("\n", " ", rna_layout_layout_labels(p$labels$subtitle)), if (roomy) 70L else 58L)
  if (is.character(p$labels$caption)) p$labels$caption <- rna_layout_layout_caption(p$labels$caption)
  if (is.character(p$labels$y)) p$labels$y <- rna_layout_layout_wrap(p$labels$y, 22L)
  if (is.character(p$labels$x)) p$labels$x <- rna_layout_layout_wrap(p$labels$x, 35L)
  p <- p + ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.1),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
    legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
    legend.spacing.y = grid::unit(3, "mm"),
    plot.margin = ggplot2::margin(12, 25, 12, 12),
    strip.text = ggplot2::element_text(margin = ggplot2::margin(5, 5, 5, 5)))

  if (!inherits(p$facet, "FacetNull")) {
    # ggproto evaluates its parent lazily. Capture it in an immutable local
    # binding before assigning p$facet, otherwise the parent becomes itself.
    facet <- local({original_facet <- p$facet; ggplot2::ggproto(NULL, original_facet)})
    facet$params <- as.list(p$facet$params)
    facet$params$labeller <- ggplot2::as_labeller(rna_layout_layout_labels)
    p$facet <- facet
  }
  if (family == "gain_loss_floor") p <- p + ggplot2::labs(
    y = "High-confidence isoforms\n(full depth)", x = "Minimum gene full-length reads")
  if (family == "gain_loss_paired") p <- p + ggplot2::labs(
    title = "Paired isoform gain and loss", y = "Isoforms", color = NULL,
    caption = "Filled: common-depth mean and 95% range. Open: full depth.")
  if (family == "paired_qc") {
    key <- sub("^paired_", "", sub("__.*$", "", record$stem))
    short <- c(ism_3prime_fragment_FL_fraction = "3-prime ISM reads (%)",
      ism_5prime_fragment_FL_fraction = "5-prime ISM reads (%)",
      ism_internal_fragment_FL_fraction = "Internal ISM reads (%)",
      log10_total_full_length_reads = "Full-length reads (log10)",
      mitochondrial_chromosome_FL_fraction = "Mitochondrial reads (%)",
      tss_within_threshold_FL_fraction = "Reads near TSS (%)",
      tts_within_threshold_FL_fraction = "Reads near TTS (%)")
    if (key %in% names(short)) p <- p + ggplot2::labs(y = unname(short[[key]]))
  }
  if (family %in% c("sample_depth", "structural_composition", "annotation_heatmap")) {
    # Keep explicit source ordering; only change how the IDs are printed.
    compact_sample_ids <- !roomy || inherits(p$coordinates, "CoordFlip")
    p <- rna_layout_discrete_labels(p, "x", function(z) rna_layout_sample_labels(z, compact = compact_sample_ids), add = TRUE)
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (family == "annotation_heatmap") 35 else 0,
      hjust = if (family == "annotation_heatmap") 1 else .5))
    is_sample_plot <- grepl("full_length_reads_by_sample|sample_depth_and_common_target|technical_annotation_coverage|RNA_structural_composition_(?!condition)", record$stem, perl = TRUE)
    if (is_sample_plot && compact_sample_ids) {
      caption <- p$labels$caption
      if (is.null(caption)) caption <- character()
      p <- p + ggplot2::labs(caption = paste(c(caption, "Sample labels: P = primary; R = recurrent."), collapse = "\n"))
    }
  }
  if (family == "structural_composition") {
    p <- rna_layout_discrete_labels(p, "fill", rna_layout_layout_labels)
    caption <- p$labels$caption
    if (is.null(caption)) caption <- character()
    p <- p + ggplot2::labs(fill = "Structural category",
      caption = paste(c(caption, "FSM: full splice match; ISM: incomplete splice match.",
        "NIC: novel in catalog; NNC: novel not in catalog."), collapse = "\n"))
  }
  if (family %in% c("technical_heatmap", "annotation_heatmap")) {
    p <- rna_layout_discrete_labels(p, "y", function(z) rna_layout_layout_wrap(rna_layout_layout_labels(z), 28L), add = TRUE)
    if (family == "technical_heatmap") p <- rna_layout_discrete_labels(p, "x", function(z) rna_layout_layout_wrap(rna_layout_layout_labels(z), 22L), add = TRUE)
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (roomy) 35 else 55, hjust = 1))
  }
  if (family == "gene_example") {
    standalone <- isTRUE(record$standalone) || grepl("__panel_|__isoform_", record$stem)
    p <- p + ggplot2::labs(y = "Isoform fraction", color = "Patient") +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = if (standalone) 2L else 1L,
        byrow = TRUE, title.position = if (standalone) "top" else "left"))
    p <- rna_layout_discrete_labels(p, "colour", rna_layout_layout_labels)
  }
  if (family == "dtu_structure") {
    p <- rna_layout_discrete_labels(p, "x", rna_layout_layout_labels, add = TRUE)
    p <- rna_layout_discrete_labels(p, "y", rna_layout_layout_labels, add = TRUE)
    p <- p + ggplot2::labs(x = "Primary structure", y = "Recurrent structure",
      caption = "FSM: full splice match; ISM: incomplete splice match.\nNIC: novel in catalog; NNC: novel not in catalog.")
  }
  if (family == "dtu_driver_heatmap") p <- p + ggplot2::labs(fill = "Change in isoform fraction")
  if (family == "null_heatmap") p <- rna_layout_heatmap_text_contrast(p)
  if (family == "enrichment") {
    p <- rna_layout_discrete_labels(p, "y", function(z) rna_layout_layout_wrap(gsub("_", " ", z), if (roomy) 48L else 40L), add = TRUE)
    d <- p$data
    columns <- intersect(c("selected_fraction", "enrichment_fraction"), names(d))
    if (is.data.frame(d) && length(columns)) {
      values <- d[[columns[1L]]]; values <- values[is.finite(values)]
      if (length(values) && all(values >= 0 & values <= 1)) {
        upper <- min(1, max(.01, max(values) * 1.12))
        p <- rna_layout_scale_change(p, "x", function(s) {
          s$limits <- c(0, upper); s$expand <- c(0, 0)
          s$breaks <- function(x) {z <- pretty(c(0, upper), n = 4L); z[z >= 0 & z <= upper]}
          s$labels <- scales::label_percent(accuracy = if (upper < .05) .1 else 1)
          s
        })
      }
    }
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(lineheight = 1.05),
      panel.grid.major.y = ggplot2::element_blank(), legend.box = "vertical")
  }
  p
}

# Rebuild only the axis geometry of an already selected enrichment plot.
# The original rows, counts, P values and FDR columns remain untouched.
rna_enrichment_layout <- function(plot, cfg, width, height) {
  d <- plot$data
  fraction <- intersect(c("selected_fraction", "enrichment_fraction"), names(d))
  count <- intersect(c("plotted_genes", "Count", "overlap_genes", "observed_genes"), names(d))
  term <- intersect(c("pathway", "Description", "pathway_label", "label"), names(d))
  if (!is.data.frame(d) || !nrow(d) || !length(fraction) || !length(count) || !length(term))
    return(list(plot = plot, width = width, height = height))
  values <- d[[fraction[1L]]]
  counts <- d[[count[1L]]]
  if (any(!is.finite(values) | values < 0 | values > 1) ||
      any(!is.finite(counts) | counts < 0 | counts != round(counts)))
    stop("Enrichment display requires saved fractions in [0,1] and integer counts.", call. = FALSE)
  labels <- as.character(d[[term[1L]]])
  labels <- sub("^(HALLMARK_|GOBP_|KEGG_MEDICUS_|KEGG_)", "", labels)
  labels <- gsub("_", " ", labels, fixed = TRUE)
  labels <- gsub("[\r\n]+", " ", labels)
  if (anyNA(labels) || any(!nzchar(labels))) stop("Unavailable enrichment labels.", call. = FALSE)
  wrap_width <- as.integer(cfg$reporting$enrichment_label_width %||% 40L)
  if (!is.finite(wrap_width) || wrap_width < 10L) stop("Invalid enrichment label width.", call. = FALSE)
  d$.layout_label <- rna_layout_layout_wrap(labels, wrap_width)
  d$.layout_lines <- vapply(strsplit(d$.layout_label, "\n", fixed = TRUE), length, integer(1))
  old_y <- tryCatch(rlang::eval_tidy(plot$mapping$y, d), error = function(e) NULL)
  order_top <- if (is.factor(old_y)) order(-as.integer(old_y), -seq_len(nrow(d))) else
    order(-values, -seq_len(nrow(d)))
  row_space <- d$.layout_lines[order_top] + .8
  total <- sum(row_space)
  d$.layout_y <- numeric(nrow(d))
  d$.layout_y[order_top] <- total - cumsum(row_space) + row_space / 2
  p <- plot
  p$data <- d
  p$mapping$y <- ggplot2::aes(y = .layout_y)$y
  p$layers <- lapply(p$layers, function(layer) {
    out <- local({ original <- layer; ggplot2::ggproto(NULL, original) })
    if (!is.null(layer$mapping$y)) {
      out$mapping <- layer$mapping
      out$mapping$y <- ggplot2::aes(y = .layout_y)$y
    }
    out
  })
  p$scales <- p$scales$clone()
  p$scales$scales <- Filter(function(s) !any(c("y", "size", "x") %in% s$aesthetics), p$scales$scales)
  upper <- min(1, max(.025, ceiling(max(values) * 40) / 40))
  ticks <- pretty(c(0, upper), n = 3L)
  ticks <- ticks[ticks >= 0 & ticks <= upper]
  point_minimum <- as.numeric(cfg$reporting$minimum_enrichment_point_size %||% 3)
  if (!is.finite(point_minimum) || point_minimum <= 0 || point_minimum > 7)
    stop("minimum_enrichment_point_size must be positive and at most 7.", call. = FALSE)
  p <- p + ggplot2::scale_y_continuous(limits = c(0, total), breaks = d$.layout_y,
      labels = d$.layout_label, expand = ggplot2::expansion(0)) +
    ggplot2::scale_x_continuous(limits = c(0, upper), breaks = ticks,
      labels = scales::label_percent(accuracy = if (upper < .05) .1 else 1),
      expand = ggplot2::expansion(mult = c(0, if (upper < 1) .04 else 0))) +
    ggplot2::scale_size_continuous(range = c(point_minimum, 7),
      breaks = rna_integer_legend_breaks(counts)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(lineheight = 1.05),
      axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0,
      barwidth = grid::unit(3, "in"), barheight = grid::unit(.14, "in"), order = 1),
      size = ggplot2::guide_legend(title.position = "top", nrow = 1L, order = 2))
  list(plot = p, width = if (max(d$.layout_lines) >= 3L) 10.5 else 9.2,
       height = max(5.2, 2.8 + .25 * total))
}

# Entry point called by the shared exporter after applying the package theme.
rna_prepare_plot <- function(plot, stem, cfg, width = 8, height = 6) {
  if (isTRUE(cfg$.layout_prepared) || identical(attr(plot, "IsoRemodel.layout_version"), "0.2.0"))
    return(list(plot = plot, width = width, height = height))
  if (inherits(plot, "patchwork")) {
    # Child layouts are prepared by their own constructors. Keep the composite
    # geometry intact and do not export spacers as standalone plots.
    attr(plot, "IsoRemodel.layout_version") <- "0.2.0"
    return(list(plot = plot, width = width, height = height))
  }
  if (!inherits(plot, "ggplot")) return(list(plot = plot, width = width, height = height))
  family <- rna_layout_family(stem)
  record <- list(stem = stem, family = family, width = width, height = height, plot = plot)
  dimensions <- rna_layout_dimensions(record, list(id = "compact", layout = "compact"))
  p <- rna_layout_repair_plot(plot, record, list(id = "compact", layout = "compact"))
  if (family == "enrichment") {
    output <- rna_enrichment_layout(p, cfg, dimensions[[1L]], dimensions[[2L]])
    p <- output$plot; dimensions <- c(output$width, output$height)
  }
  if (family == "filter_flow" && all(c("stage", "unit", "retained") %in% names(p$data))) {
    output <- rna_flow_plot(p$data, cfg, p$labels$title %||% "DTU testing-family construction")
    p <- output$plot; dimensions <- c(output$width, output$height)
  }
  # A fraction's valid range is known for these paired isoform examples. This
  # affects display coordinates only, preserving the original observations.
  if (family == "gene_example") {
    y <- tryCatch(rlang::eval_tidy(p$mapping$y, p$data), error = function(e) NULL)
    if (is.numeric(y) && length(y) && all(is.na(y) | is.finite(y) & y >= 0 & y <= 1)) {
      p$coordinates <- ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE)
      p <- rna_layout_scale_change(p, "y", function(s) { s$breaks <- seq(0, 1, .25); s })
    }
    # Retain room for more than two isoform panels in a full example.
    if (!inherits(p$facet, "FacetNull")) dimensions[1L] <- max(dimensions[1L], min(width, 14))
  }
  known_fraction <- family == "structural_composition" ||
    grepl("RNA_detected_isoforms_per_gene|depth_robust_structural_category_composition", stem)
  known_percentage <- family == "paired_qc" &&
    (grepl("fraction|percent", stem, ignore.case = TRUE) ||
     is.character(p$labels$y) && any(grepl("%", p$labels$y, fixed = TRUE)))
  if (known_fraction || known_percentage) {
    y <- tryCatch(rlang::eval_tidy(p$mapping$y, p$data), error = function(e) NULL)
    bound <- if (known_percentage) 100 else 1
    if (is.numeric(y) && length(y) && any(is.finite(y)) &&
        all(is.na(y) | is.finite(y) & y >= 0 & y <= bound)) {
      # Stacked bars have a visible total greater than a single row. Use their
      # built geometry to retain that full total, without changing observations.
      built <- tryCatch(ggplot2::ggplot_build(p), error = function(e) NULL)
      values <- y
      if (!is.null(built)) for (layer in built$data) {
        for (key in intersect(c("y", "ymin", "ymax"), names(layer)))
          if (is.numeric(layer[[key]])) values <- c(values, layer[[key]])
      }
      values <- values[is.finite(values)]
      if (length(values) && all(values >= -1e-8 & values <= bound + 1e-8)) {
        observed <- range(pmax(0, pmin(bound, values)))
        pad <- max(diff(observed) * .06, bound * .005)
        limits <- c(max(0, observed[1L] - pad), min(bound, observed[2L] + pad))
        p$coordinates <- if (inherits(p$coordinates, "CoordFlip"))
          ggplot2::coord_flip(ylim = limits, expand = FALSE) else
          ggplot2::coord_cartesian(ylim = limits, expand = FALSE)
      }
    }
  }
  attr(p, "IsoRemodel.layout_version") <- "0.2.0"
  list(plot = p, width = dimensions[[1L]], height = dimensions[[2L]])
}

rna_upset_plot <- function(tab, sets, cfg, title, subtitle = NULL, max_intersections = 40L) {
  if (!is.data.frame(tab) || !length(sets) || anyDuplicated(sets) || !all(sets %in% names(tab)))
    stop("UpSet requires a data frame and unique membership columns.", call. = FALSE)
  membership <- as.matrix(tab[, sets, drop = FALSE])
  if (anyNA(membership) || any(!membership %in% c(FALSE, TRUE, 0, 1)))
    stop("UpSet memberships must be observed logical values.", call. = FALSE)
  storage.mode(membership) <- "logical"
  if (!nrow(membership) || !any(rowSums(membership) > 0L))
    return(list(plot = NULL, width = 0, height = 0, counts = data.frame(),
                sizes = data.frame(set = sets, genes = rep(0L, length(sets))), components = list()))
  active <- membership[rowSums(membership) > 0, , drop = FALSE]
  keys <- apply(active, 1L, function(z) paste(as.integer(z), collapse = ""))
  all_keys <- sort(unique(keys))
  sizes <- data.frame(set = sets, genes = as.integer(colSums(active)), stringsAsFactors = FALSE)
  patterns <- do.call(rbind, lapply(all_keys, function(z) strsplit(z, "", fixed = TRUE)[[1L]] == "1"))
  totals <- as.integer(tabulate(match(keys, all_keys), nbins = length(all_keys)))
  counts <- data.frame(membership_key = all_keys,
    intersection = apply(patterns, 1L, function(z) paste(sets[z], collapse = " + ")),
    genes = totals, stringsAsFactors = FALSE)
  ordering <- order(-counts$genes, counts$membership_key)
  counts <- counts[ordering, , drop = FALSE]; patterns <- patterns[ordering, , drop = FALSE]
  counts$intersection_id <- seq_len(nrow(counts)); rownames(counts) <- NULL
  if (!identical(as.numeric(colSums(patterns * counts$genes)), as.numeric(sizes$genes)))
    stop("UpSet exact counts do not reconstruct full set totals.", call. = FALSE)
  max_intersections <- as.integer(cfg$reporting$maximum_upset_intersections %||% max_intersections)
  if (!is.finite(max_intersections) || max_intersections < 1L) stop("Invalid maximum_upset_intersections.")
  number <- min(nrow(counts), max_intersections)
  displayed <- counts[seq_len(number), , drop = FALSE]
  pats <- patterns[seq_len(number), , drop = FALSE]
  dots <- do.call(rbind, lapply(seq_len(number), function(i)
    data.frame(intersection_id = i, set = sets, member = pats[i, ], stringsAsFactors = FALSE)))
  dots$set <- factor(dots$set, levels = rev(sets))
  segments <- do.call(rbind, lapply(seq_len(number), function(i) {
    positions <- match(sets[pats[i, ]], rev(sets))
    data.frame(intersection_id = i, lower = min(positions), upper = max(positions))
  }))
  body <- cfg$reporting$font_size_pt %||% 14
  font <- cfg$reporting$font_family %||% "DejaVu Sans"
  labels <- rna_layout_layout_wrap(rna_layout_layout_labels(sets), 30L)
  x_limits <- c(.4, number + .6)
  bars <- ggplot2::ggplot(displayed, ggplot2::aes(intersection_id, genes)) +
    ggplot2::geom_col(fill = "#595959", width = .72) +
    ggplot2::geom_text(ggplot2::aes(label = genes), vjust = -.45,
      family = font, size = body / (72.27 / 25.4)) +
    ggplot2::scale_x_continuous(breaks = displayed$intersection_id, limits = x_limits,
      expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(breaks = rna_integer_legend_breaks,
      expand = ggplot2::expansion(mult = c(0, .18))) +
    ggplot2::labs(x = "Intersection", y = "Genes") + rna_theme(cfg)
  matrix <- ggplot2::ggplot(dots, ggplot2::aes(intersection_id, set)) +
    ggplot2::geom_point(colour = "#D7D7D7", size = 3.7) +
    ggplot2::geom_segment(data = segments, ggplot2::aes(x = intersection_id, xend = intersection_id,
      y = lower, yend = upper), inherit.aes = FALSE, colour = "#222222", linewidth = .7) +
    ggplot2::geom_point(data = dots[dots$member, , drop = FALSE], colour = "#222222", size = 3.7) +
    ggplot2::scale_x_continuous(breaks = displayed$intersection_id, limits = x_limits,
      expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_discrete(labels = stats::setNames(labels, sets), drop = FALSE,
      expand = ggplot2::expansion(add = .6)) +
    ggplot2::labs(x = "Intersection", y = NULL) + rna_theme(cfg) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  cap <- paste0("Each column is one exact intersection. Grey dots indicate absence.\n",
    if (number < nrow(counts)) paste0("Largest ", number, " of ", nrow(counts), " intersections shown; all counts saved. ") else "",
    "Genes in none of these sets are omitted.")
  combined <- patchwork::wrap_plots(bars + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()),
    matrix, ncol = 1L, heights = c(1.4, max(1, length(sets) / 4))) +
    patchwork::plot_annotation(title = title, subtitle = subtitle, caption = cap, theme = rna_theme(cfg))
  components <- list(count_bars = bars + ggplot2::labs(title = "Exact intersection sizes", subtitle = subtitle),
    membership = matrix + ggplot2::labs(title = "Intersection membership", subtitle = subtitle))
  attr(combined, "IsoRemodel.layout_version") <- "0.2.0"
  for (i in seq_along(components)) attr(components[[i]], "IsoRemodel.layout_version") <- "0.2.0"
  list(plot = combined, width = max(9.5, 3.5 + .6 * number),
       height = max(6.4, 4.2 + .55 * length(sets)), counts = counts, sizes = sizes,
       components = components)
}

rna_flow_plot <- function(data, cfg, title = "DTU testing-family construction") {
  if (!all(c("stage", "unit", "retained") %in% names(data)) ||
      anyNA(data[c("stage", "unit", "retained")]) || any(!is.finite(data$retained) | data$retained < 0))
    stop("Filtering-flow plot requires saved stage,unit,retained columns.", call. = FALSE)
  d <- data
  labels <- as.character(d$stage)
  labels <- gsub("\n", " ", labels, fixed = TRUE)
  labels[grepl("universe|Comparable", labels, ignore.case = TRUE)] <- "Comparable universe"
  labels[grepl("observed isoforms", labels, ignore.case = TRUE)] <- "Observed isoforms"
  labels[grepl("samples;|proportion", labels, ignore.case = TRUE)] <- "Proportion filter"
  labels[grepl("count filter", labels, ignore.case = TRUE)] <- "Count filter"
  labels[grepl("pilot", labels, ignore.case = TRUE)] <- "Pilot subset"
  d$.layout_stage <- factor(labels, levels = rev(unique(labels)))
  units <- as.character(d$unit)
  units[tolower(units) == "genes"] <- "Genes"; units[tolower(units) == "isoforms"] <- "Isoforms"
  if (any(!units %in% c("Genes", "Isoforms"))) stop("Unknown filtering-flow count unit.")
  d$.layout_unit <- factor(units, levels = c("Genes", "Isoforms"))
  d$.layout_x <- d$retained + .035 * stats::ave(d$retained, d$.layout_unit, FUN = max)
  p <- ggplot2::ggplot(d, ggplot2::aes(retained, .layout_stage, group = 1)) +
    ggplot2::geom_path(colour = "#4C78A8", linewidth = .9) +
    ggplot2::geom_point(colour = "#B2182B", size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(x = .layout_x, label = scales::comma(retained)),
      hjust = 0, size = (cfg$reporting$font_size_pt %||% 14) / (72.27 / 25.4),
      family = cfg$reporting$font_family %||% "DejaVu Sans") +
    ggplot2::facet_wrap(~.layout_unit, scales = "free_x", ncol = 2L) +
    ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(.02, .45))) +
    ggplot2::coord_cartesian(clip = "off") + rna_theme(cfg) +
    ggplot2::labs(title = title, x = "Retained features", y = NULL) +
    ggplot2::theme(panel.spacing.x = grid::unit(24, "pt"), plot.margin = ggplot2::margin(14, 30, 14, 16))
  components <- list()
  for (unit in unique(units)) {
    one <- p; one$data <- d[units == unit, , drop = FALSE]
    one <- one + ggplot2::facet_null() + ggplot2::labs(title = paste("DTU filtering:", tolower(unit)))
    attr(one, "IsoRemodel.layout_version") <- "0.2.0"
    components[[tolower(unit)]] <- one
  }
  attr(p, "IsoRemodel.layout_version") <- "0.2.0"
  list(plot = p, width = 11.2, height = max(5.2, 2 + .8 * length(unique(labels))), components = components)
}
