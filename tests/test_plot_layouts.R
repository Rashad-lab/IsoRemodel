# Approved layouts must preserve stored values and remain usable at small sizes.
local({
  if (!all(vapply(c("ggplot2", "scales", "rlang", "patchwork"), requireNamespace,
                  logical(1), quietly = TRUE))) {
    message("SKIP layout tests: plotting dependencies unavailable")
  } else {
    cfg <- list(reporting = list(font_family = "DejaVu Sans", font_size_pt = 14,
      enrichment_label_width = 40, minimum_enrichment_point_size = 3))
    d <- data.frame(pathway = c("KEGG_MEDICUS_A_VERY_LONG_PATHWAY_NAME_WITH_MULTIPLE_COMPONENTS_AND_EXTENDED_BIOLOGICAL_DESCRIPTION", "KEGG_SMALL"),
      selected_fraction = c(.08, .1), plotted_genes = c(4L, 5L), plot_score = c(2, 3),
      plotted_p = c(.01, .001), plotted_bh = c(.3, .3), display_tier = "raw_P_below_0p05_fallback")
    p <- ggplot2::ggplot(d, ggplot2::aes(selected_fraction, pathway, size = plotted_genes, colour = plot_score)) +
      ggplot2::geom_point() + ggplot2::labs(caption = "Nominal P < 0.05; none meet pathway FDR < 0.1")
    a <- rna_prepare_plot(p, "matched_ORA_FDR05_KEGG", cfg, 8, 5)
    stopifnot(identical(d, a$plot$data[names(d)]), a$height > 5,
      all(a$plot$scales$get_scales("size")$breaks == c(4, 5)),
      identical(a$plot$scales$get_scales("x")$limits[1L], 0),
      min(ggplot2::ggplot_build(a$plot)$data[[1L]]$size) >= 3,
      grepl("Nominal", a$plot$labels$caption, fixed = TRUE))
    again <- rna_prepare_plot(a$plot, "matched_ORA_FDR05_KEGG", cfg, a$width, a$height)
    stopifnot(identical(again, a))
    # The remodeling schema uses different column names; it shares the repairs.
    d2 <- data.frame(Description = c("HALLMARK_ALPHA", "HALLMARK_BETA"),
      enrichment_fraction = c(.08, .1), Count = c(4L, 5L), score = c(2, 3), FDR = c(.3, .3))
    p2 <- ggplot2::ggplot(d2, ggplot2::aes(enrichment_fraction, Description, size = Count, colour = score)) + ggplot2::geom_point()
    b <- rna_prepare_plot(p2, "ORA_MSigDB_Hallmark", cfg, 8, 5)
    stopifnot(identical(d2, b$plot$data[names(d2)]),
      identical(b$plot$scales$get_scales("x")$limits[1L], 0),
      identical(as.numeric(b$plot$scales$get_scales("size")$breaks), c(4, 5)))
    tab <- data.frame(gene_id = letters[1:5], check.names = FALSE)
    tab[["A + B"]] <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
    tab[["C"]] <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
    tab[["Zero set"]] <- FALSE
    upset <- rna_upset_plot(tab, c("A + B", "C", "Zero set"), cfg,
      "Test exact intersections", max_intersections = 2L)
    stopifnot(sum(upset$counts$genes) == 3L, nrow(upset$counts) == 3L,
      identical(upset$sizes$genes, c(2L, 2L, 0L)),
      all(vapply(upset$components, rna_plot_meaningful, logical(1))),
      grepl("Largest 2 of 3", upset$plot$patches$annotation$caption, fixed = TRUE))
    empty <- tab; empty[["A + B"]] <- FALSE; empty$C <- FALSE
    stopifnot(is.null(rna_upset_plot(empty, c("A + B", "C"), cfg, "Empty")$plot))
    invalid <- tab; invalid$C[1L] <- NA
    stopifnot(inherits(try(rna_upset_plot(invalid, c("A + B", "C"), cfg, "Invalid"), silent = TRUE), "try-error"))
    flow <- data.frame(stage = rep(c("Comparable universe", "Count filter"), 2),
      unit = rep(c("Genes", "Isoforms"), each = 2), retained = c(20606L, 6826L, 583610L, 43674L))
    f <- rna_flow_plot(flow, cfg)
    stopifnot(identical(flow, f$plot$data[names(flow)]), length(f$components) == 2,
      f$width >= 11, all(vapply(f$components, rna_plot_meaningful, logical(1))))
    stopifnot(!rna_plot_meaningful(patchwork::plot_spacer()),
      !rna_plot_meaningful(ggplot2::ggplot() + ggplot2::geom_blank()),
      rna_plot_meaningful(patchwork::wrap_plots(patchwork::plot_spacer(), p)))
    pct <- data.frame(condition = c("primary", "recurrent"), plotted_value = c(0, .02))
    pp <- ggplot2::ggplot(pct, ggplot2::aes(condition, plotted_value)) + ggplot2::geom_point() +
      ggplot2::labs(y = "Annotated reads (%)")
    fixed <- rna_prepare_plot(pp, "paired_tss_within_threshold_FL_fraction", cfg, 6, 4.5)
    stopifnot(identical(fixed$plot$data, pct), fixed$plot$coordinates$limits$y[1L] >= 0,
      fixed$plot$coordinates$limits$y[2L] <= 100)
    comp <- data.frame(sample = "ptn1_primary", structure = c("FSM", "ISM"), fraction = c(.7, .3))
    cp <- ggplot2::ggplot(comp, ggplot2::aes(sample, fraction, fill = structure)) + ggplot2::geom_col() +
      ggplot2::scale_y_continuous(labels = scales::percent)
    fixed <- rna_prepare_plot(cp, "RNA_structural_composition_support_filtered_FL_fraction", cfg, 8, 6)
    stopifnot(identical(fixed$plot$data, comp),
      identical(fixed$plot$coordinates$limits$y, c(0, 1)))
  }
})
message("Plot layout data-preservation and geometry checks passed.")
