# Shared utilities for the standalone RNA workflow.
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

rna_stop <- function(...) stop(..., call. = FALSE)
rna_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) rna_stop("Missing R packages: ", paste(missing, collapse = ", "),
    ". Run scripts/install_dependencies.R first.")
  invisible(TRUE)
}
rna_package_version <- function(root = getOption("IsoRemodel.code_root")) {
  file <- if (!is.null(root)) file.path(root, "DESCRIPTION") else ""
  if (nzchar(file) && file.exists(file)) {
    value <- tryCatch(read.dcf(file)[1L, "Version"], error = function(e) NULL)
    if (!is.null(value)) return(as.character(value))
  }
  tryCatch(as.character(utils::packageVersion("IsoRemodel")), error = function(e) "0.2.1")
}

rna_banner <- function(log = FALSE, startup = FALSE) {
  lines <- c(paste0("IsoRemodel v", rna_package_version()),
    "An R workflow for paired and unpaired long-read isoform analysis",
    "Pigeon import | QC | Depth-controlled gain/loss | Remodeling | DTU")
  if (isTRUE(log)) for (line in lines) rna_log(line)
  else if (isTRUE(startup)) packageStartupMessage(paste(lines, collapse = "\n"))
  else message(paste(lines, collapse = "\n"))
  invisible(lines)
}

rna_log <- function(..., level = "INFO") {
  value <- paste0(..., collapse = "")
  value <- gsub("[\r\n]+", " | ", value)
  level <- toupper(as.character(level)[1L])
  line <- paste0("[IsoRemodel] [", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
                 "] [", level, "] ", value)
  message(line)
  file <- getOption("IsoRemodel.log_file")
  if (!is.null(file)) cat(line, "\n", file = file, append = TRUE)
  invisible(line)
}

# Outer-stage progress counts stages that passed checks, not time remaining.
# Detailed resample/gene counters continue to be emitted by the numerical engines.
rna_progress <- function(current, total, stage, status = "running", detail = "", stage_index = NULL) {
  stopifnot(length(current) == 1L, length(total) == 1L,
    is.finite(current), is.finite(total), current >= 0, total >= 1, current <= total)
  started <- getOption("IsoRemodel.run_start")
  elapsed <- if (!is.null(started)) max(0, as.numeric(difftime(Sys.time(), started, units = "secs"))) else 0
  index <- stage_index %||% if (status %in% c("starting", "running", "failed", "reused")) min(total, current + 1L) else current
  entry <- data.frame(timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    run_id = getOption("IsoRemodel.run_id", ""), stage = as.character(stage),
    stage_index = as.integer(index), stage_total = as.integer(total),
    stages_passed = as.integer(current), stages_percent = round(100 * current / total, 1),
    status = as.character(status), elapsed_seconds = round(elapsed, 2),
    detail = gsub("[\r\n]+", " | ", as.character(detail)), stringsAsFactors = FALSE)
  file <- getOption("IsoRemodel.progress_file")
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    append <- file.exists(file) && file.info(file)$size > 0
    utils::write.table(entry, file, sep = "\t", row.names = FALSE, quote = TRUE,
                       col.names = !append, append = append, na = "NA")
  }
  latest <- getOption("IsoRemodel.progress_latest_file")
  if (!is.null(latest)) {
    dir.create(dirname(latest), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(entry, latest, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")
  }
  rna_log("Progress: ", current, "/", total, " stages passed (", entry$stages_percent,
    "%) | ", stage, " | ", status, " | elapsed ", round(elapsed), "s",
    if (nzchar(detail)) paste0(" | ", detail) else "",
    level = if (identical(status, "failed")) "ERROR" else "INFO")
  invisible(entry)
}

rna_dirs <- function(stage_dir) {
  ans <- list(root = stage_dir)
  for (name in c("tables", "figures", "qc", "objects", "cache")) ans[[name]] <- file.path(stage_dir, name)
  invisible(lapply(ans, dir.create, recursive = TRUE, showWarnings = FALSE))
  ans
}
rna_write_tsv <- function(x, path) {
  rna_require("readr")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(x), path, na = "NA")
  invisible(path)
}
rna_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(paste0(basename(path), ".tmp_"), dirname(path))
  backup <- paste0(path, ".previous")
  on.exit(unlink(temp), add = TRUE)
  saveRDS(object, temp, compress = "gzip")
  if (file.exists(backup)) unlink(backup)
  previous <- file.exists(path)
  if (previous && !file.rename(path, backup)) rna_stop("Cannot protect previous checkpoint: ", path)
  if (!file.rename(temp, path)) {
    if (previous && file.exists(backup)) file.rename(backup, path)
    rna_stop("Cannot finalize checkpoint: ", path)
  }
  if (file.exists(backup)) unlink(backup)
  invisible(path)
}
rna_hash <- function(object) {
  rna_require("digest")
  digest::digest(object, algo = "sha256")
}
rna_stage_seed <- function(master_seed, ...) {
  rna_require("digest")
  key <- paste(c(as.character(master_seed), unlist(list(...))), collapse = "|")
  as.integer((strtoi(substr(digest::digest(key, algo = "xxhash32"), 1L, 7L), 16L) %% 2147483646) + 1)
}
rna_scientific_config <- function(cfg) {
  out <- cfg[!startsWith(names(cfg), ".")]
  out$reporting <- out$output <- NULL
  out$execution$workers <- out$execution$resume <- out$execution$force_refit <- NULL
  out$execution$checkpoint_every_genes <- NULL
  for (key in c("run_upset", "write_supplement_workbook", "plot_all_dtu_genes",
                "maximum_gene_plots", "examples_per_class", "isoforms_per_gene_example",
                "priority_example_genes", "maximum_heatmap_genes", "force_refit")) out$dtu[[key]] <- NULL
  out$gain_loss$checkpoint_every <- NULL
  out
}
rna_font_directory <- function(cfg = NULL) {
  root <- cfg$.code_root %||% getOption("IsoRemodel.code_root") %||% getwd()
  directory <- file.path(root, "inst", "fonts")
  if (!file.exists(file.path(directory, "DejaVuSans.ttf")))
    directory <- system.file("fonts", package = "IsoRemodel")
  directory
}

rna_register_fonts <- function(cfg = NULL) {
  rna_require(c("systemfonts", "sysfonts"))
  family <- cfg$reporting$font_family %||% "DejaVu Sans"
  if (!identical(family, "DejaVu Sans")) {
    if (!family %in% systemfonts::system_fonts()$family)
      rna_stop("Configured font is not installed: ", family, ". Use bundled DejaVu Sans or install that font.")
    match <- systemfonts::match_fonts(family, weight = c("normal", "bold"))
    regular <- match$path[1L]; bold <- match$path[2L]
  } else {
    directory <- rna_font_directory(cfg)
    regular <- file.path(directory, "DejaVuSans.ttf")
    bold <- file.path(directory, "DejaVuSans-Bold.ttf")
    if (!all(file.exists(c(regular, bold)))) rna_stop("Bundled DejaVu Sans font files are missing.")
  }
  # A system-installed family cannot be registered again as a custom family.
  # Keep font setup idempotent across hundreds of exports and repeated calls.
  registered <- if ("registry_fonts" %in% getNamespaceExports("systemfonts"))
    systemfonts::registry_fonts()$family else character()
  available <- unique(c(systemfonts::system_fonts()$family, registered))
  if (!family %in% available)
    systemfonts::register_font(family, plain = regular, bold = bold)
  sysfonts::font_add(family, regular = regular, bold = bold)
  invisible(family)
}

rna_theme <- function(cfg = NULL) {
  rna_require("ggplot2")
  settings <- cfg$reporting %||% list()
  family <- settings$font_family %||% "DejaVu Sans"
  body <- settings$font_size_pt %||% 14
  heading <- settings$title_size_pt %||% 16
  ggplot2::theme_classic(base_size = body, base_family = family) +
    ggplot2::theme(plot.title.position = "plot",
      plot.title = ggplot2::element_text(size = heading, face = "bold", hjust = 0,
                                        margin = ggplot2::margin(b = 7)),
      plot.subtitle = ggplot2::element_text(size = settings$subtitle_size_pt %||% body, hjust = 0,
                                           margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(size = settings$axis_title_size_pt %||% 16),
      axis.text = ggplot2::element_text(size = settings$axis_text_size_pt %||% body, color = "black"),
      strip.text = ggplot2::element_text(size = body, face = "bold"),
      legend.title = ggplot2::element_text(size = body, face = "bold"),
      legend.text = ggplot2::element_text(size = settings$legend_text_size_pt %||% body),
      plot.caption = ggplot2::element_text(size = settings$caption_size_pt %||% body, hjust = 0),
      legend.position = "bottom", legend.box = "vertical",
      plot.margin = ggplot2::margin(10, 14, 10, 10),
      panel.grid = ggplot2::element_blank())
}

rna_clean_plot_text <- function(x) {
  if (!is.character(x)) return(x)
  x <- gsub("[Mm]odule[[:space:]]+0*[0-9]+[[:space:]]*[-:–—]?[[:space:]]*", "", x)
  x <- gsub("[Ee]xploratory[[:space:]]*", "", x)
  x <- gsub("\\bBH05\\b", "FDR <0.05", x)
  x <- gsub("\\bBH10\\b", "FDR <0.10", x)
  gsub("\\bBH\\b", "FDR", x)
}

rna_wrap_plot_text <- function(x, width) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(x)
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1L]], function(line)
    paste(strwrap(line, width = width), collapse = "\n"), character(1)), collapse = "\n")
}

rna_export_device <- function(plot, path, format, width, height, cfg) {
  dpi <- cfg$reporting$dpi %||% 600
  family <- cfg$reporting$font_family %||% "DejaVu Sans"
  outlined <- format %in% c("pdf", "svg") && !identical(cfg$reporting$outline_vector_fonts, FALSE)
  if (format == "png") {
    rna_require("ragg"); ragg::agg_png(path, width = width, height = height, units = "in", res = dpi, background = "white")
  } else if (format %in% c("jpeg", "jpg")) {
    rna_require("ragg"); ragg::agg_jpeg(path, width = width, height = height, units = "in", res = dpi, quality = 95, background = "white")
  } else if (format %in% c("svg", "editable.svg")) {
    rna_require("svglite"); svglite::svglite(path, width = width, height = height, bg = "white")
  } else if (format == "pdf") {
    if (outlined) grDevices::pdf(path, width = width, height = height, family = "sans", useDingbats = FALSE)
    else if (capabilities("cairo")) grDevices::cairo_pdf(path, width = width, height = height, family = family)
    else rna_stop("Text-based PDF requires Cairo; enable outline_vector_fonts for portable PDF export.")
  } else rna_stop("Unsupported figure format: ", format)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (outlined) {
    rna_require("showtext"); showtext::showtext_begin()
    on.exit(showtext::showtext_end(), add = TRUE, after = FALSE)
  }
  if (inherits(plot, c("ggplot", "patchwork"))) print(plot)
  else { grid::grid.newpage(); grid::grid.draw(plot) }
  invisible(path)
}

rna_style_plot <- function(plot, cfg, width) {
  if (inherits(plot, "patchwork")) {
    columns <- plot$patches$layout$ncol %||% ceiling(sqrt(length(plot)))
    columns <- max(1, as.integer(columns))
    for (i in seq_len(length(plot))) plot[[i]] <- rna_style_plot(plot[[i]], cfg, width / columns)
    plot <- plot & rna_theme(cfg)
  } else if (inherits(plot, "ggplot")) {
    for (name in names(plot$labels)) plot$labels[name] <- list(rna_clean_plot_text(plot$labels[[name]]))
    plot$labels$title <- rna_wrap_plot_text(plot$labels$title, max(25, floor(width * 7.5)))
    plot$labels$subtitle <- rna_wrap_plot_text(plot$labels$subtitle, max(30, floor(width * 8.5)))
    plot$layers <- lapply(plot$layers, function(layer) {
      if (!any(grepl("Text|Label", class(layer$geom)))) return(layer)
      revised <- ggplot2::ggproto(NULL, layer)
      revised$aes_params <- as.list(layer$aes_params)
      revised$aes_params$family <- cfg$reporting$font_family %||% "DejaVu Sans"
      if (is.null(layer$mapping$size) && is.null(plot$mapping$size))
        revised$aes_params$size <- (cfg$reporting$font_size_pt %||% 14) / (72.27 / 25.4)
      revised
    })
    plot <- plot + rna_theme(cfg)
  }
  plot
}

rna_save_plot <- function(plot, stem, dirs, cfg, width = 8, height = 6) {
  if (identical(cfg$reporting$enabled, FALSE)) return(invisible(character()))
  rna_require("ggplot2")
  rna_register_fonts(cfg)
  plot <- rna_style_plot(plot, cfg, width)
  if (!isTRUE(cfg$.layout_prepared)) {
    prepared <- rna_prepare_plot(plot, stem, cfg, width, height)
    plot <- prepared$plot; width <- prepared$width; height <- prepared$height
  }
  folder <- if (is.list(dirs)) dirs$figures else dirs
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  safe_stem <- gsub("[^A-Za-z0-9._-]", "_", stem)
  prefix <- file.path(folder, safe_stem)
  if (!identical(cfg$reporting$save_plot_objects, FALSE))
    rna_atomic_rds(list(plot = plot, width = width, height = height), paste0(prefix, ".plot.rds"))
  formats <- cfg$reporting$formats %||% c("pdf", "svg", "png", "jpeg")
  if (isTRUE(cfg$reporting$editable_svg)) formats <- unique(c(formats, "editable.svg"))
  paths <- vapply(formats, function(format)
    rna_export_device(plot, paste0(prefix, ".", format), format, width, height, cfg), character(1))
  if (inherits(plot, "patchwork") && isTRUE(cfg$reporting$separate_panels_for_composites)) {
    for (i in seq_len(length(plot))) {
      if (!rna_plot_meaningful(plot[[i]])) next
      subcfg <- cfg; subcfg$.layout_prepared <- TRUE
      subcfg$reporting$separate_panels_for_composites <- FALSE
      rna_save_plot(plot[[i]], paste0(stem, "__component_", i), dirs, subcfg,
                    min(width, 8), min(height, 7))
    }
  }
  # Every ordinary facet also has an independently exportable figure. Clone
  # only the layer data; never install recursive wrappers around ggproto methods.
  if (inherits(plot, "ggplot") && !inherits(plot, "patchwork") &&
      isTRUE(cfg$reporting$separate_panels_for_composites) && !inherits(plot$facet, "FacetNull")) {
    vars <- unique(c(names(plot$facet$params$facets), names(plot$facet$params$rows), names(plot$facet$params$cols)))
    if (length(vars) && is.data.frame(plot$data) && all(vars %in% names(plot$data))) {
      panels <- unique(plot$data[, vars, drop = FALSE])
      for (i in seq_len(nrow(panels))) {
        keep_rows <- function(data) {
          if (!is.data.frame(data)) return(data)
          present <- intersect(vars, names(data))
          if (!length(present)) return(data)
          keep <- rep(TRUE, nrow(data))
          for (v in present) keep <- keep & if (is.na(panels[[v]][i])) is.na(data[[v]]) else
            !is.na(data[[v]]) & as.character(data[[v]]) == as.character(panels[[v]][i])
          data[keep, , drop = FALSE]
        }
        piece <- plot; piece$data <- keep_rows(piece$data)
        piece$layers <- lapply(piece$layers, function(layer) {
          new <- ggplot2::ggproto(NULL, layer); new$data <- keep_rows(layer$data); new
        })
        piece <- piece + ggplot2::facet_null()
        description <- vapply(vars, function(v) paste(v, as.character(panels[[v]][i]), sep = ": "), character(1))
        piece$labels$subtitle <- rna_clean_plot_text(paste(description, collapse = "; "))
        subcfg <- cfg; subcfg$reporting$separate_panels_for_composites <- FALSE
        subcfg$.layout_prepared <- TRUE
        rna_save_plot(piece, paste0(stem, "__panel_", i), dirs, subcfg, min(width, 8), min(height, 7))
      }
    }
  }
  invisible(paths)
}

# Small bounded graphics check, useful before starting an expensive run.
rna_plot_smoke_test <- function(output_root, cfg = list()) {
  rna_require("ggplot2")
  cfg$reporting$enabled <- TRUE
  cfg$reporting$editable_svg <- TRUE
  cfg$reporting$separate_panels_for_composites <- TRUE
  cfg$reporting$formats <- c("pdf", "svg", "png", "jpeg")
  palette <- cfg$reporting$palettes$timepoints %||% cfg$reporting$palettes$timepoint %||% c(primary = "#4C78A8", recurrent = "#E45756")
  d <- expand.grid(subject = seq_len(4), condition = c("primary", "recurrent"), panel = c("Counts", "Fractions"))
  d$value <- seq_len(nrow(d)) / 3
  p <- ggplot2::ggplot(d, ggplot2::aes(condition, value, color = condition, group = subject)) +
    ggplot2::geom_line(color = "grey65", linewidth = .5) + ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(label = subject), vjust = -1.0, show.legend = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.1, .18))) +
    ggplot2::scale_color_manual(values = palette) + ggplot2::facet_wrap(~panel) +
    ggplot2::labs(title = "RNA figure typography and export check", x = NULL, y = "Value", color = "Condition")
  rna_save_plot(p, "typography_check", output_root, cfg, width = 7.2, height = 4.8)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    one <- p; one$data <- one$data[one$data$panel == "Counts", , drop = FALSE]
    one <- one + ggplot2::facet_null()
    two <- p; two$data <- two$data[two$data$panel == "Fractions", , drop = FALSE]
    two <- two + ggplot2::facet_null()
    rna_save_plot(patchwork::wrap_plots(one, two), "composite_check", output_root, cfg, 10, 4.8)
  }
  invisible(output_root)
}

rna_engine_env <- function(packages, parent = parent.frame()) {
  rna_require(packages)
  env <- new.env(parent = parent)
  for (pkg in packages) {
    exports <- getNamespaceExports(pkg)
    for (name in exports) if (!exists(name, envir = env, inherits = FALSE))
      assign(name, getExportedValue(pkg, name), envir = env)
  }
  env
}
