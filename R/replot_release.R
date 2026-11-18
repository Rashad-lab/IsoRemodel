# Release replotting reads saved statistics only. It never calls run_dtu(),
# run_dtu_enrichment(), stageR, or another analytical-stage entry point.
rna_replot_read_tsv <- function(path) {
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, check.names = FALSE, quote = "\"", comment.char = "",
                    stringsAsFactors = FALSE)
}

rna_replot_validate_dtu <- function(dtu_source_dir) {
  path <- file.path(dtu_source_dir, "objects", "DTU_complete.rds")
  if (!file.exists(path)) rna_stop("Completed DTU results are missing: ", path)
  source <- readRDS(path)
  if (!is.list(source) || !isTRUE(source$complete) || !is.list(source$cohort_results) ||
      !length(source$cohort_results) || is.null(names(source$cohort_results)) ||
      anyNA(names(source$cohort_results)) || any(!nzchar(names(source$cohort_results))) ||
      anyDuplicated(names(source$cohort_results)))
    rna_stop("DTU source is not a completed, named cohort-result handoff.")
  # Saved results produced before the independent-group branch are paired.
  design <- source$analysis_design %||% "paired"
  if (length(design) != 1L || is.na(design) || !design %in% c("paired", "unpaired"))
    rna_stop("DTU source has an invalid analysis design.")
  source$analysis_design <- design
  contract <- source$rng_contract
  required <- c("cohort", "seed", "kind", "normal_kind", "sample_kind", "contract_version")
  if (!is.data.frame(contract) || !all(required %in% names(contract)) ||
      !nrow(contract) || anyNA(contract[required]) || anyDuplicated(contract$cohort) ||
      !is.numeric(contract$seed) || any(!is.finite(contract$seed) | contract$seed < 1 |
                                      contract$seed != round(contract$seed)) ||
      !setequal(contract$cohort, names(source$cohort_results)) ||
      any(contract$kind != "L'Ecuyer-CMRG") || any(contract$normal_kind != "Inversion") ||
      any(contract$sample_kind != "Rejection") || any(contract$contract_version != "dtu-rng-v1"))
    rna_stop("DTU source lacks the explicit corrected L'Ecuyer-CMRG RNG contract.")
  status_path <- file.path(dirname(dtu_source_dir), "refit_status.rds")
  manifest_path <- file.path(dirname(dtu_source_dir), "latest_run_manifest.rds")
  if (file.exists(status_path)) {
    if (!identical(design, "paired"))
      rna_stop("A paired confirmation-refit status cannot validate unpaired DTU results.")
    status <- readRDS(status_path)
    if (!isTRUE(status$completed) || !isTRUE(status$original_results_reproduced) ||
        !identical(status$status, "original_results_reproduced"))
      rna_stop("The supplied DTU confirmation refit did not reproduce the original results.")
    expected <- dtu_rng_contract(list(kind = "L'Ecuyer-CMRG", normal_kind = "Inversion",
      sample_kind = "Rejection"), 260605L, names(source$cohort_results))
    ordering <- match(expected$cohort, contract$cohort)
    if (!identical(as.integer(contract$seed[ordering]), as.integer(expected$seed)))
      rna_stop("Confirmed DTU source does not have the original cohort seeds.")
    provenance <- status_path
  } else if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    version_ok <- tryCatch(utils::compareVersion(manifest$package_version, "0.2.0") >= 0,
                           error = function(e) FALSE)
    if (!isTRUE(manifest$completed) || !isTRUE(version_ok) ||
        !identical(manifest$status, "complete") ||
        !isTRUE(manifest$stage_records$dtu$complete))
      rna_stop("Native DTU replotting requires a completed IsoRemodel v0.2.0-or-newer run.")
    manifest_design <- manifest$analysis_design %||% manifest$design
    if (!is.null(manifest_design) && !identical(manifest_design, design))
      rna_stop("The native run manifest and DTU result have different analysis designs.")
    provenance <- manifest_path
  } else rna_stop("DTU source has neither a verified refit status nor a completed v0.2 run manifest.")
  if (!is.data.frame(source$gates) || !nrow(source$gates) ||
      !"pass" %in% names(source$gates) || !all(source$gates$pass %in% TRUE))
    rna_stop("DTU completion gates are missing or did not all pass.")
  list(result = source, paths = c(path, provenance), contract = contract, analysis_design = design)
}

rna_replot_review <- function(index, output_dir, volume_size = 60L) {
  valid <- index[index$status == "complete" & file.exists(index$preview_png), , drop = FALSE]
  if (!nrow(valid)) return(character())
  folder <- file.path(output_dir, "review_volumes")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  groups <- split(seq_len(nrow(valid)), ceiling(seq_len(nrow(valid)) / volume_size))
  paths <- character()
  for (v in seq_along(groups)) {
    path <- file.path(folder, sprintf("IsoRemodel_figures_review_volume_%03d.pdf", v))
    grDevices::pdf(path, width = 13, height = 10, useDingbats = FALSE)
    tryCatch({
      for (i in groups[[v]]) {
        row <- valid[i, , drop = FALSE]
        grid::grid.newpage()
        header <- paste0(row$figure_id, " | ", row$relative_stem)
        if (nchar(header) > 100L) {
          starts <- seq.int(1L, nchar(header), by = 100L)
          header <- paste(substring(header, starts, pmin(starts + 99L, nchar(header))), collapse = "\n")
        }
        grid::grid.text(header,
          x = .025, y = .975, just = c("left", "top"),
          gp = grid::gpar(fontsize = 11, fontface = "bold"))
        im <- png::readPNG(row$preview_png)
        ar <- dim(im)[2] / dim(im)[1]
        area <- (13 * .94) / (10 * .85)
        if (ar > area) { w <- .94; h <- .85 * area / ar }
        else { h <- .85; w <- .94 * ar / area }
        grid::grid.raster(im, x = .5, y = .49, width = grid::unit(w, "npc"),
                          height = grid::unit(h, "npc"), interpolate = FALSE)
        grid::grid.text(paste0("IsoRemodel v0.2.1 | saved results | page ", i,
          "/", nrow(valid)), x = .025, y = .015, just = "left", gp = grid::gpar(fontsize = 9))
      }
    }, finally = grDevices::dev.off())
    paths <- c(paths, path)
  }
  paths
}

rna_replot_make_preview <- function(path, dest, dpi) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  # nativeRaster stores each RGBA pixel in one integer. A floating-point RGBA
  # array can take eight times as much RAM for the same 600-dpi figure.
  im <- png::readPNG(path, native = TRUE)
  width <- dim(im)[2] / dpi; height <- dim(im)[1] / dpi
  ny <- min(dim(im)[1], max(1L, ceiling(height * 120)))
  nx <- min(dim(im)[2], max(1L, ceiling(width * 120)))
  yy <- unique(round(seq(1, dim(im)[1], length.out = ny)))
  xx <- unique(round(seq(1, dim(im)[2], length.out = nx)))
  small <- im[yy, xx, drop = FALSE]
  class(small) <- class(im)
  png::writePNG(small, dest)
  invisible(c(width = width, height = height))
}

#' Replot saved RNA results using corrected DTU and the v0.2 layouts.
replot_isoremodel <- function(config_file = NULL,
    source_dir = "F:/IsoRemodel/outputs",
    dtu_source_dir = "F:/IsoRemodel/DTU_refit/20261026T043514_45632/M06_dtu",
    output_dir = NULL, formats = NULL, dpi = NULL, strict = TRUE) {
  rna_require(c("ggplot2", "patchwork", "scales", "yaml", "digest", "readr", "dplyr",
    "tidyr", "tibble", "purrr", "forcats", "stringr", "png", "zip", "ragg", "svglite",
    "showtext", "sysfonts", "systemfonts"))
  if (is.null(config_file)) {
    root <- getOption("IsoRemodel.code_root")
    config_file <- if (!is.null(root)) file.path(root, "config", "analysis.yaml") else
      system.file("config", "analysis.yaml", package = "IsoRemodel")
  }
  cfg <- read_rna_config(config_file, validate_inputs = FALSE)
  source_dir <- rna_run_absolute(source_dir)
  dtu_source_dir <- rna_run_absolute(dtu_source_dir)
  if (!dir.exists(source_dir) || !dir.exists(dtu_source_dir))
    rna_stop("Both source_dir and dtu_source_dir must exist.")
  inside <- function(a, b) identical(tolower(a), tolower(b)) ||
    startsWith(tolower(a), paste0(tolower(b), "/"))
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "_", Sys.getpid())
  if (is.null(output_dir)) output_dir <- file.path(dirname(source_dir), "replot_v0.2.1", run_id)
  output_dir <- rna_run_absolute(output_dir)
  if (inside(output_dir, source_dir) || inside(output_dir, dtu_source_dir) ||
      inside(source_dir, output_dir) || inside(dtu_source_dir, output_dir))
    rna_stop("Replot output must be separate from the source result trees.")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    rna_stop("Replot output directory must be new or empty: ", output_dir)
  dtu <- rna_replot_validate_dtu(dtu_source_dir)
  paired_design <- isTRUE(cfg$comparison$paired)
  if (!identical(if (paired_design) "paired" else "unpaired", dtu$analysis_design))
    rna_stop("The replot configuration and saved DTU result have different analysis designs.")
  cfg$reporting$enabled <- TRUE
  cfg$reporting$save_plot_objects <- FALSE
  cfg$reporting$separate_panels_for_composites <- TRUE
  if (!is.null(formats)) cfg$reporting$formats <- as.character(formats)
  if (!is.null(dpi)) cfg$reporting$dpi <- as.numeric(dpi)
  if (any(!cfg$reporting$formats %in% c("pdf", "svg", "png", "jpeg", "jpg")))
    rna_stop("Replot formats must be PDF, SVG, PNG or JPEG.")
  if (length(cfg$reporting$dpi) != 1L || !is.finite(cfg$reporting$dpi) || cfg$reporting$dpi < 72)
    rna_stop("Replot dpi must be at least 72.")
  # A PNG is also retained for every figure's visual review, regardless of the
  # requested vector formats. No statistical cache or fitted object is written.
  cfg$reporting$formats <- unique(c(as.character(cfg$reporting$formats), "png"))
  cfg$dtu$write_supplement_workbook <- FALSE
  option_names <- c("IsoRemodel.log_file", "IsoRemodel.progress_file",
    "IsoRemodel.progress_latest_file", "IsoRemodel.run_id", "IsoRemodel.run_start")
  old_options <- stats::setNames(lapply(option_names, getOption), option_names)
  old_rng <- dtu_rng_snapshot()
  on.exit({options(old_options); dtu_rng_restore(old_rng)}, add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  options(IsoRemodel.log_file = file.path(output_dir, "run.log"),
    IsoRemodel.progress_file = file.path(output_dir, "progress.tsv"),
    IsoRemodel.progress_latest_file = file.path(output_dir, "progress_latest.tsv"),
    IsoRemodel.run_id = run_id, IsoRemodel.run_start = Sys.time())
  rna_banner(log = TRUE)
  rna_log("Replot only; source: ", source_dir, "; corrected DTU: ", dtu_source_dir)
  rna_log("No model fitting, stageR adjustment, enrichment testing or other analytical stage will run.")
  saved_paths <- sort(list.files(source_dir, pattern = "[.]plot[.]rds$", recursive = TRUE, full.names = TRUE))
  relative <- substring(saved_paths, nchar(source_dir) + 2L)
  keep <- !grepl("(^|/)M06_dtu(/|$)", relative, ignore.case = TRUE)
  stale <- saved_paths[!keep]; saved_paths <- saved_paths[keep]
  # Parents are rebuilt with meaningful standalone panels. Do not reimport old
  # children (including blank UpSet spacers) when their saved parent exists.
  stem_paths <- sub("[.]plot[.]rds$", "", saved_paths)
  parent <- sub("__(component|panel)_[0-9]+$", "", stem_paths)
  keep <- parent == stem_paths | !parent %in% stem_paths
  skipped_children <- saved_paths[!keep]; saved_paths <- saved_paths[keep]
  remodeling_path <- file.path(source_dir, "M04_remodeling", "cache", "stage_result.rds")
  if (!file.exists(remodeling_path)) rna_stop("Full saved remodeling cache is required for DTU overlap: ", remodeling_path)
  remodeling_cache <- readRDS(remodeling_path)
  remodeling <- remodeling_cache$result
  if (!isTRUE(remodeling_cache$complete) || !is.data.frame(remodeling$gene_table))
    rna_stop("The saved remodeling handoff is not complete.")
  ora_paths <- sort(list.files(dtu_source_dir, pattern = "_ORA_FDR[0-9]+_.*[.]tsv[.]gz$",
    recursive = TRUE, full.names = TRUE))
  sources <- unique(c(cfg$.config_file, saved_paths, dtu$paths, remodeling_path, ora_paths))
  before <- rna_run_file_manifest(sources)
  rna_write_tsv(before, file.path(output_dir, "source_files_before.tsv"))
  rna_write_tsv(dtu$contract, file.path(output_dir, "DTU_rng_contract.tsv"))
  yaml::write_yaml(cfg$reporting, file.path(output_dir, "reporting_settings.yaml"))
  yaml::write_yaml(cfg[!startsWith(names(cfg), ".")], file.path(output_dir, "resolved_replot_config.yaml"))
  rna_write_tsv(data.frame(source = c(stale, skipped_children), reason =
    c(rep("old_DTU_replaced_by_corrected_saved_results", length(stale)),
      rep("standalone_panels_regenerated_from_parent", length(skipped_children)))),
    file.path(output_dir, "excluded_saved_plots.tsv"))
  audit <- list(); context <- list(source = "", kind = "saved_plot")
  raw_save <- rna_save_plot
  capture_save <- function(plot, stem, dirs, cfg, width = 8, height = 6) {
    folder <- if (is.list(dirs)) dirs$figures else dirs
    rel <- substring(file.path(folder, gsub("[^A-Za-z0-9._-]", "_", stem)), nchar(output_dir) + 2L)
    warnings <- character()
    result <- tryCatch(withCallingHandlers({
      if (exists("rna_plot_meaningful", mode = "function") && !rna_plot_meaningful(plot))
        return(invisible(character()))
      paths <- raw_save(plot, stem, dirs, cfg, width, height)
      if (is.data.frame(plot$data) && ncol(plot$data)) {
        data <- plot$data
        scalar <- vapply(data, function(x) is.atomic(x) && is.null(dim(x)), logical(1))
        if (any(scalar)) rna_write_tsv(data[, scalar, drop = FALSE],
          file.path(output_dir, "plot_data", paste0(rel, ".tsv.gz")))
      }
      list(status = "complete", error = "", paths = paths)
    }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) list(status = "failed", error = conditionMessage(e), paths = character()))
    audit[[length(audit) + 1L]] <<- data.frame(relative_stem = rel, source = context$source,
      source_kind = context$kind, status = result$status, error = result$error,
      warnings = paste(unique(warnings), collapse = " | "), outputs = paste(result$paths, collapse = ";"),
      stringsAsFactors = FALSE)
    rna_write_tsv(do.call(rbind, audit), file.path(output_dir, "export_audit.tsv"))
    if (result$status == "failed") rna_log("Export failed: ", rel, ": ", result$error, level = "ERROR")
    invisible(result$paths)
  }
  rna_progress(0, 4, "saved_figures", status = "starting", detail = paste(length(saved_paths), "saved non-DTU parent figures"))
  for (i in seq_along(saved_paths)) {
    path <- saved_paths[i]; context <- list(source = path, kind = "saved_plot")
    relative <- substring(path, nchar(source_dir) + 2L)
    destination <- file.path(output_dir, dirname(relative))
    stem <- sub("[.]plot[.]rds$", "", basename(path))
    rna_log("Saved plot ", i, "/", length(saved_paths), ": ", relative)
    saved <- tryCatch(readRDS(path), error = function(e) e)
    if (inherits(saved, "error") || !is.list(saved) || is.null(saved$plot)) {
      audit[[length(audit) + 1L]] <- data.frame(relative_stem = sub("[.]plot[.]rds$", "", relative),
        source = path, source_kind = "saved_plot", status = "failed", error = "Invalid saved plot object.", warnings = "", outputs = "")
    } else capture_save(saved$plot, stem, destination, cfg, saved$width %||% 8, saved$height %||% 6)
    rm(saved)
  }
  rna_progress(1, 4, "DTU_views", status = "starting", detail = "Rendering corrected saved DTU measurements and tests")
  # Copy only pure reporting functions into an environment with a local exporter.
  # Their source data remain the verified result object loaded above.
  e <- new.env(parent = environment(dtu_plot_outputs))
  functions <- c("dtu_filter_flow_plots", "dtu_plot_outputs", "dtu_gene_examples",
    "dtu_unpaired_plot_outputs", "dtu_unpaired_gene_examples", "dtu_broad_structure_outputs",
    "dtu_upset_components", "dtu_collective_reports", "dtu_sensitivity_outputs")
  for (name in functions) { fn <- get(name, mode = "function"); environment(fn) <- e; assign(name, fn, envir = e) }
  e$rna_save_plot <- capture_save
  results <- dtu$result$cohort_results
  alphas <- sort(unique(as.numeric(cfg$multiplicity$fdr_thresholds)))
  context <- list(source = file.path(dtu_source_dir, "objects", "DTU_complete.rds"), kind = "corrected_DTU")
  dirs <- rna_dirs(file.path(output_dir, "M06_dtu"))
  report_error <- ""
  tryCatch({
    for (cohort in names(results)) {
      result <- results[[cohort]]; cd <- rna_dirs(file.path(dirs$root, cohort))
      for (alpha in alphas) {
        if (paired_design) e$dtu_plot_outputs(result, cfg, cd, alpha)
        else e$dtu_unpaired_plot_outputs(result, cfg, cd, alpha)
      }
      if (paired_design) e$dtu_gene_examples(result, cfg, cd)
      else e$dtu_unpaired_gene_examples(result, cfg, cd)
      if (paired_design && is.data.frame(result$broad_structure) && nrow(result$broad_structure))
        e$dtu_broad_structure_outputs(result$broad_structure, result$gene, cfg, cd, alphas)
      flow <- as.data.frame(result$prepared$flow)
      if (nrow(flow)) {flow$cohort <- cohort; e$dtu_filter_flow_plots(flow, cfg, dirs)}
    }
    e$dtu_sensitivity_outputs(results, cfg, dirs, alphas)
    if (paired_design) e$dtu_collective_reports(results, cfg, dirs, alphas, dtu$result$gates, remodeling)
  }, error = function(err) {report_error <<- conditionMessage(err); rna_log(report_error, level = "ERROR")})
  rna_progress(2, 4, "saved_pathways", status = "starting", detail = "Reading already-computed ORA tables; no enrichment tests")
  helper <- file.path(cfg$.code_root, "inst", "plotting", "saved_ora_plots.R")
  if (!file.exists(helper)) helper <- system.file("plotting", "saved_ora_plots.R", package = "IsoRemodel")
  ora <- new.env(parent = environment())
  sys.source(helper, envir = ora)
  for (path in ora_paths) {
    relative <- substring(path, nchar(dtu_source_dir) + 2L)
    cohort <- strsplit(relative, "/", fixed = TRUE)[[1]][1]
    if (!cohort %in% names(results)) next
    tag <- regmatches(basename(path), regexpr("FDR[0-9]+", basename(path)))
    alpha <- switch(tag, FDR05 = .05, FDR10 = .10, NA_real_)
    if (!is.finite(alpha) || !alpha %in% alphas) next
    matched <- startsWith(basename(path), "matched_")
    tab <- ora$rna_saved_dtu_read_tsv(path)
    selected <- ora$rna_saved_dtu_select_ora(tab, matched, cfg)
    if (!nrow(selected)) next
    patient <- if (grepl("/patient_[^/]+/", relative))
      sub("^patient_", "", sub("^.*?/(patient_[^/]+)/.*$", "\\1", relative)) else NULL
    p <- ora$rna_saved_dtu_ora_plot(selected, alpha, if (matched) "matched" else "standard", patient, cfg)
    context <- list(source = path, kind = "corrected_DTU_saved_ORA")
    # dirname(relative) ends in tables. Keep patient-specific folders intact.
    destination <- file.path(dirs$root, dirname(dirname(relative)), "figures")
    capture_save(p, sub("[.]tsv[.]gz$", "", basename(path)), destination, cfg, 9.2, 6)
  }
  # Saved balance diagnostics are already computed. Plot only their stored SMDs.
  balance <- dtu$result$enrichment$balance
  if (is.data.frame(balance) && nrow(balance) && all(c("cohort", "threshold", "covariate", "standardized_mean_difference") %in% names(balance))) {
    for (cohort in unique(balance$cohort)) for (alpha in alphas) {
      x <- balance[balance$cohort == cohort & balance$threshold == dtu_fdr_tag(alpha), , drop = FALSE]
      if (!nrow(x)) next
      x$display_covariate <- gsub("_", " ", sub("^log_", "", x$covariate))
      p <- ggplot2::ggplot(x, ggplot2::aes(standardized_mean_difference, display_covariate)) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
        ggplot2::geom_point(size = 3, colour = "#4E79A7") +
        ggplot2::labs(title = "Matched-background balance", subtitle = paste0("FDR < ", alpha),
          x = "Standardized mean difference", y = NULL)
      context <- list(source = file.path(dtu_source_dir, "objects", "DTU_complete.rds"), kind = "corrected_DTU")
      capture_save(p, paste0("matched_background_balance_", cohort, "_", dtu_fdr_tag(alpha)),
                   file.path(dirs$root, "enrichment", "figures"), cfg, 7, 4.8)
    }
  }
  rna_progress(3, 4, "review", status = "starting", detail = "Building the complete figure index and review PDFs")
  png_paths <- sort(list.files(output_dir, pattern = "[.]png$", recursive = TRUE, full.names = TRUE))
  png_paths <- png_paths[!grepl("/(previews|review_volumes)/", png_paths)]
  index <- vector("list", length(png_paths))
  for (i in seq_along(png_paths)) {
    rel <- substring(png_paths[i], nchar(output_dir) + 2L)
    prefix <- sub("[.]png$", "", rel)
    preview <- file.path(output_dir, "previews", rel)
    # Read only one raster at a time; high-resolution exports stay out of ZIPs.
    dimensions <- rna_replot_make_preview(png_paths[i], preview, cfg$reporting$dpi)
    width <- unname(dimensions[["width"]]); height <- unname(dimensions[["height"]])
    index[[i]] <- data.frame(figure_id = sprintf("IR02-F%05d", i), relative_stem = prefix,
      width_in = width, height_in = height, dpi = cfg$reporting$dpi, status = "complete",
      preview_png = preview, pdf = if (file.exists(file.path(output_dir, paste0(prefix, ".pdf"))))
        file.path(output_dir, paste0(prefix, ".pdf")) else NA_character_, stringsAsFactors = FALSE)
  }
  index <- if (length(index)) do.call(rbind, index) else data.frame(figure_id = character(),
    relative_stem = character(), width_in = numeric(), height_in = numeric(), dpi = numeric(),
    status = character(), preview_png = character(), pdf = character())
  rna_write_tsv(index, file.path(output_dir, "figure_index.tsv"))
  review_paths <- character(); review_error <- ""
  tryCatch({review_paths <- rna_replot_review(index, output_dir)}, error = function(err) review_error <<- conditionMessage(err))
  export_audit <- if (length(audit)) do.call(rbind, audit) else data.frame(status = character())
  rna_write_tsv(export_audit, file.path(output_dir, "export_audit.tsv"))
  after <- rna_run_file_manifest(sources)
  unchanged <- identical(before, after)
  rna_write_tsv(after, file.path(output_dir, "source_files_after.tsv"))
  complete <- nrow(index) > 0 && !any(export_audit$status == "failed") && unchanged &&
    !nzchar(report_error) && !nzchar(review_error)
  summary <- data.frame(package_version = rna_package_version(), run_id = run_id,
    analysis_design = dtu$analysis_design,
    status = if (complete) "complete" else "requires_review", figures = nrow(index),
    failed_exports = sum(export_audit$status == "failed"), source_files_unchanged = unchanged,
    analytical_modules_run = FALSE, DTU_source = dtu_source_dir, report_error = report_error,
    review_error = review_error, stringsAsFactors = FALSE)
  rna_write_tsv(summary, file.path(output_dir, "run_summary.tsv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  rna_progress(4, 4, "replot", status = if (complete) "complete" else "requires_review",
    detail = paste(nrow(index), "figures; source files unchanged:", unchanged))
  share <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
  share <- share[grepl("(^review_volumes/|[.]tsv([.]gz)?$|[.]yaml$|[.]log$|[.]txt$)", share)]
  bytes <- file.info(file.path(output_dir, share))$size
  included <- startsWith(share, "review_volumes/") | bytes < 15 * 1024^2
  rna_write_tsv(data.frame(relative_path = share, bytes = bytes, included = included,
    reason = ifelse(included, "included", "large_plot_data_retained_locally")),
    file.path(output_dir, "bundle_files.tsv"))
  share <- c(share[included], "bundle_files.tsv")
  review_zip <- file.path(output_dir, "IsoRemodel_v0.2.1_replot_review.zip")
  zip::zipr(review_zip, files = share, root = output_dir, mode = "mirror", include_directories = FALSE)
  if (!setequal(zip::zip_list(review_zip)$filename, share)) rna_stop("Review ZIP paths were not preserved.")
  rna_log("Replot finished. Output: ", output_dir, "; review: ", review_zip)
  if (isTRUE(strict) && !complete)
    rna_stop("Some plotting checks need review. Successful exports were retained. Upload ", review_zip)
  invisible(list(output_directory = output_dir, figure_index = index,
    review_pdfs = review_paths, review_zip = review_zip, summary = summary))
}
