# Dependency setup for IsoRemodel.
# source this file to define the functions, or execute with Rscript.

rna_dependency_specification <- function(include_optional = TRUE) {
  cran <- c("data.table", "digest", "dplyr", "forcats", "ggplot2", "jsonlite",
            "knitr", "MASS", "Matrix", "patchwork", "purrr", "readr", "ragg", "scales",
            "showtext", "stringr", "svglite", "sysfonts", "systemfonts",
            "tibble", "tidyr", "vegan", "yaml", "rlang", "png", "zip")
  bioc <- c("BiocParallel", "DRIMSeq", "stageR")
  if (isTRUE(include_optional)) {
    cran <- c(cran, "ComplexUpset", "ggalluvial", "ggrepel", "msigdbr",
              "openxlsx", "renv", "rmarkdown")
    bioc <- c(bioc, "AnnotationDbi", "clusterProfiler", "GO.db", "org.Hs.eg.db", "rtracklayer")
  }
  # Author-reported package versions from the completed RNA analyses.
  # New configuration/font dependencies were not part of that version report.
  known <- c(data.table = "1.17.8", digest = "0.6.39", dplyr = "1.1.4",
             forcats = "1.0.1", ggplot2 = "4.0.3", jsonlite = "2.0.0",
             knitr = "1.51", MASS = "7.3.65", Matrix = "1.7.4", patchwork = "1.3.2",
             purrr = "1.2.0", readr = "2.1.6", scales = "1.4.0",
             stringr = "1.6.0", svglite = "2.2.2", tibble = "3.3.0",
             tidyr = "1.3.2", vegan = "2.7.6", ComplexUpset = "1.3.3",
             ggalluvial = "0.12.5", ggrepel = "0.9.6", msigdbr = "25.1.1",
             openxlsx = "4.2.8.1", rmarkdown = "2.30", BiocParallel = "1.44.0",
             DRIMSeq = "1.38.0", stageR = "1.32.0", AnnotationDbi = "1.72.0",
             clusterProfiler = "4.18.4", org.Hs.eg.db = "3.22.0", rtracklayer = "1.70.1")
  data.frame(package = c(cran, bioc), source = c(rep("CRAN", length(cran)), rep("Bioconductor", length(bioc))),
             paper_version = unname(known[c(cran, bioc)]), stringsAsFactors = FALSE)
}

check_rna_dependencies <- function(include_optional = TRUE) {
  x <- rna_dependency_specification(include_optional)
  x$installed <- vapply(x$package, function(p) length(find.package(p, quiet = TRUE)) > 0L, logical(1))
  x$installed_version <- vapply(x$package, function(p) {
    if (!length(find.package(p, quiet = TRUE))) return(NA_character_)
    as.character(utils::packageVersion(p))
  }, character(1))
  x$matches_reported_paper_version <- ifelse(is.na(x$paper_version), NA,
                                           x$installed_version == x$paper_version)
  x$adapter_compatible <- ifelse(x$package == "DRIMSeq",
                                 x$installed_version == "1.38.0", TRUE)
  x
}

install_rna_dependencies <- function(upgrade = FALSE, pin_known_versions = TRUE,
                                     include_optional = TRUE, lib = .libPaths()[1L],
                                     cran_repo = "https://cloud.r-project.org") {
  message("[IsoRemodel] Checking and installing workflow dependencies.")
  if (getRversion() < "4.5.0" || getRversion() >= "4.6.0") {
    stop("The paper environment uses R 4.5.x and Bioconductor 3.22. ",
         "Use R 4.5.x for this installer. A different R/Bioconductor profile ",
         "must be validated before claiming paper reproduction.", call. = FALSE)
  }
  lib <- path.expand(lib)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (file.access(lib, 2L) != 0L) stop("Package library is not writable: ", lib, call. = FALSE)
  .libPaths(unique(c(normalizePath(lib, winslash = "/", mustWork = TRUE), .libPaths())))
  old_repos <- getOption("repos")
  on.exit(options(repos = old_repos), add = TRUE)
  options(repos = c(CRAN = cran_repo))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager", repos = cran_repo, lib = lib)
  if (!identical(as.character(BiocManager::version()), "3.22")) {
    BiocManager::install(version = "3.22", ask = FALSE, update = FALSE, lib = lib)
  }
  spec <- rna_dependency_specification(include_optional)
  if (isTRUE(pin_known_versions) && !requireNamespace("remotes", quietly = TRUE))
    utils::install.packages("remotes", repos = cran_repo, lib = lib)

  installed_version <- function(p) {
    if (!length(find.package(p, quiet = TRUE))) return(NA_character_)
    as.character(utils::packageVersion(p))
  }
  for (i in seq_len(nrow(spec))) {
    p <- spec$package[[i]]
    wanted <- spec$paper_version[[i]]
    current <- installed_version(p)
    exact <- isTRUE(pin_known_versions) && !is.na(wanted)
    # The audited adapter is always pinned, including an explicit upgrade run.
    if (p == "DRIMSeq") { wanted <- "1.38.0"; exact <- TRUE }
    if (!is.na(current) && (!exact || identical(current, wanted)) && !isTRUE(upgrade)) next
    if (!is.na(current) && exact && identical(current, wanted)) next
    message("[IsoRemodel] Installing ", p, if (exact) paste0(" ", wanted) else " (compatible repository version)")
    if (spec$source[[i]] == "CRAN") {
      if (exact) {
        # packageVersion() prints hyphenated upstream versions in normalized
        # form (e.g. Matrix 1.7-4 -> 1.7.4). remotes compares package_version
        # objects and resolves the actual repository/archive tarball name.
        remotes::install_version(p, version = wanted, repos = cran_repo,
                                 lib = lib, upgrade = "never", dependencies = NA,
                                 quiet = FALSE)
      } else {
        utils::install.packages(p, repos = cran_repo, lib = lib, dependencies = NA)
      }
    } else {
      BiocManager::install(p, version = "3.22", lib = lib,
                           ask = FALSE, update = FALSE, force = isTRUE(upgrade))
      if (exact && !identical(installed_version(p), wanted)) {
        # Historical Bioconductor versions remain in the release's source
        # repository/Archive. Try both official locations explicitly.
        repos <- BiocManager::repositories(version = "3.22")
        repos <- repos[grepl("^BioC", names(repos))]
        bases <- paste0(sub("/$", "", unname(repos)), "/src/contrib")
        urls <- unique(c(sprintf("%s/%s_%s.tar.gz", bases, p, wanted),
                         sprintf("%s/Archive/%s/%s_%s.tar.gz", bases, p, p, wanted)))
        tmp <- tempfile(fileext = ".tar.gz")
        on.exit(unlink(tmp), add = TRUE)
        got <- FALSE
        for (u in urls) {
          ok <- tryCatch({ utils::download.file(u, tmp, mode = "wb", quiet = FALSE); TRUE },
                         error = function(e) FALSE, warning = function(w) FALSE)
          if (ok && file.exists(tmp) && file.info(tmp)$size > 0) {
            utils::install.packages(tmp, repos = NULL, type = "source", lib = lib)
            got <- identical(installed_version(p), wanted)
            if (got) break
          }
        }
        if (!got) stop("Could not install the requested paper version ", p, " ", wanted,
                       ". Preserve the installation log and inspect official release archives.", call. = FALSE)
      }
    }
  }
  result <- check_rna_dependencies(include_optional)
  print(result, row.names = FALSE)
  if (any(!result$installed)) stop("Some dependencies remain unavailable.", call. = FALSE)
  if (!identical(installed_version("DRIMSeq"), "1.38.0"))
    stop("DRIMSeq must be 1.38.0 for the bundled numerical adapter.", call. = FALSE)
  if (isTRUE(pin_known_versions) && any(result$matches_reported_paper_version %in% FALSE))
    stop("Installed direct dependencies do not match all requested paper versions.", call. = FALSE)
  message("[IsoRemodel] Dependency setup completed. Restart R before running the analysis or capture_environment.R.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  lib <- .libPaths()[1L]
  hit <- grep("^--lib=", args, value = TRUE)
  if (length(hit)) lib <- sub("^--lib=", "", hit[[1L]])
  unknown <- args[!args %in% c("--upgrade", "--no-pin-known-versions", "--core-only") & !grepl("^--lib=", args)]
  if (length(unknown)) stop("Unknown installer argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  install_rna_dependencies(upgrade = "--upgrade" %in% args,
                           pin_known_versions = !"--no-pin-known-versions" %in% args,
                           include_optional = !"--core-only" %in% args, lib = lib)
}
