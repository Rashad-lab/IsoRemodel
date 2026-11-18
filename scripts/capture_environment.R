# Capture the actual machine environment. This does not install or upgrade.
capture_rna_environment <- function(output_dir = "environment", project_root = ".",
                                    write_lockfile = TRUE) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  writeLines(capture.output(Sys.info()), file.path(output_dir, "systemInfo.txt"))
  writeLines(capture.output(.libPaths()), file.path(output_dir, "library_paths.txt"))
  installer <- file.path(project_root, "scripts", "install_dependencies.R")
  env <- new.env(parent = baseenv())
  sys.source(installer, envir = env)
  versions <- env$check_rna_dependencies(include_optional = TRUE)
  utils::write.table(versions, file.path(output_dir, "package_versions.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    writeLines(as.character(BiocManager::version()), file.path(output_dir, "Bioconductor_version.txt"))
  }
  if (isTRUE(write_lockfile)) {
    if (!requireNamespace("renv", quietly = TRUE)) {
      warning("renv is not installed. Version tables were written; no lockfile was created.")
    } else {
      pkgs <- intersect(versions$package, rownames(utils::installed.packages()))
      renv::snapshot(project = project_root, packages = pkgs,
                     lockfile = file.path(output_dir, "renv.lock"),
                     prompt = FALSE, force = FALSE)
    }
  }
  message("Captured the actual environment in ", normalizePath(output_dir, winslash = "/"))
  invisible(versions)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root_arg <- grep("^--project=", args, value = TRUE)
  out_arg <- grep("^--output=", args, value = TRUE)
  root <- if (length(root_arg)) sub("^--project=", "", root_arg[[1L]]) else "."
  out <- if (length(out_arg)) sub("^--output=", "", out_arg[[1L]]) else file.path(root, "environment")
  capture_rna_environment(output_dir = out, project_root = root,
                          write_lockfile = !"--no-lockfile" %in% args)
}
