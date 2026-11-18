# Run bounded correctness tests without performing the full study analysis.
run_rna_tests <- function(repo_root = ".", pattern = "^test_.*[.]R$") {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(repo_root, "DESCRIPTION")))
    stop("repo_root must be the IsoRemodel repository directory.")
  env <- new.env(parent = globalenv())
  env$.repo_root <- repo_root
  # Explicit priority ensures helpers exist before source-derived engines.
  rfiles <- list.files(file.path(repo_root, "R"), pattern = "[.]R$", full.names = TRUE)
  priority <- c("utils.R", "config.R", "import.R", "gain_loss.R", "qc.R",
                "remodeling.R", "dtu.R", "enrichment.R", "run.R")
  rank <- match(basename(rfiles), priority)
  rank[is.na(rank)] <- length(priority) + seq_len(sum(is.na(rank)))
  for (path in rfiles[order(rank, basename(rfiles))]) sys.source(path, envir = env)
  fixture <- file.path(repo_root, "tests", "make_fixture.R")
  if (file.exists(fixture)) sys.source(fixture, envir = env)
  tests <- sort(list.files(file.path(repo_root, "tests"), pattern = pattern, full.names = TRUE))
  if (!length(tests)) stop("No test files matched: ", pattern)
  failures <- character()
  results <- vector("list", length(tests))
  for (i in seq_along(tests)) {
    path <- tests[[i]]
    message("Testing ", basename(path))
    t0 <- proc.time()[["elapsed"]]
    outcome <- tryCatch({ sys.source(path, envir = env); NULL }, error = function(e) conditionMessage(e))
    results[[i]] <- data.frame(test = basename(path), passed = is.null(outcome),
                               seconds = proc.time()[["elapsed"]] - t0,
                               error = if (is.null(outcome)) "" else outcome,
                               stringsAsFactors = FALSE)
    if (!is.null(outcome)) {
      failures <- c(failures, paste0(basename(path), ": ", outcome))
      message("FAILED: ", outcome)
    }
  }
  result <- do.call(rbind, results)
  print(result, row.names = FALSE)
  if (length(failures)) stop("RNA tests failed:\n", paste(failures, collapse = "\n"), call. = FALSE)
  message("All ", nrow(result), " bounded test files passed.")
  invisible(result)
}

if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  root <- if (length(arg)) dirname(dirname(normalizePath(sub("^--file=", "", arg[[1L]]), winslash = "/"))) else "."
  run_rna_tests(root)
}
