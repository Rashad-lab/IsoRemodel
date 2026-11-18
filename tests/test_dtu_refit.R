# Targeted, dependency-free regression checks for the DTU refit comparison.
# Loaded by tests/run_all.R after the package's R source files.
.refit_expect_error <- function(expr) {
  value <- tryCatch({ force(expr); NULL }, error = function(e) e)
  if (!inherits(value, "error")) stop("Expected an error, but the expression succeeded.")
  invisible(value)
}
.refit_expect_all_pass <- function(x) {
  stopifnot(is.data.frame(x), nrow(x) > 0L, "pass" %in% names(x),
            all(!is.na(x$pass)), all(x$pass))
}

# Table comparisons join by identity, independently of row order. Missing tests
# stay unavailable and must never silently become conservative P = 1 tests.
.refit_reference <- data.frame(gene_id = paste0("g", seq_len(6)),
  p = c(NA_real_, .01, 1, 0, Inf, -Inf),
  status = c("unavailable", rep("finite", 3), "upper", "lower"),
  stringsAsFactors = FALSE)
.refit_new <- .refit_reference[c(6, 2, 5, 1, 4, 3), , drop = FALSE]
.refit_equal <- dtu_refit_compare_table(.refit_new, .refit_reference,
  keys = "gene_id", numeric_fields = "p", categorical_fields = "status")
stopifnot(isTRUE(.refit_equal$match), nrow(.refit_equal$detail) == 6L)
.refit_expect_all_pass(.refit_equal$summary)

.refit_changed <- .refit_new
.refit_changed$p[.refit_changed$gene_id == "g1"] <- 1
.refit_not_equal <- dtu_refit_compare_table(.refit_changed, .refit_reference,
  keys = "gene_id", numeric_fields = "p", categorical_fields = "status")
stopifnot(identical(.refit_not_equal$match, FALSE),
  !.refit_not_equal$detail$p_match[.refit_not_equal$detail$gene_id == "g1"])
.refit_changed <- .refit_new
.refit_changed$p[.refit_changed$gene_id == "g5"] <- -Inf
stopifnot(!dtu_refit_compare_table(.refit_changed, .refit_reference,
  keys = "gene_id", numeric_fields = "p")$match)

# Numeric tolerance is small and explicit; categorical calls remain exact.
.refit_small <- data.frame(gene_id = "tiny", p = 0, significant = TRUE)
.refit_near <- .refit_small; .refit_near$p <- 5e-13
stopifnot(dtu_refit_compare_table(.refit_near, .refit_small,
  keys = "gene_id", numeric_fields = "p", categorical_fields = "significant")$match)
.refit_near$p <- 2e-12
stopifnot(!dtu_refit_compare_table(.refit_near, .refit_small,
  keys = "gene_id", numeric_fields = "p")$match)
.refit_near$p <- 0; .refit_near$significant <- FALSE
stopifnot(!dtu_refit_compare_table(.refit_near, .refit_small,
  keys = "gene_id", numeric_fields = "p", categorical_fields = "significant")$match)

# Full outer key coverage is required. A matching intersection is insufficient.
.refit_missing <- .refit_new[-1L, , drop = FALSE]
stopifnot(!dtu_refit_compare_table(.refit_missing, .refit_reference,
  keys = "gene_id", numeric_fields = "p")$match)
.refit_extra <- rbind(.refit_new, data.frame(gene_id = "extra", p = .5, status = "finite"))
stopifnot(!dtu_refit_compare_table(.refit_extra, .refit_reference,
  keys = "gene_id", numeric_fields = "p")$match)
.refit_duplicate <- rbind(.refit_new, .refit_new[1L, , drop = FALSE])
.refit_expect_error(dtu_refit_compare_table(.refit_duplicate, .refit_reference,
  keys = "gene_id", numeric_fields = "p"))
.refit_bad_key <- .refit_new; .refit_bad_key$gene_id[1L] <- NA_character_
.refit_expect_error(dtu_refit_compare_table(.refit_bad_key, .refit_reference,
  keys = "gene_id", numeric_fields = "p"))
.refit_expect_error(dtu_refit_compare_table(.refit_new, .refit_reference[, "gene_id", drop = FALSE],
  keys = "gene_id", numeric_fields = "p"))
.refit_expect_error(dtu_refit_compare_table(.refit_new[FALSE, , drop = FALSE],
  .refit_reference[FALSE, , drop = FALSE], keys = "gene_id", numeric_fields = "p"))

# Count-storage type does not alter input equality, but row/sample/design order
# does. These checks run before fitting and do not reorder cached data silently.
.refit_samples <- data.frame(sample_id = c("p1", "r1", "p2", "r2"),
  patient = factor(rep(c("a", "b"), each = 2)),
  timepoint = factor(rep(c("primary", "recurrent"), 2),
    levels = c("primary", "recurrent")))
.refit_counts <- data.frame(gene_id = c("g1", "g1", "g2", "g2"),
  feature_id = paste0("tx", 1:4), p1 = c(4L, 5L, 7L, 8L),
  r1 = c(5L, 4L, 8L, 7L), p2 = c(6L, 3L, 9L, 6L),
  r2 = c(3L, 6L, 6L, 9L), check.names = FALSE)
.refit_design <- stats::model.matrix(~ patient + timepoint, .refit_samples)
rownames(.refit_design) <- .refit_samples$sample_id
.refit_input <- list(counts = .refit_counts, samples = .refit_samples,
                     gene_order = c("g1", "g2"), design = .refit_design)
.refit_expect_all_pass(dtu_refit_validate_inputs(.refit_input, .refit_input, "fixture"))
.refit_double <- .refit_input
for (.refit_column in .refit_samples$sample_id)
  .refit_double$counts[[.refit_column]] <- as.double(.refit_double$counts[[.refit_column]])
.refit_expect_all_pass(dtu_refit_validate_inputs(.refit_input, .refit_double, "storage"))
.refit_bad <- .refit_input; .refit_bad$counts$p1[1] <- .refit_bad$counts$p1[1] + 1
stopifnot(any(!dtu_refit_validate_inputs(.refit_input, .refit_bad, "count")$pass))
.refit_bad <- .refit_input; .refit_bad$counts <- .refit_bad$counts[4:1, , drop = FALSE]
stopifnot(any(!dtu_refit_validate_inputs(.refit_input, .refit_bad, "features")$pass))
.refit_bad <- .refit_input; .refit_bad$samples <- .refit_bad$samples[4:1, , drop = FALSE]
stopifnot(any(!dtu_refit_validate_inputs(.refit_input, .refit_bad, "samples")$pass))
.refit_bad <- .refit_input; .refit_bad$gene_order <- rev(.refit_bad$gene_order)
stopifnot(any(!dtu_refit_validate_inputs(.refit_input, .refit_bad, "genes")$pass))
.refit_bad <- .refit_input; .refit_bad$design[1L, 1L] <- 2
stopifnot(any(!dtu_refit_validate_inputs(.refit_input, .refit_bad, "design")$pass))

# A completed production cache must identify its owning stage. Pilot, partial,
# wrong-stage, and old unwrapped cache payloads are rejected.
.refit_cache_dir <- tempfile("dtu-refit-cache-test-")
dir.create(.refit_cache_dir)
.refit_cache_file <- file.path(.refit_cache_dir, "stage_result.rds")
.refit_cache <- list(schema_version = "whole-stage-cache-v1", stage = "import",
  complete = TRUE, profile = "production", result = list(value = 42))
saveRDS(.refit_cache, .refit_cache_file)
.refit_loaded <- dtu_refit_load_stage(.refit_cache_file, "import")
stopifnot(is.list(.refit_loaded))
for (.refit_bad_case in c("incomplete", "pilot", "wrong_stage", "old_schema", "no_wrapper")) {
  .refit_bad_cache <- .refit_cache
  if (.refit_bad_case == "incomplete") .refit_bad_cache$complete <- FALSE
  if (.refit_bad_case == "pilot") .refit_bad_cache$profile <- "pilot"
  if (.refit_bad_case == "wrong_stage") .refit_bad_cache$stage <- "remodeling"
  if (.refit_bad_case == "old_schema") .refit_bad_cache$schema_version <- "unknown"
  if (.refit_bad_case == "no_wrapper") .refit_bad_cache <- .refit_bad_cache$result
  saveRDS(.refit_bad_cache, .refit_cache_file)
  .refit_expect_error(dtu_refit_load_stage(.refit_cache_file, "import"))
}
unlink(.refit_cache_dir, recursive = TRUE)
message("DTU targeted-refit comparison and input tests passed.")
