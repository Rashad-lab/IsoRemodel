# Run after sourcing R/utils.R, R/import.R and tests/make_fixture.R.
if (!exists("make_rna_fixture", mode = "function")) {
  source(file.path(.repo_root, "tests", "make_fixture.R"), local = environment())
}
test_import <- function() {
  rna_require(c("data.table", "digest"))
  root <- tempfile("rna_import_test_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- make_rna_fixture(file.path(root, "inputs"), n_pairs = 3L)
  imported <- import_pigeon(fixture$cfg, fixture$samples, file.path(root, "import"))
  stopifnot(nrow(imported$local_transcripts) == 60L,
              ncol(imported$counts_matrix) == 6L,
              all(colSums(imported$counts_matrix) == 381L),
              all(imported$count_conservation$difference == 0),
              !anyDuplicated(rownames(imported$counts_matrix)),
              identical(rownames(imported$counts_matrix), imported$isoform_annotation$isoform_harmonized_id))
  local <- imported$local_transcripts
  chain7 <- local$isoform_harmonized_id[local$sample == "ptn1_primary" & local$isoform_local == "PB.1.P.7"]
  chain8 <- local$isoform_harmonized_id[local$sample == "ptn1_primary" & local$isoform_local == "PB.1.P.8"]
  stopifnot(identical(chain7, chain8), imported$counts_matrix[chain7, "ptn1_primary"] == 60L)
  expected <- paste0("JUNC_", substr(digest::digest("chr4|+|701:800:GT-AG", algo = "xxhash64"), 1L, 16L))
  stopifnot(identical(chain7, expected))
  a <- imported$isoform_annotation[match(chain7, imported$isoform_annotation$isoform_harmonized_id), ]
  stopifnot(a$associated_gene == "ATRX", a$n_gene_candidates == 2L,
              a$mapping_resolution == "resolved_by_FL_and_sample_support")
  mono <- local[local$exons == 1L, ]
  stopifnot(nrow(mono) == 6L, length(unique(mono$isoform_harmonized_id)) == 6L,
              all(startsWith(mono$isoform_harmonized_id, "PBONLY_")))
  # Missing optional annotations remain missing, not false/zero.
  fixture2 <- make_rna_fixture(file.path(root, "missing_optional"), n_pairs = 2L)
  first <- fixture2$samples$classification_file[1L]
  raw <- utils::read.delim(first, check.names = FALSE, stringsAsFactors = FALSE)
  raw$within_CAGE_peak <- NULL
  raw$diff_to_gene_TSS <- NULL
  utils::write.table(raw, first, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  missing <- import_pigeon(fixture2$cfg, fixture2$samples, file.path(root, "missing_import"))
  first_local <- missing$local_transcripts[missing$local_transcripts$sample == "ptn1_primary", ]
  stopifnot(all(is.na(first_local$within_CAGE_peak)), all(is.na(first_local$diff_to_gene_TSS)))
  # Shuffling sample order must not change the canonical count matrix.
  shuffled <- import_pigeon(fixture$cfg, fixture$samples[rev(seq_len(nrow(fixture$samples))), ],
                             file.path(root, "shuffled"))
  stopifnot(identical(imported$counts_matrix, shuffled$counts_matrix))
  expect_error <- function(expr, pattern) {
    error <- tryCatch({force(expr); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
  }
  # Do not silently deduplicate source transcript records.
  fixture3 <- make_rna_fixture(file.path(root, "duplicated"), n_pairs = 2L)
  dup_file <- fixture3$samples$classification_file[1L]
  raw <- utils::read.delim(dup_file, check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(rbind(raw, raw[1L, ]), dup_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  expect_error(import_pigeon(fixture3$cfg, fixture3$samples, file.path(root, "duplicate_import")), "duplicated local")
  # Integer FL, not TPM or rounded continuous abundance, is required.
  fixture4 <- make_rna_fixture(file.path(root, "fractional"), n_pairs = 2L)
  fractional_file <- fixture4$samples$classification_file[1L]
  raw <- utils::read.delim(fractional_file, check.names = FALSE, stringsAsFactors = FALSE)
  raw[["FL.sample"]][1L] <- 90.5
  utils::write.table(raw, fractional_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  expect_error(import_pigeon(fixture4$cfg, fixture4$samples, file.path(root, "fractional_import")), "integer")
  message("Import tests passed: identity, aggregation, gene annotation, missingness and invalid-count/key checks.")
  invisible(TRUE)
}

test_import()
