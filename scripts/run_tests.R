args <- commandArgs(trailingOnly = TRUE)
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]), winslash = "/"))) else "."
root_arg <- grep("^--project=", args, value = TRUE)
if (length(root_arg)) root <- sub("^--project=", "", root_arg[[1L]])
source(file.path(root, "tests", "run_all.R"), local = TRUE)
run_rna_tests(root)
