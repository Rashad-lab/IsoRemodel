# Load IsoRemodel. Nothing is exported until replot_isoremodel() is called.
.ir_entry <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
if (is.null(.ir_entry)) stop("Use source() with the absolute path to run_replot.R.")
.ir_root <- dirname(normalizePath(.ir_entry, winslash = "/", mustWork = TRUE))
source(file.path(.ir_root, "load.R"), local = TRUE)
rm(.ir_entry, .ir_root)
message("Run replot_isoremodel() to export saved results with IsoRemodel v0.2.1 layouts.")
