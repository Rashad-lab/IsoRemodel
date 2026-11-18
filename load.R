# Load the repository's R functions without changing the working directory.
# source('/absolute/path/to/IsoRemodel/load.R')
.rna_entry_file <- tryCatch(sys.frame(1L)$ofile,error=function(e) NULL)
if(is.null(.rna_entry_file) || !nzchar(.rna_entry_file)) {
  .rna_entry_arg <- grep("^--file=",commandArgs(),value=TRUE)
  if(length(.rna_entry_arg)) .rna_entry_file <- sub("^--file=","",.rna_entry_arg[[1L]])
}
if(is.null(.rna_entry_file) || !nzchar(.rna_entry_file)) stop("Cannot locate load.R; use source() with its path.")
.rna_entry_root <- dirname(normalizePath(.rna_entry_file,winslash="/",mustWork=TRUE))
source(file.path(.rna_entry_root,"R","load.R"),local=TRUE)
load_rna_code(.rna_entry_root,envir=environment())
rm(list=intersect(c(".rna_entry_file",".rna_entry_root",".rna_entry_arg"),ls(all.names=TRUE)))
