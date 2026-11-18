# Source-mode loading for a repository checkout; package installation is optional.
# This file has no side effects when sourced by a test harness or package build.
load_rna_code <- function(repo_root=".",envir=parent.frame(),startup=interactive()) {
  root <- normalizePath(repo_root,winslash="/",mustWork=TRUE)
  if(!file.exists(file.path(root,"DESCRIPTION")) || !dir.exists(file.path(root,"R")))
    stop("repo_root must point to the IsoRemodel repository.",call.=FALSE)
  priority <- c("utils.R","config.R","import.R","gain_loss.R","qc.R","remodeling.R","dtu.R","enrichment.R","run.R")
  files <- list.files(file.path(root,"R"),pattern="[.]R$",full.names=TRUE)
  files <- files[basename(files)!="load.R"]
  rank <- match(basename(files),priority);rank[is.na(rank)] <- length(priority)+seq_len(sum(is.na(rank)))
  for(file in files[order(rank,basename(files))]) sys.source(file,envir=envir)
  options(IsoRemodel.code_root=root)
  if (isTRUE(startup) && !isTRUE(getOption("IsoRemodel.quiet", FALSE))) envir$rna_banner()
  invisible(root)
}
rna_load <- load_rna_code
