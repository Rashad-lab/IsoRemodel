# Branded console/file logging and accurate progress independent of model results.
local({
  stopifnot(identical(unname(read.dcf(file.path(.repo_root, "DESCRIPTION"))[1L,"Package"]), "IsoRemodel"),
            identical(rna_package_version(.repo_root), "0.2.1"))
  folder <- tempfile("isoremodel_logging_"); dir.create(folder)
  on.exit(unlink(folder, recursive=TRUE), add=TRUE)
  names <- c("IsoRemodel.log_file","IsoRemodel.progress_file","IsoRemodel.progress_latest_file",
             "IsoRemodel.run_id","IsoRemodel.run_start")
  previous <- stats::setNames(lapply(names, getOption), names)
  on.exit(options(previous), add=TRUE)
  options(IsoRemodel.log_file=file.path(folder,"run.log"),
    IsoRemodel.progress_file=file.path(folder,"progress.tsv"),
    IsoRemodel.progress_latest_file=file.path(folder,"progress_latest.tsv"),
    IsoRemodel.run_id="fixture", IsoRemodel.run_start=Sys.time())
  rna_banner(log=TRUE)
  logged <- rna_log("message with\nline",level="WARNING")
  stopifnot(grepl("[IsoRemodel]",logged,fixed=TRUE),grepl("[WARNING]",logged,fixed=TRUE),
    grepl("message with | line",logged,fixed=TRUE))
  row <- rna_progress(2,5,"remodeling",status="starting",detail="row 1\nrow 2",stage_index=4L)
  stopifnot(row$stage_index==4L,row$stages_passed==2L,row$stages_percent==40,row$elapsed_seconds>=0)
  rna_progress(2,5,"remodeling",status="failed",stage_index=4L)
  history <- utils::read.delim(file.path(folder,"progress.tsv"),stringsAsFactors=FALSE)
  latest <- utils::read.delim(file.path(folder,"progress_latest.tsv"),stringsAsFactors=FALSE)
  lines <- readLines(file.path(folder,"run.log"),warn=FALSE)
  stopifnot(nrow(history)==2L,nrow(latest)==1L,latest$status=="failed",
    all(history$stage_index==4L),all(history$stages_percent==40),
    all(startsWith(lines,"[IsoRemodel]")),any(grepl("IsoRemodel v0.2.1",lines,fixed=TRUE)))
  # The installed package's standard attachment hook resolves the branded banner.
  text <- capture.output(suppressPackageStartupMessages(.onAttach(NULL,"IsoRemodel")),type="message")
  stopifnot(length(text)==0L)
})
message("IsoRemodel banner, logging and progress tests passed.")
