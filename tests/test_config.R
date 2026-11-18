# Configuration and metadata validation never invokes the scientific engines.
test_config <- function() {
  rna_require("yaml")
  root <- if(exists(".repo_root",inherits=TRUE)) get(".repo_root",inherits=TRUE) else getwd()
  path <- file.path(root,"config","analysis.yaml")
  paper <- read_rna_config(path,validate_inputs=FALSE)
  defaults <- read_rna_metadata(paper,validate_inputs=FALSE)
  stopifnot(nrow(defaults)==12L,identical(defaults$sample,defaults$sample_id),
    identical(paper$harmonization$canonical_gene_policy,"paper_support_rank_preserve_all_counts"),
    paper$reporting$font_family=="DejaVu Sans",paper$reporting$title_size_pt==16,
    paper$reporting$font_size_pt==14,all(c("pdf","svg","png","jpeg") %in% paper$reporting$formats),
    basename(defaults$classification_file[1L])=="ptn1_primary_classification.filtered_lite_classification.txt")
  expected <- c("all_6_pairs","excluding_ptn2_5_pairs")
  stopifnot(identical(names(rna_resolve_cohorts(defaults,paper)),expected))
  allcfg <- paper;allcfg$sensitivity$mode <- "all"
  allcohorts <- rna_resolve_cohorts(defaults,allcfg)
  stopifnot(length(allcohorts)==7L,all(lengths(allcohorts[-1L])==5L),
    identical(allcohorts$excluding_ptn2_5_pairs,c("ptn1","ptn3","ptn4","ptn5","ptn6")))
  none <- paper;none$sensitivity$mode <- "none"
  stopifnot(length(rna_resolve_cohorts(defaults,none))==1L)
  error <- function(expr,pattern) {
    msg <- tryCatch({force(expr);NULL},error=function(e) conditionMessage(e))
    if(is.null(msg) || !grepl(pattern,msg,ignore.case=TRUE)) stop("Expected error matching ",pattern,"; got: ",msg)
  }
  bad <- paper;bad$gain_loss$depth_robust_probability <- 1.2
  error(rna_config_validate(bad),"depth_robust_probability")
  bad <- paper;bad$gain_loss$n_depth_curve_replicates <- bad$gain_loss$n_downsample_replicates+1
  error(rna_config_validate(bad),"depth-curve")
  bad <- paper;bad$sensitivity$patients <- "unknown_subject"
  error(rna_resolve_cohorts(defaults,bad),"unknown subject")
  four <- defaults[defaults$patient %in% paste0("ptn",1:4),,drop=FALSE]
  error(rna_resolve_cohorts(four,paper),"absolute requirement")
  folder <- tempfile("config_fixture_");dir.create(folder)
  on.exit(unlink(folder,recursive=TRUE),add=TRUE)
  fixture <- make_rna_fixture(file.path(folder,"inputs"))
  raw <- utils::read.csv(fixture$metadata_file,stringsAsFactors=FALSE,check.names=FALSE)
  # Labels may differ from the internal role names, without changing filenames.
  raw$condition <- ifelse(raw$condition=="primary","initial","relapse")
  raw$classification_file <- basename(raw$classification_file)
  raw$junctions_file <- basename(raw$junctions_file)
  meta <- file.path(folder,"samples.csv");utils::write.csv(raw,meta,row.names=FALSE,na="")
  y <- yaml::read_yaml(path,eval.expr=FALSE)
  y$inputs <- list(root_directory="inputs",sample_metadata="samples.csv")
  y$comparison$reference <- "initial";y$comparison$comparison <- "relapse"
  y$output$directory <- "new_output_does_not_exist"
  y$references$annotation_gtf <- fixture$gtf_file
  y$enrichment$reference_cache_directory <- "new_pathway_cache"
  config_path <- file.path(folder,"analysis.yaml");yaml::write_yaml(y,config_path)
  loaded <- read_rna_config(config_path,validate_inputs=TRUE)
  samples <- read_rna_metadata(loaded,validate_inputs=TRUE)
  stopifnot(all(file.exists(samples$classification_file)),
    all(is.na(samples$pigeon_summary_file)),all(is.na(samples$filter_reasons_file)),
    identical(samples$timepoint,ifelse(samples$condition=="initial","primary","recurrent")),
    !dir.exists(loaded$output$directory),identical(loaded$.config_dir,normalizePath(folder,winslash="/")))
  duplicate <- raw;duplicate$sample_id[2L] <- duplicate$sample_id[1L]
  utils::write.csv(duplicate,meta,row.names=FALSE,na="")
  error(read_rna_metadata(loaded,FALSE),"unique")
  utils::write.csv(raw[-1L,],meta,row.names=FALSE,na="")
  error(read_rna_metadata(loaded,FALSE),"exactly one sample")
  unavailable <- raw;unavailable$classification_file[1L] <- "missing_classification.txt"
  utils::write.csv(unavailable,meta,row.names=FALSE,na="")
  stopifnot(nrow(read_rna_metadata(loaded,FALSE))==12L)
  error(read_rna_metadata(loaded,TRUE),"does not exist")
  message("Configuration and metadata tests passed.")
}
test_config()
