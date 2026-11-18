# Metadata and exclusion design contracts; no statistical engine is called.
(function(){
 root <- .repo_root
 cfg <- read_rna_config(file.path(root,'config','analysis.yaml'),FALSE)
 paired <- read_rna_metadata(cfg,FALSE)
 stopifnot(identical(names(rna_resolve_cohorts(paired,cfg)),c('all_6_pairs','excluding_ptn2_5_pairs')))
 none <- cfg;none$sensitivity$mode <- 'none'
 stopifnot(length(rna_resolve_cohorts(paired,none))==1L)
 all <- cfg;all$sensitivity$mode <- 'all'
 stopifnot(length(rna_resolve_cohorts(paired,all))==7L)
 wrong_unit <- cfg;wrong_unit$cohorts <- c(cfg$cohorts,list(list(id='wrong_unit',exclude_samples='ptn1_primary')))
 stopifnot(inherits(try(rna_resolve_cohorts(paired,wrong_unit),silent=TRUE),'try-error'))
 temp <- tempfile('unpaired_metadata_');dir.create(temp);on.exit(unlink(temp,recursive=TRUE),add=TRUE)
 d <- data.frame(sample_id=paste0('S',1:8),condition=rep(c('control','treated'),each=4),
    classification_file=paste0('class',1:8,'.txt'),junctions_file=paste0('junction',1:8,'.txt'))
 utils::write.csv(d,file.path(temp,'samples.csv'),row.names=FALSE)
 cfg$comparison$paired <- FALSE;cfg$comparison$reference <- 'control';cfg$comparison$comparison <- 'treated'
 cfg$.metadata_path <- file.path(temp,'samples.csv');cfg$inputs$root_directory <- temp
 cfg$cohorts <- list(list(id='all_samples',exclude_samples=character()))
 cfg$sensitivity$mode <- 'none'; cfg$sensitivity$samples <- NULL
 x <- read_rna_metadata(cfg,FALSE)
 stopifnot(identical(x$subject_id,x$sample_id),length(rna_resolve_cohorts(x,cfg))==1L)
 cfg$sensitivity$mode <- 'selected';cfg$sensitivity$samples <- c('S1','S8')
 co <- rna_resolve_cohorts(x,cfg);stopifnot(length(co)==3L,all(lengths(co[-1])==7L),! 'S1'%in%co[[2]],! 'S8'%in%co[[3]])
 cfg$sensitivity$mode <- 'all';co <- rna_resolve_cohorts(x,cfg)
 stopifnot(length(co)==9L,all(lengths(co[-1])==7L))
 cfg$cohorts <- c(cfg$cohorts,list(list(id='exclude_one_per_group',exclude_samples=c('S1','S5'))))
 stopifnot(length(rna_resolve_cohorts(x,cfg)$exclude_one_per_group)==6L)
 fails <- function(expr){stopifnot(inherits(try(force(expr),silent=TRUE),'try-error'))}
 cfg$cohorts <- list(list(id='all_samples',exclude_samples=character()))
 fails(rna_resolve_cohorts(x[-1,],cfg)) # LOO would leave2control replicates
 cfg$sensitivity$mode <- 'selected';cfg$sensitivity$samples <- 'missing';fails(rna_resolve_cohorts(x,cfg))
 d$subject_id <- rep(c('same','other'),each=4);utils::write.csv(d,file.path(temp,'samples.csv'),row.names=FALSE)
 fails(read_rna_metadata(cfg,FALSE))
 cfg$comparison$paired <- 'false';fails(rna_config_validate(cfg))
})()
message('Paired/unpaired metadata and none/selected/all exclusion contracts passed.')
