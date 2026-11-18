# Pure edge cases plus a tiny two-replicate full-stage/resume check.
test_gain_loss <- function() {
  required <- c("Matrix","vegan","dplyr","tidyr","ggplot2","digest")
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Gain/loss tests require: ",paste(missing,collapse=", "))
  subjects <- paste0("case",LETTERS[1:4])
  metadata <- do.call(rbind,lapply(subjects,function(p) data.frame(
    sample=paste(p,c("baseline","followup"),sep="_"),patient=p,
    timepoint=c("primary","recurrent"),stringsAsFactors=FALSE)))
  ids <- c("JUNC_target","JUNC_other","JUNC_ineligible","JUNC_unavailable","PBONLY_example")
  x <- matrix(0L,5,8,dimnames=list(ids,metadata$sample))
  x[1,seq(2,8,2)] <- 10L
  x[2,seq(1,8,2)] <- 90L;x[2,seq(2,8,2)] <- 80L
  x[3,] <- 10L
  x[4,seq(1,8,2)] <- 19L;x[4,seq(2,8,2)] <- 100L
  x[5,seq(2,8,2)] <- 10L
  a <- data.frame(isoform_harmonized_id=ids,
    associated_gene=c("gene1","gene1","gene1","gene2","gene3"),
    structural_category=c("full-splice_match","full-splice_match","intergenic","full-splice_match","full-splice_match"),
    mapping_resolution="one_gene",cross_sample_comparable=c(TRUE,TRUE,TRUE,TRUE,FALSE),
    gene_mapped=TRUE,analyzed_for_gain_loss=c(TRUE,TRUE,FALSE,TRUE,FALSE),stringsAsFactors=FALSE)
  inp <- list(counts_matrix=x,isoform_annotation=a,sample_manifest=metadata,fingerprint="synthetic_gain_loss_fixture")
  cfg <- list(gain_loss=list(isoform_presence_min_fl=5L,isoform_presence_min_fraction=.1,
    gene_evaluable_min_fl=20L,primary_patient_support=3L,
    keep_structural_categories=c("full-splice_match"),require_cross_sample_comparable_id=TRUE,
    n_downsample_replicates=2L,n_depth_curve_replicates=2L,checkpoint_every=1L,
    patient_support_thresholds=c(2L,3L),gene_floor_sensitivity=c(10L,20L),
    depth_curve_fractions=c(.5,1),depth_robust_probability=.8,
    depth_robust_probability_sensitivity=c(.7,.8,.9),opposite_direction_maximum=0L),
    execution=list(seeds=list(gain_loss=20260903L),resume=TRUE,force_refit=FALSE),
    reporting=list(enabled=FALSE))
  cohorts <- list(all_pairs=subjects,omit_last=subjects[1:3])
  ctx <- rna_gain_loss_context(inp,cfg,cohorts)
  state <- rna_gain_loss_state(x,ctx,20L,TRUE)
  stopifnot(all(state$gain["JUNC_target",]),!any(state$evaluable["JUNC_unavailable",]),
    !any(state$gain["JUNC_unavailable",]),!any(state$gain["PBONLY_example",]))
  su <- rna_gain_loss_support(state,subjects)
  stopifnot(rna_gain_loss_call(su$gain,su$loss,3L,"recurrent_gained")[1])
  # A single opposing event vetoes an otherwise four-patient gain rule.
  stopifnot(!rna_gain_loss_call(4L,1L,3L,"recurrent_gained"))
  # Full gene denominator includes an isoform outside the candidate categories.
  y <- x;y[1,2] <- 5L;y[2,2] <- 35L;y[3,2] <- 60L
  st <- rna_gain_loss_state(y,ctx,20L)
  stopifnot(!st$gain["JUNC_target",1],!st$present_recurrent["JUNC_target",1])
  # Presence requires BOTH count and fraction, not their logical OR.
  y <- x;y[1,2] <- 4L;y[2,2] <- 26L;y[3,2] <- 0L
  stopifnot(!rna_gain_loss_state(y,ctx,20L)$gain["JUNC_target",1])
  # Exact target conservation and deterministic same-seed sampling.
  set.seed(41);r1 <- rna_gain_loss_rarefy(x,100L)
  set.seed(41);r2 <- rna_gain_loss_rarefy(x,100L)
  stopifnot(identical(r1,r2),all(colSums(r1)==100L),all(r1<=x))
  # Exercise full implementation, no-event handling, shared depth and resume.
  d <- tempfile("gain_loss_test_");dir.create(d)
  on.exit(unlink(d,recursive=TRUE),add=TRUE)
  one <- run_gain_loss(inp,cfg,cohorts,d)
  two <- run_gain_loss(inp,cfg,cohorts,d)
  stopifnot(one$complete,two$complete,one$rarefaction$n_completed==2L,
    identical(one$rarefaction$gain_counts,two$rarefaction$gain_counts),
    all(one$tables$normalization$normalization_target_used==one$rarefaction_target),
    all(one$depth_robust_gain_loss$isoform_harmonized_id %in% one$depth_robust_candidates$isoform_harmonized_id))
  changed <- cfg;changed$execution$seeds$gain_loss <- 20260904L
  three <- run_gain_loss(inp,changed,cohorts,file.path(d,"changed_seed"))
  stopifnot(!identical(one$fingerprint,three$fingerprint))
  none <- cfg;none$gain_loss$isoform_presence_min_fl <- 200L
  none$gain_loss$gene_evaluable_min_fl <- 200L
  empty <- run_gain_loss(inp,none,cohorts,file.path(d,"no_events"))
  stopifnot(empty$complete,nrow(empty$depth_robust_candidates)==0L,
    nrow(empty$depth_robust_gain_loss)==0L,
    all(empty$tables$robust_summary$depth_robust==0L))
  message("Gain/loss tests passed.")
}
test_gain_loss()
