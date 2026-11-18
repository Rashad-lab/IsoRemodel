# Depth-controlled paired isoform gain/loss.
# Adapted from CancerCommunications_revision_Module02 (workflow 0.1.1).
# All full-length counts are rarefied before restricting candidate isoforms.

rna_gain_loss_context <- function(inputs, cfg, cohorts) {
  x <- as.matrix(inputs$counts_matrix)
  a <- as.data.frame(inputs$isoform_annotation, stringsAsFactors = FALSE)
  s <- as.data.frame(inputs$sample_manifest, stringsAsFactors = FALSE)
  if (!identical(rownames(x), as.character(a$isoform_harmonized_id)))
    stop("Gain/loss count rows and annotation IDs are not aligned.")
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(abs(x-round(x)) > 1e-8))
    stop("Gain/loss requires finite nonnegative integer FL counts.")
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      anyDuplicated(s$sample) || !setequal(colnames(x), s$sample))
    stop("Gain/loss requires unique, matching count and metadata sample IDs.")
  s <- s[match(colnames(x), s$sample), , drop = FALSE]
  patients <- unique(as.character(s$patient))
  if (!length(cohorts) || is.null(names(cohorts)) || any(!nzchar(names(cohorts))) ||
      anyDuplicated(names(cohorts))) stop("Named, unique cohorts are required.")
  if (!setequal(cohorts[[1L]], patients))
    stop("The first gain/loss cohort must contain all imported subjects.")
  pairs <- lapply(patients, function(p) {
    z <- s[s$patient == p, , drop=FALSE]
    if (nrow(z) != 2L || !setequal(z$timepoint, c("primary", "recurrent")))
      stop("Gain/loss requires one reference/comparison sample for subject ", p)
    c(primary=match(z$sample[z$timepoint == "primary"],colnames(x)),
      recurrent=match(z$sample[z$timepoint == "recurrent"],colnames(x)))
  })
  names(pairs) <- patients
  for (id in names(cohorts)) {
    p <- as.character(cohorts[[id]])
    if (!length(p) || anyDuplicated(p) || !all(p %in% patients))
      stop("Invalid patient membership for cohort ", id)
    if (length(p) < as.integer(cfg$gain_loss$primary_patient_support))
      stop("Gain/loss cohort ", id, " has fewer subjects than its configured support requirement.")
  }
  gene_mapped <- !is.na(a$associated_gene) & nzchar(a$associated_gene) & a$associated_gene != "unknown"
  if ("gene_mapped" %in% names(a)) gene_mapped <- gene_mapped & a$gene_mapped %in% TRUE
  categories <- as.character(unlist(cfg$gain_loss$keep_structural_categories))
  retained <- !is.na(a$structural_category) & a$structural_category %in% categories
  comparable <- if ("cross_sample_comparable" %in% names(a)) a$cross_sample_comparable %in% TRUE else startsWith(a$isoform_harmonized_id,"JUNC_")
  legacy <- gene_mapped & retained
  eligible <- legacy & (!isTRUE(cfg$gain_loss$require_cross_sample_comparable_id) | comparable)
  a$gene_mapped <- gene_mapped
  a$cross_sample_comparable <- comparable
  a$structural_category_retained <- retained
  a$legacy_inclusive_gain_loss_eligible <- legacy
  a$analyzed_for_gain_loss <- eligible
  genes <- sort(unique(as.character(a$associated_gene[gene_mapped])))
  gi <- match(a$associated_gene, genes)
  valid <- which(gene_mapped)
  projection <- Matrix::sparseMatrix(i=gi[valid], j=valid, x=1,
    dims=c(length(genes), nrow(x)), dimnames=list(genes,rownames(x)))
  list(counts=x, annotation=a, samples=s, patients=patients, pairs=pairs,
       genes=genes, gene_index=gi, projection=projection,
       eligible=which(eligible), legacy_eligible=which(legacy),
       cohorts=cohorts, settings=cfg$gain_loss)
}

rna_gain_loss_state <- function(counts, context, gene_floor, return_events=FALSE,
                                eligible_index=context$eligible) {
  a <- context$annotation; ids <- rownames(counts); np <- length(context$patients)
  gene_total <- as.matrix(context$projection %*% counts)
  blank <- matrix(FALSE, nrow(counts), np, dimnames=list(ids,context$patients))
  gain <- loss <- evaluable <- present_primary <- present_recurrent <- blank
  summaries <- vector("list",np); events <- if (return_events) vector("list",np) else NULL
  for (j in seq_along(context$patients)) {
    p <- context$patients[j]; pair <- context$pairs[[p]]
    ip <- unname(pair["primary"]); ir <- unname(pair["recurrent"])
    eg <- gene_total[,ip] >= gene_floor & gene_total[,ir] >= gene_floor
    ii <- eligible_index[eg[context$gene_index[eligible_index]]]
    gp <- gene_total[context$gene_index[ii],ip]; gr <- gene_total[context$gene_index[ii],ir]
    cp <- counts[ii,ip]; cr <- counts[ii,ir]
    fp <- cp/gp; fr <- cr/gr
    pp <- cp >= context$settings$isoform_presence_min_fl & fp >= context$settings$isoform_presence_min_fraction
    pr <- cr >= context$settings$isoform_presence_min_fl & fr >= context$settings$isoform_presence_min_fraction
    ga <- pr & !pp; lo <- pp & !pr
    evaluable[ii,j] <- TRUE; present_primary[ii,j] <- pp; present_recurrent[ii,j] <- pr
    gain[ii,j] <- ga; loss[ii,j] <- lo
    summaries[[j]] <- data.frame(patient=p, gene_evaluable_min_fl=gene_floor,
      n_evaluable_genes=sum(eg), n_evaluable_isoforms=length(ii),
      n_gained=sum(ga), n_lost=sum(lo), n_shared_present=sum(pp & pr),
      n_absent_both=sum(!pp & !pr), effective_depth_primary=sum(gene_total[eg,ip]),
      effective_depth_recurrent=sum(gene_total[eg,ir]),
      raw_depth_primary=sum(counts[,ip]),raw_depth_recurrent=sum(counts[,ir]))
    if (return_events) events[[j]] <- data.frame(
      isoform_harmonized_id=ids[ii],associated_gene=a$associated_gene[ii],
      structural_category=a$structural_category[ii],patient=rep(p,length(ii)),
      primary_sample=rep(colnames(counts)[ip],length(ii)),
      recurrent_sample=rep(colnames(counts)[ir],length(ii)),gene_evaluable=rep(TRUE,length(ii)),
      FL_primary=as.integer(cp),FL_recurrent=as.integer(cr),
      gene_total_FL_primary=as.numeric(gp),gene_total_FL_recurrent=as.numeric(gr),
      isoform_fraction_primary=as.numeric(fp),isoform_fraction_recurrent=as.numeric(fr),
      present_primary=pp,present_recurrent=pr,
      event=ifelse(ga,"recurrent_gained",ifelse(lo,"primary_lost",ifelse(pp & pr,"shared_present","absent_both"))),
      stringsAsFactors=FALSE)
  }
  list(gain=gain,loss=loss,evaluable=evaluable,present_primary=present_primary,
       present_recurrent=present_recurrent,patient_summary=dplyr::bind_rows(summaries),
       events=if(return_events) dplyr::bind_rows(events) else NULL,gene_totals=gene_total)
}

rna_gain_loss_support <- function(state, patients) {
  if (!all(patients %in% colnames(state$gain))) stop("Unknown subject in gain/loss cohort.")
  lapply(state[c("gain","loss","evaluable","present_primary","present_recurrent")],
         function(x) as.integer(rowSums(x[,patients,drop=FALSE])))
}

rna_gain_loss_call <- function(gain,loss,threshold,direction,opposed_maximum=0L) {
  if (identical(direction,"recurrent_gained")) return(gain >= threshold & loss <= opposed_maximum)
  if (identical(direction,"primary_lost")) return(loss >= threshold & gain <= opposed_maximum)
  stop("Unknown gain/loss direction: ",direction)
}

rna_gain_loss_rarefy <- function(counts,target) {
  if (!is.finite(target) || target < 1 || target > .Machine$integer.max || target != floor(target))
    stop("Rarefaction target must be a positive representable integer.")
  if (any(colSums(counts) < target)) stop("A library is shallower than the requested common depth.")
  z <- t(vegan::rrarefy(t(counts), sample=as.integer(target)))
  dimnames(z) <- dimnames(counts);storage.mode(z) <- "integer"
  if (any(colSums(z) != target)) stop("Exact rarefaction count-conservation assertion failed.")
  z
}

run_gain_loss <- function(inputs,cfg,cohorts,stage_dir) {
  rna_require(c("Matrix","vegan","dplyr","tidyr","tibble","ggplot2","digest"))
  dirs <- rna_dirs(stage_dir)
  ctx <- rna_gain_loss_context(inputs,cfg,cohorts)
  g <- cfg$gain_loss
  g$n_downsample_replicates <- as.integer(g$n_downsample_replicates)
  g$n_depth_curve_replicates <- as.integer(g$n_depth_curve_replicates)
  g$checkpoint_every <- as.integer(g$checkpoint_every)
  g$patient_support_thresholds <- sort(unique(c(as.integer(unlist(g$patient_support_thresholds)),as.integer(g$primary_patient_support))))
  g$gene_floor_sensitivity <- sort(unique(c(as.numeric(unlist(g$gene_floor_sensitivity)),g$gene_evaluable_min_fl)))
  g$depth_curve_fractions <- sort(unique(c(as.numeric(unlist(g$depth_curve_fractions)),1)))
  g$depth_robust_probability_sensitivity <- sort(unique(c(as.numeric(unlist(g$depth_robust_probability_sensitivity)),g$depth_robust_probability)))
  opposed <- as.integer(g$opposite_direction_maximum %||% 0L)
  if (g$n_downsample_replicates < 2L || g$n_depth_curve_replicates < 1L ||
      g$n_depth_curve_replicates > g$n_downsample_replicates || g$checkpoint_every < 1L ||
      length(g$isoform_presence_min_fl)!=1L || !is.finite(g$isoform_presence_min_fl) ||
      g$isoform_presence_min_fl < 1L ||
      length(g$isoform_presence_min_fraction)!=1L || !is.finite(g$isoform_presence_min_fraction) ||
      g$isoform_presence_min_fraction <= 0 || g$isoform_presence_min_fraction > 1 ||
      length(g$gene_evaluable_min_fl)!=1L || !is.finite(g$gene_evaluable_min_fl) ||
      g$gene_evaluable_min_fl < g$isoform_presence_min_fl ||
      any(g$patient_support_thresholds < 1L) || any(g$patient_support_thresholds > length(ctx$patients)) ||
      !is.finite(g$depth_robust_probability) || g$depth_robust_probability <= 0 ||
      g$depth_robust_probability > 1 || any(g$depth_curve_fractions <= 0 | g$depth_curve_fractions > 1))
    stop("Invalid gain/loss simulation/support settings.")
  ctx$settings <- g
  x <- ctx$counts; ids <- rownames(x); a <- ctx$annotation; np <- length(ctx$patients)
  nr <- nrow(x);seed <- as.integer(cfg$execution$seeds$gain_loss %||% 20260903L)
  primary_samples <- ctx$samples$sample[ctx$samples$patient %in% cohorts[[1L]]]
  target <- min(colSums(x[,primary_samples,drop=FALSE]))
  if (!is.finite(target) || target < 1) stop("At least one primary-cohort library is empty.")
  code_identity <- rna_hash(list(context=body(rna_gain_loss_context),state=body(rna_gain_loss_state),
    support=body(rna_gain_loss_support),call=body(rna_gain_loss_call),
    rarefy=body(rna_gain_loss_rarefy),runner=body(run_gain_loss)))
  fingerprint <- rna_hash(list(engine="gain_loss_0.1.0",code=code_identity,input=inputs$fingerprint,
    counts=x,annotation=a,cohorts=cohorts,settings=g,seed=seed,
    rng=c("L'Ecuyer-CMRG","Inversion","Rejection"),r_version=as.character(getRversion()),
    target=target,vegan=as.character(utils::packageVersion("vegan"))))
  write <- function(z,name,qc=FALSE) rna_write_tsv(z,file.path(if(qc) dirs$qc else dirs$tables,name))
  call <- function(s,k,d) rna_gain_loss_call(s$gain,s$loss,k,d,opposed)
  directions <- c("recurrent_gained","primary_lost")
  support_rows <- function(state,depth,replicate=NA_integer_,floor=g$gene_evaluable_min_fl) {
    rows <- list();n <- 0L
    for (cn in names(cohorts)) {
      su <- rna_gain_loss_support(state,cohorts[[cn]])
      for (k in g$patient_support_thresholds) for (d in directions) {
        n <- n+1L; rows[[n]] <- data.frame(cohort=cn,n_pairs=length(cohorts[[cn]]),
          analysis_depth=depth,replicate=replicate,gene_evaluable_min_fl=floor,
          support_threshold=k,direction=d,n_high_confidence_isoforms=sum(call(su,k,d)),
          opposite_direction_maximum=opposed,stringsAsFactors=FALSE)
      }
    }
    dplyr::bind_rows(rows)
  }
  detection <- function(counts,r,target_depth,fraction,depth) {
    gt <- as.matrix(ctx$projection %*% counts)
    z <- data.frame(replicate=r,analysis_depth=depth,target_depth=target_depth,
      depth_fraction_of_common_target=fraction,sample=colnames(counts),
      observed_total_FL=as.numeric(colSums(counts)),
      detected_isoforms_FLge1=as.integer(colSums(counts>=1)),
      retained_category_isoforms_FLge1=as.integer(colSums(counts[ctx$eligible,,drop=FALSE]>=1)),
      detected_genes_FLge1=as.integer(colSums(gt>=1)),
      genes_FLge_primary_floor=as.integer(colSums(gt>=g$gene_evaluable_min_fl)))
    dplyr::left_join(z,ctx$samples,by="sample")
  }
  full <- rna_gain_loss_state(x,ctx,g$gene_evaluable_min_fl,TRUE)
  write(full$events,"paired_isoform_events_full_depth.tsv.gz")
  full$events <- NULL;full$gene_totals <- NULL
  write(full$patient_summary,"patient_gain_loss_full_depth.tsv")
  depth <- detection(x,NA_integer_,target,1,"full_depth")
  depth$total_full_length_reads <- depth$observed_total_FL
  depth$rarefaction_target_FL <- target
  depth$retained_depth_fraction <- target/depth$observed_total_FL
  write(depth,"sample_depth_detection_summary.tsv")
  normalization <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    ss <- ctx$samples$sample[ctx$samples$patient %in% cohorts[[cn]]]
    data.frame(cohort=cn,n_pairs=length(cohorts[[cn]]),
      shallowest_library_within_cohort=min(colSums(x[,ss,drop=FALSE])),
      normalization_target_used=target,reference_cohort=names(cohorts)[1L])
  }))
  write(normalization,"normalization_target_audit.tsv",TRUE)
  full_counts <- support_rows(full,"full_depth")
  write(full_counts,"support_threshold_counts_full_depth.tsv")
  full_support <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    su <- rna_gain_loss_support(full,cohorts[[cn]])
    z <- data.frame(isoform_harmonized_id=ids,cohort=cn,n_pairs=length(cohorts[[cn]]),
      n_gain_patients=su$gain,n_loss_patients=su$loss,n_evaluable_patients=su$evaluable,
      n_present_primary_patients=su$present_primary,n_present_recurrent_patients=su$present_recurrent)
    dplyr::left_join(z,a,by="isoform_harmonized_id")
  }))
  write(full_support,"isoform_cross_patient_support_full_depth.tsv.gz")
  full_sets <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    su <- rna_gain_loss_support(full,cohorts[[cn]])
    dplyr::bind_rows(lapply(g$patient_support_thresholds,function(k)
      dplyr::bind_rows(lapply(directions,function(d) {
        ii <- which(call(su,k,d));z <- data.frame(isoform_harmonized_id=ids[ii],
          cohort=rep(cn,length(ii)),n_pairs=rep(length(cohorts[[cn]]),length(ii)),
          support_threshold=rep(k,length(ii)),direction=rep(d,length(ii)),
          n_gain_patients=su$gain[ii],n_loss_patients=su$loss[ii],n_evaluable_patients=su$evaluable[ii])
        dplyr::left_join(z,a,by="isoform_harmonized_id")
      }))))
  }))
  write(full_sets,"high_confidence_gain_loss_full_depth.tsv.gz")
  floors_patient <- floors_cohort <- list()
  for (v in g$gene_floor_sensitivity) {
    st <- if(v == g$gene_evaluable_min_fl) full else rna_gain_loss_state(x,ctx,v)
    floors_patient[[as.character(v)]] <- st$patient_summary
    floors_cohort[[as.character(v)]] <- support_rows(st,"full_depth",floor=v)
  }
  floor_patient <- dplyr::bind_rows(floors_patient);floor_cohort <- dplyr::bind_rows(floors_cohort)
  write(floor_patient,"gene_floor_patient_sensitivity.tsv")
  write(floor_cohort,"gene_floor_support_sensitivity.tsv")
  legacy <- rna_gain_loss_state(x,ctx,g$gene_evaluable_min_fl,eligible_index=ctx$legacy_eligible)
  identifier_patient <- dplyr::bind_rows(
    dplyr::mutate(full$patient_summary,analysis_universe="comparable_junction_chains"),
    dplyr::mutate(legacy$patient_summary,analysis_universe="including_sample_specific_ids"))
  identifier_support <- dplyr::bind_rows(support_rows(full,"comparable_junction_chains"),
                                        support_rows(legacy,"including_sample_specific_ids"))
  write(identifier_patient,"sample_specific_fallback_sensitivity_by_patient.tsv")
  write(identifier_support,"sample_specific_fallback_support_sensitivity.tsv")
  rm(legacy,floors_patient,floors_cohort)

  cache <- file.path(dirs$cache,paste0("rarefaction_",substr(fingerprint,1,20),".rds"))
  new_state <- function() list(fingerprint=fingerprint,n_completed=0L,
    gain_counts=matrix(0L,nr,np,dimnames=list(NULL,ctx$patients)),
    loss_counts=matrix(0L,nr,np,dimnames=list(NULL,ctx$patients)),
    evaluable_counts=matrix(0L,nr,np,dimnames=list(NULL,ctx$patients)),
    high_confidence_counts=setNames(lapply(names(cohorts),function(z)
      setNames(lapply(directions,function(d) integer(nr)),directions)),names(cohorts)),
    support_threshold_counts=list(),patient_event_counts=list(),
    common_target_detection=list(),depth_curve_event_counts=list(),
    depth_curve_detection=list(),depth_audit=list())
  acc <- if(isTRUE(cfg$execution$resume) && !isTRUE(cfg$execution$force_refit) && file.exists(cache)) readRDS(cache) else new_state()
  if (!identical(acc$fingerprint,fingerprint) || acc$n_completed > g$n_downsample_replicates)
    stop("Gain/loss cache fingerprint or completed-replicate count is invalid.")
  oldkind <- RNGkind();hadseed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadseed) oldseed <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit({do.call(RNGkind,as.list(oldkind));if(hadseed) assign(".Random.seed",oldseed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("L'Ecuyer-CMRG",normal.kind="Inversion",sample.kind="Rejection")
  first <- acc$n_completed+1L
  if(first <= g$n_downsample_replicates) for(r in seq.int(first,g$n_downsample_replicates)) {
    set.seed(seed+r)
    rare <- rna_gain_loss_rarefy(x,target)
    st <- rna_gain_loss_state(rare,ctx,g$gene_evaluable_min_fl)
    acc$gain_counts <- acc$gain_counts+st$gain
    acc$loss_counts <- acc$loss_counts+st$loss
    acc$evaluable_counts <- acc$evaluable_counts+st$evaluable
    for(cn in names(cohorts)) {
      su <- rna_gain_loss_support(st,cohorts[[cn]])
      for(d in directions) acc$high_confidence_counts[[cn]][[d]] <-
        acc$high_confidence_counts[[cn]][[d]]+as.integer(call(su,g$primary_patient_support,d))
    }
    acc$support_threshold_counts[[r]] <- support_rows(st,"common_depth_rarefaction",r)
    acc$patient_event_counts[[r]] <- dplyr::mutate(st$patient_summary,replicate=r,analysis_depth="common_depth_rarefaction")
    acc$common_target_detection[[r]] <- detection(rare,r,target,1,"common_depth_rarefaction")
    acc$depth_audit[[r]] <- data.frame(replicate=r,minimum_sample_depth=min(colSums(rare)),
      maximum_sample_depth=max(colSums(rare)),target_depth=target,exact_target_pass=all(colSums(rare)==target))
    if(r <= g$n_depth_curve_replicates) for(fi in seq_along(g$depth_curve_fractions)) {
      f <- g$depth_curve_fractions[fi];ct <- max(1L,as.integer(floor(target*f)))
      if(f==1) {cx <- rare;cs <- st} else {
        set.seed(seed+1000000L+r*100L+fi)
        cx <- rna_gain_loss_rarefy(rare,ct);cs <- rna_gain_loss_state(cx,ctx,g$gene_evaluable_min_fl)
      }
      nm <- paste(r,fi,sep="_")
      acc$depth_curve_event_counts[[nm]] <- dplyr::mutate(cs$patient_summary,
        replicate=r,target_depth=ct,depth_fraction_of_common_target=f)
      acc$depth_curve_detection[[nm]] <- detection(cx,r,ct,f,"depth_curve")
    }
    acc$n_completed <- r
    if(r %% g$checkpoint_every == 0L || r == g$n_downsample_replicates) {
      rna_atomic_rds(acc,cache);rna_log("Gain/loss rarefaction ",r,"/",g$n_downsample_replicates)
    }
  }
  if(acc$n_completed != g$n_downsample_replicates) stop("Gain/loss rarefactions are incomplete.")
  B <- acc$n_completed
  for(nm in c("support_threshold_counts","patient_event_counts","common_target_detection",
               "depth_curve_event_counts","depth_curve_detection","depth_audit")) {
    if(is.list(acc[[nm]]) && !is.data.frame(acc[[nm]])) acc[[nm]] <- dplyr::bind_rows(acc[[nm]])
  }
  write(acc$depth_audit,"rarefaction_exact_depth_audit.tsv",TRUE)
  write(acc$patient_event_counts,"patient_gain_loss_common_depth_replicates.tsv.gz")
  write(acc$support_threshold_counts,"support_threshold_counts_common_depth_replicates.tsv.gz")
  write(acc$common_target_detection,"sample_detection_common_depth_replicates.tsv.gz")
  write(acc$depth_curve_event_counts,"gain_loss_depth_curve_replicates.tsv.gz")
  write(acc$depth_curve_detection,"sample_detection_depth_curve_replicates.tsv.gz")

  candidates <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    su <- rna_gain_loss_support(full,cohorts[[cn]])
    dplyr::bind_rows(lapply(directions,function(d) {
      ii <- which(call(su,g$primary_patient_support,d));n <- length(ii)
      probability <- acc$high_confidence_counts[[cn]][[d]][ii]/B
      z <- data.frame(isoform_harmonized_id=ids[ii],cohort=rep(cn,n),n_pairs=rep(length(cohorts[[cn]]),n),
        direction=rep(d,n),support_threshold=rep(g$primary_patient_support,n),
        full_n_gain_patients=su$gain[ii],full_n_loss_patients=su$loss[ii],
        full_n_evaluable_patients=su$evaluable[ii],rarefaction_high_confidence_probability=probability,
        depth_robust_probability_threshold=rep(g$depth_robust_probability,n),
        depth_robust=probability >= g$depth_robust_probability)
      dplyr::left_join(z,a,by="isoform_harmonized_id")
    }))
  }))
  robust <- candidates[candidates$depth_robust,,drop=FALSE]
  write(candidates,"depth_robust_gain_loss_candidates.tsv.gz")
  write(robust,"depth_robust_gain_loss.tsv.gz")
  stability <- dplyr::bind_rows(lapply(seq_along(ctx$patients),function(j) {
    ii <- ctx$eligible;n <- length(ii);ev <- acc$evaluable_counts[ii,j]
    data.frame(isoform_harmonized_id=ids[ii],patient=rep(ctx$patients[j],n),
      p_evaluable=ev/B,p_gain=acc$gain_counts[ii,j]/B,p_loss=acc$loss_counts[ii,j]/B,
      p_gain_given_evaluable=ifelse(ev>0,acc$gain_counts[ii,j]/pmax(ev,1),NA_real_),
      p_loss_given_evaluable=ifelse(ev>0,acc$loss_counts[ii,j]/pmax(ev,1),NA_real_))
  }))
  write(stability,"isoform_patient_rarefaction_stability.tsv.gz")
  trace <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    rr <- robust[robust$cohort==cn,c("isoform_harmonized_id","associated_gene","structural_category","cohort","direction"),drop=FALSE]
    if(!nrow(rr)) return(NULL)
    zz <- merge(rr,data.frame(patient=cohorts[[cn]]),by=NULL)
    ix <- cbind(match(zz$isoform_harmonized_id,ids),match(zz$patient,ctx$patients))
    zz$full_depth_patient_event <- ifelse(full$gain[ix],"recurrent_gained",
      ifelse(full$loss[ix],"primary_lost",ifelse(full$evaluable[ix],"evaluable_no_change","not_evaluable")))
    dplyr::left_join(zz,stability,by=c("isoform_harmonized_id","patient"))
  }))
  if(!ncol(trace)) trace <- data.frame(isoform_harmonized_id=character(),associated_gene=character(),
    structural_category=character(),cohort=character(),direction=character(),patient=character(),
    full_depth_patient_event=character(),p_evaluable=double(),p_gain=double(),p_loss=double(),
    p_gain_given_evaluable=double(),p_loss_given_evaluable=double())
  write(trace,"depth_robust_gain_loss_patient_trace.tsv.gz")
  robust_summary <- dplyr::bind_rows(lapply(names(cohorts),function(cn)
    dplyr::bind_rows(lapply(directions,function(d) {
      cc <- candidates[candidates$cohort==cn & candidates$direction==d,,drop=FALSE]
      data.frame(cohort=cn,n_pairs=length(cohorts[[cn]]),direction=d,
        full_depth_high_confidence=nrow(cc),depth_robust=sum(cc$depth_robust),
        retained_fraction=if(nrow(cc)) mean(cc$depth_robust) else NA_real_)
    }))))
  write(robust_summary,"depth_robust_gain_loss_summary.tsv")
  probability_sensitivity <- dplyr::bind_rows(lapply(names(cohorts),function(cn)
    dplyr::bind_rows(lapply(directions,function(d)
      dplyr::bind_rows(lapply(g$depth_robust_probability_sensitivity,function(thr) {
        z <- candidates[candidates$cohort==cn & candidates$direction==d,,drop=FALSE]
        data.frame(cohort=cn,n_pairs=length(cohorts[[cn]]),direction=d,stability_threshold=thr,
          n_full_depth_high_confidence=nrow(z),n_depth_robust=sum(z$rarefaction_high_confidence_probability>=thr),
          retained_fraction=if(nrow(z)) mean(z$rarefaction_high_confidence_probability>=thr) else NA_real_)
      }))))))
  write(probability_sensitivity,"depth_robust_probability_sensitivity.tsv")
  support_summary <- acc$support_threshold_counts |>
    dplyr::group_by(.data$cohort,.data$n_pairs,.data$support_threshold,.data$direction) |>
    dplyr::summarise(n_replicates=dplyr::n_distinct(.data$replicate),
      mean_high_confidence_isoforms=mean(.data$n_high_confidence_isoforms),
      median_high_confidence_isoforms=stats::median(.data$n_high_confidence_isoforms),
      lower_95=as.numeric(stats::quantile(.data$n_high_confidence_isoforms,.025,names=FALSE)),
      upper_95=as.numeric(stats::quantile(.data$n_high_confidence_isoforms,.975,names=FALSE)),.groups="drop") |>
    dplyr::left_join(full_counts[,c("cohort","n_pairs","support_threshold","direction","n_high_confidence_isoforms")],
      by=c("cohort","n_pairs","support_threshold","direction"))
  write(support_summary,"support_threshold_rarefaction_summary.tsv")
  patient_summary <- acc$patient_event_counts |>
    tidyr::pivot_longer(cols=c("n_gained","n_lost"),names_to="direction",values_to="n_isoform_events") |>
    dplyr::mutate(direction=ifelse(.data$direction=="n_gained","recurrent_gained","primary_lost")) |>
    dplyr::group_by(.data$patient,.data$direction) |>
    dplyr::summarise(n_replicates=dplyr::n_distinct(.data$replicate),mean_events=mean(.data$n_isoform_events),
      median_events=stats::median(.data$n_isoform_events),
      lower_95=as.numeric(stats::quantile(.data$n_isoform_events,.025,names=FALSE)),
      upper_95=as.numeric(stats::quantile(.data$n_isoform_events,.975,names=FALSE)),
      mean_evaluable_genes=mean(.data$n_evaluable_genes),
      mean_evaluable_isoforms=mean(.data$n_evaluable_isoforms),
      mean_effective_depth_primary=mean(.data$effective_depth_primary),
      mean_effective_depth_recurrent=mean(.data$effective_depth_recurrent),.groups="drop")
  write(patient_summary,"patient_gain_loss_rarefaction_summary.tsv")
  cohort_burden_reps <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    z <- acc$patient_event_counts[acc$patient_event_counts$patient %in% cohorts[[cn]],,drop=FALSE]
    z <- dplyr::summarise(dplyr::group_by(z,.data$replicate),n_gained=sum(.data$n_gained),n_lost=sum(.data$n_lost),.groups="drop")
    z$cohort <- cn;z$n_pairs <- length(cohorts[[cn]])
    tidyr::pivot_longer(z,cols=c("n_gained","n_lost"),names_to="direction",values_to="common_depth_events") |>
      dplyr::mutate(direction=ifelse(.data$direction=="n_gained","recurrent_gained","primary_lost"))
  }))
  cohort_summary <- cohort_burden_reps |>
    dplyr::group_by(.data$cohort,.data$n_pairs,.data$direction) |>
    dplyr::summarise(n_replicates=dplyr::n_distinct(.data$replicate),
      mean_common_depth_events=mean(.data$common_depth_events),
      median_common_depth_events=stats::median(.data$common_depth_events),
      lower_95=as.numeric(stats::quantile(.data$common_depth_events,.025,names=FALSE)),
      upper_95=as.numeric(stats::quantile(.data$common_depth_events,.975,names=FALSE)),.groups="drop")
  cohort_summary$full_depth_events <- vapply(seq_len(nrow(cohort_summary)),function(i) {
    cn <- cohort_summary$cohort[i];d <- cohort_summary$direction[i]
    ps <- full$patient_summary[full$patient_summary$patient %in% cohorts[[cn]],,drop=FALSE]
    sum(ps[[if(d=="recurrent_gained") "n_gained" else "n_lost"]])
  },numeric(1))
  write(cohort_burden_reps,"cohort_gain_loss_common_depth_replicates.tsv.gz")
  write(cohort_summary,"cohort_gain_loss_full_vs_common_depth.tsv")

  # Fisher tests use robust events versus non-events in the pair-evaluable
  # eligible-isoform background, with a separate category family per direction.
  enrichment <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    bg <- a$analyzed_for_gain_loss %in% TRUE & rowSums(full$evaluable[,cohorts[[cn]],drop=FALSE])>0
    bgids <- ids[bg]; bgcats <- a$structural_category[bg]
    dplyr::bind_rows(lapply(directions,function(d) {
      ev <- bgids %in% robust$isoform_harmonized_id[robust$cohort==cn & robust$direction==d]
      out <- dplyr::bind_rows(lapply(as.character(unlist(g$keep_structural_categories)),function(cat) {
        ca <- bgcats==cat;aa <- sum(ev & ca);bb <- sum(ev & !ca);cc <- sum(!ev & ca);dd <- sum(!ev & !ca)
        ft <- if((aa+bb)>0 && (cc+dd)>0) tryCatch(stats::fisher.test(matrix(c(aa,bb,cc,dd),2,byrow=TRUE)),error=function(e) NULL) else NULL
        n <- sum(ev);nb <- length(ev);nc <- sum(ca)
        data.frame(cohort=cn,n_pairs=length(cohorts[[cn]]),direction=d,structural_category=cat,
          n_event_category=aa,n_event_total=n,n_background_category=nc,n_background_total=nb,
          event_fraction=if(n>0) aa/n else NA_real_,background_fraction=if(nb>0) nc/nb else NA_real_,
          fold_enrichment=if(n>0 && nc>0) (aa/n)/(nc/nb) else NA_real_,
          odds_ratio=if(!is.null(ft)) unname(ft$estimate) else NA_real_,
          p_value=if(!is.null(ft)) ft$p.value else NA_real_)
      }))
      out$FDR <- stats::p.adjust(out$p_value,"BH");out
    }))
  }))
  write(enrichment,"depth_robust_structural_category_enrichment.tsv")
  full_patient_direction <- full$patient_summary |>
    tidyr::pivot_longer(cols=c("n_gained","n_lost"),names_to="direction",values_to="n_isoform_events") |>
    dplyr::mutate(direction=ifelse(.data$direction=="n_gained","recurrent_gained","primary_lost"),
      directional_raw_depth=ifelse(.data$direction=="recurrent_gained",.data$raw_depth_recurrent,.data$raw_depth_primary),
      directional_effective_depth=ifelse(.data$direction=="recurrent_gained",.data$effective_depth_recurrent,.data$effective_depth_primary),
      minimum_pair_raw_depth=pmin(.data$raw_depth_primary,.data$raw_depth_recurrent),
      minimum_pair_effective_depth=pmin(.data$effective_depth_primary,.data$effective_depth_recurrent))
  associations <- dplyr::bind_rows(lapply(names(cohorts),function(cn)
    dplyr::bind_rows(lapply(directions,function(d)
      dplyr::bind_rows(lapply(c("directional_raw_depth","directional_effective_depth","minimum_pair_raw_depth","minimum_pair_effective_depth"),function(v)
        dplyr::bind_rows(lapply(c("pearson","spearman"),function(method) {
          z <- full_patient_direction[full_patient_direction$patient %in% cohorts[[cn]] & full_patient_direction$direction==d,,drop=FALSE]
          ok <- is.finite(z[[v]]) & is.finite(z$n_isoform_events);n <- sum(ok)
          ct <- if(n>=4L && length(unique(z[[v]][ok]))>1L && length(unique(z$n_isoform_events[ok]))>1L)
            tryCatch(suppressWarnings(stats::cor.test(z[[v]][ok],z$n_isoform_events[ok],method=method,exact=FALSE)),error=function(e) NULL) else NULL
          data.frame(cohort=cn,direction=d,depth_metric=v,correlation=method,n=n,
            estimate=if(!is.null(ct)) unname(ct$estimate) else NA_real_,
            p_value=if(!is.null(ct)) ct$p.value else NA_real_,status=if(is.null(ct)) "unavailable" else "tested")
        })))))))) |>
    dplyr::group_by(.data$cohort,.data$correlation) |>
    dplyr::mutate(FDR=stats::p.adjust(.data$p_value,"BH")) |>
    dplyr::ungroup()
  write(associations,"gain_loss_depth_associations.tsv")
  asymmetry <- dplyr::bind_rows(lapply(names(cohorts),function(cn) {
    z <- full$patient_summary[full$patient_summary$patient %in% cohorts[[cn]],,drop=FALSE]
    w <- tryCatch(suppressWarnings(stats::wilcox.test(z$n_gained,z$n_lost,paired=TRUE,exact=FALSE)),error=function(e) NULL)
    data.frame(cohort=cn,n_pairs=nrow(z),full_depth_gained=sum(z$n_gained),full_depth_lost=sum(z$n_lost),
      paired_wilcoxon_p=if(is.null(w)) NA_real_ else w$p.value,
      n_rarefactions=B,n_rarefactions_losses_exceed_gains=sum(
        cohort_burden_reps$common_depth_events[cohort_burden_reps$cohort==cn & cohort_burden_reps$direction=="primary_lost"] >
        cohort_burden_reps$common_depth_events[cohort_burden_reps$cohort==cn & cohort_burden_reps$direction=="recurrent_gained"]),
      interpretation="Rarefactions reuse biological specimens; replicate asymmetry is descriptive.")
  }))
  write(asymmetry,"gain_loss_asymmetry_summary.tsv")

  tables <- list(full_support_counts=full_counts,full_support=full_support,full_sets=full_sets,
    sample_depth=depth,normalization=normalization,gene_floor_patient=floor_patient,
    gene_floor_support=floor_cohort,identifier_patient=identifier_patient,
    identifier_support=identifier_support,robust_summary=robust_summary,
    probability_sensitivity=probability_sensitivity,structural_enrichment=enrichment,
    full_patient_direction=full_patient_direction,depth_associations=associations,
    asymmetry=asymmetry,cohort_burden_replicates=cohort_burden_reps)
  result <- list(engine_version="0.1.0",fingerprint=fingerprint,complete=TRUE,
    cohorts=cohorts,rarefaction_target=target,full_state=full,
    rarefaction=acc,depth_robust_gain_loss=robust,depth_robust_candidates=candidates,
    depth_robust_patient_trace=trace,isoform_patient_stability=stability,
    cohort_summary=cohort_summary,support_summary=support_summary,
    patient_summary=patient_summary,tables=tables,paths=dirs,settings=g)
  rna_atomic_rds(result,file.path(dirs$objects,"gain_loss_results.rds"))
  plot_gain_loss_results(result,cfg,dirs)
  gates <- data.frame(gate=c("All rarefaction replicates completed","All rarefactions have exact target depth","Robust set is subset of full-depth candidates"),
    pass=c(B==g$n_downsample_replicates,all(acc$depth_audit$exact_target_pass),all(robust$isoform_harmonized_id %in% candidates$isoform_harmonized_id)))
  write(gates,"completion_gates.tsv",TRUE)
  if(any(!gates$pass)) stop("Gain/loss completion checks failed.")
  result
}

plot_gain_loss_results <- function(result,cfg,dirs=result$paths) {
  if(identical(cfg$reporting$enabled,FALSE)) return(invisible(NULL))
  rna_require(c("ggplot2","dplyr","tidyr"))
  event_palette <- unlist(cfg$reporting$palettes$events %||%
    c(recurrent_gained="#D55E00",primary_lost="#0072B2"),use.names=TRUE)
  timepoint_palette <- unlist(cfg$reporting$palettes$timepoints %||%
    cfg$reporting$palettes$timepoint %||%
    c(primary="#4C78A8",recurrent="#E45756"),use.names=TRUE)
  direction_labels <- c(recurrent_gained="Gained",primary_lost="Lost")
  timepoint_labels <- c(primary=tools::toTitleCase(cfg$comparison$reference %||% "Primary"),
    recurrent=tools::toTitleCase(cfg$comparison$comparison %||% "Recurrent"))
  save <- function(p,name,w=8,h=6) {
    for(aesthetic in intersect(c("colour","fill"),names(p$mapping))) {
      variable <- paste(deparse(p$mapping[[aesthetic]]),collapse="")
      if(grepl("direction",variable,fixed=TRUE)) {
        scale <- if(aesthetic=="colour") ggplot2::scale_colour_manual else ggplot2::scale_fill_manual
        p <- p+scale(values=event_palette,labels=direction_labels)
      } else if(grepl("timepoint",variable,fixed=TRUE)) {
        scale <- if(aesthetic=="colour") ggplot2::scale_colour_manual else ggplot2::scale_fill_manual
        p <- p+scale(values=timepoint_palette,labels=timepoint_labels)
      }
    }
    xvar <- paste(deparse(p$mapping$x),collapse="")
    if(grepl("direction",xvar,fixed=TRUE)) p <- p+ggplot2::scale_x_discrete(labels=direction_labels)
    if(grepl("structural_category",xvar,fixed=TRUE)) p <- p+ggplot2::scale_x_discrete(labels=c(
      "full-splice_match"="Full splice match","incomplete-splice_match"="Incomplete splice match",
      "novel_in_catalog"="Novel in catalog","novel_not_in_catalog"="Novel not in catalog"))
    rna_save_plot(p,name,dirs,cfg,width=w,height=h)
  }
  theme <- rna_theme(cfg)
  summary <- result$patient_summary
  full <- result$full_state$patient_summary |>
    tidyr::pivot_longer(cols=c("n_gained","n_lost"),names_to="direction",values_to="full_depth_events") |>
    dplyr::mutate(direction=ifelse(.data$direction=="n_gained","recurrent_gained","primary_lost"))
  combined <- dplyr::left_join(summary,full[,c("patient","direction","full_depth_events")],by=c("patient","direction"))
  p <- ggplot2::ggplot(combined,ggplot2::aes(x=.data$patient,y=.data$mean_events,color=.data$direction))+
    ggplot2::geom_linerange(ggplot2::aes(ymin=.data$lower_95,ymax=.data$upper_95),position=ggplot2::position_dodge(.35))+
    ggplot2::geom_point(position=ggplot2::position_dodge(.35),size=3)+
    ggplot2::geom_point(ggplot2::aes(y=.data$full_depth_events),shape=1,size=3,position=ggplot2::position_dodge(.35))+
    ggplot2::labs(title="Paired isoform gain and loss",x="Patient",y="Isoform events",color="Direction",
      caption="Filled: common-depth mean and 95% range. Open: full depth.")+theme
  save(p,"paired_gain_loss_full_vs_common_depth",8,5)
  p <- ggplot2::ggplot(result$support_summary,ggplot2::aes(x=.data$support_threshold,y=.data$mean_high_confidence_isoforms,color=.data$direction))+
    ggplot2::geom_line()+ggplot2::geom_point(size=2.5)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$lower_95,ymax=.data$upper_95),width=.1)+
    ggplot2::facet_wrap(~cohort,scales="free_y")+
    ggplot2::labs(title="Cross-patient support after depth control",x="Supporting patients",y="Isoforms",color="Direction")+theme
  save(p,"gain_loss_support_threshold_sensitivity",9,5)
  p <- ggplot2::ggplot(result$tables$robust_summary,ggplot2::aes(x=.data$direction,y=.data$depth_robust,fill=.data$direction))+
    ggplot2::geom_col(width=.65)+ggplot2::geom_text(ggplot2::aes(label=.data$depth_robust),vjust=-.3)+
    ggplot2::facet_wrap(~cohort)+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.15)))+
    ggplot2::labs(title="Depth-robust isoform gain and loss",x=NULL,y="Isoforms",fill="Direction")+theme
  save(p,"depth_robust_gain_loss_counts",8,5)
  p <- ggplot2::ggplot(result$tables$probability_sensitivity,ggplot2::aes(x=.data$stability_threshold,y=.data$n_depth_robust,color=.data$direction))+
    ggplot2::geom_line()+ggplot2::geom_point(size=2.5)+ggplot2::facet_wrap(~cohort)+
    ggplot2::labs(title="Rarefaction robustness threshold",x="Required probability",y="Retained isoforms",color="Direction")+theme
  save(p,"depth_robust_probability_sensitivity",8,5)
  curves <- result$rarefaction$depth_curve_event_counts |>
    tidyr::pivot_longer(cols=c("n_gained","n_lost"),names_to="direction",values_to="events") |>
    dplyr::mutate(direction=ifelse(.data$direction=="n_gained","recurrent_gained","primary_lost")) |>
    dplyr::group_by(.data$patient,.data$direction,.data$depth_fraction_of_common_target) |>
    dplyr::summarise(mean_events=mean(.data$events),.groups="drop")
  p <- ggplot2::ggplot(curves,ggplot2::aes(x=.data$depth_fraction_of_common_target,y=.data$mean_events,color=.data$direction))+
    ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::facet_wrap(~patient,scales="free_y")+
    ggplot2::labs(title="Gain and loss across sequencing depths",x="Fraction of common depth",y="Mean isoform events",color="Direction")+theme
  save(p,"gain_loss_depth_curves",10,6)
  floor <- result$tables$gene_floor_support
  floor <- floor[floor$support_threshold==result$settings$primary_patient_support,,drop=FALSE]
  p <- ggplot2::ggplot(floor,ggplot2::aes(x=.data$gene_evaluable_min_fl,y=.data$n_high_confidence_isoforms,color=.data$direction))+
    ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::facet_wrap(~cohort)+
    ggplot2::labs(title="Gene read-support requirement",x="Minimum gene FL in both samples",y="Full-depth high-confidence isoforms",color="Direction")+theme
  save(p,"gain_loss_gene_floor_sensitivity",8,5)
  en <- result$tables$structural_enrichment
  p <- ggplot2::ggplot(en,ggplot2::aes(x=.data$structural_category,y=.data$event_fraction,fill=.data$direction))+
    ggplot2::geom_col(position="dodge",width=.7)+ggplot2::facet_wrap(~cohort)+
    ggplot2::coord_flip()+ggplot2::labs(title="Structural composition of depth-robust events",x=NULL,y="Fraction of robust isoforms",fill="Direction")+theme
  save(p,"depth_robust_structural_category_composition",9,6)
  p <- ggplot2::ggplot(en,ggplot2::aes(x=.data$structural_category,y=.data$fold_enrichment,color=.data$direction))+
    ggplot2::geom_hline(yintercept=1,linetype=2,color="grey60")+ggplot2::geom_point(size=3,na.rm=TRUE)+
    ggplot2::facet_wrap(~cohort)+ggplot2::coord_flip()+
    ggplot2::labs(title="Structural enrichment of depth-robust events",x=NULL,y="Fold enrichment",color="Direction")+theme
  save(p,"depth_robust_structural_category_enrichment",9,6)
  p <- ggplot2::ggplot(result$tables$sample_depth,ggplot2::aes(x=.data$sample,y=.data$observed_total_FL,fill=.data$timepoint))+
    ggplot2::geom_col()+ggplot2::geom_hline(yintercept=result$rarefaction_target,linetype=2)+
    ggplot2::coord_flip()+ggplot2::labs(title="Full-length RNA library depth",x=NULL,y="FL counts",fill="Condition")+theme
  save(p,"sample_depth_and_common_target",8,6)
  invisible(NULL)
}
