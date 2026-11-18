# YAML and sample metadata are the public interface to this workflow.
# Relative configuration paths use the YAML directory. Relative Pigeon paths
# use inputs.root_directory; metadata sample IDs never derive from filenames.
rna_config_absolute <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}
rna_config_path <- function(path, base) {
  if (is.null(path) || !length(path)) return(NULL)
  if (length(path) != 1L || is.na(path)) stop("A configured path must be one nonmissing string.", call.=FALSE)
  path <- trimws(as.character(path))
  if (!nzchar(path)) return(NULL)
  path <- gsub("\\\\", "/", path.expand(path))
  if (!rna_config_absolute(path)) path <- file.path(base,path)
  normalizePath(path,winslash="/",mustWork=FALSE)
}
rna_config_get <- function(cfg,key) {
  out <- cfg
  for (part in strsplit(key,".",fixed=TRUE)[[1L]]) out <- out[[part]]
  out
}
rna_config_number <- function(cfg,key,lower=0,upper=Inf,integer=FALSE,strict_lower=FALSE,required=TRUE) {
  x <- unlist(rna_config_get(cfg,key),use.names=FALSE)
  if (!length(x)) {
    if(required) stop("Missing numeric setting: ",key,call.=FALSE)
    return(invisible(NULL))
  }
  if(!is.numeric(x) || any(!is.finite(x)) || any(if(strict_lower) x<=lower else x<lower) ||
     any(x>upper) || (integer && any(x!=floor(x))))
    stop("Invalid numeric setting ",key,": expected ",if(integer) "integer " else "numeric ",
      "values ",if(strict_lower) "> " else ">= ",lower," and <= ",upper,".",call.=FALSE)
  invisible(x)
}
rna_config_validate <- function(cfg) {
  required <- c("inputs","comparison","harmonization","gain_loss","remodeling","dtu",
                "multiplicity","execution","reporting","sensitivity")
  absent <- setdiff(required,names(cfg))
  if(length(absent)) stop("Configuration lacks sections: ",paste(absent,collapse=", "),call.=FALSE)
  reference <- cfg$comparison$reference; comparison <- cfg$comparison$comparison
  if(length(reference)!=1L || length(comparison)!=1L || is.na(reference) || is.na(comparison) ||
     !nzchar(reference) || !nzchar(comparison) || identical(reference,comparison))
    stop("comparison.reference and comparison.comparison must be different nonempty condition labels.",call.=FALSE)
  if(length(cfg$comparison$paired)!=1L || !is.logical(cfg$comparison$paired) || is.na(cfg$comparison$paired))
    stop("comparison.paired must be true or false.",call.=FALSE)
  if (!isTRUE(cfg$comparison$paired)) {
    for (key in c("min_samples_per_condition", "n_bootstrap", "n_permutations"))
      rna_config_number(cfg,paste0("unpaired.",key),if(key=="min_samples_per_condition") 3 else 2,integer=TRUE)
    rna_config_number(cfg,"unpaired.minimum_prevalence_difference",0,1,strict_lower=TRUE)
    rna_config_number(cfg,"unpaired.minimum_evaluable_per_condition",3,integer=TRUE)
    if(cfg$unpaired$minimum_evaluable_per_condition < cfg$unpaired$min_samples_per_condition)
      stop("unpaired.minimum_evaluable_per_condition must be at least min_samples_per_condition.",call.=FALSE)
  }
  for(key in c("incomplete_pairs","technical_replicates"))
    if(!identical(cfg$comparison[[key]],"error")) stop("comparison.",key," must be error; automatic dropping or pooling is not implemented.",call.=FALSE)
  if(!identical(cfg$harmonization$identity_version,"junction_chain_paper_v1") ||
     !identical(cfg$harmonization$canonical_gene_policy,"paper_support_rank_preserve_all_counts"))
    stop("Use harmonization identity_version junction_chain_paper_v1 and canonical_gene_policy paper_support_rank_preserve_all_counts.",call.=FALSE)
  if(!cfg$sensitivity$mode %in% c("selected","all","none")) stop("sensitivity.mode must be selected, all, or none.",call.=FALSE)
  if(!cfg$execution$profile %in% c("production","pilot","smoke")) stop("execution.profile must be production, pilot, or smoke.",call.=FALSE)
  if(!identical(cfg$multiplicity$method,"BH")) stop("Only the paper's BH FDR adjustment is implemented.",call.=FALSE)
  policies <- c("harmonization.no_junction_policy"="sample_specific_not_cross_sample_comparable",
    "import.counts_unit"="integer_full_length_counts","import.duplicate_local_keys"="error",
    "gain_loss.rarefaction_target_policy"="minimum_primary_cohort_library",
    "sensitivity.common_depth_reference"="primary_cohort",
    "sensitivity.remodeling_isoform_retention_reference"="primary_cohort",
    "dtu.count_source"="full_depth_integer_counts")
  for(key in names(policies)) if(!identical(rna_config_get(cfg,key),unname(policies[[key]])))
    stop(key," must be ",policies[[key]]," for this implementation.",call.=FALSE)
  if(!isTRUE(cfg$gain_loss$sensitivity_reuses_primary_depth) ||
     !isTRUE(cfg$sensitivity$dtu_filter_and_refit_each_cohort))
    stop("Sensitivity analyses must reuse the primary common depth and refilter/refit DTU in each cohort.",call.=FALSE)
  for(key in c("execution.workers","execution.checkpoint_every_genes","gain_loss.n_downsample_replicates",
    "gain_loss.n_depth_curve_replicates","gain_loss.checkpoint_every","gain_loss.isoform_presence_min_fl",
    "gain_loss.gene_evaluable_min_fl","gain_loss.primary_patient_support","gain_loss.gene_floor_sensitivity",
    "gain_loss.patient_support_thresholds","remodeling.definition.min_gene_sample_count",
    "remodeling.definition.min_isoform_total_fl","remodeling.definition.min_evaluable_pairs",
    "remodeling.definition.dominance_support_n","remodeling.definition.patient_support_thresholds",
    "remodeling.definition.primary_patient_support","remodeling.bootstrap.n_bootstrap_replicates",
    "remodeling.split_null.n_split_null_replicates","remodeling.split_null.checkpoint_every_null_replicates",
    "remodeling.geometry.n_geometry_per_stratum","remodeling.geometry.geometry_isoform_numbers",
    "remodeling.geometry.geometry_depths","remodeling.cohort_test.gene_null_replicates",
    "remodeling.cohort_test.gene_null_max_isoforms","remodeling.cohort_test.covariance_initial_replicates",
    "remodeling.cohort_test.covariance_maximum_replicates","remodeling.cohort_test.profile_grid_points",
    "remodeling.patient_test.patient_test_min_condition_fl","remodeling.patient_test.patient_test_min_pooled_component_fl",
    "remodeling.patient_test.patient_test_max_isoforms","remodeling.validation.validation_replicates_per_scenario",
    "remodeling.validation.validation_null_replicates","remodeling.refinement.refinement_replicates",
    "dtu.minimum_samples_isoform_fraction","dtu.minimum_gene_count","dtu.minimum_samples_gene_count",
    "dtu.minimum_isoform_count","dtu.minimum_samples_isoform_count","dtu.minimum_isoforms_per_gene"))
    rna_config_number(cfg,key,1,integer=TRUE)
  for(key in c("gain_loss.opposite_direction_maximum","dtu.maximum_gene_plots","dtu.pilot_max_genes",
    "remodeling.cohort_test.gene_null_max_genes","remodeling.refinement.refinement_max_genes"))
    rna_config_number(cfg,key,0,integer=TRUE)
  for(key in c("multiplicity.fdr_thresholds","gain_loss.isoform_presence_min_fraction",
    "gain_loss.depth_robust_probability","gain_loss.depth_robust_probability_sensitivity","gain_loss.depth_curve_fractions",
    "remodeling.definition.jsd_threshold","remodeling.definition.max_delta_if_threshold",
    "remodeling.definition.entropy_direction_threshold","remodeling.definition.direction_delta_if",
    "remodeling.definition.min_retained_mass","remodeling.definition.dominance_mean_if",
    "remodeling.definition.dominance_margin","remodeling.definition.dominance_sample_min_if",
    "remodeling.definition.secondary_min_mass","remodeling.bootstrap.stability_probability",
    "remodeling.bootstrap.bootstrap_confidence","remodeling.bootstrap.bootstrap_mcse_target",
    "remodeling.bootstrap.minimum_valid_fraction","remodeling.split_null.empirical_null_overall_target",
    "remodeling.split_null.empirical_null_stratum_target","remodeling.split_null.empirical_null_convergence_tolerance",
    "remodeling.cohort_test.gene_model_minimum_success_fraction",
    "remodeling.cohort_test.gene_model_minimum_stratum_success_fraction",
    "remodeling.patient_test.patient_test_max_fraction_expected_below_5",
    "remodeling.validation.validation_minimum_production_availability",
    "remodeling.validation.validation_minimum_stratum_availability","remodeling.validation.validation_maximum_null_call_rate",
    "remodeling.refinement.refinement_trigger_p","dtu.minimum_isoform_fraction","dtu.patient_support_min_abs_delta_if"))
    rna_config_number(cfg,key,0,1,strict_lower=TRUE)
  for(key in c("import","gain_loss","qc","remodeling","dtu"))
    if(length(cfg$stages[[key]])!=1L || !is.logical(cfg$stages[[key]]) || is.na(cfg$stages[[key]]))
      stop("stages.",key," must be true or false.",call.=FALSE)
  if(!isTRUE(cfg$stages$import)) stop("stages.import must remain true: downstream stages need the raw imported inputs.",call.=FALSE)
  if(cfg$gain_loss$n_downsample_replicates < 2L || cfg$gain_loss$n_depth_curve_replicates > cfg$gain_loss$n_downsample_replicates)
    stop("Gain/loss needs at least two rarefactions and depth-curve replicates cannot exceed rarefaction replicates.",call.=FALSE)
  if(cfg$gain_loss$gene_evaluable_min_fl < cfg$gain_loss$isoform_presence_min_fl)
    stop("gain_loss.gene_evaluable_min_fl cannot be smaller than isoform_presence_min_fl.",call.=FALSE)
  if(cfg$remodeling$cohort_test$covariance_initial_replicates > cfg$remodeling$cohort_test$covariance_maximum_replicates)
    stop("Initial covariance replicates exceed the configured maximum.",call.=FALSE)
  if(cfg$dtu$minimum_isoforms_per_gene < 2L) stop("DTU requires at least two isoforms per gene.",call.=FALSE)
  dtu_rng_spec(cfg$dtu$rng)
  formats <- unlist(cfg$reporting$formats,use.names=FALSE)
  if(!length(formats) || any(!formats %in% c("pdf","svg","png","jpeg","jpg","tiff")))
    stop("Unsupported or empty reporting.formats.",call.=FALSE)
  rna_config_number(cfg,"reporting.dpi",72,integer=TRUE)
  rna_config_number(cfg,"reporting.font_size_pt",1)
  seeds <- unlist(cfg$execution$seeds,use.names=FALSE)
  if(!is.numeric(seeds) || !length(seeds) || any(!is.finite(seeds)) || any(seeds<1 | seeds>2000000000 | seeds!=floor(seeds)))
    stop("execution.seeds must contain positive integers <= 2000000000.",call.=FALSE)
  invisible(cfg)
}

read_rna_config <- function(config_file,validate_inputs=TRUE) {
  rna_require("yaml")
  if(length(config_file)!=1L || !file.exists(config_file)) stop("Configuration file does not exist: ",config_file,call.=FALSE)
  file <- normalizePath(config_file,winslash="/",mustWork=TRUE)
  cfg <- yaml::read_yaml(file,eval.expr=FALSE)
  if(!is.list(cfg)) stop("The YAML root must be a mapping.",call.=FALSE)
  cfg$.config_file <- file;cfg$.config_dir <- dirname(file)
  candidate_root <- dirname(cfg$.config_dir)
  cfg$.code_root <- getOption("IsoRemodel.code_root") %||%
    if(file.exists(file.path(candidate_root,"DESCRIPTION"))) candidate_root else getwd()
  cfg$.code_root <- normalizePath(cfg$.code_root,winslash="/",mustWork=TRUE)
  cfg$stages <- utils::modifyList(setNames(as.list(rep(TRUE,5L)),c("import","gain_loss","qc","remodeling","dtu")),cfg$stages %||% list())
  cfg$reporting <- utils::modifyList(list(enabled=TRUE,formats=c("pdf","svg","png","jpeg"),dpi=600,
    font_family="DejaVu Sans",font_size_pt=14,title_size_pt=16,axis_title_size_pt=16,
    axis_text_size_pt=14,legend_text_size_pt=14,subtitle_size_pt=14,caption_size_pt=14,
    scientific_labels_only=TRUE,separate_panels_for_composites=TRUE,replot_from_saved_results=TRUE),cfg$reporting %||% list())
  cfg$execution$force_refit <- cfg$execution$force_refit %||% FALSE
  cfg$dtu$rng <- dtu_rng_spec(cfg$dtu$rng)
  if (is.null(cfg$comparison$paired)) cfg$comparison$paired <- TRUE
  cfg$unpaired <- utils::modifyList(list(min_samples_per_condition=3L,
    n_bootstrap=1000L,n_permutations=1000L,minimum_prevalence_difference=.5),cfg$unpaired %||% list())
  cfg$unpaired$minimum_evaluable_per_condition <- cfg$unpaired$minimum_evaluable_per_condition %||% cfg$unpaired$min_samples_per_condition
  cfg$inputs$root_directory <- rna_config_path(cfg$inputs$root_directory %||% ".",cfg$.config_dir)
  cfg$.metadata_path <- rna_config_path(cfg$inputs$sample_metadata,cfg$.config_dir)
  if(is.null(cfg$.metadata_path)) stop("inputs.sample_metadata is required.",call.=FALSE)
  cfg$inputs$sample_metadata <- cfg$.metadata_path
  cfg$output$directory <- rna_config_path(cfg$output$directory %||% "../outputs/RNA_public_workflow",cfg$.config_dir)
  for(key in c("annotation_gtf","gene_id_map","pathway_term2gene","pathway_term_names"))
    cfg$references[[key]] <- rna_config_path(cfg$references[[key]],cfg$.config_dir)
  cfg$enrichment$reference_cache_directory <- rna_config_path(cfg$enrichment$reference_cache_directory %||%
    file.path(cfg$output$directory,"references","pathways"),cfg$.config_dir)
  rna_config_validate(cfg)
  if(isTRUE(validate_inputs)) {
    samples <- read_rna_metadata(cfg,validate_inputs=TRUE)
    rna_resolve_cohorts(samples,cfg)
    for(key in c("annotation_gtf","gene_id_map","pathway_term2gene","pathway_term_names")) {
      p <- cfg$references[[key]]
      if(!is.null(p) && !file.exists(p)) stop("Configured reference file does not exist (references.",key,"): ",p,call.=FALSE)
    }
    if(isTRUE(cfg$stages$dtu) && isTRUE(cfg$enrichment$enabled) && isTRUE(cfg$enrichment$matched_background) &&
       is.null(cfg$references$annotation_gtf)) stop("Matched DTU enrichment requires references.annotation_gtf.",call.=FALSE)
  }
  cfg
}

read_rna_metadata <- function(cfg,validate_inputs=TRUE) {
  path <- cfg$.metadata_path %||% cfg$inputs$sample_metadata
  if(is.null(path) || !file.exists(path)) stop("Sample metadata CSV does not exist: ",path,call.=FALSE)
  x <- utils::read.csv(path,stringsAsFactors=FALSE,check.names=FALSE,colClasses="character",na.strings=c("","NA"))
  if(anyDuplicated(names(x))) stop("Sample metadata has duplicated column names.",call.=FALSE)
  required <- c("sample_id","condition","classification_file","junctions_file")
  if (isTRUE(cfg$comparison$paired)) required <- c(required,"subject_id")
  if (!isTRUE(cfg$comparison$paired) && "sample_id" %in% names(x)) {
    if(!"subject_id" %in% names(x)) x$subject_id <- x$sample_id
    blank <- is.na(x$subject_id) | !nzchar(trimws(x$subject_id))
    x$subject_id[blank] <- x$sample_id[blank]
  }
  missing <- setdiff(required,names(x))
  if(length(missing)) stop("Sample metadata lacks columns: ",paste(missing,collapse=", "),call.=FALSE)
  for(key in names(x)) x[[key]] <- trimws(x[[key]])
  for(key in c("fl_count_column","pigeon_summary_file","filter_reasons_file","barcode"))
    if(!key %in% names(x)) x[[key]] <- NA_character_
  if(!"include" %in% names(x)) x$include <- "TRUE"
  include <- tolower(x$include)
  if(any(is.na(include) | !include %in% c("true","false","1","0","yes","no")))
    stop("Metadata include must be TRUE/FALSE, 1/0, or yes/no for every row.",call.=FALSE)
  x$include <- include %in% c("true","1","yes")
  x <- x[x$include,,drop=FALSE]
  if(!nrow(x)) stop("No samples have include=TRUE.",call.=FALSE)
  for(key in required) if(anyNA(x[[key]]) || any(!nzchar(x[[key]])))
    stop("Included metadata rows have missing ",key," values.",call.=FALSE)
  if(anyDuplicated(x$sample_id)) stop("Included sample_id values must be unique.",call.=FALSE)
  conditions <- c(cfg$comparison$reference,cfg$comparison$comparison)
  if(any(!x$condition %in% conditions)) stop("Metadata conditions must match the two configured comparison labels.",call.=FALSE)
  if (isTRUE(cfg$comparison$paired)) {
    pairs <- table(x$subject_id,factor(x$condition,levels=conditions))
    if(any(pairs!=1L)) stop("Each included subject must have exactly one sample in each configured condition; incomplete pairs and technical replicates are not silently dropped.",call.=FALSE)
    if(nrow(pairs)<2L) stop("At least two complete subject pairs are required.",call.=FALSE)
  } else {
    if(anyNA(x$subject_id) || any(!nzchar(x$subject_id)) || anyDuplicated(x$subject_id))
      stop("Unpaired mode requires independent samples with unique subject_id values. Repeated subjects require paired mode; do not treat technical replicates as independent samples.",call.=FALSE)
    totals <- table(factor(x$condition,levels=conditions))
    if(any(totals < cfg$unpaired$min_samples_per_condition))
      stop("Unpaired mode requires at least ",cfg$unpaired$min_samples_per_condition," independent samples per condition.",call.=FALSE)
  }
  x$sample <- x$sample_id;x$patient <- x$subject_id
  x$timepoint <- ifelse(x$condition==cfg$comparison$reference,"primary","recurrent")
  missing_count <- is.na(x$fl_count_column) | !nzchar(x$fl_count_column)
  x$fl_count_column[missing_count] <- cfg$import$count_column_default %||% "FL.sample"
  for(key in c("classification_file","junctions_file","pigeon_summary_file","filter_reasons_file")) {
    x[[key]] <- vapply(x[[key]],function(p) {
      if(is.na(p) || !nzchar(p)) return(NA_character_)
      rna_config_path(p,cfg$inputs$root_directory %||% dirname(path))
    },character(1))
    if(isTRUE(validate_inputs)) {
      present <- !is.na(x[[key]]) & nzchar(x[[key]])
      if(any(present & !file.exists(x[[key]]))) stop("Metadata ",key," file does not exist: ",
        paste(x[[key]][present & !file.exists(x[[key]])],collapse=", "),call.=FALSE)
    }
  }
  for(key in c("classification_file","junctions_file")) if(anyDuplicated(x[[key]]))
    stop("Included samples must not reuse the same ",key," path.",call.=FALSE)
  rownames(x) <- NULL
  x
}

rna_resolve_cohorts <- function(samples,cfg) {
  if (!isTRUE(cfg$comparison$paired)) return(rna_resolve_unpaired_cohorts(samples,cfg))
  patients <- unique(as.character(samples$patient %||% samples$subject_id))
  if(length(patients)<2L || anyNA(patients) || any(!nzchar(patients))) stop("Invalid subject IDs for cohort resolution.",call.=FALSE)
  specs <- cfg$cohorts %||% list(list(id=paste0("all_",length(patients),"_pairs"),exclude_subjects=character()))
  if(!length(specs) || !is.list(specs[[1L]])) stop("cohorts must be a YAML sequence of cohort definitions.",call.=FALSE)
  cohorts <- list()
  add <- function(id,exclude) {
    if(length(id)!=1L || is.na(id) || !nzchar(id) || !grepl("^[A-Za-z0-9_.-]+$",id)) stop("Invalid cohort ID: ",id,call.=FALSE)
    exclude <- as.character(unlist(exclude,use.names=FALSE))
    if(anyDuplicated(exclude) || any(!exclude %in% patients)) stop("Cohort ",id," excludes duplicated or unknown subjects.",call.=FALSE)
    members <- patients[!patients %in% exclude]
    if(id %in% names(cohorts)) {
      if(!identical(cohorts[[id]],members)) stop("Conflicting definitions for cohort ",id,call.=FALSE)
      return(invisible(NULL))
    }
    cohorts[[id]] <<- members
  }
  for(spec in specs) {
    if(length(unlist(spec$exclude_samples)))
      stop("Paired cohorts use exclude_subjects, not exclude_samples.",call.=FALSE)
    add(spec$id,spec$exclude_subjects)
  }
  if(!setequal(cohorts[[1L]],patients)) stop("The first cohort must include all imported subjects.",call.=FALSE)
  mode <- cfg$sensitivity$mode %||% "selected"
  omit <- if(mode=="all") patients else if(mode=="none") character() else as.character(unlist(cfg$sensitivity$patients,use.names=FALSE))
  if(mode=="selected" && !length(omit)) stop("sensitivity.patients must name at least one subject in selected mode.",call.=FALSE)
  if(anyDuplicated(omit) || any(!omit %in% patients)) stop("Sensitivity patients include duplicated or unknown subject IDs: ",paste(omit,collapse=", "),call.=FALSE)
  pattern <- cfg$sensitivity$cohort_id_pattern %||% "excluding_{patient}_{remaining_pairs}_pairs"
  for(p in omit) {
    id <- gsub("{patient}",p,pattern,fixed=TRUE)
    id <- gsub("{remaining_pairs}",as.character(length(patients)-1L),id,fixed=TRUE)
    add(id,p)
  }
  for(id in names(cohorts)) {
    n <- length(cohorts[[id]])
    if(n<2L) stop("Cohort ",id," has fewer than two complete pairs.",call.=FALSE)
    requirements <- c()
    if(isTRUE(cfg$stages$gain_loss)) requirements <- c(requirements,gain_loss=cfg$gain_loss$primary_patient_support)
    if(isTRUE(cfg$stages$remodeling)) requirements <- c(requirements,
      remodeling_evaluable=cfg$remodeling$definition$min_evaluable_pairs,
      remodeling_support=cfg$remodeling$definition$primary_patient_support,
      dominance_support=cfg$remodeling$definition$dominance_support_n)
    if(length(requirements) && any(requirements>n)) stop("Cohort ",id," has ",n,
      " pairs, below the configured absolute requirement(s): ",
      paste(names(requirements)[requirements>n],requirements[requirements>n],sep="=",collapse=", "),
      ". Change the cohort or explicitly change the scientific support setting.",call.=FALSE)
    if(isTRUE(cfg$stages$dtu)) for(key in c("minimum_samples_isoform_fraction","minimum_samples_gene_count","minimum_samples_isoform_count"))
      if(cfg$dtu[[key]]>2L*n) stop("Cohort ",id," has fewer samples than dtu.",key," requires.",call.=FALSE)
  }
  cohorts
}


# Unpaired cohorts contain sample IDs; paired cohorts above retain subject IDs.
rna_resolve_unpaired_cohorts <- function(samples,cfg) {
  ids <- as.character(samples$sample %||% samples$sample_id)
  conditions <- c(cfg$comparison$reference,cfg$comparison$comparison)
  if(anyNA(ids)||anyDuplicated(ids)||any(!nzchar(ids))) stop("Invalid independent sample IDs.",call.=FALSE)
  specs <- cfg$cohorts %||% list(list(id="all_samples",exclude_samples=character()))
  cohorts <- list()
  add <- function(id,exclude) {
    if(length(id)!=1 || is.na(id)||!grepl("^[A-Za-z0-9_.-]+$",id))stop("Invalid cohort ID.",call.=FALSE)
    exclude <- as.character(unlist(exclude,use.names=FALSE))
    if(anyDuplicated(exclude)||any(!exclude %in% ids))stop("Cohort ",id," excludes duplicate or unknown sample IDs.",call.=FALSE)
    members <- ids[!ids %in% exclude]
    if(id %in% names(cohorts)) {
      if(!identical(cohorts[[id]],members))stop("Conflicting cohort definitions: ",id,call.=FALSE)
      return(invisible(NULL))
    }
    cohorts[[id]] <<- members
  }
  for(spec in specs) {
    if(length(unlist(spec$exclude_subjects)))stop("Unpaired cohorts use exclude_samples, not exclude_subjects.",call.=FALSE)
    add(spec$id,spec$exclude_samples)
  }
  if(!length(cohorts)||!setequal(cohorts[[1]],ids))stop("The first unpaired cohort must contain all imported samples.",call.=FALSE)
  mode <- cfg$sensitivity$mode %||% "none"
  omit <- if(mode=="all") ids else if(mode=="selected") as.character(unlist(cfg$sensitivity$samples,use.names=FALSE)) else character()
  if(mode=="selected"&&!length(omit))stop("sensitivity.samples must name at least one sample in unpaired selected mode.",call.=FALSE)
  if(anyDuplicated(omit)||any(!omit %in% ids))stop("Sensitivity contains duplicate or unknown sample IDs.",call.=FALSE)
  pattern <- cfg$sensitivity$sample_cohort_id_pattern %||% "excluding_{sample}_{remaining_samples}_samples"
  safe <- make.unique(gsub("[^A-Za-z0-9_.-]","_",ids),sep="_")
  for(id in omit) {
    label <- gsub("{sample}",safe[match(id,ids)],pattern,fixed=TRUE)
    label <- gsub("{remaining_samples}",as.character(length(ids)-1),label,fixed=TRUE)
    add(label,id)
  }
  minimum <- as.integer(cfg$unpaired$min_samples_per_condition %||% 3L)
  for(id in names(cohorts)) {
    take <- samples[match(cohorts[[id]],ids),,drop=FALSE]
    totals <- table(factor(take$condition,levels=conditions))
    if(any(totals<minimum))stop("Cohort ",id," leaves fewer than ",minimum," independent samples in a condition; change exclusions or the design. Minimum sizes are not reduced automatically.",call.=FALSE)
    if(isTRUE(cfg$stages$dtu))for(key in c("minimum_samples_isoform_fraction","minimum_samples_gene_count","minimum_samples_isoform_count"))
      if(cfg$dtu[[key]]>nrow(take))stop("Cohort ",id," has fewer samples than dtu.",key," requires.",call.=FALSE)
  }
  cohorts
}
