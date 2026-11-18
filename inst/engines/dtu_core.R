# Source-derived DTU engine. Extracted from CancerCommunications_revision_Module06
# workflow/model contract 1.3.0. Paired design, complete correction families and
# numerical adapter are retained. Private DRIMSeq helper adaptations target
# DRIMSeq 1.38.0 (Malgorzata Nowicka; GPL >=3); see NOTICE.
# Portability changes: explicit metadata cohorts, serial numerical audits,
# stable cohort seeds, empty-family returns and alpha-specific stageR tags.
# This file is evaluated in a private environment by run_dtu().

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  x
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

finite_mean <- function(x) {
  x <- safe_numeric(x)
  if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}

finite_median <- function(x) {
  x <- safe_numeric(x)
  if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
}

finite_max <- function(x) {
  x <- safe_numeric(x)
  if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
}

safe_name <- function(x) {
  output <- str_replace_all(as.character(x), "[^A-Za-z0-9._-]+", "_")
  str_replace_all(output, "^_+|_+$", "")
}

validate_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " lacks required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

pick_column <- function(data, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) {
    if (required) {
      stop("Could not identify ", label, ". Tried: ", paste(candidates, collapse = ", "))
    }
    return(NULL)
  }
  hit[[1]]
}

strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", clean_character(x))

display_structural_category <- function(x) {
  value <- str_to_lower(coalesce(clean_character(x), "unclassified"))
  case_when(
    value %in% c("fsm", "full-splice_match", "full_splice_match", "full splice match") ~ "FSM",
    value %in% c("ism", "incomplete-splice_match", "incomplete_splice_match", "incomplete splice match") ~ "ISM",
    value %in% c("nic", "novel_in_catalog", "novel-in-catalog", "novel in catalog") ~ "NIC",
    value %in% c("nnc", "novel_not_in_catalog", "novel-not-in-catalog", "novel not in catalog") ~ "NNC",
    str_detect(value, "antisense") ~ "Antisense",
    str_detect(value, "fusion") ~ "Fusion",
    str_detect(value, "intergenic") ~ "Intergenic",
    str_detect(value, "genic") ~ "Genic",
    TRUE ~ str_to_sentence(str_replace_all(value, "[_-]+", " "))
  )
}

# The independent design has no subject nuisance coefficient. DRIMSeq removes
# the tested timepoint coefficient to construct its intercept-only null.
m06_design_matrix <- function(samples, paired = TRUE) {
  if (isTRUE(paired)) model.matrix(~ patient + timepoint, data = samples)
  else model.matrix(~ timepoint, data = samples)
}

prepare_dtu_family <- function(cohort_name) {
  paired_design <- !identical(cfg$paired, FALSE)
  cohort_samples <- sample_manifest %>%
    filter(if (paired_design) patient %in% analysis_cohorts[[cohort_name]] else sample %in% analysis_cohorts[[cohort_name]]) %>%
    mutate(
      sample_id = sample,
      patient = factor(patient),
      timepoint = factor(timepoint, levels = c("primary", "recurrent"))
    ) %>%
    select(sample_id, patient, timepoint)

  detected_in_cohort <- rowSums(
    counts_matrix[, cohort_samples$sample_id, drop = FALSE], na.rm = TRUE
  ) > 0
  detected_isoform_complexity <- isoform_annotation %>%
    transmute(
      gene_id = clean_character(associated_gene),
      gene_mapped = coalesce(as.logical(gene_mapped), FALSE),
      cross_sample_comparable = coalesce(as.logical(cross_sample_comparable), FALSE),
      detected_in_cohort = detected_in_cohort
    ) %>%
    filter(
      gene_mapped, cross_sample_comparable, detected_in_cohort,
      !is.na(gene_id)
    ) %>%
    count(gene_id, name = "detected_isoforms")

  included_isoforms <- isoform_annotation %>%
    mutate(
      isoform_harmonized_id = as.character(isoform_harmonized_id),
      gene_id = clean_character(associated_gene),
      structural_category = coalesce(clean_character(structural_category), "Unclassified"),
      analyzed_for_gain_loss = coalesce(as.logical(analyzed_for_gain_loss), FALSE)
    ) %>%
    filter(analyzed_for_gain_loss, !is.na(gene_id)) %>%
    distinct(isoform_harmonized_id, .keep_all = TRUE)

  count_index <- match(included_isoforms$isoform_harmonized_id, rownames(counts_matrix))
  if (anyNA(count_index)) stop("Internal count/annotation mismatch while preparing ", cohort_name)
  count_mat <- counts_matrix[count_index, cohort_samples$sample_id, drop = FALSE]
  storage.mode(count_mat) <- "numeric"
  nonzero <- rowSums(count_mat) > 0
  count_mat <- count_mat[nonzero, , drop = FALSE]
  included_isoforms <- included_isoforms[nonzero, , drop = FALSE]
  if (!nrow(count_mat)) return(empty_dtu_prepared(cohort_name, cohort_samples,
    tibble(cohort = cohort_name, stage = "Comparable isoform universe", genes = 0L, isoforms = 0L)))

  # Expression matching must not be conditioned on which isoforms later survive
  # the DTU proportion/count filters. Pool the complete comparable Module 02
  # gain/loss universe first, then carry that gene-level depth into the frozen
  # testing family as a separate covariate.
  pooled_library_depth <- tibble(
    sample_id = cohort_samples$sample_id,
    pooled_comparable_library_FL = as.numeric(colSums(count_mat))
  )
  if (any(pooled_library_depth$pooled_comparable_library_FL <= 0)) {
    stop("A comparable-isoform library has zero depth in ", cohort_name, ".")
  }

  pooled_gene_expression <- as_tibble(
    rowsum(count_mat, included_isoforms$gene_id, reorder = FALSE),
    rownames = "gene_id"
  ) %>%
    pivot_longer(
      all_of(cohort_samples$sample_id),
      names_to = "sample_id", values_to = "pooled_gene_count"
    ) %>%
    left_join(pooled_library_depth, by = "sample_id") %>%
    mutate(
      pooled_gene_CPM =
        1e6 * pooled_gene_count / pooled_comparable_library_FL
    ) %>%
    group_by(gene_id) %>%
    summarise(
      pooled_mean_gene_count = mean(pooled_gene_count),
      pooled_median_gene_count = median(pooled_gene_count),
      pooled_minimum_gene_count = min(pooled_gene_count),
      pooled_mean_gene_CPM = mean(pooled_gene_CPM),
      pooled_median_gene_CPM = median(pooled_gene_CPM),
      .groups = "drop"
    )

  flow <- tibble(
    cohort = cohort_name,
    stage = "Comparable isoform universe",
    genes = n_distinct(included_isoforms$gene_id),
    isoforms = nrow(included_isoforms)
  )

  at_least_two <- included_isoforms %>%
    count(gene_id, name = "n_isoforms") %>%
    filter(n_isoforms >= cfg$filters$minimum_isoforms_per_gene) %>%
    pull(gene_id)
  keep <- included_isoforms$gene_id %in% at_least_two
  count_mat <- count_mat[keep, , drop = FALSE]
  included_isoforms <- included_isoforms[keep, , drop = FALSE]
  flow <- bind_rows(flow, tibble(
    cohort = cohort_name,
    stage = paste0(">=", cfg$filters$minimum_isoforms_per_gene, " observed isoforms"),
    genes = n_distinct(included_isoforms$gene_id),
    isoforms = nrow(included_isoforms)
  ))
  if (!nrow(count_mat)) return(empty_dtu_prepared(cohort_name, cohort_samples, flow))

  gene_total_mat <- rowsum(count_mat, included_isoforms$gene_id, reorder = FALSE)
  gene_row_index <- match(included_isoforms$gene_id, rownames(gene_total_mat))
  gene_total_per_isoform <- gene_total_mat[gene_row_index, , drop = FALSE]
  within_gene_fraction <- count_mat / pmax(gene_total_per_isoform, 1)
  proportion_keep <- rowSums(
    within_gene_fraction >= cfg$filters$minimum_isoform_fraction
  ) >= cfg$filters$minimum_samples_isoform_fraction
  count_mat <- count_mat[proportion_keep, , drop = FALSE]
  included_isoforms <- included_isoforms[proportion_keep, , drop = FALSE]

  after_prop_genes <- included_isoforms %>%
    count(gene_id, name = "n_isoforms") %>%
    filter(n_isoforms >= cfg$filters$minimum_isoforms_per_gene) %>%
    pull(gene_id)
  keep <- included_isoforms$gene_id %in% after_prop_genes
  count_mat <- count_mat[keep, , drop = FALSE]
  included_isoforms <- included_isoforms[keep, , drop = FALSE]
  flow <- bind_rows(flow, tibble(
    cohort = cohort_name,
    stage = paste0(">=", 100 * cfg$filters$minimum_isoform_fraction, "% in >=", cfg$filters$minimum_samples_isoform_fraction, " samples; >=", cfg$filters$minimum_isoforms_per_gene, " retained"),
    genes = n_distinct(included_isoforms$gene_id),
    isoforms = nrow(included_isoforms)
  ))

  counts_wide <- bind_cols(
    included_isoforms %>% transmute(
      gene_id,
      feature_id = paste(gene_id, isoform_harmonized_id, sep = "|"),
      isoform_harmonized_id,
      structural_category
    ),
    as_tibble(count_mat, .name_repair = "minimal")
  )
  names(counts_wide)[seq.int(5L, ncol(counts_wide))] <- cohort_samples$sample_id
  feature_map <- counts_wide %>%
    select(gene_id, feature_id, isoform_harmonized_id, structural_category)
  dm_counts <- counts_wide %>%
    select(gene_id, feature_id, all_of(cohort_samples$sample_id)) %>%
    as.data.frame(stringsAsFactors = FALSE)

  if (anyDuplicated(dm_counts$feature_id)) stop("Duplicate feature IDs in ", cohort_name)
  if (!nrow(dm_counts)) return(empty_dtu_prepared(cohort_name, cohort_samples, flow))
  dm <- DRIMSeq::dmDSdata(counts = dm_counts, samples = as.data.frame(cohort_samples))
  dm <- DRIMSeq::dmFilter(
    dm,
    min_samps_gene_expr = cfg$filters$minimum_samples_gene_count,
    min_gene_expr = cfg$filters$minimum_gene_count,
    min_samps_feature_expr = cfg$filters$minimum_samples_isoform_count,
    min_feature_expr = cfg$filters$minimum_isoform_count
  )

  filtered_counts <- as_tibble(DRIMSeq::counts(dm))
  retained_gene_ids <- filtered_counts %>% count(gene_id, name = "retained_isoforms") %>%
    filter(retained_isoforms >= cfg$filters$minimum_isoforms_per_gene) %>% pull(gene_id)
  filtered_counts <- filtered_counts %>% filter(gene_id %in% retained_gene_ids)
  if (!nrow(filtered_counts) || n_distinct(filtered_counts$gene_id) == 0L) {
    return(empty_dtu_prepared(cohort_name, cohort_samples, flow))
  }
  dm <- DRIMSeq::dmDSdata(counts = as.data.frame(filtered_counts), samples = as.data.frame(cohort_samples))
  feature_map <- feature_map %>% semi_join(filtered_counts, by = c("gene_id", "feature_id"))
  flow <- bind_rows(flow, tibble(
    cohort = cohort_name,
    stage = "DRIMSeq count filter",
    genes = n_distinct(filtered_counts$gene_id),
    isoforms = nrow(filtered_counts)
  ))

  gene_depth <- filtered_counts %>%
    pivot_longer(all_of(cohort_samples$sample_id), names_to = "sample_id", values_to = "count") %>%
    group_by(gene_id, sample_id) %>%
    summarise(gene_count = sum(count), .groups = "drop") %>%
    group_by(gene_id) %>%
    summarise(
      tested_family_mean_gene_count = mean(gene_count),
      tested_family_median_gene_count = median(gene_count),
      tested_family_minimum_gene_count = min(gene_count),
      .groups = "drop"
    )

  if (cfg$analysis_profile == "pilot" && n_distinct(filtered_counts$gene_id) > cfg$pilot_max_genes) {
    selected_genes <- gene_depth %>%
      arrange(desc(tested_family_mean_gene_count), gene_id) %>%
      slice_head(n = cfg$pilot_max_genes) %>%
      pull(gene_id)
    filtered_counts <- filtered_counts %>% filter(gene_id %in% selected_genes)
    feature_map <- feature_map %>% filter(gene_id %in% selected_genes)
    dm <- DRIMSeq::dmDSdata(
      counts = as.data.frame(filtered_counts),
      samples = as.data.frame(cohort_samples)
    )
    gene_depth <- gene_depth %>% filter(gene_id %in% selected_genes)
    flow <- bind_rows(flow, tibble(
      cohort = cohort_name,
      stage = "Pilot cap by mean gene count",
      genes = n_distinct(filtered_counts$gene_id),
      isoforms = nrow(filtered_counts)
    ))
  }

  tested_isoforms <- filtered_counts %>% count(gene_id, name = "tested_isoforms")
  testing_family <- gene_depth %>%
    inner_join(tested_isoforms, by = "gene_id") %>%
    inner_join(pooled_gene_expression, by = "gene_id") %>%
    left_join(detected_isoform_complexity, by = "gene_id") %>%
    mutate(
      gencode_lookup_id = case_when(
        str_detect(gene_id, "ENSG[0-9]+") ~
          strip_ensembl_version(str_extract(gene_id, "ENSG[0-9]+(?:\\.[0-9]+)?")),
        gene_id %in% gencode_complexity$gene_id ~ gene_id,
        str_remove(gene_id, "_.*$") %in% gencode_complexity$gene_id ~
          str_remove(gene_id, "_.*$"),
        TRUE ~ gene_id
      )
    ) %>%
    left_join(
      gencode_complexity %>% rename(gencode_lookup_id = gene_id),
      by = "gencode_lookup_id"
    ) %>%
    mutate(
      detected_isoforms = as.integer(detected_isoforms),
      gencode_annotated_transcripts = as.integer(gencode_annotated_transcripts),
      gencode_complexity_available = !is.na(gencode_annotated_transcripts)
    ) %>%
    select(-gencode_lookup_id) %>%
    arrange(gene_id)
  if (any(!is.finite(testing_family$detected_isoforms)) ||
      any(testing_family$detected_isoforms < testing_family$tested_isoforms)) {
    stop(
      "Detected-isoform complexity is missing or smaller than the tested ",
      "isoform family for ", cohort_name, "."
    )
  }

  design <- m06_design_matrix(DRIMSeq::samples(dm), paired = paired_design)
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) stop(if (paired_design) "Paired" else "Unpaired", " DTU design is not estimable for ", cohort_name)
  coef_candidates <- colnames(design)[grepl("^timepoint", colnames(design))]
  if (length(coef_candidates) != 1L) {
    stop("Expected one recurrent-timepoint coefficient for ", cohort_name, "; found: ",
         paste(coef_candidates, collapse = ", "))
  }

  list(
    cohort = cohort_name,
    dm = dm,
    counts = filtered_counts,
    samples = as_tibble(DRIMSeq::samples(dm)),
    design = design,
    coefficient = coef_candidates[[1]],
    feature_map = feature_map,
    testing_family = testing_family,
    flow = flow
  )
}

make_m06_drimseq_engine <- function() {
  # This local adapter targets the exact upstream implementation audited here.
  # It does not change the installed package or its namespace.
  if (as.character(utils::packageVersion("DRIMSeq")) != "1.38.0") {
    stop("DTU numerical adapter requires DRIMSeq 1.38.0; review the adapter before using another version.")
  }
  ns <- asNamespace('DRIMSeq')
  e <- new.env(parent = ns)
  helpers <- c('dmDS_estimateCommonPrecision','dmDS_profileLikCommon',
    'dmDS_estimateTagwisePrecision','dmDS_profileLik','dmDS_CRadjustment',
    'dmDS_CRadjustmentRegression_gene','dmDS_fit','dmDS_fitRegression_gene',
    'dm_fitRegression')
  for (nm in helpers) {
    f <- get(nm, envir = ns, inherits = FALSE)
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  # Correct the diagonal trigamma difference; the original plus sign does
  # not equal the derivative of DRIMSeq's score. Keep the remaining expression.
  h <- get("dm_Hessian_regG_prop", envir = ns, inherits = FALSE)
  h_text <- paste(deparse(body(h), width.cutoff = 500L), collapse = "\n")
  old_term <- "ltgp_y[, jp] + ltgp[, jp]"
  if (!grepl(old_term, h_text, fixed = TRUE)) {
    stop("DRIMSeq Hessian differs from the audited 1.38.0 implementation.")
  }
  body(h) <- parse(text = gsub(
    old_term, "ltgp_y[, jp] - ltgp[, jp]", h_text, fixed = TRUE
  ))[[1L]]
  environment(h) <- e
  e$dm_Hessian_regG_prop <- h
  # Same BFGS objective and tolerance, with additional iterations for convergence.
  e$optim <- function(...) {
    args <- list(...)
    args$control$maxit <- 1000L
    do.call(stats::optim, args)
  }
  e$dm_CRadjustmentRegression <- function(y, x, prec, prop) {
    if (any(!is.finite(prop))) return(NA_real_)
    prop <- t(prop)
    n <- nrow(prop); q <- ncol(prop)
    log_adjustment <- function(H) {
      if (any(!is.finite(H))) return(NA_real_)
      de <- determinant(n * (-H), logarithm = TRUE)
      if (de$sign != 1 || !is.finite(as.numeric(de$modulus))) return(NA_real_)
      as.numeric(de$modulus) / 2
    }
    H <- dm_Hessian_regG_prop(y = t(y), prec = prec, prop = prop, x = x)
    adj <- log_adjustment(H)
    if (!is.finite(adj)) {
      logit_prop <- log(prop / prop[, q])
      par <- c(MASS::ginv(x) %*% logit_prop[, -q, drop = FALSE])
      if (any(!is.finite(par))) return(NA_real_)
      H <- stats::optimHess(par = par, fn = dm_lik_regG, gr = dm_score_regG,
                           x = x, prec = prec, y = y)
      adj <- log_adjustment(H)
    }
    adj
  }
  environment(e$dm_CRadjustmentRegression) <- e
  # One-way null fits use patient groups for the paired design and one
  # intercept-only group for the independent design. A group with only one
  # observed isoform has a valid boundary maximum (log likelihood zero).
  # Remove zero-total categories only within that group fit and restore their
  # exact zero probabilities; no counts or features leave the testing family.
  for (nm in c("dmDS_fitManyGroups_gene", "dm_fitManyGroups",
               "bbDS_fit", "bbDS_fitManyGroups_gene", "bb_fitManyGroups",
               "bbDS_fitRegression_gene", "bb_fitRegression")) {
    f <- get(nm, envir = ns, inherits = FALSE)
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  e$original_one_group <- get("dm_fitOneGroup", envir = ns)
  e$dm_fitOneGroup <- function(y, prec, prop_mode = "constrOptim", prop_tol = 1e-12) {
    q <- nrow(y)
    keep <- which(rowSums(y) > 0)
    if (!is.finite(prec) || q < 2L) {
      return(original_one_group(y, prec, prop_mode, prop_tol))
    }
    if (!length(keep)) return(list(prop = rep(1 / q, q), lik = 0))
    prop <- rep(0, q)
    if (length(keep) == 1L) {
      prop[keep] <- 1
      return(list(prop = prop, lik = 0))
    }
    fitted <- original_one_group(y[keep, , drop = FALSE], prec, prop_mode, prop_tol)
    prop[keep] <- fitted$prop
    list(prop = prop, lik = fitted$lik)
  }
  environment(e$dm_fitOneGroup) <- e

  # Algebraically identical DM likelihood and score, evaluated without
  # subtracting nearly equal probabilities or gamma/digamma infinities.
  # y=0 contributes exactly zero to log Gamma(y+a)-log Gamma(a).
  e$stable_dm_parts <- function(b, x, prec, y) {
    eta <- cbind(x %*% matrix(b, ncol(x), nrow(y) - 1L), 0)
    eta <- eta - apply(eta, 1L, max)
    log_prop <- eta - log(rowSums(exp(eta)))
    prop <- exp(log_prop)
    counts <- t(y)
    a <- prec * prop
    positive <- counts > 0
    terms <- matrix(0, nrow(counts), ncol(counts))
    score_mass <- terms
    terms[positive] <- lgamma(counts[positive] + a[positive]) -
      lgamma(1 + a[positive]) + log(prec) + log_prop[positive]
    score_mass[positive] <- 1 + a[positive] * (
      digamma(counts[positive] + a[positive]) - digamma(1 + a[positive])
    )
    list(
      lik = sum(lgamma(prec) - lgamma(prec + rowSums(counts))) + sum(terms),
      score = c(crossprod(x, (score_mass - prop * rowSums(score_mass))[
        , -ncol(counts), drop = FALSE
      ])),
      prop = prop
    )
  }
  e$stable_dm_lik <- function(b, x, prec, y) stable_dm_parts(b, x, prec, y)$lik
  e$stable_dm_score <- function(b, x, prec, y) stable_dm_parts(b, x, prec, y)$score
  for (nm in c("stable_dm_lik", "stable_dm_score")) environment(e[[nm]]) <- e
  e$original_regression <- e$dm_fitRegression
  e$dm_fitRegression <- function(y, design, prec, coef_mode = "optim", coef_tol = 1e-12) {
    original <- tryCatch(
      original_regression(y, design, prec, coef_mode, coef_tol),
      error = function(error) NULL
    )
    if (!is.finite(prec) || (!is.null(original) && is.finite(original$lik) &&
                            all(is.finite(original$fit)))) return(original)
    p <- ncol(design)
    q <- nrow(y)
    if (is.null(original)) {
      original <- list(b = matrix(NA_real_, q, p, dimnames = list(rownames(y), NULL)),
                       lik = NA_real_, fit = matrix(NA_real_, q, ncol(y), dimnames = dimnames(y)))
    }
    initial <- t(y) + 1 # Starting values only; the objective uses original counts.
    initial <- initial / rowSums(initial)
    start <- c(MASS::ginv(design) %*% log(initial[, -q, drop = FALSE] / initial[, q]))
    fitted <- tryCatch(stats::optim(
      start, fn = stable_dm_lik, gr = stable_dm_score,
      x = design, prec = prec, y = y, method = "BFGS",
      control = list(fnscale = -1, reltol = coef_tol, maxit = 10000L)
    ), error = function(error) NULL)
    if (is.null(fitted) || fitted$convergence != 0L || !is.finite(fitted$value)) {
      return(original)
    }
    evaluated <- stable_dm_parts(fitted$par, design, prec, y)
    coefficients <- rbind(t(matrix(fitted$par, p, q - 1L)), rep(0, p))
    proportions <- t(evaluated$prop)
    rownames(coefficients) <- rownames(y)
    rownames(proportions) <- rownames(y)
    list(b = coefficients, lik = evaluated$lik, fit = proportions)
  }
  environment(e$dm_fitRegression) <- e

  # stageR feature evidence uses the BB marginals of the same fitted DM
  # probabilities. Evaluate the p=0 and p=1 boundaries by their exact limits.
  e$stable_bb_likelihood <- function(y, prec, prop) {
    if (!is.finite(prec) || any(!is.finite(prop))) return(rep(NA_real_, nrow(y)))
    if (is.null(dim(prop))) prop <- matrix(prop, nrow(y), ncol(y))
    total <- colSums(y)
    remainder <- matrix(total, nrow(y), ncol(y), byrow = TRUE) - y
    log_rising <- function(count, a) {
      out <- matrix(0, nrow(count), ncol(count))
      positive <- count > 0
      out[positive] <- lgamma(count[positive] + a[positive]) -
        lgamma(1 + a[positive]) + log(a[positive])
      out
    }
    rowSums(log_rising(y, prec * prop) + log_rising(remainder, prec * (1 - prop))) +
      sum(lgamma(prec) - lgamma(prec + total))
  }
  e$bb_fitOneGroup <- function(y, prec, prop) {
    list(prop = prop, lik = stable_bb_likelihood(y, prec, prop))
  }
  e$bb_lik_regG_prop <- function(y, prec, prop) {
    stable_bb_likelihood(t(y), prec, t(prop))
  }
  for (nm in c("bb_fitOneGroup", "bb_lik_regG_prop")) environment(e[[nm]]) <- e

  # Audit every tagwise grid point without changing its likelihood or precision.
  e$record_tagwise <- FALSE
  e$profile_audit <- list()
  e$tagwise_implementation <- e$dmDS_estimateTagwisePrecision
  e$dmDS_estimateTagwisePrecision <- function(...) {
    record_tagwise <<- TRUE
    on.exit(record_tagwise <<- FALSE, add = TRUE)
    tagwise_implementation(...)
  }
  environment(e$dmDS_estimateTagwisePrecision) <- e
  profile <- e$dmDS_profileLik
  profile_body <- body(profile)
  audit_expression <- quote({
    if (record_tagwise) {
      profile_audit[[length(profile_audit) + 1L]] <<- data.frame(
        gene_id = names(counts), grid_precision = prec,
        grid_log_likelihood = as.numeric(fit$lik),
        grid_CR_adjustment = as.numeric(adj),
        grid_adjusted_log_likelihood = as.numeric(adj_lik),
        stringsAsFactors = FALSE
      )
    }
  })
  if (!identical(profile_body[[length(profile_body)]], quote(return(adj_lik)))) {
    stop("Unexpected DRIMSeq profile-likelihood return expression.")
  }
  body(profile) <- as.call(c(
    as.list(profile_body)[-length(profile_body)],
    list(audit_expression, quote(return(adj_lik)))
  ))
  e$dmDS_profileLik <- profile
  specs <- c(dmPrecision='dmDSdata', dmFit='dmDSprecision', dmTest='dmDSfit')
  out <- lapply(names(specs), function(nm) {
    method <- methods::selectMethod(get(nm, envir = ns), specs[[nm]])
    # S4 dispatch wrappers embed .local as an already-created closure. Changing
    # the outer wrapper's environment alone leaves .local in the namespace.
    method_body <- body(method)
    if (length(method_body) < 2L || !is.call(method_body[[2L]]) ||
        !identical(method_body[[2L]][[2L]], as.name(".local")) ||
        !is.function(method_body[[2L]][[3L]])) {
      stop("Unexpected DRIMSeq S4 method wrapper: ", nm)
    }
    f <- method_body[[2L]][[3L]]
    environment(f) <- e
    f
  })
  names(out) <- names(specs)
  out$profile_audit <- function() {
    if (!length(e$profile_audit)) return(data.frame())
    do.call(rbind, e$profile_audit)
  }
  out
}

m06_reference_order <- function(counts, samples, paired = TRUE) {
  if (!isTRUE(paired)) return(m06_unpaired_reference_order(counts, samples))
  sample_ids <- as.character(samples$sample_id)
  patients <- unique(as.character(samples$patient))
  mat <- as.matrix(counts[,sample_ids,drop=FALSE])
  pair_sums <- vapply(patients, function(pt) rowSums(mat[,as.character(samples$patient)==pt,drop=FALSE]), numeric(nrow(mat)))
  min_pair <- apply(pair_sums,1,min)
  gene_indices <- split(seq_len(nrow(counts)),factor(counts$gene_id,levels=unique(counts$gene_id)))
  audit <- lapply(gene_indices,function(idx){
    current <- tail(idx,1L)
    best <- idx[order(-min_pair[idx],-rowSums(mat[idx,,drop=FALSE]),counts$feature_id[idx])][1]
    selected <- if(min_pair[current]>0) current else best
    data.frame(gene_id=counts$gene_id[current],original_reference=counts$feature_id[current],
      selected_reference=counts$feature_id[selected],minimum_pair_reference_count=min_pair[selected],
      reference_changed=selected!=current, reference_supported_in_all_pairs=min_pair[selected]>0,
      reference_is_best_available=min_pair[selected]>0 || min_pair[selected]==max(min_pair[idx]),
      patient_pairs_without_gene_counts=sum(colSums(pair_sums[idx,,drop=FALSE])==0),
      samples_without_gene_counts=sum(colSums(mat[idx,,drop=FALSE])==0),
      design_rank_with_gene_counts=qr(model.matrix(~patient+timepoint,data=samples)[colSums(mat[idx,,drop=FALSE])>0,,drop=FALSE])$rank,
      isoforms_absent_in_primary=sum(rowSums(mat[idx,as.character(samples$timepoint)=="primary",drop=FALSE])==0),
      isoforms_absent_in_recurrent=sum(rowSums(mat[idx,as.character(samples$timepoint)=="recurrent",drop=FALSE])==0),
      isoforms_absent_in_at_least_one_pair=sum(rowSums(pair_sums[idx,,drop=FALSE]==0)>0),
      isoform_patient_pairs_with_zero_counts=sum(pair_sums[idx,,drop=FALSE]==0))
  })
  audit <- do.call(rbind,audit);rownames(audit)<-NULL
  reordered <- unlist(Map(function(idx,ref)c(idx[counts$feature_id[idx]!=ref],idx[counts$feature_id[idx]==ref]),gene_indices,audit$selected_reference),use.names=FALSE)
  list(counts=counts[reordered,,drop=FALSE],audit=audit)
}

m06_unpaired_reference_order <- function(counts, samples) {
  sample_ids <- as.character(samples$sample_id)
  conditions <- c("primary", "recurrent")
  mat <- as.matrix(counts[, sample_ids, drop = FALSE])
  pooled <- vapply(conditions, function(condition)
    rowSums(mat[, as.character(samples$timepoint) == condition, drop = FALSE]), numeric(nrow(mat)))
  pooled <- matrix(pooled, nrow = nrow(mat), ncol = length(conditions))
  minimum <- apply(pooled, 1L, min)
  indices <- split(seq_len(nrow(counts)), factor(counts$gene_id, levels = unique(counts$gene_id)))
  design <- m06_design_matrix(samples, paired = FALSE)
  audit <- lapply(indices, function(idx) {
    current <- tail(idx, 1L)
    best <- idx[order(-minimum[idx], -rowSums(mat[idx, , drop = FALSE]), counts$feature_id[idx])][1L]
    selected <- if (minimum[current] > 0) current else best
    data.frame(gene_id = counts$gene_id[current], original_reference = counts$feature_id[current],
      selected_reference = counts$feature_id[selected],
      minimum_pair_reference_count = NA_real_, minimum_condition_reference_count = minimum[selected],
      reference_changed = selected != current, reference_supported_in_all_pairs = NA,
      reference_supported_in_both_conditions = minimum[selected] > 0,
      reference_is_best_available = minimum[selected] > 0 || minimum[selected] == max(minimum[idx]),
      patient_pairs_without_gene_counts = NA_integer_,
      samples_without_gene_counts = sum(colSums(mat[idx, , drop = FALSE]) == 0),
      design_rank_with_gene_counts = qr(design[colSums(mat[idx, , drop = FALSE]) > 0, , drop = FALSE])$rank,
      isoforms_absent_in_primary = sum(pooled[idx, 1L] == 0),
      isoforms_absent_in_recurrent = sum(pooled[idx, 2L] == 0),
      isoforms_absent_in_at_least_one_pair = NA_integer_,
      isoform_patient_pairs_with_zero_counts = NA_integer_, analysis_design = "unpaired")
  })
  audit <- do.call(rbind, audit); rownames(audit) <- NULL
  reordered <- unlist(Map(function(idx, ref) c(idx[counts$feature_id[idx] != ref],
    idx[counts$feature_id[idx] == ref]), indices, audit$selected_reference), use.names = FALSE)
  list(counts = counts[reordered, , drop = FALSE], audit = audit)
}

fit_drimseq_model <- function(prepared) {
  previous_rng <- dtu_rng_snapshot()
  on.exit(dtu_rng_restore(previous_rng), add = TRUE)
  cohort_name <- prepared$cohort
  rng <- dtu_rng_spec(cfg$rng)
  cohort_seed <- dtu_cohort_seed(cfg$seed, cohort_name)
  rng_contract <- dtu_rng_contract(rng, cfg$seed, cohort_name)
  cache_path <- file.path(m06$cache, paste0("drimseq_fit_", cohort_name, ".rds"))
  fit_hash <- digest::digest(
    list(
      input_fingerprint = input_fingerprint,
      cohort = cohort_name,
      testing_genes = prepared$testing_family$gene_id,
      feature_ids = prepared$feature_map$feature_id,
      design = prepared$design,
      coefficient = prepared$coefficient,
      DRIMSeq_version = as.character(utils::packageVersion("DRIMSeq")),
      model_contract_version = cfg$model_contract_version,
      rng_contract = rng_contract
    ),
    algo = "sha256"
  )

  cached <- if (!cfg$force_refit && file.exists(cache_path)) readRDS(cache_path) else NULL
  if (!is.null(cached) && identical(cached$hash, fit_hash) &&
      !is.null(cached$result$precision_grid_audit) && !is.null(cached$result$precision_boundary_summary)) {
    log_msg("Reusing ", if (identical(cfg$paired, FALSE)) "unpaired" else "paired", " DRIMSeq cache for ", cohort_name)
    write_table(cached$result$reference_audit, file.path(m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_reference_isoform_audit.tsv"))
    write_table(cached$result$precision_grid_audit, file.path(m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_precision_grid_audit.tsv.gz"))
    write_table(cached$result$precision_boundary_summary, file.path(m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_precision_boundary_summary.tsv"))
    return(cached$result)
  }

  # Keep the pre-correction cache for provenance; the model-contract hash forces
  # one complete refit rather than mixing old and corrected precision estimates.
  if (file.exists(cache_path)) {
    archive_dir <- file.path(m06$cache, "archive")
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    archive_path <- file.path(archive_dir, paste0(
      tools::file_path_sans_ext(basename(cache_path)), "_",
      format(Sys.time(), "%Y%m%dT%H%M%S"), ".rds"
    ))
    if (!file.copy(cache_path, archive_path, overwrite = FALSE)) {
      stop("Could not preserve the previous DRIMSeq cache: ", cache_path)
    }
  }
  ordered <- m06_reference_order(
    as.data.frame(prepared$counts), as.data.frame(prepared$samples), paired = !identical(cfg$paired, FALSE)
  )
  dm_for_fit <- DRIMSeq::dmDSdata(
    counts = ordered$counts, samples = as.data.frame(prepared$samples)
  )
  write_table(ordered$audit, file.path(
    m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_reference_isoform_audit.tsv"
  ))
  engine <- make_m06_drimseq_engine()
  bpparam <- make_bpparam(cfg$n_workers)
  fit_once <- function(dm_object) {
    # Constructing SerialParam above precedes the scoped RNG selection, exactly
    # as in the audited original fit. No installed DRIMSeq namespace is changed.
    dtu_with_rng(rng, cohort_seed, {
      fitted <- engine$dmPrecision(dm_object, design = prepared$design, BPPARAM = bpparam)
      fitted <- engine$dmFit(fitted, design = prepared$design, BPPARAM = bpparam)
      engine$dmTest(fitted, coef = prepared$coefficient, BPPARAM = bpparam)
    })
  }

  initial_error <- NA_character_
  failed_genes <- character()
  fitted_dm <- tryCatch(fit_once(dm_for_fit), error = function(e) {
    stop(if (identical(cfg$paired, FALSE)) "Unpaired" else "Paired", " DRIMSeq fit failed for ", cohort_name,
         "; the frozen family was retained and no genes were removed: ",
         conditionMessage(e), call. = FALSE)
  })

  profile_detail <- as_tibble(engine$profile_audit())
  validate_columns(profile_detail, c(
    "gene_id", "grid_precision", "grid_log_likelihood", "grid_CR_adjustment",
    "grid_adjusted_log_likelihood"
  ), "DRIMSeq numerical profile audit")
  profile_summary <- profile_detail %>%
    group_by(gene_id) %>%
    summarise(
      precision_grid_points = n(),
      nonfinite_grid_likelihoods = sum(!is.finite(grid_log_likelihood)),
      nonfinite_grid_CR_adjustments = sum(!is.finite(grid_CR_adjustment)),
      nonfinite_grid_profiles = sum(!is.finite(grid_adjusted_log_likelihood)),
      .groups = "drop"
    ) %>%
    left_join(ordered$audit, by = "gene_id")
  if (!setequal(profile_summary$gene_id, prepared$testing_family$gene_id) ||
      any(profile_summary$precision_grid_points != 21L)) {
    stop("The complete 21-point DRIMSeq precision grid was not audited for ", cohort_name)
  }
  write_table(profile_detail, file.path(
    m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_precision_grid_audit.tsv.gz"
  ))
  write_table(profile_summary, file.path(
    m06_root, cohort_folder[[cohort_name]], "qc", "DRIMSeq_precision_boundary_summary.tsv"
  ))
  gene_success_raw <- as_tibble(DRIMSeq::results(fitted_dm, level = "gene"))
  for (optional_column in c("lr", "df")) {
    if (!optional_column %in% names(gene_success_raw)) gene_success_raw[[optional_column]] <- NA_real_
  }
  precision_table <- as_tibble(DRIMSeq::genewise_precision(fitted_dm)) %>%
    transmute(gene_id = as.character(gene_id), precision = safe_numeric(genewise_precision))
  likelihood_table <- tibble(
    gene_id = names(methods::slot(fitted_dm, "lik_full")),
    full_model_log_likelihood = as.numeric(methods::slot(fitted_dm, "lik_full")),
    null_model_log_likelihood = as.numeric(methods::slot(fitted_dm, "lik_null"))
  )
  gene_success <- gene_success_raw %>%
    select(-any_of("precision")) %>%
    left_join(precision_table, by = "gene_id") %>%
    left_join(likelihood_table, by = "gene_id") %>%
    left_join(profile_summary, by = "gene_id") %>%
    mutate(
      CR_boundary_associated_P1 = !is.finite(precision) &
        precision_grid_points == 21L & nonfinite_grid_likelihoods == 0L &
        nonfinite_grid_CR_adjustments > 0L &
        ((if (identical(cfg$paired, FALSE)) FALSE else isoforms_absent_in_at_least_one_pair > 0L) |
           isoforms_absent_in_primary > 0L | isoforms_absent_in_recurrent > 0L) &
        reference_is_best_available,
      numerical_boundary_reason = case_when(
        !CR_boundary_associated_P1 ~ NA_character_,
        !identical(cfg$paired, FALSE) & patient_pairs_without_gene_counts > 0L ~ "No gene counts in one or more patient pairs; paired nuisance effects lack information",
        isoforms_absent_in_primary > 0L | isoforms_absent_in_recurrent > 0L ~ "Condition-specific zero isoform counts; boundary precision profile",
        TRUE ~ if (identical(cfg$paired, FALSE)) "Independent-group boundary precision profile" else "Patient-specific zero isoform counts; boundary precision profile"
      )
    ) %>%
    transmute(
      gene_id = as.character(gene_id),
      precision = safe_numeric(precision),
      full_model_log_likelihood, null_model_log_likelihood,
      precision_grid_points, nonfinite_grid_likelihoods,
      nonfinite_grid_CR_adjustments, nonfinite_grid_profiles,
      isoforms_absent_in_at_least_one_pair,
      isoform_patient_pairs_with_zero_counts, CR_boundary_associated_P1,
      numerical_boundary_reason, reference_is_best_available,
      patient_pairs_without_gene_counts, samples_without_gene_counts,
      design_rank_with_gene_counts, isoforms_absent_in_primary, isoforms_absent_in_recurrent,
      numerical_failure_stage = case_when(
        CR_boundary_associated_P1 ~ "precision_CR_boundary_associated",
        !is.finite(precision) ~ "precision",
        !is.finite(full_model_log_likelihood) ~ "full_model",
        !is.finite(null_model_log_likelihood) ~ "null_model",
        !is.finite(pvalue) ~ "likelihood_ratio_test",
        TRUE ~ "none"
      ),
      likelihood_ratio = safe_numeric(lr),
      df = safe_numeric(df),
      dtu_p_value = safe_numeric(pvalue)
    )
  feature_success <- as_tibble(DRIMSeq::results(fitted_dm, level = "feature")) %>%
    transmute(
      gene_id = as.character(gene_id),
      feature_id = as.character(feature_id),
      isoform_dtu_p_value = safe_numeric(pvalue)
    )

  primary_samples <- as.character(prepared$samples$sample_id[
    as.character(prepared$samples$timepoint) == "primary"
  ])
  recurrent_samples <- as.character(prepared$samples$sample_id[
    as.character(prepared$samples$timepoint) == "recurrent"
  ])
  if (!length(primary_samples) || !length(recurrent_samples)) {
    stop("Could not identify both primary and recurrent samples for ", cohort_name)
  }
  count_diagnostics <- prepared$counts
  count_diagnostics$feature_total_primary <- rowSums(
    as.matrix(count_diagnostics[, primary_samples, drop = FALSE]), na.rm = TRUE
  )
  count_diagnostics$feature_total_recurrent <- rowSums(
    as.matrix(count_diagnostics[, recurrent_samples, drop = FALSE]), na.rm = TRUE
  )
  count_diagnostics <- count_diagnostics %>%
    group_by(gene_id) %>%
    summarise(
      tested_isoforms_in_boundary_diagnostic = n(),
      primary_zero_recurrent_positive_isoforms = sum(
        feature_total_primary == 0 & feature_total_recurrent > 0
      ),
      recurrent_zero_primary_positive_isoforms = sum(
        feature_total_recurrent == 0 & feature_total_primary > 0
      ),
      condition_specific_zero_boundary =
        primary_zero_recurrent_positive_isoforms > 0 |
        recurrent_zero_primary_positive_isoforms > 0,
      .groups = "drop"
    )

  family <- prepared$testing_family %>%
    select(gene_id) %>%
    left_join(gene_success, by = "gene_id") %>%
    left_join(count_diagnostics, by = "gene_id") %>%
    mutate(
      raw_dtu_p_value = dtu_p_value,
      fit_status = case_when(
        is.finite(raw_dtu_p_value) ~ "finite_test",
        gene_id %in% failed_genes ~ "fit_error_p1",
        CR_boundary_associated_P1 %in% TRUE ~ "boundary_nonfinite_p1",
        TRUE ~ "other_nonfinite_p1"
      ),
      dtu_p_value = if_else(is.finite(dtu_p_value), pmin(pmax(dtu_p_value, 0), 1), 1),
      dtu_bh = p.adjust(dtu_p_value, method = "BH"),
      dtu_primary_bh05 = dtu_bh < cfg$fdr$primary,
      dtu_exploratory_bh10 = dtu_bh < cfg$fdr$exploratory,
      dtu_exploratory_only_bh10 = dtu_bh >= cfg$fdr$primary & dtu_bh < cfg$fdr$exploratory,
      dtu_fdr = dtu_bh, dtu_fdr05 = dtu_primary_bh05,
      dtu_fdr10 = dtu_exploratory_bh10, dtu_fdr05_to10 = dtu_exploratory_only_bh10
    )

  failure_table <- family %>%
    filter(fit_status != "finite_test") %>%
    transmute(
      cohort = cohort_name,
      gene_id,
      fit_status, numerical_failure_stage, precision,
      full_model_log_likelihood, null_model_log_likelihood,
      precision_grid_points, nonfinite_grid_likelihoods,
      nonfinite_grid_CR_adjustments, nonfinite_grid_profiles,
      isoforms_absent_in_at_least_one_pair,
      isoform_patient_pairs_with_zero_counts, CR_boundary_associated_P1,
      numerical_boundary_reason, reference_is_best_available,
      patient_pairs_without_gene_counts, samples_without_gene_counts,
      design_rank_with_gene_counts, isoforms_absent_in_primary, isoforms_absent_in_recurrent,
      tested_isoforms_in_boundary_diagnostic,
      primary_zero_recurrent_positive_isoforms,
      recurrent_zero_primary_positive_isoforms,
      condition_specific_zero_boundary,
      initial_full_fit_error = initial_error
    )

  result <- list(
    cohort = cohort_name,
    rng_contract = rng_contract,
    model_contract_version = cfg$model_contract_version,
    R_version = R.version.string,
    DRIMSeq_version = as.character(utils::packageVersion("DRIMSeq")),
    fitted_dm = fitted_dm,
    gene = family,
    feature = feature_success,
    failures = failure_table,
    initial_full_fit_error = initial_error,
    reference_audit = ordered$audit,
    precision_grid_audit = profile_detail,
    precision_boundary_summary = profile_summary,
    numerical_adapter = "DRIMSeq 1.38.0: corrected Hessian, supported reference, exact one-way/BB boundaries, stable DM retry with convergence check"
  )
  atomic_save_rds(list(hash = fit_hash, result = result), cache_path)
  result
}

extract_adjusted_component <- function(component, id_name, value_name, fallback_ids = NULL) {
  if (is.null(component)) return(tibble(!!id_name := character(), !!value_name := numeric()))
  if (is.matrix(component) || is.data.frame(component)) {
    values <- as.numeric(component[, 1])
    identifiers <- rownames(component)
  } else {
    values <- as.numeric(component)
    identifiers <- names(component)
  }
  if (is.null(identifiers) || length(identifiers) != length(values)) identifiers <- fallback_ids
  if (is.null(identifiers) || length(identifiers) != length(values)) {
    stop("stageR returned adjusted values without recoverable identifiers.")
  }
  tibble(!!id_name := as.character(identifiers), !!value_name := values)
}

normalize_stageR_adjusted <- function(adjusted, tx2gene, p_screen,
                                      p_confirmation_matrix, alpha_tag) {
  gene_value_name <- paste0("stageR_gene_adjusted_", alpha_tag)
  feature_value_name <- paste0("stageR_isoform_adjusted_", alpha_tag)

  # Current stageR returns one transcript-indexed table with columns
  # geneID, txID, gene, and transcript. Retain the legacy list branch so
  # older Bioconductor installations fail safely rather than misalign IDs.
  if (is.matrix(adjusted) || is.data.frame(adjusted)) {
    adjusted_table <- as_tibble(as.data.frame(adjusted, stringsAsFactors = FALSE))
    validate_columns(
      adjusted_table,
      c("geneID", "txID", "gene", "transcript"),
      "stageR adjusted transcript table"
    )
    feature_out <- adjusted_table %>%
      transmute(
        feature_id = as.character(txID),
        gene_id = as.character(geneID),
        !!feature_value_name := safe_numeric(transcript)
      ) %>%
      distinct(feature_id, .keep_all = TRUE)
    gene_out <- adjusted_table %>%
      transmute(
        gene_id = as.character(geneID),
        !!gene_value_name := safe_numeric(gene)
      ) %>%
      group_by(gene_id) %>%
      summarise(
        !!gene_value_name := {
          values <- unique(.data[[gene_value_name]][is.finite(.data[[gene_value_name]])])
          if (length(values)) values[[1]] else NA_real_
        },
        .groups = "drop"
      )
    return(list(gene = gene_out, feature = feature_out, return_shape = "transcript_table"))
  }

  if (is.list(adjusted) && all(c("gene", "transcript") %in% names(adjusted))) {
    gene_out <- extract_adjusted_component(
      adjusted$gene, "gene_id", gene_value_name,
      fallback_ids = names(p_screen)
    ) %>%
      distinct(gene_id, .keep_all = TRUE)
    feature_out <- extract_adjusted_component(
      adjusted$transcript, "feature_id", feature_value_name,
      fallback_ids = rownames(p_confirmation_matrix)
    ) %>%
      left_join(
        rename(tx2gene, feature_id = transcript, gene_id = gene),
        by = "feature_id"
      ) %>%
      distinct(feature_id, .keep_all = TRUE)
    return(list(gene = gene_out, feature = feature_out, return_shape = "legacy_list"))
  }

  stop(
    "Unsupported stageR getAdjustedPValues() return structure. Class: ",
    paste(class(adjusted), collapse = "/"),
    ". Expected the documented transcript table or a legacy gene/transcript list."
  )
}

run_stageR_alpha <- function(cohort_name, alpha) {
  fit <- dtu_fits[[cohort_name]]
  scheduled_features <- dtu_prepared[[cohort_name]]$feature_map %>%
    select(gene_id, feature_id) %>% distinct()
  if (anyDuplicated(scheduled_features$feature_id) ||
      anyDuplicated(fit$feature$feature_id) || anyDuplicated(fit$gene$gene_id) ||
      !setequal(scheduled_features$gene_id, fit$gene$gene_id)) {
    stop("stageR frozen gene/feature family is inconsistent for ", cohort_name)
  }
  feature <- scheduled_features %>%
    left_join(fit$feature, by = c("gene_id", "feature_id")) %>%
    mutate(
      isoform_dtu_p_value = if_else(
        is.finite(isoform_dtu_p_value),
        pmin(pmax(isoform_dtu_p_value, 0), 1),
        1
      )
    ) %>%
    distinct(feature_id, .keep_all = TRUE)

  if (!nrow(feature)) {
    return(list(
      available = FALSE,
      reason = "no_successful_isoform_tests",
      object = NULL,
      gene = tibble(gene_id = character()),
      feature = tibble(feature_id = character(), gene_id = character())
    ))
  }

  tx2gene <- feature %>%
    transmute(transcript = feature_id, gene = gene_id) %>%
    distinct(transcript, .keep_all = TRUE)
  gene_screen <- fit$gene %>%
    filter(gene_id %in% tx2gene$gene) %>%
    distinct(gene_id, .keep_all = TRUE)

  # Retain every scheduled gene AND its real transcripts. stageR uses the
  # number of supplied genes in its confirmation scaling as well as the screen.
  # Supplying full-family BH values after dropping failed genes is insufficient.
  stopifnot(nrow(gene_screen) == nrow(fit$gene))
  p_screen <- gene_screen$dtu_bh
  names(p_screen) <- gene_screen$gene_id
  p_confirmation <- feature$isoform_dtu_p_value
  names(p_confirmation) <- feature$feature_id
  p_confirmation <- p_confirmation[tx2gene$transcript]
  p_confirmation_matrix <- matrix(
    as.numeric(p_confirmation),
    ncol = 1L,
    dimnames = list(tx2gene$transcript, "recurrent_vs_primary")
  )

  stage_object <- stageR::stageRTx(
    pScreen = p_screen,
    pConfirmation = p_confirmation_matrix,
    tx2gene = as.data.frame(tx2gene),
    pScreenAdjusted = TRUE
  )
  stage_object <- stageR::stageWiseAdjustment(
    stage_object,
    method = "dtu",
    alpha = alpha
  )
  adjusted <- stageR::getAdjustedPValues(
    stage_object,
    onlySignificantGenes = FALSE,
    order = FALSE
  )

  alpha_tag <- dtu_alpha_tag(alpha)
  normalized <- normalize_stageR_adjusted(
    adjusted, tx2gene, p_screen, p_confirmation_matrix, alpha_tag
  )

  list(
    available = TRUE,
    reason = "success",
    object = stage_object,
    gene = normalized$gene,
    feature = normalized$feature,
    return_shape = normalized$return_shape,
    scheduled_genes = nrow(fit$gene),
    screening_genes = length(p_screen),
    confirmation_family_genes = n_distinct(tx2gene$gene),
    screening_discoveries = sum(p_screen <= alpha),
    confirmation_alpha = alpha * sum(p_screen <= alpha) / length(p_screen),
    conservative_feature_P1 = sum(!is.finite(
      scheduled_features %>% left_join(fit$feature, by = c("gene_id", "feature_id")) %>%
        pull(isoform_dtu_p_value)
    ))
  )
}

build_broad_structural_category_means <- function(cohort_name) {
  samples <- dtu_prepared[[cohort_name]]$samples %>%
    transmute(
      sample_id = as.character(sample_id),
      timepoint = as.character(timepoint)
    )
  annotation <- isoform_annotation %>%
    transmute(
      matrix_row = row_number(),
      gene_id = clean_character(associated_gene),
      structural_category_display = display_structural_category(structural_category),
      gene_mapped = coalesce(as.logical(gene_mapped), FALSE),
      cross_sample_comparable = coalesce(as.logical(cross_sample_comparable), FALSE),
      analyzed_for_gain_loss = coalesce(as.logical(analyzed_for_gain_loss), FALSE)
    ) %>%
    filter(
      gene_mapped, cross_sample_comparable, analyzed_for_gain_loss,
      !is.na(gene_id),
      gene_id %in% dtu_prepared[[cohort_name]]$testing_family$gene_id
    )
  count_mat <- counts_matrix[
    annotation$matrix_row,
    samples$sample_id,
    drop = FALSE
  ]
  detected <- rowSums(count_mat, na.rm = TRUE) > 0
  annotation <- annotation[detected, , drop = FALSE]
  count_mat <- count_mat[detected, , drop = FALSE]
  if (!nrow(annotation)) {
    return(tibble(
      gene_id = character(), structural_category_display = character(),
      mean_primary_category_fraction = numeric(),
      mean_recurrent_category_fraction = numeric()
    ))
  }

  annotation <- annotation %>%
    mutate(
      category_key = paste(
        gene_id, structural_category_display,
        sep = "___M06_STRUCTURAL_CATEGORY___"
      )
    )
  category_lookup <- annotation %>%
    distinct(category_key, gene_id, structural_category_display)
  category_counts <- rowsum(count_mat, annotation$category_key, reorder = FALSE)

  category_long <- as_tibble(category_counts, rownames = "category_key") %>%
    pivot_longer(
      all_of(samples$sample_id),
      names_to = "sample_id",
      values_to = "category_count"
    ) %>%
    left_join(category_lookup, by = "category_key") %>%
    group_by(gene_id, sample_id) %>%
    mutate(
      gene_count = sum(category_count),
      category_fraction = ifelse(gene_count > 0, category_count / gene_count, NA_real_)
    ) %>%
    ungroup() %>%
    left_join(samples, by = "sample_id")

  output <- category_long %>%
    group_by(gene_id, structural_category_display, timepoint) %>%
    summarise(
      mean_category_fraction = finite_mean(category_fraction),
      evaluable_samples = sum(is.finite(category_fraction)),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = timepoint,
      values_from = c(mean_category_fraction, evaluable_samples),
      names_sep = "_"
    )
  for (column_name in c(
    "mean_category_fraction_primary", "mean_category_fraction_recurrent"
  )) {
    if (!column_name %in% names(output)) output[[column_name]] <- NA_real_
  }
  output %>%
    rename(
      mean_primary_category_fraction = mean_category_fraction_primary,
      mean_recurrent_category_fraction = mean_category_fraction_recurrent
    )
}
