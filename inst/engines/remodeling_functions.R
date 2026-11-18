# Source-derived statistical functions from Module 04 v0.6.0.
# Jensen-Shannon/dominance/bootstrap/additive REML/profile mixture/GPD implementation.
# Exact published engine functions are retained; orchestration is public and metadata-driven.
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

parse_numeric_csv <- function(x, integer = FALSE) {
  values <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  output <- suppressWarnings(as.numeric(values))
  if (!length(output) || any(!is.finite(output))) {
    stop("Could not parse numeric parameter list: ", x)
  }
  if (integer) output <- as.integer(round(output))
  sort(unique(output))
}

safe_norm <- function(path) {
  ifelse(
    is.na(path) | !nzchar(path),
    NA_character_,
    normalizePath(path, winslash = "/", mustWork = FALSE)
  )
}

atomic_save_rds <- function(object, path) {
  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  previous_path <- paste0(path, ".previous")
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, compress = "gzip")
  invisible(readRDS(temporary_path))
  if (file.exists(previous_path)) unlink(previous_path)
  had_previous <- file.exists(path)
  if (had_previous && !file.rename(path, previous_path)) {
    stop("Could not protect the previous cache: ", path)
  }
  if (!file.rename(temporary_path, path)) {
    if (had_previous && file.exists(previous_path)) file.rename(previous_path, path)
    stop("Could not finalize cache: ", path)
  }
  if (file.exists(previous_path)) unlink(previous_path)
  invisible(path)
}

deterministic_seed <- function(master_seed, ...) {
  key <- paste(c("M04", as.character(master_seed), as.character(list(...))), collapse = "|")
  hash <- digest::digest(key, algo = "sha256", serialize = FALSE)
  high <- strtoi(substr(hash, 1L, 6L), base = 16L)
  low <- strtoi(substr(hash, 7L, 12L), base = 16L)
  folded <- as.double(high) * 16777216 + as.double(low)
  as.integer((folded %% (.Machine$integer.max - 1L)) + 1L)
}

with_seed <- function(seed, expression) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(expression)
}


probability_vector <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0) || sum(x) <= 0) {
    return(NULL)
  }
  x / sum(x)
}

entropy_bits <- function(p) {
  p <- probability_vector(p)
  if (is.null(p)) return(NA_real_)
  keep <- p > 0
  -sum(p[keep] * log2(p[keep]))
}

jsd_base2 <- function(p, q) {
  p <- probability_vector(p)
  q <- probability_vector(q)
  if (is.null(p) || is.null(q) || length(p) != length(q)) return(NA_real_)
  midpoint <- 0.5 * (p + q)
  kl2 <- function(a, b) {
    keep <- a > 0
    if (!any(keep)) return(0)
    sum(a[keep] * log2(a[keep] / b[keep]))
  }
  min(1, max(0, 0.5 * kl2(p, midpoint) + 0.5 * kl2(q, midpoint)))
}

jsd_natural <- function(p, q) {
  p <- probability_vector(p)
  q <- probability_vector(q)
  if (is.null(p) || is.null(q) || length(p) != length(q)) return(NA_real_)
  midpoint <- 0.5 * (p + q)
  kln <- function(a, b) {
    keep <- a > 0
    if (!any(keep)) return(0)
    sum(a[keep] * log(a[keep] / b[keep]))
  }
  0.5 * kln(p, midpoint) + 0.5 * kln(q, midpoint)
}

jsd_contributions <- function(p, q) {
  p <- probability_vector(p)
  q <- probability_vector(q)
  if (is.null(p) || is.null(q) || length(p) != length(q)) {
    return(rep(NA_real_, max(length(p), length(q))))
  }
  midpoint <- 0.5 * (p + q)
  contribution_p <- contribution_q <- rep(0, length(p))
  keep_p <- p > 0
  keep_q <- q > 0
  contribution_p[keep_p] <- p[keep_p] * log2(p[keep_p] / midpoint[keep_p])
  contribution_q[keep_q] <- q[keep_q] * log2(q[keep_q] / midpoint[keep_q])
  pmax(0, 0.5 * contribution_p + 0.5 * contribution_q)
}

total_variation <- function(p, q) {
  p <- probability_vector(p)
  q <- probability_vector(q)
  if (is.null(p) || is.null(q) || length(p) != length(q)) return(NA_real_)
  0.5 * sum(abs(p - q))
}

joint_score <- function(jsd, max_delta_if) {
  pmax(jsd / cfg$thresholds$jsd, max_delta_if / cfg$thresholds$delta_if)
}

top_unique_index <- function(p, tolerance = 1e-12) {
  if (!length(p) || anyNA(p)) return(NA_integer_)
  top_value <- max(p)
  positions <- which(abs(p - top_value) <= tolerance)
  if (length(positions) != 1L) return(NA_integer_)
  positions
}

dominance_summary <- function(group_mean, sample_compositions, transcript_ids) {
  if (!length(group_mean) || anyNA(group_mean)) {
    return(list(
      has_dominant = FALSE, top_id = NA_character_, dominant_id = NA_character_,
      dominant_mean_if = NA_real_, runner_up_mean_if = NA_real_,
      dominance_margin = NA_real_, support_n = 0L,
      evaluable_n = nrow(sample_compositions), mean_if_pass = FALSE,
      margin_pass = FALSE, support_pass = FALSE,
      reason = "dominance_summary_unavailable"
    ))
  }
  order_index <- order(group_mean, decreasing = TRUE, transcript_ids)
  top_index <- order_index[[1]]
  runner_index <- if (length(order_index) >= 2L) order_index[[2]] else NA_integer_
  top_if <- group_mean[[top_index]]
  runner_if <- if (is.na(runner_index)) 0 else group_mean[[runner_index]]
  margin <- top_if - runner_if
  support <- if (nrow(sample_compositions)) {
    apply(sample_compositions, 1L, function(composition) {
      index <- top_unique_index(composition)
      !is.na(index) && index == top_index &&
        composition[[index]] >= cfg$thresholds$dominance_sample_min_if
    })
  } else logical()
  support_n <- sum(support %in% TRUE)
  mean_if_pass <- top_if >= cfg$thresholds$dominance_mean_if
  margin_pass <- margin >= cfg$thresholds$dominance_margin
  support_pass <- support_n >= cfg$thresholds$dominance_support_n
  has_dominant <- mean_if_pass && margin_pass && support_pass
  failed_n <- sum(!c(mean_if_pass, margin_pass, support_pass))
  reason <- if (has_dominant) {
    "clear_dominant"
  } else if (failed_n > 1L) {
    "multiple_criteria_failed"
  } else if (!mean_if_pass) {
    "top_mean_if_below_threshold"
  } else if (!margin_pass) {
    "dominance_margin_below_threshold"
  } else {
    "patient_support_below_threshold"
  }
  list(
    has_dominant = has_dominant,
    top_id = transcript_ids[[top_index]],
    dominant_id = if (has_dominant) transcript_ids[[top_index]] else NA_character_,
    dominant_mean_if = top_if,
    runner_up_mean_if = runner_if,
    dominance_margin = margin,
    support_n = as.integer(support_n),
    evaluable_n = as.integer(nrow(sample_compositions)),
    mean_if_pass = mean_if_pass,
    margin_pass = margin_pass,
    support_pass = support_pass,
    reason = reason
  )
}

secondary_remodeling <- function(
    primary_mean, recurrent_mean, primary_dominant, recurrent_dominant,
    primary_has_dominant, recurrent_has_dominant, transcript_ids
) {
  unavailable <- function(primary_mass = NA_real_, recurrent_mass = NA_real_) {
    list(
      secondary_jsd = NA_real_, secondary_max_abs_delta_if = NA_real_,
      secondary_mass_primary = primary_mass,
      secondary_mass_recurrent = recurrent_mass,
      secondary_driver = NA_character_, secondary_driver_delta_if = NA_real_,
      secondary_mass_pass = FALSE, secondary_effect_pass = FALSE
    )
  }
  if (!isTRUE(primary_has_dominant) || !isTRUE(recurrent_has_dominant) ||
      is.na(primary_dominant) || is.na(recurrent_dominant) ||
      !identical(primary_dominant, recurrent_dominant)) return(unavailable())
  keep <- transcript_ids != primary_dominant
  primary_mass <- sum(primary_mean[keep])
  recurrent_mass <- sum(recurrent_mean[keep])
  if (!any(keep) || primary_mass <= 0 || recurrent_mass <= 0) {
    return(unavailable(primary_mass, recurrent_mass))
  }
  secondary_delta <- recurrent_mean[keep] - primary_mean[keep]
  driver_index <- which.max(abs(secondary_delta))[[1]]
  secondary_jsd <- jsd_base2(
    primary_mean[keep] / primary_mass,
    recurrent_mean[keep] / recurrent_mass
  )
  secondary_max <- max(abs(secondary_delta))
  mass_pass <- min(primary_mass, recurrent_mass) >= cfg$thresholds$secondary_min_mass
  effect_pass <- mass_pass &&
    (secondary_jsd >= cfg$thresholds$jsd ||
       secondary_max >= cfg$thresholds$delta_if)
  list(
    secondary_jsd = secondary_jsd,
    secondary_max_abs_delta_if = secondary_max,
    secondary_mass_primary = primary_mass,
    secondary_mass_recurrent = recurrent_mass,
    secondary_driver = transcript_ids[keep][[driver_index]],
    secondary_driver_delta_if = secondary_delta[[driver_index]],
    secondary_mass_pass = mass_pass,
    secondary_effect_pass = effect_pass
  )
}

classify_architecture <- function(
    quantifiable, effect_pass, primary_dominance, recurrent_dominance,
    secondary
) {
  if (!isTRUE(quantifiable)) return("NQ")
  if (!isTRUE(effect_pass)) return("Below_effect")
  primary_has <- isTRUE(primary_dominance$has_dominant)
  recurrent_has <- isTRUE(recurrent_dominance$has_dominant)
  if (primary_has && recurrent_has) {
    if (!is.na(primary_dominance$dominant_id) &&
        !is.na(recurrent_dominance$dominant_id) &&
        identical(primary_dominance$dominant_id, recurrent_dominance$dominant_id)) {
      return(if (isTRUE(secondary$secondary_effect_pass)) "B" else "Unresolved_effect")
    }
    return("A")
  }
  if (primary_has && !recurrent_has) return("C")
  if (!primary_has && recurrent_has) return("D")
  "E"
}

cohort_gene_view <- function(object, included_patients) {
  pair_table <- sample_by_patient %>% filter(patient %in% included_patients)
  valid <- object$pair_quantifiable[pair_table$patient] %in% TRUE
  pair_table <- pair_table[valid, , drop = FALSE]
  if (object$n_isoforms_retained < 2L ||
      nrow(pair_table) < cfg$thresholds$min_evaluable_pairs) {
    return(list(quantifiable = FALSE, pair_table = pair_table))
  }
  primary_counts <- object$counts[pair_table$primary, , drop = FALSE]
  recurrent_counts <- object$counts[pair_table$recurrent, , drop = FALSE]
  primary_composition <- primary_counts / rowSums(primary_counts)
  recurrent_composition <- recurrent_counts / rowSums(recurrent_counts)
  list(
    quantifiable = TRUE,
    pair_table = pair_table,
    primary_counts = primary_counts,
    recurrent_counts = recurrent_counts,
    primary_composition = primary_composition,
    recurrent_composition = recurrent_composition
  )
}

analyze_gene_observed <- function(object, included_patients, cohort) {
  view <- cohort_gene_view(object, included_patients)
  base <- tibble(
    cohort = cohort,
    gene_id = object$gene_id,
    quantifiable = isTRUE(view$quantifiable),
    n_isoforms_total = object$n_isoforms_total,
    n_isoforms_retained = object$n_isoforms_retained,
    pooled_gene_FL = object$pooled_gene_count,
    pooled_retained_FL = object$pooled_retained_count,
    n_evaluable_pairs = nrow(view$pair_table)
  )
  if (!isTRUE(view$quantifiable)) {
    return(list(
      gene = base %>% mutate(
        quantifiability_reason = case_when(
          n_isoforms_retained < 2L ~ "fewer_than_two_retained_isoforms",
          TRUE ~ "insufficient_complete_quantifiable_pairs"
        ),
        jsd = NA_real_, jsd_natural_log = NA_real_, max_abs_delta_if = NA_real_,
        total_variation_distance = NA_real_, joint_effect_score = NA_real_,
        entropy_primary = NA_real_, entropy_recurrent = NA_real_,
        delta_entropy = NA_real_, median_paired_delta_entropy = NA_real_,
        effective_isoform_ratio = NA_real_, driver_isoform = NA_character_,
        driver_delta_if = NA_real_, driver_patient_support_n = NA_integer_,
        has_dominant_primary = NA, dominant_primary = NA_character_,
        dominant_primary_mean_if = NA_real_, dominant_primary_margin = NA_real_,
        dominant_primary_support_n = NA_integer_, dominance_primary_reason = NA_character_,
        has_dominant_recurrent = NA, dominant_recurrent = NA_character_,
        dominant_recurrent_mean_if = NA_real_, dominant_recurrent_margin = NA_real_,
        dominant_recurrent_support_n = NA_integer_, dominance_recurrent_reason = NA_character_,
        secondary_jsd = NA_real_, secondary_max_abs_delta_if = NA_real_,
        secondary_mass_primary = NA_real_, secondary_mass_recurrent = NA_real_,
        secondary_driver = NA_character_, secondary_driver_delta_if = NA_real_,
        secondary_mass_pass = NA, secondary_effect_pass = NA,
        effect_pass = FALSE, dual_effect_pass = FALSE,
        observed_architecture = "NQ", entropy_direction = NA_character_
      ),
      isoform = tibble(),
      patient = tibble()
    ))
  }
  primary_mean <- colMeans(view$primary_composition)
  recurrent_mean <- colMeans(view$recurrent_composition)
  names(primary_mean) <- names(recurrent_mean) <- object$transcript_ids
  delta_if <- recurrent_mean - primary_mean
  absolute_delta <- abs(delta_if)
  max_delta <- max(absolute_delta)
  driver_candidates <- which(abs(absolute_delta - max_delta) <= 1e-12)
  driver_index <- driver_candidates[[1]]
  driver <- object$transcript_ids[[driver_index]]
  driver_patient_delta <- view$recurrent_composition[, driver_index] -
    view$primary_composition[, driver_index]
  expected_sign <- sign(delta_if[[driver_index]])
  driver_support <- sign(driver_patient_delta) == expected_sign &
    abs(driver_patient_delta) >= cfg$thresholds$direction_delta_if
  jsd_value <- jsd_base2(primary_mean, recurrent_mean)
  natural_value <- jsd_natural(primary_mean, recurrent_mean)
  tvd_value <- total_variation(primary_mean, recurrent_mean)
  primary_entropy_by_patient <- apply(view$primary_composition, 1L, entropy_bits)
  recurrent_entropy_by_patient <- apply(view$recurrent_composition, 1L, entropy_bits)
  paired_delta_entropy <- recurrent_entropy_by_patient - primary_entropy_by_patient
  delta_entropy <- mean(paired_delta_entropy)
  primary_dominance <- dominance_summary(
    primary_mean, view$primary_composition, object$transcript_ids
  )
  recurrent_dominance <- dominance_summary(
    recurrent_mean, view$recurrent_composition, object$transcript_ids
  )
  secondary <- secondary_remodeling(
    primary_mean, recurrent_mean,
    primary_dominance$dominant_id, recurrent_dominance$dominant_id,
    primary_dominance$has_dominant, recurrent_dominance$has_dominant,
    object$transcript_ids
  )
  effect_pass <- jsd_value >= cfg$thresholds$jsd ||
    max_delta >= cfg$thresholds$delta_if
  observed_architecture <- classify_architecture(
    TRUE, effect_pass, primary_dominance, recurrent_dominance, secondary
  )
  gene <- base %>% mutate(
    quantifiability_reason = "quantifiable",
    jsd = jsd_value,
    jsd_natural_log = natural_value,
    max_abs_delta_if = max_delta,
    total_variation_distance = tvd_value,
    joint_effect_score = joint_score(jsd_value, max_delta),
    entropy_primary = mean(primary_entropy_by_patient),
    entropy_recurrent = mean(recurrent_entropy_by_patient),
    delta_entropy = delta_entropy,
    median_paired_delta_entropy = median(paired_delta_entropy),
    effective_isoform_ratio = 2^delta_entropy,
    driver_isoform = driver,
    driver_delta_if = delta_if[[driver_index]],
    driver_tie_n = length(driver_candidates),
    driver_patient_support_n = sum(driver_support %in% TRUE),
    has_dominant_primary = primary_dominance$has_dominant,
    dominant_primary = primary_dominance$dominant_id,
    dominant_primary_mean_if = primary_dominance$dominant_mean_if,
    dominant_primary_margin = primary_dominance$dominance_margin,
    dominant_primary_support_n = primary_dominance$support_n,
    dominance_primary_reason = primary_dominance$reason,
    has_dominant_recurrent = recurrent_dominance$has_dominant,
    dominant_recurrent = recurrent_dominance$dominant_id,
    dominant_recurrent_mean_if = recurrent_dominance$dominant_mean_if,
    dominant_recurrent_margin = recurrent_dominance$dominance_margin,
    dominant_recurrent_support_n = recurrent_dominance$support_n,
    dominance_recurrent_reason = recurrent_dominance$reason,
    secondary_jsd = secondary$secondary_jsd,
    secondary_max_abs_delta_if = secondary$secondary_max_abs_delta_if,
    secondary_mass_primary = secondary$secondary_mass_primary,
    secondary_mass_recurrent = secondary$secondary_mass_recurrent,
    secondary_driver = secondary$secondary_driver,
    secondary_driver_delta_if = secondary$secondary_driver_delta_if,
    secondary_mass_pass = secondary$secondary_mass_pass,
    secondary_effect_pass = secondary$secondary_effect_pass,
    effect_pass = effect_pass,
    dual_effect_pass = jsd_value >= cfg$thresholds$jsd &&
      max_delta >= cfg$thresholds$delta_if,
    observed_architecture = observed_architecture,
    entropy_direction = case_when(
      delta_entropy >= cfg$thresholds$entropy_delta ~ "more_diffuse_in_recurrence",
      delta_entropy <= -cfg$thresholds$entropy_delta ~ "more_focused_in_recurrence",
      TRUE ~ "entropy_stable"
    )
  )
  contributions <- jsd_contributions(primary_mean, recurrent_mean)
  isoform <- tibble(
    cohort = cohort,
    gene_id = object$gene_id,
    isoform_harmonized_id = object$transcript_ids,
    mean_if_primary = primary_mean,
    mean_if_recurrent = recurrent_mean,
    delta_if = delta_if,
    abs_delta_if = absolute_delta,
    jsd_contribution = contributions,
    jsd_contribution_fraction = if (jsd_value > 0) contributions / jsd_value else 0,
    is_driver = object$transcript_ids == driver
  )
  patient <- tibble(
    cohort = cohort,
    gene_id = object$gene_id,
    patient = view$pair_table$patient,
    driver_isoform = driver,
    driver_delta_if = as.numeric(driver_patient_delta),
    supports_group_direction = as.logical(driver_support),
    entropy_primary = primary_entropy_by_patient,
    entropy_recurrent = recurrent_entropy_by_patient,
    delta_entropy = paired_delta_entropy
  )
  list(gene = gene, isoform = isoform, patient = patient)
}
jsd_rows <- function(primary, recurrent) {
  primary <- as.matrix(primary)
  recurrent <- as.matrix(recurrent)
  midpoint <- 0.5 * (primary + recurrent)
  primary_term <- recurrent_term <- matrix(0, nrow(primary), ncol(primary))
  keep_primary <- primary > 0
  keep_recurrent <- recurrent > 0
  primary_term[keep_primary] <- primary[keep_primary] *
    log2(primary[keep_primary] / midpoint[keep_primary])
  recurrent_term[keep_recurrent] <- recurrent[keep_recurrent] *
    log2(recurrent[keep_recurrent] / midpoint[keep_recurrent])
  pmin(1, pmax(0, 0.5 * rowSums(primary_term) + 0.5 * rowSums(recurrent_term)))
}

row_means_from_draws <- function(values, draws) {
  values <- as.numeric(values)
  matrix(values[draws], nrow = nrow(draws), ncol = ncol(draws)) %>% rowMeans()
}

dominance_state_matches <- function(candidate, observed_has, observed_id) {
  if (!identical(isTRUE(candidate$has_dominant), isTRUE(observed_has))) return(FALSE)
  if (!isTRUE(observed_has)) return(TRUE)
  identical(as.character(candidate$dominant_id), as.character(observed_id))
}

bootstrap_dominance_matrix <- function(
    group_mean_matrix, sample_compositions, draws, transcript_ids
) {
  group_mean_matrix <- as.matrix(group_mean_matrix)
  b <- nrow(group_mean_matrix)
  k <- ncol(group_mean_matrix)
  row_id <- seq_len(b)
  top_index <- max.col(group_mean_matrix, ties.method = "first")
  top_value <- group_mean_matrix[cbind(row_id, top_index)]
  runner_matrix <- group_mean_matrix
  runner_matrix[cbind(row_id, top_index)] <- -Inf
  runner_index <- max.col(runner_matrix, ties.method = "first")
  runner_value <- if (k > 1L) {
    runner_matrix[cbind(row_id, runner_index)]
  } else rep(0, b)
  margin <- top_value - runner_value
  support_n <- integer(b)
  for (draw_column in seq_len(ncol(draws))) {
    composition <- sample_compositions[draws[, draw_column], , drop = FALSE]
    sample_top_index <- max.col(composition, ties.method = "first")
    sample_top_value <- composition[cbind(row_id, sample_top_index)]
    unique_top <- rowSums(
      abs(composition - sample_top_value) <= 1e-12
    ) == 1L
    support_n <- support_n + as.integer(
      unique_top & sample_top_index == top_index &
        composition[cbind(row_id, top_index)] >=
          cfg$thresholds$dominance_sample_min_if
    )
  }
  mean_if_pass <- top_value >= cfg$thresholds$dominance_mean_if
  margin_pass <- margin >= cfg$thresholds$dominance_margin
  support_pass <- support_n >= cfg$thresholds$dominance_support_n
  has_dominant <- mean_if_pass & margin_pass & support_pass
  top_id <- transcript_ids[top_index]
  list(
    has_dominant = has_dominant,
    top_index = top_index,
    top_id = top_id,
    dominant_id = ifelse(has_dominant, top_id, NA_character_),
    dominant_mean_if = top_value,
    runner_up_mean_if = runner_value,
    dominance_margin = margin,
    support_n = support_n,
    mean_if_pass = mean_if_pass,
    margin_pass = margin_pass,
    support_pass = support_pass
  )
}

bootstrap_gene <- function(object, observed_row, included_patients, cohort_name) {
  view <- cohort_gene_view(object, included_patients)
  if (!isTRUE(view$quantifiable)) {
    return(list(
      summary = tibble(
        cohort = cohort_name, gene_id = object$gene_id,
        bootstrap_replicates = cfg$bootstrap$replicates,
        bootstrap_valid_replicates = 0L, bootstrap_valid_fraction = 0,
        jsd_bootstrap_lower = NA_real_, jsd_bootstrap_upper = NA_real_,
        max_abs_delta_if_bootstrap_lower = NA_real_,
        max_abs_delta_if_bootstrap_upper = NA_real_,
        delta_entropy_bootstrap_lower = NA_real_,
        delta_entropy_bootstrap_upper = NA_real_,
        effect_pass_bootstrap_probability = NA_real_,
        observed_class_bootstrap_probability = NA_real_,
        observed_driver_bootstrap_probability = NA_real_,
        primary_dominance_state_probability = NA_real_,
        recurrent_dominance_state_probability = NA_real_
      ),
      class_probability = tibble()
    ))
  }
  b <- cfg$bootstrap$replicates
  n_pair <- nrow(view$pair_table)
  k <- object$n_isoforms_retained
  seed <- deterministic_seed(cfg$seed, "paired_bootstrap", cohort_name, object$gene_id)
  result <- with_seed(seed, {
    draws <- matrix(sample.int(n_pair, b * n_pair, replace = TRUE), nrow = b, ncol = n_pair)
    primary_boot <- matrix(0, nrow = b, ncol = k)
    recurrent_boot <- matrix(0, nrow = b, ncol = k)
    for (isoform_index in seq_len(k)) {
      primary_boot[, isoform_index] <- row_means_from_draws(
        view$primary_composition[, isoform_index], draws
      )
      recurrent_boot[, isoform_index] <- row_means_from_draws(
        view$recurrent_composition[, isoform_index], draws
      )
    }
    jsd_value <- jsd_rows(primary_boot, recurrent_boot)
    delta_if <- recurrent_boot - primary_boot
    abs_delta <- abs(delta_if)
    max_delta <- apply(abs_delta, 1L, max)
    driver_index <- max.col(abs_delta, ties.method = "first")
    observed_driver_index <- match(observed_row$driver_isoform, object$transcript_ids)
    driver_match <- driver_index == observed_driver_index
    primary_entropy <- apply(view$primary_composition, 1L, entropy_bits)
    recurrent_entropy <- apply(view$recurrent_composition, 1L, entropy_bits)
    paired_delta_entropy <- recurrent_entropy - primary_entropy
    delta_entropy <- row_means_from_draws(paired_delta_entropy, draws)
    effect_pass <- jsd_value >= cfg$thresholds$jsd |
      max_delta >= cfg$thresholds$delta_if
    primary_dominance <- bootstrap_dominance_matrix(
      primary_boot, view$primary_composition, draws, object$transcript_ids
    )
    recurrent_dominance <- bootstrap_dominance_matrix(
      recurrent_boot, view$recurrent_composition, draws, object$transcript_ids
    )
    both_dominant <- primary_dominance$has_dominant &
      recurrent_dominance$has_dominant
    same_dominant <- both_dominant &
      primary_dominance$top_index == recurrent_dominance$top_index

    row_id <- seq_len(b)
    secondary_primary <- primary_boot
    secondary_recurrent <- recurrent_boot
    shared_index <- primary_dominance$top_index
    secondary_primary[cbind(row_id, shared_index)] <- 0
    secondary_recurrent[cbind(row_id, shared_index)] <- 0
    secondary_mass_primary <- rowSums(secondary_primary)
    secondary_mass_recurrent <- rowSums(secondary_recurrent)
    secondary_valid <- same_dominant & secondary_mass_primary > 0 &
      secondary_mass_recurrent > 0
    secondary_jsd <- rep(NA_real_, b)
    if (any(secondary_valid)) {
      normalized_primary <- secondary_primary[secondary_valid, , drop = FALSE] /
        secondary_mass_primary[secondary_valid]
      normalized_recurrent <- secondary_recurrent[secondary_valid, , drop = FALSE] /
        secondary_mass_recurrent[secondary_valid]
      secondary_jsd[secondary_valid] <- jsd_rows(
        normalized_primary, normalized_recurrent
      )
    }
    secondary_delta <- recurrent_boot - primary_boot
    secondary_delta[cbind(row_id, shared_index)] <- 0
    secondary_abs_delta <- abs(secondary_delta)
    secondary_driver_index <- max.col(secondary_abs_delta, ties.method = "first")
    secondary_max_delta <- secondary_abs_delta[
      cbind(row_id, secondary_driver_index)
    ]
    secondary_mass_pass <- pmin(
      secondary_mass_primary, secondary_mass_recurrent
    ) >= cfg$thresholds$secondary_min_mass
    secondary_effect_pass <- secondary_valid & secondary_mass_pass &
      (secondary_jsd >= cfg$thresholds$jsd |
         secondary_max_delta >= cfg$thresholds$delta_if)

    classes <- ifelse(effect_pass, "E", "Below_effect")
    classes[effect_pass & both_dominant & !same_dominant] <- "A"
    classes[effect_pass & same_dominant & secondary_effect_pass] <- "B"
    classes[effect_pass & same_dominant & !secondary_effect_pass] <-
      "Unresolved_effect"
    classes[effect_pass & primary_dominance$has_dominant &
              !recurrent_dominance$has_dominant] <- "C"
    classes[effect_pass & !primary_dominance$has_dominant &
              recurrent_dominance$has_dominant] <- "D"

    primary_dom_match <- if (isTRUE(observed_row$has_dominant_primary)) {
      primary_dominance$has_dominant &
        primary_dominance$dominant_id == observed_row$dominant_primary
    } else {
      !primary_dominance$has_dominant
    }
    recurrent_dom_match <- if (isTRUE(observed_row$has_dominant_recurrent)) {
      recurrent_dominance$has_dominant &
        recurrent_dominance$dominant_id == observed_row$dominant_recurrent
    } else {
      !recurrent_dominance$has_dominant
    }
    list(
      jsd = jsd_value,
      max_delta = max_delta,
      delta_entropy = delta_entropy,
      effect_pass = effect_pass,
      classes = classes,
      driver_match = driver_match,
      primary_dom_match = primary_dom_match,
      recurrent_dom_match = recurrent_dom_match
    )
  })
  valid <- is.finite(result$jsd) & is.finite(result$max_delta) &
    is.finite(result$delta_entropy)
  n_valid <- sum(valid)
  alpha <- 1 - cfg$bootstrap$confidence
  interval <- function(x) {
    if (!any(valid)) return(c(NA_real_, NA_real_))
    as.numeric(quantile(x[valid], c(alpha / 2, 1 - alpha / 2),
                        names = FALSE, type = cfg$bootstrap$quantile_type))
  }
  jsd_interval <- interval(result$jsd)
  delta_interval <- interval(result$max_delta)
  entropy_interval <- interval(result$delta_entropy)
  class_prob <- table(factor(result$classes[valid], levels = null_architecture_levels)) /
    max(1L, n_valid)
  observed_class <- as.character(observed_row$observed_architecture)
  observed_class_probability <- if (observed_class %in% null_architecture_levels) {
    unname(class_prob[[observed_class]])
  } else NA_real_
  list(
    summary = tibble(
      cohort = cohort_name,
      gene_id = object$gene_id,
      bootstrap_replicates = b,
      bootstrap_valid_replicates = n_valid,
      bootstrap_valid_fraction = n_valid / b,
      jsd_bootstrap_lower = jsd_interval[[1]],
      jsd_bootstrap_upper = jsd_interval[[2]],
      max_abs_delta_if_bootstrap_lower = delta_interval[[1]],
      max_abs_delta_if_bootstrap_upper = delta_interval[[2]],
      delta_entropy_bootstrap_lower = entropy_interval[[1]],
      delta_entropy_bootstrap_upper = entropy_interval[[2]],
      effect_pass_bootstrap_probability = if (any(valid)) mean(result$effect_pass[valid]) else NA_real_,
      observed_class_bootstrap_probability = observed_class_probability,
      observed_driver_bootstrap_probability = if (any(valid)) mean(result$driver_match[valid]) else NA_real_,
      primary_dominance_state_probability = if (any(valid)) mean(result$primary_dom_match[valid]) else NA_real_,
      recurrent_dominance_state_probability = if (any(valid)) mean(result$recurrent_dom_match[valid]) else NA_real_,
      random_seed = seed
    ),
    class_probability = tibble(
      cohort = cohort_name,
      gene_id = object$gene_id,
      bootstrap_architecture = null_architecture_levels,
      probability = as.numeric(class_prob),
      n_valid_replicates = n_valid
    )
  )
}

run_bootstrap_cohort <- function(cohort_name, included_patients) {
  observed_lookup <- observed_gene %>%
    filter(cohort == cohort_name) %>%
    arrange(gene_id)
  target_ids <- observed_lookup$gene_id[observed_lookup$quantifiable]
  cache_path <- file.path(m04$cache, paste0("bootstrap_", cohort_name, ".rds"))
  configuration_hash <- digest::digest(
    list(
      input_fingerprint = cache_input_fingerprint,
      cohort = cohort_name,
      included_patients = included_patients,
      bootstrap = cfg$bootstrap,
      thresholds = cache_thresholds,
      genes = target_ids
    ), algo = "sha256"
  )
  legacy_configuration_hashes <- vapply(
    legacy_input_fingerprints_v0_5,
    function(legacy_fingerprint) {
      digest::digest(
        list(
          input_fingerprint = legacy_fingerprint,
          cohort = cohort_name,
          included_patients = included_patients,
          bootstrap = cfg$bootstrap,
          thresholds = legacy_thresholds_v0_5,
          genes = target_ids
        ),
        algo = "sha256"
      )
    },
    character(1)
  )
  state <- list(configuration_hash = configuration_hash, completed = list())
  if (cfg$resume_cached_results && file.exists(cache_path)) {
    candidate <- readRDS(cache_path)
    if (candidate$configuration_hash %in%
        c(configuration_hash, legacy_configuration_hashes)) {
      state <- candidate
      state$configuration_hash <- configuration_hash
      log_msg("Loaded compatible bootstrap cache for ", cohort_name, ".")
    }
  }
  pending <- setdiff(target_ids, names(state$completed))
  if (length(pending)) {
    log_msg("Bootstrap ", cohort_name, ": ", length(pending), " genes pending.")
  }
  for (position in seq_along(pending)) {
    gene_id <- pending[[position]]
    observed_row <- observed_lookup[match(gene_id, observed_lookup$gene_id), , drop = FALSE]
    state$completed[[gene_id]] <- bootstrap_gene(
      gene_objects[[gene_id]], observed_row, included_patients, cohort_name
    )
    if (position %% cfg$checkpoint_every_genes == 0L || position == length(pending)) {
      atomic_save_rds(state, cache_path)
      log_msg(
        "Bootstrap ", cohort_name, ": ", length(state$completed), " / ",
        length(target_ids), " genes cached."
      )
    }
  }
  list(
    summary = rbindlist(lapply(state$completed, `[[`, "summary"), fill = TRUE) %>% as_tibble(),
    class_probability = rbindlist(lapply(state$completed, `[[`, "class_probability"), fill = TRUE) %>% as_tibble()
  )
}
ilr_basis <- function(k) {
  if (!is.finite(k) || k < 2L) stop("ILR basis requires at least two parts.")
  basis <- matrix(0, nrow = k, ncol = k - 1L)
  for (column in seq_len(k - 1L)) {
    basis[seq_len(column), column] <- 1 / sqrt(column * (column + 1))
    basis[column + 1L, column] <- -column / sqrt(column * (column + 1))
  }
  stopifnot(
    max(abs(colSums(basis))) < 1e-10,
    max(abs(crossprod(basis) - diag(k - 1L))) < 1e-10
  )
  basis
}

ilr_transform <- function(composition, basis) {
  composition <- as.matrix(composition)
  if (any(!is.finite(composition)) || any(composition <= 0)) {
    stop("ILR transformation requires finite, strictly positive compositions.")
  }
  log(composition) %*% basis
}

ilr_inverse <- function(coordinates, basis) {
  coordinates <- as.matrix(coordinates)
  log_parts <- coordinates %*% t(basis)
  log_parts <- sweep(log_parts, 1L, apply(log_parts, 1L, max), "-")
  parts <- exp(log_parts)
  parts / rowSums(parts)
}

smoothed_composition <- function(counts, pseudocount = 0.5) {
  counts <- as.matrix(counts)
  totals <- rowSums(counts)
  composition <- sweep(
    counts + pseudocount,
    1L,
    totals + ncol(counts) * pseudocount,
    "/"
  )
  list(composition = composition, totals = totals)
}

multinomial_composition_rows <- function(probability, totals) {
  probability <- as.matrix(probability)
  totals <- as.integer(round(totals))
  output <- matrix(0, nrow = nrow(probability), ncol = ncol(probability))
  for (row_index in seq_len(nrow(probability))) {
    count <- as.numeric(rmultinom(1, totals[[row_index]], probability[row_index, ]))
    output[row_index, ] <- count / max(1, totals[[row_index]])
  }
  output
}

covariance_split_qc <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  n <- nrow(coordinates)
  dimension <- ncol(coordinates)
  unavailable <- function(reason) list(available = FALSE, stable = FALSE, reason = reason)
  if (n < 4L || dimension < 1L || any(!is.finite(coordinates))) {
    return(unavailable("invalid_covariance_draws"))
  }
  half <- floor(n / 2L)
  first <- matrix(cov(coordinates[seq_len(half), , drop = FALSE]),
                  nrow = dimension, ncol = dimension)
  second <- matrix(cov(coordinates[(half + 1L):n, , drop = FALSE]),
                   nrow = dimension, ncol = dimension)
  full <- matrix(cov(coordinates), nrow = dimension, ncol = dimension)
  full <- 0.5 * (full + t(full))
  trace_full <- sum(diag(full))
  trace_discrepancy <- abs(sum(diag(first)) - sum(diag(second))) /
    max(abs(trace_full), .Machine$double.eps)
  frobenius_scale <- sqrt(sum(full^2))
  frobenius_discrepancy <- sqrt(sum((first - second)^2)) /
    max(frobenius_scale, .Machine$double.eps)
  eigenvalues <- tryCatch(
    eigen(full, symmetric = TRUE, only.values = TRUE)$values,
    error = function(error) numeric()
  )
  if (length(eigenvalues) != dimension || any(!is.finite(eigenvalues))) {
    return(unavailable("covariance_eigendecomposition_failed"))
  }
  machine_tolerance <- 100 * .Machine$double.eps * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -machine_tolerance) {
    return(unavailable("covariance_materially_non_psd"))
  }
  if (min(eigenvalues) <= 0) {
    return(unavailable("covariance_singular"))
  }
  stable <- trace_discrepancy <= cfg$gene_null$covariance_trace_tolerance &&
    frobenius_discrepancy <= cfg$gene_null$covariance_frobenius_tolerance
  list(
    available = TRUE,
    stable = stable,
    reason = if (stable) "available" else "split_half_unresolved",
    covariance = full,
    trace = trace_full,
    trace_relative_discrepancy = trace_discrepancy,
    frobenius_relative_discrepancy = frobenius_discrepancy,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    condition_number = max(eigenvalues) / min(eigenvalues)
  )
}

finite_count_covariance <- function(counts, basis, seed) {
  counts <- as.numeric(counts)
  total <- sum(counts)
  k <- length(counts)
  invalid <- function(reason, replicates = 0L, qc = NULL) {
    list(available = FALSE, reason = reason, covariance = matrix(numeric()),
         replicates = as.integer(replicates), qc = qc)
  }
  if (k < 2L || any(!is.finite(counts)) || any(counts < 0) ||
      total < 1 || any(counts != round(counts))) {
    return(invalid("finite_count_covariance_input_invalid"))
  }
  probability <- (counts + cfg$gene_null$prior_strength) /
    (total + k * cfg$gene_null$prior_strength)
  result <- with_seed(seed, {
    coordinate_draws <- matrix(numeric(), nrow = 0L, ncol = ncol(basis))
    target <- cfg$gene_null$covariance_initial_replicates
    repeat {
      needed <- target - nrow(coordinate_draws)
      if (needed > 0L) {
        count_draw <- t(rmultinom(needed, size = as.integer(total), prob = probability))
        smoothed <- sweep(
          count_draw + cfg$gene_null$prior_strength,
          1L,
          total + k * cfg$gene_null$prior_strength,
          "/"
        )
        coordinate_draws <- rbind(coordinate_draws, ilr_transform(smoothed, basis))
      }
      qc <- covariance_split_qc(coordinate_draws)
      if (!isTRUE(qc$available) || isTRUE(qc$stable) ||
          target >= cfg$gene_null$covariance_maximum_replicates) break
      target <- min(2L * target, cfg$gene_null$covariance_maximum_replicates)
    }
    list(coordinates = coordinate_draws, qc = qc)
  })
  qc <- result$qc
  n_draw <- nrow(result$coordinates)
  if (!isTRUE(qc$available)) return(invalid(qc$reason, n_draw, qc))
  if (!isTRUE(qc$stable)) {
    return(invalid("finite_count_covariance_split_half_unresolved", n_draw, qc))
  }
  list(
    available = TRUE, reason = "available", covariance = qc$covariance,
    replicates = as.integer(n_draw), qc = qc, estimated_probability = probability
  )
}

valid_reml_inputs <- function(difference, sampling_covariance) {
  difference <- as.matrix(difference)
  n <- nrow(difference)
  dimension <- ncol(difference)
  if (n < 2L || dimension < 1L || any(!is.finite(difference)) ||
      length(sampling_covariance) != n) return(FALSE)
  all(vapply(sampling_covariance, function(value) {
    value <- as.matrix(value)
    identical(dim(value), c(dimension, dimension)) &&
      all(is.finite(value)) && max(abs(value - t(value))) <= 1e-10 &&
      !is.null(tryCatch(chol(value), error = function(error) NULL))
  }, logical(1)))
}

reml_evaluate <- function(variance, difference, sampling_covariance, information = FALSE) {
  variance <- as.numeric(variance)[[1]]
  difference <- as.matrix(difference)
  n <- nrow(difference)
  dimension <- ncol(difference)
  invalid <- function() list(objective = Inf, information = NA_real_, mean = numeric())
  if (!is.finite(variance) || variance < 0 || n < 2L || dimension < 1L ||
      length(sampling_covariance) != n) return(invalid())
  identity <- diag(dimension)
  precision <- vector("list", n)
  log_determinant <- numeric(n)
  for (index in seq_len(n)) {
    covariance <- 0.5 * (sampling_covariance[[index]] +
                           t(sampling_covariance[[index]])) + variance * identity
    factor <- tryCatch(chol(covariance), error = function(error) NULL)
    if (is.null(factor)) return(invalid())
    precision[[index]] <- chol2inv(factor)
    log_determinant[[index]] <- 2 * sum(log(diag(factor)))
  }
  summed_precision <- Reduce("+", precision)
  summed_factor <- tryCatch(chol(summed_precision), error = function(error) NULL)
  if (is.null(summed_factor)) return(invalid())
  right_hand_side <- Reduce("+", lapply(seq_len(n), function(index) {
    precision[[index]] %*% difference[index, ]
  }))
  mean_hat <- as.numeric(chol2inv(summed_factor) %*% right_hand_side)
  quadratic <- sum(vapply(seq_len(n), function(index) {
    residual <- difference[index, ] - mean_hat
    as.numeric(crossprod(residual, precision[[index]] %*% residual))
  }, numeric(1)))
  objective <- sum(log_determinant) +
    2 * sum(log(diag(summed_factor))) + quadratic
  reml_information <- NA_real_
  if (isTRUE(information)) {
    summed_inverse <- chol2inv(summed_factor)
    squared <- lapply(precision, function(value) value %*% value)
    cubed <- lapply(seq_len(n), function(index) squared[[index]] %*% precision[[index]])
    sum_squared <- Reduce("+", squared)
    sum_cubed <- Reduce("+", cubed)
    trace_projection_squared <-
      sum(vapply(squared, function(value) sum(diag(value)), numeric(1))) -
      2 * sum(diag(summed_inverse %*% sum_cubed)) +
      sum(diag(summed_inverse %*% sum_squared %*%
                 summed_inverse %*% sum_squared))
    if (is.finite(trace_projection_squared)) {
      reml_information <- 0.5 * max(0, trace_projection_squared)
    }
  }
  list(objective = as.numeric(objective), information = reml_information, mean = mean_hat)
}

fit_additive_reml <- function(difference, sampling_covariance) {
  invalid <- function(reason) list(
    estimate = NA_real_, objective_hat = NA_real_, objective_zero = NA_real_,
    upper = NA_real_, boundary = NA, converged = FALSE,
    expansion_count = 0L, curvature = NA_real_, information = NA_real_,
    mean = numeric(), reason = reason
  )
  if (!valid_reml_inputs(difference, sampling_covariance)) {
    return(invalid("invalid_reml_inputs"))
  }
  difference <- as.matrix(difference)
  scale <- mean(apply(difference, 2L, var))
  if (!is.finite(scale) || scale <= 0) scale <- 1e-4
  upper <- max(1e-4, 8 * scale)
  objective <- function(value) reml_evaluate(
    value, difference, sampling_covariance, FALSE
  )$objective
  objective_zero <- objective(0)
  expansion <- 0L
  repeat {
    grid <- c(0, exp(seq(
      log(max(upper * 1e-8, .Machine$double.eps)),
      log(upper), length.out = 65L
    )))
    values <- vapply(grid, objective, numeric(1))
    index <- which.min(values)
    if (index < length(grid) || expansion >= 12L) break
    upper <- upper * 4
    expansion <- expansion + 1L
  }
  if (!length(index) || !is.finite(values[[index]]) || index == length(grid)) {
    output <- invalid("reml_upper_bracket_unresolved")
    output$objective_zero <- objective_zero
    output$upper <- upper
    output$expansion_count <- expansion
    return(output)
  }
  if (index == 1L) {
    estimate <- 0
    minimum <- values[[1]]
  } else {
    optimized <- tryCatch(
      optimize(
        objective,
        interval = c(grid[max(1L, index - 1L)], grid[min(length(grid), index + 1L)]),
        tol = .Machine$double.eps^0.5
      ),
      error = function(error) NULL
    )
    if (is.null(optimized)) return(invalid("reml_optimizer_failed"))
    estimate <- optimized$minimum
    minimum <- optimized$objective
    if (is.finite(objective_zero) && objective_zero <= minimum + 1e-10) {
      estimate <- 0
      minimum <- objective_zero
    }
  }
  evaluation <- reml_evaluate(estimate, difference, sampling_covariance, TRUE)
  step <- max(1e-7, 1e-3 * max(estimate, scale))
  curvature <- if (estimate > step) {
    (objective(estimate + step) - 2 * minimum + objective(estimate - step)) / step^2
  } else {
    (objective(estimate + 2 * step) - 2 * objective(estimate + step) + minimum) / step^2
  }
  converged <- is.finite(minimum) && is.finite(evaluation$information) &&
    evaluation$information > 0 && is.finite(curvature)
  list(
    estimate = as.numeric(estimate), objective_hat = as.numeric(minimum),
    objective_zero = as.numeric(objective_zero), upper = as.numeric(upper),
    boundary = isTRUE(estimate == 0), converged = converged,
    expansion_count = as.integer(expansion), curvature = as.numeric(curvature),
    information = as.numeric(evaluation$information), mean = evaluation$mean,
    reason = if (converged) "available" else "reml_fit_unavailable"
  )
}

profile_variance_mixture <- function(difference, sampling_covariance) {
  fit <- fit_additive_reml(difference, sampling_covariance)
  invalid <- function(reason, fit, grid = numeric(), weight = numeric(),
                      upper = NA_real_, expansion_count = 0L) {
    list(
      available = FALSE, reason = reason, fit = fit, grid = grid, weight = weight,
      mean = NA_real_, sd = NA_real_, q05 = NA_real_, q50 = NA_real_, q95 = NA_real_,
      upper = upper, tail_resolved = FALSE,
      expansion_count = as.integer(expansion_count)
    )
  }
  if (!isTRUE(fit$converged) || !is.finite(fit$estimate)) {
    return(invalid(fit$reason, fit))
  }
  profile_scale <- 8 * max(fit$estimate, 1e-6)
  upper <- max(fit$upper, profile_scale)
  positive_minimum <- max(profile_scale * 1e-9, .Machine$double.eps)
  expansion <- 0L
  tail_resolved <- FALSE
  repeat {
    grid <- unique(c(0, exp(seq(
      log(positive_minimum), log(upper),
      length.out = cfg$gene_null$profile_grid_points - 1L
    ))))
    evaluated <- lapply(grid, function(value) {
      reml_evaluate(value, difference, sampling_covariance, TRUE)
    })
    objective <- vapply(evaluated, `[[`, numeric(1), "objective")
    information <- vapply(evaluated, `[[`, numeric(1), "information")
    valid <- is.finite(objective) & is.finite(information) & information > 0
    log_density <- rep(-Inf, length(grid))
    log_density[valid] <- -0.5 * objective[valid] + 0.5 * log(information[valid])
    finite <- is.finite(log_density)
    if (any(finite)) {
      maximum <- max(log_density[finite])
      mode_index <- which.max(replace(log_density, !finite, -Inf))
      tail_resolved <- mode_index < length(grid) &&
        is.finite(log_density[[length(grid)]]) &&
        log_density[[length(grid)]] <= maximum - cfg$gene_null$profile_tail_log_drop
    }
    if (tail_resolved || expansion >= 12L) break
    upper <- upper * 4
    expansion <- expansion + 1L
  }
  if (!any(is.finite(log_density))) {
    return(invalid("profile_density_unavailable", fit, grid, numeric(), upper, expansion))
  }
  if (!tail_resolved) {
    return(invalid("profile_tail_unresolved", fit, grid, numeric(), upper, expansion))
  }
  cell_width <- numeric(length(grid))
  cell_width[[1]] <- (grid[[2]] - grid[[1]]) / 2
  cell_width[[length(grid)]] <-
    (grid[[length(grid)]] - grid[[length(grid) - 1L]]) / 2
  if (length(grid) > 2L) {
    cell_width[2L:(length(grid) - 1L)] <-
      (grid[3L:length(grid)] - grid[1L:(length(grid) - 2L)]) / 2
  }
  log_mass <- log_density + log(cell_width)
  maximum_mass <- max(log_mass[is.finite(log_mass)])
  weight <- exp(log_mass - maximum_mass)
  weight[!is.finite(weight)] <- 0
  weight <- weight / sum(weight)
  discrete_quantile <- function(probability) grid[which(cumsum(weight) >= probability)[[1]]]
  mixture_mean <- sum(grid * weight)
  list(
    available = TRUE, reason = "available", fit = fit, grid = grid, weight = weight,
    mean = mixture_mean,
    sd = sqrt(sum((grid - mixture_mean)^2 * weight)),
    q05 = discrete_quantile(0.05), q50 = discrete_quantile(0.50),
    q95 = discrete_quantile(0.95), upper = upper, tail_resolved = TRUE,
    expansion_count = as.integer(expansion)
  )
}

simulate_joint_scores <- function(fit, variance, replicates, seed) {
  variance <- as.numeric(variance)
  replicates <- as.integer(replicates)
  if (length(variance) != replicates || any(!is.finite(variance)) || any(variance < 0)) {
    return(numeric())
  }
  with_seed(seed, {
    scores <- numeric(replicates)
    starts <- seq.int(1L, replicates, by = 128L)
    for (start in starts) {
      finish <- min(replicates, start + 127L)
      batch_size <- finish - start + 1L
      batch_variance <- variance[start:finish]
      n_units <- nrow(fit$midpoint)
      row_index <- rep(seq_len(n_units), times = batch_size)
      midpoint <- fit$midpoint[row_index, , drop = FALSE]
      standard_deviation <- rep(
        rep(sqrt(batch_variance), each = n_units),
        times = fit$dimension
      )
      difference <- matrix(
        rnorm(batch_size * length(fit$midpoint), sd = standard_deviation),
        nrow = batch_size * n_units, ncol = fit$dimension
      )
      primary_probability <- ilr_inverse(midpoint - 0.5 * difference, fit$basis)
      recurrent_probability <- ilr_inverse(midpoint + 0.5 * difference, fit$basis)
      primary <- multinomial_composition_rows(
        primary_probability, rep(fit$primary_totals, times = batch_size)
      )
      recurrent <- multinomial_composition_rows(
        recurrent_probability, rep(fit$recurrent_totals, times = batch_size)
      )
      primary_mean <- vapply(seq_len(ncol(primary)), function(part) {
        rowMeans(matrix(primary[, part], nrow = batch_size,
                        ncol = n_units, byrow = TRUE))
      }, numeric(batch_size))
      recurrent_mean <- vapply(seq_len(ncol(recurrent)), function(part) {
        rowMeans(matrix(recurrent[, part], nrow = batch_size,
                        ncol = n_units, byrow = TRUE))
      }, numeric(batch_size))
      primary_mean <- matrix(primary_mean, nrow = batch_size, ncol = ncol(primary))
      recurrent_mean <- matrix(recurrent_mean, nrow = batch_size, ncol = ncol(recurrent))
      scores[start:finish] <- joint_score(
        jsd_rows(primary_mean, recurrent_mean),
        apply(abs(recurrent_mean - primary_mean), 1L, max)
      )
    }
    scores
  })
}

monte_carlo_p <- function(scores, observed_score) {
  if (!length(scores) || any(!is.finite(scores)) || !is.finite(observed_score)) return(NULL)
  exceedances <- sum(scores >= observed_score - 1e-15)
  p_value <- (1 + exceedances) / (length(scores) + 1)
  quantiles <- quantile(scores, c(0.50, 0.90, 0.95, 0.99), names = FALSE, type = cfg$bootstrap$quantile_type)
  list(
    p = p_value, exceedances = as.integer(exceedances),
    replicates = as.integer(length(scores)), minimum_p = 1 / (length(scores) + 1),
    monte_carlo_se = sqrt(p_value * (1 - p_value) / (length(scores) + 1)),
    mean = mean(scores), sd = sd(scores), q50 = quantiles[[1]], q90 = quantiles[[2]],
    q95 = quantiles[[3]], q99 = quantiles[[4]]
  )
}

gpd_survival <- function(excess, scale, shape) {
  excess <- as.numeric(excess)
  if (!is.finite(scale) || scale <= 0 || !is.finite(shape) ||
      any(!is.finite(excess)) || any(excess < 0)) {
    return(rep(NA_real_, length(excess)))
  }
  if (abs(shape) < 1e-7) {
    return(exp(-excess / scale))
  }
  support <- 1 + shape * excess / scale
  output <- rep(0, length(excess))
  inside <- support > 0
  output[inside] <- support[inside]^(-1 / shape)
  pmin(1, pmax(0, output))
}

gpd_negative_loglikelihood <- function(parameter, excess) {
  scale <- exp(parameter[[1]])
  shape <- parameter[[2]]
  support <- 1 + shape * excess / scale
  if (!is.finite(scale) || scale <= 0 || !is.finite(shape) ||
      any(!is.finite(support)) || any(support <= 0)) {
    return(1e100)
  }
  if (abs(shape) < 1e-7) {
    return(length(excess) * log(scale) + sum(excess) / scale)
  }
  value <- length(excess) * log(scale) +
    (1 / shape + 1) * sum(log(support))
  if (is.finite(value)) value else 1e100
}

fit_gpd_threshold <- function(scores, observed_score, tail_fraction, branch) {
  unavailable <- function(reason, threshold = NA_real_, tail_size = 0L) {
    tibble(
      gpd_branch = branch,
      gpd_tail_fraction = tail_fraction,
      gpd_threshold = threshold,
      gpd_tail_size = as.integer(tail_size),
      gpd_scale = NA_real_,
      gpd_shape = NA_real_,
      gpd_tail_probability = NA_real_,
      gpd_conditional_survival = NA_real_,
      gpd_tail_p = NA_real_,
      gpd_ks_distance = NA_real_,
      gpd_ks_limit = NA_real_,
      gpd_gof_pass = FALSE,
      gpd_optimizer_converged = FALSE,
      gpd_available = FALSE,
      gpd_reason = reason
    )
  }
  if (!length(scores) || any(!is.finite(scores)) || !is.finite(observed_score)) {
    return(unavailable("invalid_scores_or_observed_statistic"))
  }
  threshold <- as.numeric(quantile(
    scores, probs = 1 - tail_fraction, names = FALSE, type = cfg$bootstrap$quantile_type
  ))
  excess <- scores[scores > threshold] - threshold
  excess <- excess[is.finite(excess) & excess > 0]
  if (observed_score <= threshold) {
    return(unavailable("observed_statistic_not_in_fitted_tail", threshold, length(excess)))
  }
  if (length(excess) < cfg$gpd$minimum_tail_size ||
      length(unique(signif(excess, 12))) < 10L) {
    return(unavailable("insufficient_or_degenerate_tail", threshold, length(excess)))
  }

  scale_start <- max(mean(excess), stats::sd(excess), .Machine$double.eps^0.25)
  upper_scale <- max(scale_start * 1e4, max(excess) * 1e3, 1)
  starts <- c(-0.25, 0, 0.25)
  fits <- lapply(starts, function(shape_start) {
    tryCatch(
      optim(
        par = c(log(scale_start), shape_start),
        fn = gpd_negative_loglikelihood,
        excess = excess,
        method = "L-BFGS-B",
        lower = c(log(.Machine$double.eps^0.25), -0.49),
        upper = c(log(upper_scale), 2)
      ),
      error = function(error) NULL
    )
  })
  valid_fit <- vapply(fits, function(value) {
    !is.null(value) && is.finite(value$value) && value$value < 1e99 &&
      all(is.finite(value$par))
  }, logical(1))
  if (!any(valid_fit)) {
    return(unavailable("GPD_optimizer_failed", threshold, length(excess)))
  }
  fits <- fits[valid_fit]
  selected <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  scale <- exp(selected$par[[1]])
  shape <- selected$par[[2]]
  conditional_survival <- gpd_survival(observed_score - threshold, scale, shape)
  if (length(conditional_survival) != 1L || !is.finite(conditional_survival) ||
      conditional_survival <= 0 || conditional_survival > 1) {
    return(unavailable("observed_statistic_outside_fitted_GPD_support", threshold,
                       length(excess)))
  }
  tail_probability <- (length(excess) + 1) / (length(scores) + 1)
  tail_p <- min(1, max(.Machine$double.xmin,
                       tail_probability * conditional_survival))
  ordered_excess <- sort(excess)
  fitted_cdf <- 1 - gpd_survival(ordered_excess, scale, shape)
  empirical_cdf <- (seq_along(ordered_excess) - 0.5) / length(ordered_excess)
  ks_distance <- max(abs(empirical_cdf - fitted_cdf))
  ks_limit <- cfg$gpd$ks_multiplier / sqrt(length(ordered_excess))

  tibble(
    gpd_branch = branch,
    gpd_tail_fraction = tail_fraction,
    gpd_threshold = threshold,
    gpd_tail_size = as.integer(length(excess)),
    gpd_scale = scale,
    gpd_shape = shape,
    gpd_tail_probability = tail_probability,
    gpd_conditional_survival = conditional_survival,
    gpd_tail_p = tail_p,
    gpd_ks_distance = ks_distance,
    gpd_ks_limit = ks_limit,
    gpd_gof_pass = is.finite(ks_distance) && ks_distance <= ks_limit,
    gpd_optimizer_converged = selected$convergence == 0L,
    gpd_available = selected$convergence == 0L,
    gpd_reason = if (selected$convergence == 0L) "available" else
      "GPD_optimizer_did_not_converge"
  )
}

estimate_gpd_tail <- function(scores, observed_score, empirical_result, branch) {
  empty_summary <- list(
    raw_p = NA_real_, hybrid_p = empirical_result$p,
    hybrid_method = "empirical_fallback", available = FALSE,
    threshold_stable = FALSE, selected_fraction = NA_real_,
    shape = NA_real_, scale = NA_real_, ks_distance = NA_real_,
    reason = "GPD_not_required_or_unavailable"
  )
  if (!cfg$gpd$enabled) {
    empty_summary$hybrid_method <- "GPD_disabled_empirical"
    empty_summary$reason <- "GPD_disabled"
    return(list(summary = empty_summary, diagnostics = tibble()))
  }
  if (empirical_result$exceedances >= cfg$gpd$empirical_exceedance_trigger) {
    empty_summary$hybrid_method <- "empirical_exceedances_sufficient"
    empty_summary$reason <- "empirical_resolution_sufficient"
    return(list(summary = empty_summary, diagnostics = tibble()))
  }

  diagnostics <- map_dfr(cfg$gpd$tail_fractions, function(tail_fraction) {
    fit_gpd_threshold(scores, observed_score, tail_fraction, branch)
  })
  eligible <- diagnostics %>%
    filter(gpd_available, gpd_gof_pass, is.finite(gpd_tail_p), gpd_tail_p > 0)
  if (!nrow(eligible)) {
    empty_summary$reason <- "no_GPD_threshold_passed_fit_and_GOF_checks"
    return(list(summary = empty_summary, diagnostics = diagnostics))
  }
  log10_range <- diff(range(log10(eligible$gpd_tail_p)))
  threshold_stable <- nrow(eligible) >= 2L && is.finite(log10_range) &&
    log10_range <= cfg$gpd$log10_stability_tolerance
  selected_index <- which.min(abs(eligible$gpd_tail_fraction - 0.05))
  selected <- eligible[selected_index, , drop = FALSE]
  summary <- list(
    raw_p = selected$gpd_tail_p[[1]],
    hybrid_p = if (threshold_stable) selected$gpd_tail_p[[1]] else empirical_result$p,
    hybrid_method = if (threshold_stable) "GPD_tail" else
      "empirical_fallback_unstable_GPD",
    available = TRUE,
    threshold_stable = threshold_stable,
    selected_fraction = selected$gpd_tail_fraction[[1]],
    shape = selected$gpd_shape[[1]],
    scale = selected$gpd_scale[[1]],
    ks_distance = selected$gpd_ks_distance[[1]],
    reason = if (threshold_stable) "available_and_threshold_stable" else
      "GPD_threshold_sensitivity_exceeded_tolerance"
  )
  list(summary = summary, diagnostics = diagnostics)
}
