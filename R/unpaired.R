# Independent-sample methods. These functions never construct subject pairs.
# Their results use explicit unpaired method labels rather than paired classes.
rna_unpaired_settings <- function(cfg) {
  out <- utils::modifyList(list(min_samples_per_condition = 3L,
    minimum_evaluable_per_condition = NULL, n_bootstrap = 1000L,
    n_permutations = 1000L, minimum_prevalence_difference = .5), cfg$unpaired %||% list())
  out$minimum_evaluable_per_condition <- out$minimum_evaluable_per_condition %||% out$min_samples_per_condition
  for (k in c("min_samples_per_condition", "minimum_evaluable_per_condition", "n_bootstrap", "n_permutations"))
    if (length(out[[k]]) != 1L || !is.numeric(out[[k]]) || !is.finite(out[[k]]) ||
        out[[k]] != round(out[[k]]) || out[[k]] < if (grepl("samples|evaluable", k)) 3 else 2)
      stop("Invalid unpaired setting: ", k, call. = FALSE)
  if (out$minimum_evaluable_per_condition < out$min_samples_per_condition)
    stop("unpaired.minimum_evaluable_per_condition must be at least min_samples_per_condition.", call. = FALSE)
  if (length(out$minimum_prevalence_difference) != 1L ||
      !is.finite(out$minimum_prevalence_difference) || out$minimum_prevalence_difference <= 0 ||
      out$minimum_prevalence_difference > 1) stop("Invalid unpaired.minimum_prevalence_difference.")
  out
}

rna_unpaired_context <- function(inputs, cfg, cohorts) {
  x <- as.matrix(inputs$counts_matrix); a <- as.data.frame(inputs$isoform_annotation)
  s <- as.data.frame(inputs$sample_manifest)
  if (!identical(rownames(x), as.character(a$isoform_harmonized_id)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      any(!is.finite(x) | x < 0 | x != round(x))) stop("Invalid unpaired integer count matrix.")
  if (!all(c("sample", "patient", "timepoint") %in% names(s)) || anyDuplicated(s$sample) ||
      !setequal(s$sample, colnames(x))) stop("Unpaired metadata/count identities do not match.")
  s <- s[match(colnames(x), s$sample), , drop = FALSE]
  if (any(!s$timepoint %in% c("primary", "recurrent")) || anyDuplicated(s$patient))
    stop("Unpaired mode requires independent subjects and the two configured conditions.")
  u <- rna_unpaired_settings(cfg)
  if (!length(cohorts) || is.null(names(cohorts)) || anyDuplicated(names(cohorts)) ||
      !setequal(cohorts[[1]], s$sample)) stop("First unpaired cohort must contain all imported samples.")
  for (name in names(cohorts)) {
    ids <- as.character(cohorts[[name]])
    if (anyDuplicated(ids) || any(!ids %in% s$sample)) stop("Unknown/duplicate samples in unpaired cohort ", name)
    n <- table(factor(s$timepoint[s$sample %in% ids], levels = c("primary", "recurrent")))
    if (any(n < u$min_samples_per_condition))
      stop("Unpaired cohort ", name, " has fewer than ", u$min_samples_per_condition, " independent samples per condition.")
  }
  mapped <- !is.na(a$associated_gene) & nzchar(a$associated_gene) & a$associated_gene != "unknown"
  if ("gene_mapped" %in% names(a)) mapped <- mapped & a$gene_mapped %in% TRUE
  if (!any(mapped)) stop("Unpaired analyses require at least one mapped gene.")
  comparable <- if ("cross_sample_comparable" %in% names(a)) a$cross_sample_comparable %in% TRUE else startsWith(a$isoform_harmonized_id, "JUNC_")
  gene <- as.character(a$associated_gene)
  totals <- rowsum(x[mapped, , drop = FALSE], gene[mapped], reorder = TRUE)
  list(counts = x, annotation = a, samples = s, mapped = mapped, comparable = comparable,
    gene = gene, gene_totals = totals, cohorts = cohorts, settings = u)
}

rna_unpaired_cache <- function(path, key, compute, cfg) {
  if (!isTRUE(cfg$execution$force_refit) && !identical(cfg$execution$resume, FALSE) && file.exists(path)) {
    saved <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.list(saved) && identical(saved$key, key) && isTRUE(saved$complete)) {
      rna_log("Reusing unpaired checkpoint: ", basename(path)); return(saved$result)
    }
  }
  value <- compute()
  rna_atomic_rds(list(key = key, complete = TRUE, result = value), path)
  value
}

rna_unpaired_key <- function(ctx, cfg, stage, ids, extra = NULL) rna_hash(list(
  schema = "independent-samples-v1", stage = stage, ids = ids,
  count_rows = rownames(ctx$counts), counts = ctx$counts,
  metadata = ctx$samples, annotation = ctx$annotation,
  settings = rna_scientific_config(cfg), extra = extra,
  runtime = list(R = R.version.string, RNG = dtu_rng_spec(NULL),
    packages = vapply(c("vegan", "digest"), function(pkg)
      tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "Not installed"), character(1))),
  code = list(settings = body(rna_unpaired_settings), context = body(rna_unpaired_context),
    seed = body(rna_stage_seed), with_rng = body(dtu_with_rng), rng_spec = body(dtu_rng_spec),
    prevalence = body(rna_unpaired_prevalence),
    rarefy = body(rna_unpaired_rarefy), geometry = body(rna_unpaired_geometry),
    permutation = body(rna_unpaired_permutation), gain_loss = body(run_unpaired_gain_loss),
    remodeling = body(run_unpaired_remodeling))))

rna_unpaired_rarefy <- function(counts, target) {
  if (length(target) != 1L || !is.finite(target) || target < 1 ||
      target > .Machine$integer.max || target != round(target))
    stop("Unpaired rarefaction target must be a positive representable integer.")
  # The production path is the same multivariate hypergeometric library
  # sampler used by the paired branch. The base fallback is bounded to tiny
  # fixtures rather than millions of interpreted per-isoform draws.
  if (requireNamespace("vegan", quietly = TRUE)) {
    out <- t(vegan::rrarefy(t(counts), sample = as.integer(target)))
    dimnames(out) <- dimnames(counts)
    if (any(colSums(out) != target) || any(out > counts)) stop("Invalid unpaired common-depth draw.")
    return(out)
  }
  if (nrow(counts) > 10000L) stop("Production unpaired rarefaction requires vegan.")
  # Multivariate hypergeometric sampling without expanding individual reads.
  out <- matrix(0, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  for (j in seq_len(ncol(counts))) {
    available <- sum(counts[, j]); remaining <- target
    if (available < target) stop("Common depth exceeds a source library.")
    for (i in which(counts[, j] > 0)) {
      if (remaining <= 0) break
      n <- counts[i, j]
      draw <- if (remaining == available) n else stats::rhyper(1, n, available - n, remaining)
      out[i, j] <- draw; remaining <- remaining - draw; available <- available - n
    }
    if (sum(out[, j]) != target) stop("Unpaired rarefaction did not conserve the target depth.")
  }
  out
}

rna_unpaired_prevalence <- function(counts, ctx, cfg, ids) {
  a <- ctx$annotation; s <- ctx$samples
  keep <- ctx$mapped & ctx$comparable
  categories <- as.character(unlist(cfg$gain_loss$keep_structural_categories))
  if (length(categories)) keep <- keep & a$structural_category %in% categories
  ii <- which(keep); jj <- match(ids, s$sample)
  if (!length(ii)) stop("No comparable isoforms meet the configured prevalence category criteria.")
  gt <- rowsum(counts[ctx$mapped, , drop = FALSE], ctx$gene[ctx$mapped], reorder = TRUE)
  denom <- gt[match(ctx$gene[ii], rownames(gt)), jj, drop = FALSE]
  y <- counts[ii, jj, drop = FALSE]
  evaluable <- denom >= cfg$gain_loss$gene_evaluable_min_fl
  present <- y >= cfg$gain_loss$isoform_presence_min_fl & y / pmax(denom, 1) >= cfg$gain_loss$isoform_presence_min_fraction
  ref <- s$timepoint[jj] == "primary"; cmp <- !ref
  nr <- rowSums(evaluable[, ref, drop = FALSE]); nc <- rowSums(evaluable[, cmp, drop = FALSE])
  pr <- rowSums(present[, ref, drop = FALSE] & evaluable[, ref, drop = FALSE])
  pc <- rowSums(present[, cmp, drop = FALSE] & evaluable[, cmp, drop = FALSE])
  er <- ifelse(nr > 0, pr / nr, NA_real_); ec <- ifelse(nc > 0, pc / nc, NA_real_)
  valid <- nr >= ctx$settings$minimum_evaluable_per_condition & nc >= ctx$settings$minimum_evaluable_per_condition
  delta <- ec - er; threshold <- ctx$settings$minimum_prevalence_difference
  data.frame(isoform_harmonized_id = a$isoform_harmonized_id[ii], gene_id = ctx$gene[ii],
    structural_category = a$structural_category[ii], reference_samples = sum(ref), comparison_samples = sum(cmp),
    reference_evaluable = nr, comparison_evaluable = nc, reference_present = pr, comparison_present = pc,
    reference_prevalence = er, comparison_prevalence = ec, prevalence_difference = delta,
    quantifiable = valid, prevalence_class = ifelse(!valid, "not_quantifiable",
      ifelse(delta >= threshold, "comparison_enriched", ifelse(delta <= -threshold, "comparison_depleted", "no_large_prevalence_difference"))),
    design = "unpaired", method = "independent_sample_detection_prevalence", stringsAsFactors = FALSE)
}

run_unpaired_gain_loss <- function(inputs, cfg, cohorts, stage_dir) {
  ctx <- rna_unpaired_context(inputs, cfg, cohorts); dirs <- rna_dirs(stage_dir)
  depth <- min(colSums(ctx$counts)); if (depth < 1) stop("Unpaired rarefaction requires nonempty libraries.")
  B <- as.integer(cfg$gain_loss$n_downsample_replicates)
  results <- list()
  for (cohort in names(cohorts)) {
    ids <- cohorts[[cohort]]; key <- rna_unpaired_key(ctx, cfg, "prevalence", ids, depth)
    results[[cohort]] <- rna_unpaired_cache(file.path(dirs$cache, paste0(cohort, "_unpaired_prevalence.rds")), key, function() {
      full <- rna_unpaired_prevalence(ctx$counts, ctx, cfg, ids)
      enriched <- depleted <- valid <- numeric(nrow(full))
      for (b in seq_len(B)) {
        seed <- rna_stage_seed(cfg$execution$seeds$gain_loss, "unpaired-common-depth", b)
        rare <- dtu_with_rng(NULL, seed, rna_unpaired_rarefy(ctx$counts, depth))
        t <- rna_unpaired_prevalence(rare, ctx, cfg, ids)
        valid <- valid + t$quantifiable
        enriched <- enriched + (t$prevalence_class == "comparison_enriched")
        depleted <- depleted + (t$prevalence_class == "comparison_depleted")
        if (b == 1 || b %% 10 == 0 || b == B) rna_log("Unpaired prevalence ", cohort, ": rarefaction ", b, "/", B)
      }
      full$rarefaction_evaluable_probability <- valid / B
      # Denominator is ALL requested draws, so unavailable draws never inflate confidence.
      full$comparison_enriched_probability <- enriched / B
      full$comparison_depleted_probability <- depleted / B
      full$depth_robust <- full$quantifiable & ((full$prevalence_class == "comparison_enriched" & enriched / B >= cfg$gain_loss$depth_robust_probability) |
        (full$prevalence_class == "comparison_depleted" & depleted / B >= cfg$gain_loss$depth_robust_probability))
      full$cohort <- cohort; full$common_depth <- depth
      full
    }, cfg)
  }
  tab <- do.call(rbind, results); rownames(tab) <- NULL
  rna_write_tsv(tab, file.path(dirs$tables, "unpaired_isoform_prevalence.tsv.gz"))
  sample_depth <- data.frame(sample = ctx$samples$sample, condition = ctx$samples$condition %||% ctx$samples$timepoint, total_FL = colSums(ctx$counts), common_depth = depth)
  rna_write_tsv(sample_depth, file.path(dirs$qc, "unpaired_sample_depth.tsv"))
  if (!identical(cfg$reporting$enabled, FALSE) && nrow(tab)) {
    p <- ggplot2::ggplot(tab, ggplot2::aes(reference_prevalence, comparison_prevalence, colour = prevalence_class)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70") + ggplot2::geom_point(alpha = .5, na.rm = TRUE) +
      ggplot2::facet_wrap(~cohort) + ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      ggplot2::labs(title = "Independent-sample isoform prevalence", x = paste(cfg$comparison$reference, "prevalence"),
        y = paste(cfg$comparison$comparison, "prevalence"), colour = "Prevalence class")
    rna_save_plot(p, "unpaired_isoform_prevalence", dirs, cfg, 9, 6)
  }
  list(complete = TRUE, design = "unpaired", cohort_results = results, tables = list(prevalence = tab, sample_depth = sample_depth), rarefaction_target = depth)
}

rna_unpaired_geometry <- function(a, b) {
  p <- rowMeans(a); q <- rowMeans(b); mid <- (p + q) / 2
  entropy <- function(x) -sum(x[x > 0] * log2(x[x > 0]))
  kl <- function(x) sum(ifelse(x > 0, x * log2(pmax(x, .Machine$double.xmin) / pmax(mid, .Machine$double.xmin)), 0))
  c(jsd = (kl(p) + kl(q)) / 2, total_variation = sum(abs(p - q)) / 2,
    max_abs_delta_if = max(abs(q - p)), entropy_difference = entropy(q) - entropy(p))
}

rna_unpaired_permutation <- function(a, b, repetitions) {
  x <- cbind(a, b); n <- ncol(a); N <- ncol(x)
  observed <- rna_unpaired_geometry(a, b)[["jsd"]]
  n_assignments <- choose(N, n)
  if (is.finite(n_assignments) && n_assignments <= repetitions) {
    assignments <- utils::combn(seq_len(N), n)
    value <- apply(assignments, 2, function(i) rna_unpaired_geometry(x[, i, drop = FALSE], x[, -i, drop = FALSE])[["jsd"]])
    list(p = mean(value >= observed - 1e-12), draws = length(value), method = "exact_independent_label_permutation", assignments = n_assignments)
  } else {
    value <- replicate(repetitions, {i <- sample.int(N, n); rna_unpaired_geometry(x[, i, drop = FALSE], x[, -i, drop = FALSE])[["jsd"]]})
    list(p = (1 + sum(value >= observed - 1e-12)) / (repetitions + 1), draws = repetitions,
      method = "monte_carlo_independent_label_permutation_plus_one", assignments = n_assignments)
  }
}

run_unpaired_remodeling <- function(inputs, cfg, cohorts, stage_dir) {
  ctx <- rna_unpaired_context(inputs, cfg, cohorts); dirs <- rna_dirs(stage_dir)
  u <- ctx$settings; definition <- cfg$remodeling$definition
  confidence <- cfg$remodeling$bootstrap$bootstrap_confidence %||% .95
  if (length(confidence) != 1L || !is.finite(confidence) || confidence <= 0 || confidence > 1)
    stop("Invalid remodeling.bootstrap.bootstrap_confidence.")
  primary_retained <- ctx$mapped & ctx$comparable & rowSums(ctx$counts) >= (definition$min_isoform_total_fl %||% 5)
  genes <- rownames(ctx$gene_totals); results <- list()
  for (cohort in names(cohorts)) {
    ids <- cohorts[[cohort]]; key <- rna_unpaired_key(ctx, cfg, "remodeling", ids, primary_retained)
    results[[cohort]] <- rna_unpaired_cache(file.path(dirs$cache, paste0(cohort, "_unpaired_remodeling.rds")), key, function() {
      rows <- vector("list", length(genes)); jj <- match(ids, ctx$samples$sample)
      groups <- ctx$samples$timepoint[jj]
      block_size <- as.integer(cfg$execution$checkpoint_every_genes %||% 100L)
      if (!is.finite(block_size) || block_size < 1L) stop("Invalid gene checkpoint block size.")
      block_dir <- file.path(dirs$cache, "unpaired_gene_blocks", cohort, key)
      dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
      blocks <- split(seq_along(genes), ceiling(seq_along(genes) / block_size))
      for (block_number in seq_along(blocks)) {
        positions <- blocks[[block_number]]
        block_path <- file.path(block_dir, sprintf("block_%06d.rds", block_number))
        saved <- if (!isTRUE(cfg$execution$force_refit) && !identical(cfg$execution$resume, FALSE) && file.exists(block_path))
          tryCatch(readRDS(block_path), error = function(e) NULL) else NULL
        if (is.list(saved) && identical(saved$key, key) && isTRUE(saved$complete) &&
            identical(saved$genes, genes[positions]) && length(saved$rows) == length(positions)) {
          rows[positions] <- saved$rows
          rna_log("Reusing unpaired remodeling ", cohort, ": completed gene block ", block_number, "/", length(blocks))
          next
        }
        for (g in positions) {
        gene <- genes[g]; ii <- which(primary_retained & ctx$gene == gene)
        total <- ctx$gene_totals[gene, jj]; retained <- if (length(ii)) colSums(ctx$counts[ii, jj, drop = FALSE]) else rep(0, length(jj))
        good <- total >= (definition$min_gene_sample_count %||% 20) & retained / pmax(total, 1) >= (definition$min_retained_mass %||% .8)
        nr <- sum(good & groups == "primary"); nc <- sum(good & groups == "recurrent")
        ok <- length(ii) >= 2 && nr >= u$minimum_evaluable_per_condition && nc >= u$minimum_evaluable_per_condition
        row <- data.frame(cohort = cohort, gene_id = gene, design = "unpaired", quantifiable = ok,
          reference_evaluable = nr, comparison_evaluable = nc, retained_isoforms = length(ii),
          jsd = NA_real_, total_variation = NA_real_, max_abs_delta_if = NA_real_, entropy_difference = NA_real_,
          jsd_lower = NA_real_, jsd_upper = NA_real_, max_delta_lower = NA_real_, max_delta_upper = NA_real_,
          bootstrap_confidence = confidence,
          gene_model_p = NA_real_, permutation_draws = 0L, permutation_method = "not_tested",
          measurement_universe = "retained_cross_sample_comparable_isoforms_conditional_on_retained_mass",
          unpaired_effect_class = "not_quantifiable", stringsAsFactors = FALSE)
        if (ok) {
          x <- sweep(ctx$counts[ii, jj[good], drop = FALSE], 2, retained[good], "/")
          ref <- x[, groups[good] == "primary", drop = FALSE]; cmp <- x[, groups[good] == "recurrent", drop = FALSE]
          metric <- rna_unpaired_geometry(ref, cmp)
          row[names(metric)] <- as.list(metric)
          seed <- rna_stage_seed(cfg$execution$seeds$remodeling, "unpaired-remodeling", cohort, gene)
          inferred <- dtu_with_rng(NULL, seed, {
            boot <- replicate(u$n_bootstrap, rna_unpaired_geometry(ref[, sample.int(nr, nr, replace = TRUE), drop = FALSE], cmp[, sample.int(nc, nc, replace = TRUE), drop = FALSE]))
            list(ci = apply(boot[c("jsd", "max_abs_delta_if"), , drop = FALSE], 1, stats::quantile,
                   probs = c((1 - confidence) / 2, (1 + confidence) / 2), names = FALSE),
                 test = rna_unpaired_permutation(ref, cmp, u$n_permutations))
          })
          row$jsd_lower <- inferred$ci[1, "jsd"]; row$jsd_upper <- inferred$ci[2, "jsd"]
          row$max_delta_lower <- inferred$ci[1, "max_abs_delta_if"]; row$max_delta_upper <- inferred$ci[2, "max_abs_delta_if"]
          row$gene_model_p <- inferred$test$p; row$permutation_draws <- inferred$test$draws; row$permutation_method <- inferred$test$method
          row$unpaired_effect_class <- if (metric[["jsd"]] >= (definition$jsd_threshold %||% .1) || metric[["max_abs_delta_if"]] >= (definition$max_delta_if_threshold %||% .1)) "usage_shift" else "no_large_shift"
        }
        rows[[g]] <- row
          if (g == 1 || g %% 100 == 0 || g == length(genes)) rna_log("Unpaired remodeling ", cohort, ": gene ", g, "/", length(genes))
        }
        rna_atomic_rds(list(key = key, complete = TRUE, genes = genes[positions], rows = rows[positions]), block_path)
      }
      tab <- do.call(rbind, rows)
      tab$gene_model_q <- NA_real_; tested <- which(tab$quantifiable)
      tab$gene_model_q[tested] <- stats::p.adjust(tab$gene_model_p[tested], "BH")
      tab$genes_in_correction_family <- length(tested)
      tab$minimum_attainable_permutation_p <- ifelse(tab$permutation_draws > 0,
        ifelse(tab$permutation_method == "exact_independent_label_permutation",
          ifelse(tab$reference_evaluable == tab$comparison_evaluable, 2, 1) / tab$permutation_draws,
          1 / (tab$permutation_draws + 1)), NA_real_)
      tab$downstream_primary_remodeling_set <- tab$quantifiable & tab$gene_model_q < min(as.numeric(unlist(cfg$multiplicity$fdr_thresholds))) & tab$unpaired_effect_class == "usage_shift"
      tab$downstream_primary_remodeling_set[is.na(tab$downstream_primary_remodeling_set)] <- FALSE
      tab$inference_method <- "independent_group_bootstrap_and_label_permutation"
      tab
    }, cfg)
  }
  gene_table <- do.call(rbind, results); rownames(gene_table) <- NULL
  rna_write_tsv(gene_table, file.path(dirs$tables, "unpaired_gene_remodeling.tsv.gz"))
  rna_write_tsv(data.frame(isoform_harmonized_id = ctx$annotation$isoform_harmonized_id, retained = primary_retained), file.path(dirs$qc, "primary_cohort_retained_isoform_universe.tsv.gz"))
  if (!identical(cfg$reporting$enabled, FALSE) && nrow(gene_table)) {
    p <- ggplot2::ggplot(gene_table, ggplot2::aes(max_abs_delta_if, jsd, colour = unpaired_effect_class)) + ggplot2::geom_point(alpha = .6, na.rm = TRUE) +
      ggplot2::facet_wrap(~cohort) + ggplot2::labs(title = "Independent-sample isoform remodeling", x = "Maximum change in mean isoform fraction", y = "Jensen-Shannon divergence", colour = "Effect class")
    rna_save_plot(p, "unpaired_remodeling_effects", dirs, cfg, 9, 6)
  }
  list(complete = TRUE, design = "unpaired", gene_table = gene_table, cohort_results = results,
    patient_inference = NULL, enrichment = list(status = "not_applicable_paired_architecture_sets", complete = TRUE),
    method = "independent_group_bootstrap_and_label_permutation", paired_architecture_classes = "not_applicable")
}

run_unpaired_qc <- function(inputs, cfg, cohorts, gain_loss = NULL, stage_dir) {
  ctx <- rna_unpaired_context(inputs, cfg, cohorts); dirs <- rna_dirs(stage_dir)
  x <- ctx$counts
  sample_qc <- data.frame(sample = ctx$samples$sample, subject = ctx$samples$patient,
    condition = ctx$samples$condition %||% ctx$samples$timepoint,
    total_full_length_reads = colSums(x), detected_isoforms = colSums(x > 0),
    detected_genes = colSums(ctx$gene_totals > 0), design = "unpaired", stringsAsFactors = FALSE)
  tests <- list()
  for (cohort in names(cohorts)) for (metric in c("total_full_length_reads", "detected_isoforms", "detected_genes")) {
    ii <- match(cohorts[[cohort]], sample_qc$sample); v <- sample_qc[[metric]][ii]
    ref <- ctx$samples$timepoint[ii] == "primary"; observed <- mean(v[!ref]) - mean(v[ref])
    seed <- rna_stage_seed(cfg$execution$seeds$import, "unpaired-qc", cohort, metric)
    null <- dtu_with_rng(NULL, seed, replicate(ctx$settings$n_permutations, {s <- sample.int(length(v), sum(ref)); mean(v[-s]) - mean(v[s])}))
    tests[[length(tests) + 1]] <- data.frame(cohort = cohort, metric = metric,
      reference_mean = mean(v[ref]), comparison_mean = mean(v[!ref]), mean_difference = observed,
      p_value = (1 + sum(abs(null) >= abs(observed) - 1e-12)) / (length(null) + 1),
      design = "unpaired", method = "independent_label_permutation_two_sided")
  }
  tests <- do.call(rbind, tests); tests$fdr <- NA_real_
  for (cohort in unique(tests$cohort)) {i <- tests$cohort == cohort; tests$fdr[i] <- stats::p.adjust(tests$p_value[i], "BH")}
  rna_write_tsv(sample_qc, file.path(dirs$tables, "unpaired_sample_qc.tsv"))
  rna_write_tsv(tests, file.path(dirs$tables, "unpaired_condition_qc_tests.tsv"))
  if (!identical(cfg$reporting$enabled, FALSE)) {
    p <- ggplot2::ggplot(sample_qc, ggplot2::aes(sample, total_full_length_reads, fill = condition)) + ggplot2::geom_col() +
      ggplot2::labs(title = "Full-length reads in independent libraries", x = "Sample", y = "Full-length reads", fill = "Condition")
    rna_save_plot(p, "full_length_reads_by_sample", dirs, cfg, 9, 5)
  }
  list(complete = TRUE, design = "unpaired", sample_qc = sample_qc, tests = tests,
    paired_technical_associations = "not_applicable_independent_samples")
}
