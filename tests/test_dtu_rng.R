# Base-R regression checks: DTU precision sampling must reproduce the paper's
# RNG algorithm without inheriting or changing the caller's RNG state.
(function() {
  original <- dtu_rng_snapshot()
  on.exit(dtu_rng_restore(original), add = TRUE)
  spec <- list(kind = "L'Ecuyer-CMRG", normal_kind = "Inversion", sample_kind = "Rejection")
  same_state <- function(a, b) identical(a$kind, b$kind) &&
    identical(a$seed_present, b$seed_present) && identical(a$seed, b$seed)
  stopifnot(identical(dtu_rng_spec(), spec),
    identical(dtu_cohort_seed(260605L, "all_6_pairs"), 261605L),
    identical(dtu_cohort_seed(260605L, "excluding_ptn2_5_pairs"), 262605L))

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(125L)
  before <- dtu_rng_snapshot()
  # Validate the effective configuration without consuming the caller's stream.
  stopifnot(identical(dtu_rng_spec(list(kind = "L'Ecuyer-CMRG")), spec),
            same_state(before, dtu_rng_snapshot()))
  contract_default <- dtu_rng_contract(NULL, 260605L, c("all_6_pairs", "excluding_ptn2_5_pairs"))
  after_default <- dtu_rng_snapshot()
  stopifnot(same_state(before, after_default), identical(contract_default$seed, c(261605L, 262605L)))
  # The actual common-precision subset-size convention is round(0.1 * N).
  sample_a <- dtu_with_rng(spec, 261605L, sample.int(6826L, round(.1 * 6826L)))
  stopifnot(same_state(before, dtu_rng_snapshot()))
  error <- try(dtu_with_rng(spec, 261605L, {
    runif(12L)
    stop("intentional scoped RNG error")
  }), silent = TRUE)
  stopifnot(inherits(error, "try-error"), same_state(before, dtu_rng_snapshot()))

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(261605L)
  legacy_draw <- sample.int(6826L, round(.1 * 6826L))
  stopifnot(identical(sample_a, legacy_draw))
  set.seed(777L)
  before_other <- dtu_rng_snapshot()
  sample_b <- dtu_with_rng(spec, 261605L, sample.int(6826L, round(.1 * 6826L)))
  contract_other <- dtu_rng_contract(spec, 260605L, c("all_6_pairs", "excluding_ptn2_5_pairs"))
  stopifnot(identical(sample_a, sample_b), same_state(before_other, dtu_rng_snapshot()),
            identical(serialize(contract_default, NULL), serialize(contract_other, NULL)))

  # Algorithm and seed must both participate in effective cache/provenance data.
  twister <- spec; twister$kind <- "Mersenne-Twister"
  contract_twister <- dtu_rng_contract(twister, 260605L, c("all_6_pairs", "excluding_ptn2_5_pairs"))
  contract_seed <- dtu_rng_contract(spec, 260606L, c("all_6_pairs", "excluding_ptn2_5_pairs"))
  stopifnot(!identical(serialize(contract_default, NULL), serialize(contract_twister, NULL)),
            !identical(serialize(contract_default, NULL), serialize(contract_seed, NULL)),
            !identical(sample_a, dtu_with_rng(twister, 261605L, sample.int(6826L, 683L))))

  # Missing .Random.seed is a valid caller state, including after failure.
  rm(".Random.seed", envir = .GlobalEnv)
  before_absent <- dtu_rng_snapshot()
  stopifnot(!before_absent$seed_present)
  dtu_rng_spec(spec)
  stopifnot(same_state(before_absent, dtu_rng_snapshot()))
  invisible(dtu_with_rng(spec, 262605L, runif(3L)))
  stopifnot(same_state(before_absent, dtu_rng_snapshot()))
  error <- try(dtu_with_rng(spec, 262605L, stop("intentional missing-seed error")), silent = TRUE)
  stopifnot(inherits(error, "try-error"), same_state(before_absent, dtu_rng_snapshot()))
  invalid <- list(list(kind = "invalid-RNG"), list(kind = NA_character_),
                  list(kind = c("L'Ecuyer-CMRG", "Mersenne-Twister")),
                  list(sample_kind = "invalid-sampler"), list(typo = "L'Ecuyer-CMRG"),
                  list(normal_kind = NULL), "L'Ecuyer-CMRG")
  for (bad in invalid) {
    error <- try(dtu_rng_spec(bad), silent = TRUE)
    stopifnot(inherits(error, "try-error"), same_state(before_absent, dtu_rng_snapshot()))
  }
  for (bad_seed in list(NA_real_, 0L, 1.5, numeric(), c(1L, 2L), "261605")) {
    error <- try(dtu_with_rng(spec, bad_seed, runif(1L)), silent = TRUE)
    stopifnot(inherits(error, "try-error"), same_state(before_absent, dtu_rng_snapshot()))
  }
})()
message("DTU RNG contract and caller-state regression tests passed.")
