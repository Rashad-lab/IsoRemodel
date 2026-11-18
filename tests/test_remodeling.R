# Mathematical and decision-rule tests against the actual vendored engine functions.
# These are quick tests: no production simulations or biological hit requirements.
local({
  root <- if (exists(".repo_root",inherits=TRUE)) get(".repo_root",inherits=TRUE) else getwd()
  e <- .rna_remodeling_test_environment(root)
  equal <- function(x,y,tol=1e-10) stopifnot(isTRUE(all.equal(x,y,tolerance=tol,check.attributes=FALSE)))
  equal(e$jsd_base2(c(1,0),c(0,1)),1)
  equal(e$jsd_base2(c(.1,.2,.7),c(.1,.2,.7)),0)
  p <- c(.65,.25,.10); q <- c(.15,.35,.50)
  equal(e$jsd_base2(p,q),e$jsd_natural(p,q)/log(2))
  equal(sum(e$jsd_contributions(p,q)),e$jsd_base2(p,q))
  equal(e$total_variation(p,q),.5*sum(abs(q-p)))
  stopifnot(max(abs(q-p)) <= e$total_variation(p,q)+1e-12)
  # Exact zero values are descriptive data, not imputed pseudocounts.
  equal(e$entropy_bits(c(1,0,0)),0)
  equal(e$entropy_bits(c(.5,.5,0)),1)
  # The 0.5 pseudocount is restricted to finite-count log-ratio modeling.
  smooth <- e$smoothed_composition(matrix(c(10,0,0,10),nrow=2,byrow=TRUE),.5)
  stopifnot(all(smooth$composition>0))
  equal(rowSums(smooth$composition),c(1,1))
  for (k in c(2L,3L,6L)) {
    basis <- e$ilr_basis(k)
    equal(crossprod(basis),diag(k-1L))
    equal(colSums(basis),rep(0,k-1L))
    pp <- matrix(seq_len(2*k),nrow=2); pp <- pp/rowSums(pp)
    equal(e$ilr_inverse(e$ilr_transform(pp,basis),basis),pp)
  }
  # Dominance requires a unique top transcript, mean fraction, margin and support.
  ids <- c("A","B","C")
  comp <- matrix(rep(c(.7,.2,.1),each=4),nrow=4)
  dom <- e$dominance_summary(colMeans(comp),comp,ids)
  stopifnot(dom$has_dominant,dom$dominant_id=="A",dom$support_n==4L)
  tied <- matrix(rep(c(.45,.45,.1),each=4),nrow=4)
  tie_dom <- e$dominance_summary(colMeans(tied),tied,ids)
  stopifnot(!tie_dom$has_dominant)
  secondary <- e$secondary_remodeling(c(.60,.39,.01),c(.60,.01,.39),"A","A",TRUE,TRUE,ids)
  stopifnot(secondary$secondary_effect_pass,secondary$secondary_mass_pass)
  equal(secondary$secondary_max_abs_delta_if,.38)
  same <- list(has_dominant=TRUE,dominant_id="A")
  different <- list(has_dominant=TRUE,dominant_id="B")
  absent <- list(has_dominant=FALSE,dominant_id=NA_character_)
  classify <- e$classify_architecture
  stopifnot(
    classify(TRUE,TRUE,same,different,secondary)=="A",
    classify(TRUE,TRUE,same,same,secondary)=="B",
    classify(TRUE,TRUE,same,absent,secondary)=="C",
    classify(TRUE,TRUE,absent,same,secondary)=="D",
    classify(TRUE,TRUE,absent,absent,secondary)=="E",
    classify(TRUE,FALSE,same,same,secondary)=="Below_effect",
    classify(FALSE,TRUE,same,different,secondary)=="NQ",
    classify(TRUE,TRUE,same,same,list(secondary_effect_pass=FALSE))=="Unresolved_effect"
  )
  # Monte Carlo uses the plus-one estimate, including at zero exceedances.
  mc <- e$monte_carlo_p(c(.1,.2,.3),.8)
  equal(mc$p,1/4)
  stopifnot(mc$exceedances==0L)
  # The actual additive REML fit keeps the zero biological-variance boundary.
  difference <- matrix(rep(.2,4),ncol=1)
  sampling_covariance <- rep(list(matrix(.02,nrow=1,ncol=1)),4)
  fit <- e$fit_additive_reml(difference,sampling_covariance)
  stopifnot(fit$converged,fit$boundary,fit$estimate==0,is.finite(fit$information))
  profile <- e$profile_variance_mixture(difference,sampling_covariance)
  stopifnot(profile$available,profile$tail_resolved,all(profile$weight>=0))
  equal(sum(profile$weight),1)
  stopifnot(!e$valid_reml_inputs(difference,rep(list(matrix(-.1,nrow=1)),4)))
  # Finite-count covariance is reproducible and positive definite; this is the
  # source multinomial estimator rather than an analytical variance substitute.
  e$cfg$gene_null$covariance_initial_replicates <- 128L
  e$cfg$gene_null$covariance_maximum_replicates <- 256L
  e$cfg$gene_null$covariance_trace_tolerance <- 2
  e$cfg$gene_null$covariance_frobenius_tolerance <- 2
  cv1 <- e$finite_count_covariance(c(60,40),e$ilr_basis(2),128)
  cv2 <- e$finite_count_covariance(c(60,40),e$ilr_basis(2),128)
  stopifnot(cv1$available,all(eigen(cv1$covariance,symmetric=TRUE)$values>0))
  equal(cv1$covariance,cv2$covariance)
  # Deterministic stage namespaces and the caller's RNG are preserved.
  stopifnot(e$deterministic_seed(1,"bootstrap","geneA") != e$deterministic_seed(2,"bootstrap","geneA"))
  set.seed(113); before <- .Random.seed
  first <- e$with_seed(412,runif(8))
  stopifnot(identical(before,.Random.seed))
  equal(first,e$with_seed(412,runif(8)))
  # YAML adapters accept arbitrary threshold lists and do not force six subjects.
  c2 <- e$rna_remodeling_engine_config(list(
    multiplicity=list(fdr_thresholds=c(.01,.05,.10)),
    remodeling=list(definition=list(min_evaluable_pairs=3,dominance_support_n=3),
      geometry=list(geometry_isoform_numbers=c(2,5))),
    execution=list(seeds=list(remodeling=99),profile="smoke")))
  equal(c2$inference$fdr_thresholds,c(.01,.05,.10))
  equal(c2$geometry$isoform_numbers,c(2L,5L))
  stopifnot(c2$thresholds$min_evaluable_pairs==3L,c2$seed==99L)
  message("Remodeling mathematical and classification tests passed.")
})
