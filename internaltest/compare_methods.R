## Compare bamp() results across MCMC methods ("pg" vs "taylor") and, within
## method="pg", across pg_engine ("C" vs "R") on one dataset: every variant
## gets the same model spec and MCMC length, so the number of iterations (raw
## and post-thinning), DIC, Rhat, convergence, runtime and the fitted effects
## are directly comparable.

## Internal: fit `variants` (a named list of bamp() argument overrides, e.g.
## list(pg = list(method = "pg"), taylor = list(method = "taylor"))) on the same
## data/model spec, and collect n_iter (raw) / n_kept (thinned, per chain) /
## DIC / pD / mean deviance / max Rhat (PSRF, the identified-quantity criterion
## used by checkConvergence()) / runtime.
.compare_bamp_variants <- function(cases, population, periods_per_agegroup, variants,
                                    age, period, cohort, overdisp,
                                    mcmc.options, hyperpar, parallel, verbose) {
  labels <- names(variants)
  fits <- stats::setNames(vector("list", length(labels)), labels)
  runtime_sec <- stats::setNames(numeric(length(labels)), labels)
  used_mcmc <- stats::setNames(vector("list", length(labels)), labels)

  for (nm in labels) {
    if (isTRUE(verbose)) message("Fitting ", nm, " ...")
    t0 <- Sys.time()
    args <- utils::modifyList(
      list(cases = cases, population = population, age = age, period = period,
           cohort = cohort, overdisp = overdisp, periods_per_agegroup = periods_per_agegroup,
           mcmc.options = mcmc.options, hyperpar = hyperpar, dic = TRUE,
           parallel = parallel, verbose = FALSE),
      variants[[nm]])
    used_mcmc[[nm]] <- args$mcmc.options
    fits[[nm]] <- tryCatch(do.call(bamp::bamp, args), error = function(e) {
      message("  FAILED (", nm, "): ", conditionMessage(e))
      NULL
    })
    runtime_sec[[nm]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  table <- do.call(rbind, lapply(labels, function(nm) {
    f <- fits[[nm]]
    ## "real" iterations per chain as actually requested for this variant (a
    ## variant may override mcmc.options, e.g. taylor doubling number_of_iterations
    ## below); NA if left as "auto" (never resolved back onto the fitted object,
    ## so not recoverable from it).
    n_iter <- suppressWarnings(as.numeric(used_mcmc[[nm]]$number_of_iterations))
    ## a fit can "succeed" (no R error) with every chain discarded by bamp()'s
    ## own convergence-based chain pruning (warning, not error) -- guard for
    ## that the same way checkConvergence() does, so rhat/n_kept stay scalar
    ## NA instead of empty vectors that would break the data.frame() below.
    has_samples <- !is.null(f) && length(f$samples) > 0 && !is.null(f$samples$intercept) &&
      coda::nchain(f$samples$intercept) >= 1
    rhat <- if (has_samples) {
      vals <- tryCatch(bamp:::.apc_identified_psrf(f), error = function(e) NA_real_)
      if (length(vals) == 0 || !any(is.finite(vals))) NA_real_ else max(vals, na.rm = TRUE)
    } else NA_real_
    ## "ausgedünnt": samples actually kept per chain after burn-in + thinning
    n_kept <- if (has_samples)
      tryCatch(as.integer(coda::niter(f$samples$intercept)), error = function(e) NA_integer_)
    else NA_integer_
    ## effective sample size of the deviance chain, pooled over chains: makes
    ## mixing quality (autocorrelation) visible directly, independent of n_kept
    ## -- e.g. taylor's block-MH proposals typically mix worse than pg's exact
    ## Gibbs conditionals, which inflates pD/DIC (pD = mean(D) - D(eta_bar) is
    ## sensitive to the posterior spread of D) even when posterior medians agree.
    ess <- if (has_samples)
      tryCatch(sum(coda::effectiveSize(f$samples$deviance)), error = function(e) NA_real_)
    else NA_real_
    ## deviance/DIC/pD can likewise come back as numeric(0) (not NA) when bamp()
    ## discarded every chain -- coerce anything not a single scalar to NA so a
    ## fit that "succeeded" with 0 surviving chains still yields one table row.
    scalar_or_na <- function(val, na = NA_real_) if (length(val) == 1 && !is.na(val)) val else na
    data.frame(
      variant = nm,
      n_iter = n_iter,
      n_kept = n_kept,
      dic = scalar_or_na(f$deviance$DIC),
      pD = scalar_or_na(f$deviance$pD),
      mean_deviance = scalar_or_na(f$deviance$mean.deviance),
      rhat = rhat,
      ess = ess,
      converged = if (!is.null(f)) isTRUE(is.finite(rhat) && rhat <= 1.1) else NA,
      runtime_sec = round(runtime_sec[[nm]], 1),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }))

  structure(list(table = table, models = fits), class = "bampmethodcompare")
}

compare_bamp_methods <- function(
    cases,
    population,
    periods_per_agegroup,
    age = "rw1",
    period = "rw1",
    cohort = "rw1",
    overdisp = FALSE,
    methods = c("pg", "taylor"),
    mcmc.options = list(number_of_iterations = 40000, burn_in = 20000, step = 10, tuning = 500),
    hyperpar = list(age = c(1, 0.5), period = c(1, 5e-04), cohort = c(1, 5e-04), overdisp = c(1, 0.05)),
    parallel = TRUE,
    verbose = TRUE
) {
  ## taylor (the legacy block-MH sampler) mixes more slowly than pg and needs
  ## more iterations for a fair comparison -- double it here rather than in
  ## the shared mcmc.options, so pg keeps the length the caller asked for.
  ## Only scales a numeric number_of_iterations; "auto" is left alone (bamp()
  ## resolves it internally, so there is nothing to double here).
  variants <- stats::setNames(lapply(methods, function(m) {
    v <- list(method = m)
    if (identical(m, "taylor") && is.numeric(mcmc.options$number_of_iterations))
      v$mcmc.options <- list(number_of_iterations = 2 * mcmc.options$number_of_iterations)
    v
  }), methods)
  .compare_bamp_variants(cases, population, periods_per_agegroup, variants,
                          age, period, cohort, overdisp, mcmc.options, hyperpar, parallel, verbose)
}

## Compare the native C and pure-R implementations of the pg engine (method is
## always "pg"); the two are documented to agree to numerical tolerance, so
## this is mainly a runtime/Rhat sanity check, not a model-selection tool.
compare_pg_engines <- function(
    cases,
    population,
    periods_per_agegroup,
    age = "rw1",
    period = "rw1",
    cohort = "rw1",
    overdisp = TRUE,
    engines = c("C", "R"),
    mcmc.options = list(number_of_iterations = 20000, burn_in = 10000, step = 10, tuning = 500),
    hyperpar = list(age = c(1, 0.5), period = c(1, 5e-04), cohort = c(1, 5e-04), overdisp = c(1, 0.05)),
    parallel = TRUE,
    verbose = TRUE
) {
  variants <- stats::setNames(
    lapply(engines, function(e) list(method = "pg", pg_engine = e)), engines)
  .compare_bamp_variants(cases, population, periods_per_agegroup, variants,
                          age, period, cohort, overdisp, mcmc.options, hyperpar, parallel, verbose)
}

#' @export
print.bampmethodcompare <- function(x, ...) {
  print(x$table, row.names = FALSE)
  invisible(x)
}

## Overlay the fitted age/period/cohort effect curves of every successfully
## fitted variant: posterior median (solid) plus pointwise 80% (0.1/0.9,
## dashed) and 95% (0.025/0.975, more dashed) credible bands -- the same
## widening-dash convention plot.apc() uses for multiple quantiles. Age,
## period and cohort are each their own plot (no par(mfrow)), so on an
## interactive device (RStudio Plots pane, X11, ...) they land on three
## separate pages/figures rather than stacked into one.
#' @export
plot.bampmethodcompare <- function(x, quantiles = c(0.025, 0.1, 0.5, 0.9, 0.975), ...) {
  fits <- Filter(Negate(is.null), x$models)
  if (length(fits) < 1) { message("No successful fits to plot."); return(invisible(NULL)) }
  palette <- c("blue", "green", "red", "orange")
  cols <- stats::setNames(palette[seq_along(fits)], names(fits))

  q <- length(quantiles)
  lty <- pmin(1L + round(abs(seq_len(q) - (q + 1) / 2)), 6L)  # solid median, dashed tails

  eff <- lapply(fits, function(f)
    tryCatch(stats::effects(f, quantiles = quantiles), error = function(e) NULL))

  for (part in c("age", "period", "cohort")) {
    has <- vapply(eff, function(e) !is.null(e) && !is.null(e[[part]]), logical(1))
    if (!any(has)) next
    curves <- lapply(eff[has], `[[`, part)   # each: n_levels x q matrix
    yr <- range(unlist(curves), na.rm = TRUE)
    plot(seq_len(nrow(curves[[1]])), curves[[1]][, 1], type = "n", ylim = yr,
         ylab = paste(part, "effect"), xlab = part)
    for (k in seq_along(curves))
      for (i in seq_len(q))
        lines(curves[[k]][, i], lty = lty[i], col = cols[has][k])
    graphics::legend("topleft", legend = names(fits)[has], col = cols[has], lty = 1, bty = "n")
  }
}

if (FALSE) {
  library(bamp)
  source("internaltest/compare_methods.R")
  source("internaltest/nordstat.R")  # for bamp_set_to_data()

  sets <- data.table::as.data.table(readRDS(
    "/Users/volkerschmid/projects/Apc-Full-vs-Empirical-Bayes--R.Inla/nord_data_prepare/data/nordcan_processed/nordcan_apc_sets.rds"))
  dat <- bamp_set_to_data(sets$data[[which(sets$set_id == "Sweden__male__Lip")]])
  dat <- bamp_set_to_data(sets$data[[which(sets$set_id == "Denmark__male__Lung")]])

  # method comparison: pg vs taylor
  cmp <- compare_bamp_methods(dat$cases, dat$population, dat$periods_per_agegroup)
  print(cmp)   # DIC / pD / deviance / Rhat / convergence / runtime, one row per method
  plot(cmp)    # fitted age/period/cohort effects, one colour per method

  selectModel(dat$cases,dat$population,periods_per_agegroup = dat$periods_per_agegroup, screen = list(number_of_iterations = 10000, burn_in = 5000,
                                          step = 5, tuning = 200), dic_margin = 0.1)
  bamp(dat$cases,dat$population,"rw1+het","rw1","rw1",periods_per_agegroup = dat$periods_per_agegroup, mcmc.options = list(number_of_iterations = 10000, burn_in = 5000,
                                          step = 5, tuning = 200))
  selectModel(dat$cases,dat$population,periods_per_agegroup = dat$periods_per_agegroup, screen = list(number_of_iterations = 10000, burn_in = 5000,
                                          step = 5, tuning = 200), dic_margin = 0.1, try_heterogeneity = TRUE)
  # engine comparison: pg's native C engine vs the pure-R reference
  cmpe <- compare_pg_engines(dat$cases, dat$population, dat$periods_per_agegroup)
  print(cmpe)
  plot(cmpe)
}
