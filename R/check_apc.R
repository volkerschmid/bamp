#' Check apc object, whether MCMC has converged
#'
#' @param x An apc object
#' @param info logical; print more information (including the raw per-effect
#'   diagnostic, which is affected by the age-period-cohort identifiability and
#'   should not be used on its own, see Details)
#' @param level level of check; 1 uses point estimate, 2 uses upper C.I.
#' @param auto logical; should be TRUE if called automatically from \code{\link{bamp}}
#'
#' @description
#' This function uses Gelman and Rubin's R (potential scale reduction factor) to
#' check convergence. All checked quantities should have R<1.1.
#' \code{\link{bamp}} runs at least four MCMC chains by default (more if
#' \code{parallel} is more than four).
#'
#' @details
#' In an age-period-cohort model the age, period and cohort effects are linearly
#' dependent (Clayton and Schifflers, 1987): a linear trend can be moved between
#' the three effects without changing the likelihood. The individual effect
#' chains can therefore drift along this non-identified direction even when the
#' model has fully converged, which makes a naive Gelman-R on the raw effects
#' report spurious non-convergence.
#'
#' \code{checkConvergence} therefore assesses the quantities that are actually
#' identified: the smoothing precisions and the fitted linear predictor
#' (log-odds) in every cell of the Lexis diagram, which is invariant to the
#' trend re-allocation. With \code{info=TRUE} the raw per-effect diagnostic is
#' also printed for reference.
#'
#' @import coda
#' @return logical; TRUE if check is fine.
#' @export
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' checkConvergence(model)
#' }
checkConvergence <- function(x, info = FALSE, level = 2, auto = FALSE)
{
  ## guards preserved from volkerschmid/bamp 2.2.0: without samples or with a
  ## single chain the identified-quantity PSRFs below are all NA, and
  ## !any(NA > thr, na.rm = TRUE) would spuriously report convergence -- so bail
  ## out explicitly here instead.
  if (length(x$samples) == 0 || is.null(x$samples$intercept)) {
    if (!auto) message("No MCMC samples available (all chains may have failed).")
    return(FALSE)
  }
  if (coda::nchain(x$samples$intercept) < 2) {
    if (!auto) message("Only 1 MCMC chain available: Gelman-Rubin diagnostic requires at least 2 chains.")
    return(FALSE)
  }
  j <- level
  thr <- 1.1

  ## ---- identified quantities: smoothing precisions + fitted log-odds -------
  ## shared with selectModel(), which must judge convergence identically; see
  ## .apc_identified_psrf.
  psrf <- .apc_identified_psrf(x, level = j)

  ## ---- raw per-effect diagnostic (NOT identified for full APC) ----
  maxpsrf <- function(mcl) {
    if (is.null(mcl) || length(mcl) < 2) return(NA_real_)
    d <- tryCatch(
      coda::gelman.diag(mcl, multivariate = FALSE, autoburnin = FALSE)$psrf,
      error = function(e) NULL)
    if (is.null(d)) return(NA_real_)
    max(d[, min(j, ncol(d))], na.rm = TRUE)
  }
  has_age    <- !is.null(x$model$age)    && x$model$age    != " "
  has_period <- !is.null(x$model$period) && x$model$period != " "
  has_cohort <- !is.null(x$model$cohort) && x$model$cohort != " "
  raw <- c()
  if (has_age)    raw <- c(raw, age    = maxpsrf(x$samples$age))
  if (has_period) raw <- c(raw, period = maxpsrf(x$samples$period))
  if (has_cohort) raw <- c(raw, cohort = maxpsrf(x$samples$cohort))

  ok <- !any(psrf > thr, na.rm = TRUE)

  if (!ok && !auto) {
    cat("Warning: MCMC chains did not converge (identified quantities)!\n")
  }
  if (info) {
    cat("Gelman-Rubin R (column = ",
        if (j == 1) "point estimate" else "upper C.I.", ")\n", sep = "")
    cat("  Identified quantities (used for the decision):\n")
    for (nm in names(psrf))
      cat(sprintf("    %-12s %.3f%s\n", nm, psrf[nm],
                  ifelse(!is.na(psrf[nm]) && psrf[nm] > thr, "  <-- > 1.1", "")))
    if (length(raw)) {
      cat("  Raw per-effect (drifts along the non-identified A-P-C trend; not decisive):\n")
      for (nm in names(raw))
        cat(sprintf("    %-12s %.3f\n", nm, raw[nm]))
    }
  }
  return(ok)
}

## Internal: the identified-quantity PSRFs that checkConvergence() judges
## convergence on -- the smoothing precisions of the effects present in the
## model, plus the fitted linear predictor (log-odds). The intercept and
## overdispersion precision are deliberately excluded: they are noisier and
## not decisive for judging convergence of the smooth effects. This is the
## SINGLE shared implementation -- selectModel()'s screening/refit fits call
## it too, so both judge convergence identically (the original bug this
## guards against: selectModel() judged convergence on the fitted log-odds
## alone, so it could adopt a model whose smoothing precisions had not
## actually converged).
.apc_identified_psrf <- function(x, level = 2) {
  maxpsrf <- function(mcl) {
    if (is.null(mcl) || length(mcl) < 2) return(NA_real_)
    d <- tryCatch(
      coda::gelman.diag(mcl, multivariate = FALSE, autoburnin = FALSE)$psrf,
      error = function(e) NULL)
    if (is.null(d)) return(NA_real_)
    max(d[, min(level, ncol(d))], na.rm = TRUE)
  }
  has_age    <- !is.null(x$model$age)    && x$model$age    != " "
  has_period <- !is.null(x$model$period) && x$model$period != " "
  has_cohort <- !is.null(x$model$cohort) && x$model$cohort != " "

  psrf <- c()
  if (has_age)    psrf <- c(psrf, age_prec    = maxpsrf(x$samples$age_parameter))
  if (has_period) psrf <- c(psrf, period_prec = maxpsrf(x$samples$period_parameter))
  if (has_cohort) psrf <- c(psrf, cohort_prec = maxpsrf(x$samples$cohort_parameter))
  c(psrf, fitted = .apc_eta_psrf(x, level = level))
}

## Internal: maximum Gelman-R over the fitted linear predictor eta = mu + age +
## period + cohort across all Lexis cells. eta is invariant to the A-P-C trend
## re-allocation, so it is the right thing to monitor. Returns NA if it cannot
## be computed (e.g. <2 chains).
.apc_eta_psrf <- function(x, level = 2, max_cells = 400L) {
  s <- x$samples
  nchains <- length(s$intercept)
  if (is.null(nchains) || nchains < 2) return(NA_real_)

  has_age    <- !is.null(x$model$age)    && x$model$age    != " "
  has_period <- !is.null(x$model$period) && x$model$period != " "
  has_cohort <- !is.null(x$model$cohort) && x$model$cohort != " "

  ppa <- x$data$periods_per_agegroup
  I <- if (has_age)    ncol(as.matrix(s$age[[1]]))    else 1L
  J <- if (has_period) ncol(as.matrix(s$period[[1]])) else
       if (!is.null(x$data$cases)) ncol(x$data$cases) else 1L
  if (is.null(J) || is.na(J)) J <- 1L

  ## enumerate (i,j) cells, subsample if very large
  cells <- expand.grid(i = seq_len(I), jj = seq_len(J))
  if (nrow(cells) > max_cells)
    cells <- cells[round(seq(1, nrow(cells), length.out = max_cells)), , drop = FALSE]

  eta_list <- vector("list", nchains)
  for (c in seq_len(nchains)) {
    mu <- as.matrix(s$intercept[[c]])[, 1]
    A  <- if (has_age)    as.matrix(s$age[[c]])    else NULL
    P  <- if (has_period) as.matrix(s$period[[c]]) else NULL
    C  <- if (has_cohort) as.matrix(s$cohort[[c]]) else NULL
    S  <- length(mu)
    E  <- matrix(mu, nrow = S, ncol = nrow(cells))
    for (col in seq_len(nrow(cells))) {
      i <- cells$i[col]; jj <- cells$jj[col]
      if (has_age)    E[, col] <- E[, col] + A[, i]
      if (has_period) E[, col] <- E[, col] + P[, jj]
      if (has_cohort) E[, col] <- E[, col] + C[, coh(i, jj, I, ppa)]
    }
    eta_list[[c]] <- coda::mcmc(E)
  }
  d <- tryCatch(
    coda::gelman.diag(coda::as.mcmc.list(eta_list),
                      multivariate = FALSE, autoburnin = FALSE)$psrf,
    error = function(e) NULL)
  if (is.null(d)) return(NA_real_)
  max(d[, min(level, ncol(d))], na.rm = TRUE)
}
