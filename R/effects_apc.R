#' Effects from Fitted APC Model
#'
#' @param object an apc object
#' @param mean logical. If TRUE, mean effects are computed
#' @param quantiles Scalar or vector of quantiles to compute (only if mean=FALSE)
#' @param update logical. If TRUE, the apc object including the effects is returned
#' @param convention display-layer gauge convention for the linear trend (drift)
#'   in a full age-period-cohort model; one of \code{"age"} (default),
#'   \code{"period"}, \code{"cohort"} or \code{"none"}. See Details.
#' @param combined logical. For heterogeneity models (\code{"rw1+het"} /
#'   \code{"rw2+het"}), if \code{TRUE} the returned effect is the full effect
#'   (smooth + iid heterogeneity component); if \code{FALSE} (default) only the
#'   smooth component is returned. Ignored for models without heterogeneity.
#' @param ... Additional arguments will be ignored
#'
#' @details
#' In a full age-period-cohort model the age, period and cohort effects are
#' identifiable only up to one shared linear trend (drift), because a linear
#' trend can be moved between the three effects without changing the fitted
#' rates (Clayton and Schifflers, 1987). When the raw posterior samples are
#' summarised directly, this non-identified trend makes the effect curves drift
#' between MCMC iterations and between runs, so they are not reproducible.
#'
#' \code{convention} fixes that single degree of freedom for display:
#' \code{"age"} removes the linear slope of the age effect (age is shown as
#' curvature about a zero trend, the drift is shown in the period and cohort
#' effects); \code{"period"} and \code{"cohort"} pin the corresponding effect's
#' slope to zero instead; \code{"none"} returns the raw, un-gauged effects
#' (previous behaviour, curves may differ between runs). The gauge is applied
#' per MCMC draw before the quantiles are computed, and only for full APC models
#' -- for models without all three effects there is no trend aliasing and the
#' argument is ignored.
#'
#' Fixing the gauge removes the run-to-run linear-trend (drift) component, which
#' is typically the dominant source of non-reproducibility in the effect curves;
#' the residual curvature and Monte-Carlo sampling noise are unaffected, so two
#' independent runs are made much closer but need not agree exactly. The
#' zero-slope property holds for every individual MCMC draw; because quantiles
#' are non-linear, the summarised median/quantile curve has only approximately
#' (not exactly) zero slope.
#'
#' The gauge is display-only: it never modifies the stored samples
#' (\code{object$samples}), and the fitted rates, the predictions from
#' \code{\link{predict_apc}} and the DIC are invariant to it -- only the way the
#' common linear trend is split among the three curves changes. The convention
#' actually used is recorded in \code{attr(result, "gauge_convention")}.
#'
#' @return List of age, period, cohort effects or apc object including effects (if update=TRUE)
#' @export
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' effects(model)
#' }
effects.apc<-function(object, mean=FALSE, quantiles=0.5, update=FALSE,
                      convention=c("age","period","cohort","none"), combined=FALSE, ...)
{
  convention <- match.arg(convention)
  x<-object
  key <- list(mean=mean, quantiles=quantiles, convention=convention, combined=combined)
  # reuse a cached result only if it was computed with the same settings;
  # recompute otherwise. Either way the return must honour `update` below.
  if (!is.null(x$effects) && identical(attr(x$effects,"settings"), key)) {
    effects <- x$effects
  } else {
    g <- .apc_regauge(object, convention)
    # combined=TRUE: add the iid heterogeneity component (age2/period2/cohort2)
    # to the smooth curve. The het part is sum-to-zero and not involved in the
    # trend gauge, so it adds directly. No-op for models without het.
    if (isTRUE(combined)) g <- .apc_add_het(object, g)

    summ <- function(s) {
      z <- summary(s, quantiles=quantiles)
      if (mean) z$statistics[,1] else z$quantiles
    }
    age    <- summ(g$age)
    period <- summ(g$period)
    cohort <- summ(g$cohort)

    effects<-list("age"=age, "period"=period, "cohort"=cohort)
    attr(effects,"settings") <- key
    attr(effects,"gauge_convention") <- attr(g, "gauge_convention")
  }
  if (update){x$effects<-effects; return(x)}
  return(effects)
}

## Internal: add the iid heterogeneity component to each smooth effect, per draw
## per chain. The het slots (age2/period2/cohort2) exist only for "+het" models;
## where absent the smooth effect is returned unchanged. The het part is already
## sum-to-zero and gauge-independent, so it adds directly to the (re-gauged)
## smooth curve to give the full effect.
.apc_add_het <- function(object, g) {
  addh <- function(sm, het) {
    if (is.null(sm) || is.null(het)) return(sm)
    coda::as.mcmc.list(mapply(function(a, h) coda::mcmc(as.matrix(a) + as.matrix(h)),
                              sm, het, SIMPLIFY = FALSE))
  }
  g$age    <- addh(g$age,    object$samples[["age2"]])
  g$period <- addh(g$period, object$samples[["period2"]])
  g$cohort <- addh(g$cohort, object$samples[["cohort2"]])
  g
}

## Internal: re-gauge the age/period/cohort effect samples to a fixed display
## convention. In a full age-period-cohort model the three effects share one
## non-identified linear-trend (drift) direction; the transform
##   theta_i += f*M*i, phi_j -= f*j, psi_k += f*k   (mu -= f*I*M)
## leaves the fitted log-odds unchanged. For each MCMC draw f is chosen so that
## the chosen effect has zero least-squares slope, after which each effect is
## re-centred to sum to zero. The stored samples (object$samples) are never
## modified, so fitted rates, predictions and DIC are invariant. For any model
## that is not full APC (or for convention = "none") there is no trend aliasing
## and the raw centred samples are returned unchanged. Returns a list with the
## re-gauged age/period/cohort and (unchanged) intercept as coda::mcmc.list
## objects; attr(., "gauge_convention") is the convention actually applied.
.apc_regauge <- function(object, convention = c("age","period","cohort","none")) {
  convention <- match.arg(convention)
  m <- object$model
  # an effect is "present" only if its model slot is a single non-blank string;
  # bamp encodes an absent effect as " " (and NULL pre-conversion), but guard
  # against NA / empty / length-0 so the && below can never see NA or logical(0)
  pres <- function(z) length(z) == 1L && !is.na(z) && z != " "
  full <- pres(m[["age"]]) && pres(m[["period"]]) && pres(m[["cohort"]])
  # covariate models: the period/cohort effect is stored as the absolute
  # (covariate-scaled) contribution phi_j*x_j, so the linear-trend transform
  # (which assumes an additive phi_j) is no longer eta-preserving -- skip the
  # gauge and return the raw effects (as for non-full-APC models).
  has_cov <- !is.null(object[["covariate"]])
  eff <- if (!full || convention == "none" || has_cov) "none" else convention

  # exact [[ ]] extraction: object$samples also holds age_parameter/period_parameter/
  # cohort_parameter, and $ partial-matching could silently grab those if a primary
  # slot were ever missing.
  out <- list(age = object$samples[["age"]], period = object$samples[["period"]],
              cohort = object$samples[["cohort"]], intercept = object$samples[["intercept"]])
  if (eff == "none") { attr(out, "gauge_convention") <- "none"; return(out) }

  M <- object$data$periods_per_agegroup
  nch <- length(out$age)                # out$age/period/cohort: exact raw slots
  slope <- function(Mat) {            # per-sample least-squares slope on centred index
    idx <- seq_len(ncol(Mat)) - mean(seq_len(ncol(Mat)))
    as.numeric(Mat %*% idx) / sum(idx^2)
  }
  age2 <- period2 <- cohort2 <- vector("list", nch)
  for (ch in seq_len(nch)) {
    Th <- as.matrix(out$age[[ch]])
    Ph <- as.matrix(out$period[[ch]])
    Ps <- as.matrix(out$cohort[[ch]])
    I <- ncol(Th); J <- ncol(Ph); K <- ncol(Ps)
    f <- switch(eff,
                age    = -slope(Th) / M,
                period =  slope(Ph),
                cohort = -slope(Ps))
    Th <- Th + outer(f * M, seq_len(I))
    Ph <- Ph + outer(-f,    seq_len(J))
    Ps <- Ps + outer(f,     seq_len(K))
    Th <- Th - rowMeans(Th); Ph <- Ph - rowMeans(Ph); Ps <- Ps - rowMeans(Ps)
    age2[[ch]] <- coda::mcmc(Th); period2[[ch]] <- coda::mcmc(Ph); cohort2[[ch]] <- coda::mcmc(Ps)
  }
  out$age    <- coda::as.mcmc.list(age2)
  out$period <- coda::as.mcmc.list(period2)
  out$cohort <- coda::as.mcmc.list(cohort2)
  attr(out, "gauge_convention") <- eff
  out
}
