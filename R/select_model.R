#' Automatic model selection for age-period-cohort models
#'
#' @description
#' Searches over age-period-cohort model specifications and returns the one best
#' supported by the data, by Deviance Information Criterion (DIC). It answers the
#' practical questions "is a first- or second-order random walk more appropriate
#' for each effect?", "are the data overdispersed?" and (optionally) "is extra
#' heterogeneity warranted?" without the user fitting every model by hand.
#'
#' @details
#' The search is a \strong{greedy forward selection by complexity}. It starts
#' from the simplest model -- a first-order random walk (\code{"rw1"}) for every
#' effect that is present, no overdispersion -- and at each round considers every
#' candidate that is exactly one step more complex than the current best: an
#' effect upgraded from \code{"rw1"} to \code{"rw2"}, overdispersion switched on,
#' or (if \code{try_heterogeneity = TRUE}) heterogeneity added to an effect. All
#' candidates in a round are fitted and the one with the lowest DIC is adopted,
#' but \strong{only if} it improves DIC by at least \code{dic_margin} \emph{and}
#' its chains converged; otherwise the search stops. This costs a handful of fits
#' rather than the full grid, follows an interpretable path, and -- through the
#' margin -- prefers the simpler model unless the data clearly favour the more
#' complex one. Each distinct specification is fitted at most once (results are
#' cached).
#'
#' Model comparison uses DIC (lower is better), which rewards fit and penalises
#' effective complexity (\code{pD}); see \code{\link{bamp}}. A specification that
#' did not converge is never selected, because a low DIC from a chain that has
#' not mixed is not trustworthy -- convergence is judged on the identified fitted
#' values via the same criterion as \code{\link{checkConvergence}} (the maximum
#' Gelman-Rubin statistic of the fitted log-odds across Lexis cells must be at or
#' below \code{psrf_tol}).
#'
#' For speed and fairness all candidates are fitted with the same short
#' \code{screen} MCMC settings; the selected model is then optionally refitted
#' (\code{refit = TRUE}) with the longer \code{final} settings before being
#' returned. Fitting uses \code{method = "pg"}, which is robust on the sparse,
#' rare-event data where the legacy IWLS sampler can fail to converge.
#'
#' Pin an axis to exclude it from the search by passing a fixed value: e.g.
#' \code{age = "rw2"} fixes the age effect (it is not searched), \code{age = " "}
#' removes the age effect entirely, and \code{overdispersion = FALSE} forbids
#' overdispersion. Any axis left \code{NULL} is searched.
#'
#' @param cases number of cases (matrix, periods x age groups), as in \code{\link{bamp}}.
#' @param population population number, as in \code{\link{bamp}}.
#' @param periods_per_agegroup periods per age group.
#' @param age,period,cohort optional fixed value for an effect (\code{"rw1"},
#'   \code{"rw2"}, \code{"rw1+het"}, \code{"rw2+het"} or \code{" "} for absent).
#'   If \code{NULL} (default) the effect is present and its random-walk order is
#'   searched over \code{"rw1"}/\code{"rw2"}.
#' @param overdispersion optional fixed logical. If \code{NULL} (default),
#'   whether to include overdispersion is part of the search.
#' @param try_heterogeneity logical; if \code{TRUE} the search may also add
#'   heterogeneity (\code{"+het"}) to an effect. Default \code{FALSE}.
#' @param dic_margin minimum DIC improvement required to adopt a more complex
#'   model (parsimony threshold). Default 4 (a conventional "clearly better"
#'   DIC difference).
#' @param psrf_tol convergence tolerance: a fit counts as converged if its
#'   maximum fitted-value Gelman-Rubin statistic is at or below this. Default 1.1.
#' @param screen list of MCMC settings used for the comparison fits, kept
#'   moderate for speed (default 10000 iterations, 5000 burn-in, step 5). Passed
#'   as \code{mcmc.options}. If no candidate converges at this length the search
#'   warns and selects by DIC only; increase these settings and re-run.
#' @param final list of MCMC settings used to refit the selected model, or
#'   \code{"auto"} to use the data-adaptive default of \code{\link{bamp}}.
#' @param refit logical; if \code{TRUE} (default) refit the selected model with
#'   the \code{final} settings and return it; if \code{FALSE} return the
#'   screening fit of the selected model.
#' @param hyperpar hyper-parameter list passed to \code{\link{bamp}}; defaults
#'   include heterogeneity hyper-parameters so \code{"+het"} models can be fitted.
#' @param parallel passed to \code{\link{bamp}} (chains run in parallel).
#' @param verbose logical; if \code{TRUE} (default) report progress and the
#'   running best model.
#' @param ... further arguments passed to \code{\link{bamp}} (e.g.
#'   \code{prior_scale}, \code{pg_engine}).
#'
#' @return A list (class \code{"apcselect"}) with elements \itemize{
#'   \item \code{table}: a data frame of every specification fitted, with its DIC,
#'     effective number of parameters \code{pD}, mean deviance, convergence flag,
#'     maximum fitted-value PSRF and fit time, ordered by DIC.
#'   \item \code{best}: the selected specification (a named list).
#'   \item \code{model}: the fitted \code{\link{apc}} object for the selected
#'     specification (refitted with \code{final} settings if \code{refit = TRUE}).
#'   \item \code{path}: the sequence of specifications adopted by the greedy search.}
#'
#' @seealso \code{\link{bamp}}, \code{\link{checkConvergence}}
#' @export
#' @examples
#' \dontrun{
#' data(apc)
#' sel <- selectModel(cases, population, periods_per_agegroup = 5)
#' sel$table          # ranked comparison of the models tried
#' sel$best           # the chosen specification
#' plot(sel$model)    # the refitted best model
#' }
selectModel <- function(cases, population, periods_per_agegroup,
                        age = NULL, period = NULL, cohort = NULL,
                        overdispersion = NULL,
                        try_heterogeneity = FALSE,
                        dic_margin = 4,
                        psrf_tol = 1.1,
                        screen = list(number_of_iterations = 10000, burn_in = 5000,
                                      step = 5, tuning = 200),
                        final = "auto",
                        refit = TRUE,
                        hyperpar = list(age = c(1, 0.5), period = c(1, 5e-4),
                                        cohort = c(1, 5e-4), overdisp = c(1, 0.05),
                                        age_het = c(1, 0.05), period_het = c(1, 0.05),
                                        cohort_het = c(1, 0.05)),
                        parallel = TRUE, verbose = TRUE, ...) {

  ## --- which effects are present, and which axes are free to search ---------
  ## an effect is present unless pinned to " "; its order is searched unless the
  ## user pinned a concrete value.
  eff_present <- function(x) is.null(x) || !identical(x, " ")
  has <- c(age = eff_present(age), period = eff_present(period), cohort = eff_present(cohort))
  pinned <- list(age = age, period = period, cohort = cohort)
  free_order <- vapply(names(pinned), function(nm) has[[nm]] && is.null(pinned[[nm]]),
                       logical(1))
  free_od <- is.null(overdispersion)

  ## a "spec" is list(age, period, cohort, overdisp); absent effects are " "
  base_eff <- function(nm) if (!has[[nm]]) " " else if (!is.null(pinned[[nm]])) pinned[[nm]] else "rw1"
  spec0 <- list(age = base_eff("age"), period = base_eff("period"), cohort = base_eff("cohort"),
                overdisp = if (!free_od) isTRUE(overdispersion) else FALSE)

  spec_key <- function(s) paste(s$age, s$period, s$cohort, if (s$overdisp) "od" else "no-od", sep = "|")
  spec_label <- function(s) sprintf("age=%s period=%s cohort=%s%s",
                                    s$age, s$period, s$cohort, if (s$overdisp) " +overdisp" else "")

  ## --- one fit + its summary, cached ----------------------------------------
  cache <- new.env(parent = emptyenv())
  fit_spec <- function(s, mcmc) {
    key <- spec_key(s)
    if (!is.null(cache[[key]])) return(cache[[key]])
    t0 <- Sys.time()
    m <- tryCatch(
      bamp(cases, population, age = s$age, period = s$period, cohort = s$cohort,
           overdisp = s$overdisp, periods_per_agegroup = periods_per_agegroup,
           mcmc.options = mcmc, hyperpar = hyperpar, method = "pg",
           parallel = parallel, dic = TRUE, verbose = FALSE, ...),
      error = function(e) structure(list(err = conditionMessage(e)), class = "bampfail"))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(m, "bampfail") || is.null(m$deviance$DIC)) {
      res <- list(spec = s, dic = NA_real_, pD = NA_real_, mean_dev = NA_real_,
                  psrf = NA_real_, converged = FALSE, secs = secs, model = NULL)
    } else {
      psrf <- tryCatch(max(.apc_eta_psrf(m), na.rm = TRUE), error = function(e) NA_real_)
      res <- list(spec = s, dic = m$deviance$DIC, pD = m$deviance$pD,
                  mean_dev = m$deviance$mean.deviance, psrf = psrf,
                  converged = is.finite(psrf) && psrf <= psrf_tol, secs = secs, model = m)
    }
    cache[[key]] <- res
    if (verbose)
      cat(sprintf("  %-44s DIC=%-9s %s  [%.0fs]\n", spec_label(s),
                  if (is.na(res$dic)) "NA" else formatC(res$dic, format = "f", digits = 1),
                  if (res$converged) sprintf("converged (PSRF %.3f)", res$psrf)
                  else if (is.na(res$psrf)) "FIT FAILED" else sprintf("NOT converged (PSRF %.3f)", res$psrf),
                  res$secs))
    res
  }

  ## --- neighbours: each is exactly one step more complex than s -------------
  neighbours <- function(s) {
    out <- list()
    for (nm in c("age", "period", "cohort")) {
      if (!has[[nm]] || !free_order[[nm]]) next
      cur <- s[[nm]]
      ## rw1 -> rw2 (preserving any +het suffix)
      if (grepl("^rw1", cur)) { n <- s; n[[nm]] <- sub("^rw1", "rw2", cur); out <- c(out, list(n)) }
      ## + heterogeneity (rwX -> rwX+het)
      if (try_heterogeneity && !grepl("\\+het$", cur)) {
        n <- s; n[[nm]] <- paste0(cur, "+het"); out <- c(out, list(n))
      }
    }
    if (free_od && !s$overdisp) { n <- s; n$overdisp <- TRUE; out <- c(out, list(n)) }
    out
  }

  ## --- greedy forward selection ---------------------------------------------
  if (verbose) cat("Model search (greedy forward selection by DIC):\n")
  best <- fit_spec(spec0, screen)
  path <- list(best$spec)
  if (!best$converged && verbose)
    cat("  note: the baseline model did not converge at the screening length;\n",
        "        consider longer `screen` settings.\n", sep = "")
  repeat {
    cands <- neighbours(best$spec)
    if (!length(cands)) break
    fits <- lapply(cands, fit_spec, mcmc = screen)
    ok <- Filter(function(f) f$converged && is.finite(f$dic), fits)
    if (!length(ok)) break
    cand_best <- ok[[which.min(vapply(ok, function(f) f$dic, numeric(1)))]]
    ## adopt only if it beats the current best by the parsimony margin; if the
    ## current baseline itself did not converge, accept the first converged
    ## improvement regardless of margin.
    improved <- if (best$converged) (cand_best$dic < best$dic - dic_margin) else TRUE
    if (improved) {
      if (verbose) cat(sprintf("  -> adopt %s (dDIC = %.1f)\n",
                               spec_label(cand_best$spec),
                               if (is.finite(best$dic)) cand_best$dic - best$dic else NA))
      best <- cand_best; path <- c(path, list(best$spec))
    } else break
  }

  ## If nothing converged at the screening length, the greedy path is not
  ## trustworthy. Fall back to ranking ALL fitted specs by DIC purely for
  ## information, pick the lowest-DIC one, and warn loudly -- a non-converged DIC
  ## is not a sound basis for selection, so the screen length is the real fix.
  selection_converged <- best$converged
  if (!selection_converged) {
    finite_fits <- Filter(function(f) is.finite(f$dic), as.list(cache))
    if (length(finite_fits)) {
      best <- finite_fits[[which.min(vapply(finite_fits, function(f) f$dic, numeric(1)))]]
      path <- list(best$spec)
    }
    warning("No candidate model converged at the screening MCMC length; the ",
            "selection is by DIC only and is NOT convergence-backed. Increase ",
            "`screen` (more iterations) and re-run, and check the returned model ",
            "with checkConvergence().", call. = FALSE)
  }

  ## --- comparison table (everything fitted), ordered by DIC -----------------
  all_fits <- as.list(cache)
  tab <- do.call(rbind, lapply(all_fits, function(f) data.frame(
    age = f$spec$age, period = f$spec$period, cohort = f$spec$cohort,
    overdisp = f$spec$overdisp, DIC = f$dic, pD = f$pD, mean_deviance = f$mean_dev,
    converged = f$converged, max_psrf = f$psrf, secs = round(f$secs, 1),
    stringsAsFactors = FALSE, row.names = NULL)))
  tab <- tab[order(tab$DIC, na.last = TRUE), , drop = FALSE]
  rownames(tab) <- NULL
  ## flag the selected model
  tab$selected <- mapply(function(a, p, c0, od)
    a == best$spec$age && p == best$spec$period && c0 == best$spec$cohort && od == best$spec$overdisp,
    tab$age, tab$period, tab$cohort, tab$overdisp)

  if (verbose) {
    cat(sprintf("\nSelected: %s  (DIC %.1f)\n", spec_label(best$spec), best$dic))
  }

  ## --- refit the winner at the final length ---------------------------------
  final_model <- best$model
  if (refit) {
    fmcmc <- if (identical(final, "auto"))
      list(number_of_iterations = "auto", burn_in = "auto", step = "auto", tuning = 200) else final
    if (verbose) cat("Refitting the selected model at the final MCMC length ...\n")
    final_model <- bamp(cases, population, age = best$spec$age, period = best$spec$period,
                        cohort = best$spec$cohort, overdisp = best$spec$overdisp,
                        periods_per_agegroup = periods_per_agegroup, mcmc.options = fmcmc,
                        hyperpar = hyperpar, method = "pg", parallel = parallel,
                        dic = TRUE, verbose = FALSE, ...)
  }

  structure(list(table = tab, best = best$spec, model = final_model,
                 path = path, dic_margin = dic_margin, psrf_tol = psrf_tol,
                 converged = selection_converged),
            class = "apcselect")
}

#' @export
print.apcselect <- function(x, ...) {
  cat("Automatic APC model selection (by DIC; lower is better)\n\n")
  show <- data.frame(
    age = x$table$age, period = x$table$period, cohort = x$table$cohort,
    overdisp = ifelse(x$table$overdisp, "yes", "no"),
    DIC = round(x$table$DIC, 1), pD = round(x$table$pD, 1),
    converged = ifelse(x$table$converged, "yes", "no"),
    max_psrf = round(x$table$max_psrf, 3),
    sel = ifelse(x$table$selected, "<--", ""),
    stringsAsFactors = FALSE)
  print(show, row.names = FALSE)
  b <- x$best
  cat(sprintf("\nSelected model (<--): age=%s period=%s cohort=%s%s\n",
              b$age, b$period, b$cohort, if (b$overdisp) " +overdisp" else ""))
  if (isFALSE(x$converged))
    cat("WARNING: no model converged at the screening length -- this selection is\n",
        "by DIC only and is not convergence-backed. Re-run with longer `screen`\n",
        "settings and verify with checkConvergence().\n", sep = "")
  cat("Access the fitted model with $model (e.g. plot(result$model)).\n")
  invisible(x)
}
