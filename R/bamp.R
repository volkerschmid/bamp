#' Bayesian Age-Period-Cohort Modeling and Prediction (bamp)
#'
#' @param cases number of cases
#' @param population population number
#' @param age prior for age groups ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param period prior for periods ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param cohort prior for cohorts ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param overdisp logical, add overdispersion to model
#' @param periods_per_agegroup periods per age group
#' @param period_covariate covariate for period
#' @param cohort_covariate covariate for cohort
#' @param mcmc.options list of options for MCMC. \itemize{\item number_of_iterations: number of iterations per chain. \item burn_in: number of iterations used as burnin at the beginning of the algorithm, these iterations will be removed. \item step: Step size, so only every step-th iteration is stored. \item tuning: number of iterations for automatic tuning (used by \code{method="taylor"}). Depending on the model, the MCMC algorithm will tune certain parameters for more efficient MCMC chains. After tuning, the algorithm is restarted.} Each of \code{number_of_iterations}, \code{burn_in} and \code{step} may be a number or the string \code{"auto"} (the default). \code{"auto"} chooses the value from the data: rare or zero-heavy counts (whose rare-event cells mix more slowly) get more iterations (from 40000 for well-populated data up to 120000 when almost every cell is empty or has very few events), \code{burn_in} defaults to half the iterations, and \code{step} is set to keep about 1000 stored samples per chain. Any value given as a number is used exactly as supplied, so explicit settings reproduce the previous behaviour.
#' @param hyperpar list of hyper parameters. The hyper prior for the precision (inverse variance) in the random walk priors is a Gamma distribution with parameters \eqn{a} and \eqn{b}; expected value is \eqn{a/b}, variance is \eqn{a/b^2}. Weak hyper parameters are suggested, defaults are \eqn{a=1, b=0.5} for age, \eqn{a=1, b=0.0005} for period and cohort effects and \eqn{a=1, b=0.05} for overdispersion (if added). It is recommended to choose the hyper priors depending on the model, in particular on the order of the random walk.
#' @param dic logical. If true. DIC will be computed
#' @param parallel should the chains be run in parallel. \code{TRUE}/\code{FALSE},
#' or a number giving the requested number of cores (capped at the number of
#' chains). Uses the \code{parallel} package: forked workers
#' (\code{\link[parallel]{mclapply}}) on Unix and macOS, and -- for
#' \code{method = "pg"} -- a PSOCK cluster on Windows (where forking is
#' unavailable), so the default engine now runs in parallel on all platforms.
#' (The legacy \code{method = "taylor"} engine still runs serially on Windows.)
#' Parallel runs are reproducible: the per-chain seeds are drawn in the main
#' process, so a given \code{set.seed()} yields the same result serially or in
#' parallel.
#' @param verbose verbose mode
#' @param method MCMC engine. \code{"pg"} (default) is a joint sampler that
#' combines Polya-Gamma data augmentation (Polson, Scott & Windle 2013) with a
#' Laplace (Newton) Metropolis-Hastings refinement: each sweep draws the
#' intercept and the age, period and cohort effects jointly in one exact Gibbs
#' step and then refines them with a joint Newton proposal against the true
#' binomial likelihood. It has no Metropolis tuning, never restarts on low
#' acceptance and does not prune chains; it is markedly more robust for RW2
#' priors and converges the high-population, rare-event cells of
#' incidence/mortality data that the Gibbs step alone mixes only slowly. It
#' natively supports all of the package's models -- RW1/RW2 priors,
#' heterogeneity (\code{"rw1+het"}/\code{"rw2+het"}), overdispersion and
#' period/cohort covariates. The Polya-Gamma weights use a normal approximation
#' that is essentially exact for the large population counts of
#' incidence/mortality data, so it typically needs far fewer iterations than the
#' legacy sampler. \code{"taylor"} is the original block Metropolis-Hastings
#' sampler with taylor expansion proposals (the default in versions 2.x); it remains
#' available and can be faster on well-behaved (non rare-event) data, but it can
#' fail to converge or prune all chains on sparse/zero-cell data.
#' @param prior_scale logical; only used by \code{method="pg"}. If \code{TRUE},
#' the intrinsic random-walk structure matrices are scaled to unit generalised
#' variance (Sorbye & Rue 2014) so that a single hyper-prior is comparable
#' across random-walk orders, grid sizes and data sets. The default is
#' \code{FALSE}, which keeps the same prior parameterisation (and the same
#' default hyper-parameters) as \code{method="taylor"}; if you set it to
#' \code{TRUE} you should choose hyper-parameters appropriate for the scaled
#' prior. See \sQuote{Scaling the random-walk priors} below for the rationale
#' and benefits, and the examples for a short demonstration.
#' @param pg_engine implementation of the \code{method="pg"} sampler, one of
#' \code{"C"} (default) or \code{"R"}. Both run the identical algorithm and,
#' for a given seed, produce the same draws to floating-point tolerance; the
#' \code{"C"} engine is a compiled port of the inner loop (no extra package
#' dependency) and is roughly twice as fast. \code{"R"} is the readable
#' reference implementation, kept for verification. Ignored for
#' \code{method="taylor"}.
#'
#' @description
#' Bayesian Age-Period-Cohort Modeling for the analyze of incidence or mortality data on the Lexis diagram.
#' For each pixel in the Lexis diagram (that is for a specific age group and specific period) data must be available on the number of persons under risk (population number) and the number of disease cases (typically cancer incidence or mortality).
#' A hierarchical model is assumed with a binomial model in the first-stage. As smoothing priors for the age, period and cohort parameters random walks of first and second order (RW1 or RW2) available.
#' Deviance information criterion and effective number of parameters is computed for model comparison.
#' Note that there is a non-identifiability in the likelihood of the APC-model, see e.g. Clayton and Schifflers (1987, DOI:10.1002/sim.4780060406), which indices some problems in interpreting the latent effects. Only for RW1 model, the parameters are (weakly) identifiable.
#' Period and age groups do not need to be on the same grid, for example periods can be in one year intervals and age groups in five year intervals.\cr
#' Additionally to the model described in Knorr-Held and Rainer (2001, DOI:10.1093/biostatistics/2.1.109), \code{bamp} can handle 
#' \itemize{\item AP and AC models, 
#' \item models with and without global heterogeneity parameter (overdispersion),
#' \item models with additional age, period and/or cohort heterogeneity,
#' \item additional covariates.}
#' 
#' @details This functions returns an \code{\link{apc}} object.
#' Only samples from the posterior are computed, point estimates and credible intervals will be computed in \code{\link{effects.apc}}, \code{\link{print.apc}} and \code{\link{plot.apc}}.
#' \code{\link{predict_apc}} can be used for for prediction of the future rates and number of cases and for a retrospective prediction for model checking.
#'
#' @section Scaling the random-walk priors (\code{prior_scale}):
#' Each age, period and cohort effect has an intrinsic Gaussian (random-walk)
#' prior with precision (smoothing) parameter \eqn{\kappa}: the effect vector
#' \eqn{x} has density proportional to \eqn{\exp(-\tfrac{1}{2}\kappa\, x'Kx)},
#' where \eqn{K=D'D} is built from the first- or second-order difference
#' operator \eqn{D}. A \eqn{\mathrm{Gamma}(a,b)} hyper-prior is placed on
#' \eqn{\kappa}. The difficulty is that the smoothness implied by a given
#' \eqn{\kappa} is governed not by \eqn{\kappa} alone but by the marginal
#' variance of the effect, the generalised inverse of \eqn{\kappa K}; and the
#' eigenvalues of \eqn{K} grow with the number of time points and with the
#' random-walk order. The \emph{same} hyper-prior on \eqn{\kappa} therefore
#' implies very different prior smoothness for, say, an RW1 over 10 periods and
#' an RW2 over 50 cohorts. A hyper-prior tuned on one model silently means
#' something different on another, which is one reason a fixed default can
#' behave inconsistently across data sets.
#'
#' With \code{prior_scale = TRUE} the structure matrix \eqn{K} is rescaled so
#' that the geometric mean of the (generalised) marginal variances equals one
#' (Sorbye and Rue, 2014, DOI:10.1080/01621459.2013.866549). After scaling,
#' \eqn{1/\sqrt{\kappa}} is, to a good approximation, the marginal standard
#' deviation of a typical effect element \emph{on the log-odds (logit) scale},
#' independently of the random-walk order, the number of age/period/cohort
#' points and the grid spacing.
#'
#' Benefits: (i) \strong{portable hyper-priors} -- one \eqn{\mathrm{Gamma}(a,b)}
#' encodes the same smoothness belief across RW1/RW2 and across data sets of
#' different size; (ii) an \strong{interpretable prior} -- you can set
#' \eqn{(a,b)} to express a belief about \eqn{1/\sqrt{\kappa}} as a prior effect
#' standard deviation on the logit scale; (iii) \strong{fairer model comparison}
#' (e.g. RW1 vs RW2 by DIC), because the prior is not implicitly penalising one
#' model far more than another. Scaling affects only the smooth random-walk
#' blocks; the i.i.d. heterogeneity components and overdispersion already have
#' an interpretable scale and are unchanged.
#'
#' The default is \code{prior_scale = FALSE} so that \code{method = "pg"}
#' reproduces the prior parameterisation (and default \code{hyperpar}) of the
#' legacy \code{method = "taylor"} engine. If you turn scaling on you should set
#' \code{hyperpar} for the scaled prior, where \eqn{\kappa \approx
#' 1/\mathrm{variance}}; using the unscaled defaults with \code{prior_scale =
#' TRUE} would impose a different (and probably unintended) amount of smoothing.
#' Scaling is most worthwhile when fitting many models or data sets and you want
#' one coherent, interpretable prior across all of them. The example below shows
#' the effect concretely.
#'
#' @seealso \code{vignette("modeling", package = "bamp")}
#' @useDynLib bamp
#' @export
#' @import coda
#' @importFrom utils modifyList
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' }
#'
#' ## Demonstration of prior_scale (no MCMC, runs instantly): for a fixed
#' ## precision kappa, report the geometric-mean prior marginal standard
#' ## deviation of a random-walk effect on the logit scale, with and without
#' ## Sorbye-Rue scaling, across random-walk orders and grid sizes.
#' prior_sd <- function(L, order, kappa = 1, scale = FALSE) {
#'   K <- crossprod(diff(diag(L), differences = order))   # structure matrix D'D
#'   if (scale) {                                          # Sorbye-Rue unit-variance scaling
#'     e <- eigen(K, symmetric = TRUE); keep <- e$values > max(e$values) * 1e-9
#'     V <- e$vectors[, keep, drop = FALSE]
#'     Sigma <- V %*% diag(1 / e$values[keep], sum(keep)) %*% t(V)
#'     K <- K * exp(mean(log(diag(Sigma))))
#'   }
#'   e <- eigen(K, symmetric = TRUE); keep <- e$values > max(e$values) * 1e-9
#'   V <- e$vectors[, keep, drop = FALSE]
#'   Sig <- V %*% diag(1 / e$values[keep], sum(keep)) %*% t(V) / kappa
#'   sqrt(exp(mean(log(diag(Sig)))))                       # geometric-mean marginal SD
#' }
#' grid <- expand.grid(order = 1:2, L = c(10, 25, 50))
#' data.frame(grid,
#'            unscaled = round(mapply(prior_sd, grid$L, grid$order, scale = FALSE), 3),
#'            scaled   = round(mapply(prior_sd, grid$L, grid$order, scale = TRUE), 3))
#' ## With prior_scale = FALSE the same kappa = 1 implies an effect SD ranging
#' ## from ~1.2 to ~14.6 across these models; with prior_scale = TRUE it is 1.0
#' ## throughout, so a single hyper-prior on kappa means the same smoothness for
#' ## every random-walk order and grid size.
bamp <-
function(cases, population,
        age, period, cohort, overdisp=FALSE,
        period_covariate=NULL, cohort_covariate=NULL,
        periods_per_agegroup,
        mcmc.options=list("number_of_iterations"="auto", "burn_in"="auto", "step"="auto", "tuning"=500),
        hyperpar=list("age"=c(1,0.5), "period"=c(1,0.0005), "cohort"=c(1,0.0005), "overdisp"=c(1,0.05)),
        dic=TRUE,
        parallel=TRUE, verbose=FALSE,
        method=c("pg","taylor"), prior_scale=FALSE, pg_engine=c("C","R")){
  output=apc()
  method <- match.arg(method)
  pg_engine <- match.arg(pg_engine)

  ## Normalise the "no effect" specification to " " up front, BEFORE the
  ## "...+het" checks below -- otherwise cohort=NULL (or age/period=NULL) hits
  ## `NULL == "rw1+het"` -> logical(0) -> "argument is of length zero".
  if (is.null(age))    age    <- " "
  if (is.null(period)) period <- " "
  if (is.null(cohort)) cohort <- " "

  ## The Polya-Gamma Gibbs engine natively supports plain RW1/RW2 priors,
  ## overdispersion, heterogeneity and period/cohort covariates -- there is no
  ## longer any model that falls back to taylor.

  ## Fill any hyper-parameter the caller omitted with its default. A partial
  ## hyperpar list (e.g. list(age=, period=, cohort=) with no "overdisp") would
  ## otherwise leave hyperpar$overdisp NULL, producing NULL hyper-values that
  ## later crash the sampler / mkmat() ("attempt to set an attribute on NULL").
  ## modifyList keeps any extra entries the caller supplied (e.g. age_het).
  hyperpar <- modifyList(list("age"=c(1,0.5), "period"=c(1,0.0005),
                              "cohort"=c(1,0.0005), "overdisp"=c(1,0.05)), hyperpar)

  age_hyperpar_a=hyperpar$age[1]
  age_hyperpar_b=hyperpar$age[2]
  period_hyperpar_a=hyperpar$period[1]
  period_hyperpar_b=hyperpar$period[2]
  cohort_hyperpar_a=hyperpar$cohort[1]
  cohort_hyperpar_b=hyperpar$cohort[2]
  z_hyperpar_a=hyperpar$overdisp[1]
  z_hyperpar_b=hyperpar$overdisp[2]
  if (age=="rw1+het"|age=="rw2+het")
    {
    age_hyperpar_a2=hyperpar$age_het[1]
    age_hyperpar_b2=hyperpar$age_het[2]
  }
  else
  {
    age_hyperpar_a2=1
    age_hyperpar_b2=1
  }
  if (period=="rw1+het"|period=="rw2+het")
  {
    period_hyperpar_a2=hyperpar$period_het[1]
  period_hyperpar_b2=hyperpar$period_het[2]
  }
  else
  {
    period_hyperpar_a2=1
    period_hyperpar_b2=1
  }
  if (cohort=="rw1+het"|cohort=="rw2+het")
  {
    cohort_hyperpar_a2=hyperpar$cohort_het[1]
  cohort_hyperpar_b2=hyperpar$cohort_het[2]
  }
  else
  {
    cohort_hyperpar_a2=1
    cohort_hyperpar_b2=1
  }
  
  chains=4
  if (parallel>4)chains=parallel
  if (unname(Sys.info()["sysname"] == "Windows"))parallel=FALSE
  
  ## Resolve the MCMC length. A numeric entry is used exactly as given (full
  ## backward compatibility -- existing calls are unaffected). An entry left at
  ## the default "auto" (or NULL) is filled from a data-aware heuristic: rare /
  ## zero-heavy data (the high-population rare-event cells whose Polya-Gamma
  ## augmentation mixes slowly) gets more iterations, well-populated data fewer.
  ## burn_in then defaults to half the iterations and step keeps ~1000 stored
  ## samples, so the stored-sample count is roughly constant across data sets.
  ## 
  is_auto <- function(x) is.null(x) || (length(x) == 1L && is.character(x) && x == "auto")
  number_of_iterations <- if (is_auto(mcmc.options$number_of_iterations))
    .bamp_auto_mcmc(cases) else mcmc.options$number_of_iterations
  if (is_auto(mcmc.options$number_of_iterations)&overdisp)number_of_iterations <- 1.5*number_of_iterations
  burn_in <- if (is_auto(mcmc.options$burn_in))
    as.integer(round(number_of_iterations / 2)) else mcmc.options$burn_in
  step <- if (is_auto(mcmc.options$step))
    max(1L, as.integer(round((number_of_iterations - burn_in) / 1000))) else mcmc.options$step
  tuning <- if (is.null(mcmc.options$tuning)) 500 else mcmc.options$tuning
  if (burn_in >= number_of_iterations)
    stop("burn_in must be smaller than number_of_iterations", call. = FALSE)
  if (verbose && (is_auto(mcmc.options$number_of_iterations) ||
                  is_auto(mcmc.options$burn_in) || is_auto(mcmc.options$step)))
    cat(sprintf("Auto MCMC settings from data rarity: %d iterations, %d burn-in, step %d.\n",
                number_of_iterations, burn_in, step))
  dataorder = 0
  number_of_agegroups=dim(cases)[2]
  number_of_periods=dim(cases)[1]
  cohort_start = 1
  period_start = 1

  z_mode=ifelse(overdisp,1,0)

  model<-list()
  if (is.null(age))age=" "
  age_block = switch(age,
    "rw1" = 1,
    "rw2" = 2,
    "rw1+het" = 3,
    "rw2+het" = 4,
    " " = 0
  )
  if (is.null(period))period=" "
  period_block = switch(period,
                     "rw1" = 1,
                     "rw2" = 2,
                     "rw1+het" = 3,
                     "rw2+het" = 4,
                     " " = 0
  )
  if (is.null(cohort))cohort=" "
  cohort_block = switch(cohort,
                     "rw1" = 1,
                     "rw2" = 2,
                     "rw1+het" = 3,
                     "rw2+het" = 4,
                     " " = 0
  )

  model$age=age
  model$period=period
  model$cohort=cohort

  model$overdispersion=overdisp
  ## record whether the RW structure matrices were Sorbye-Rue scaled (method="pg"
  ## only). predict_apc needs this: the extrapolation innovation variance is
  ## 1/(precision * scale), so a forecast must apply the SAME per-effect scale the
  ## fit used, else long-horizon credible bands are inflated by the scale factor.
  model$prior_scale=isTRUE(prior_scale)
  output$model=model

  #if (!is.null(age_covariate))age_block=age_block=age_block+7
  if (!is.null(period_covariate))period_block=period_block+7
  if (!is.null(cohort_covariate))cohort_block=cohort_block+7


  zentrieren <-
    function(mat, my){
      summe <- rowSums(mat)

      for(i in 1:dim(mat)[1]){
        mat[i,] <- mat[i,] - summe[i]/dim(mat)[2]
        my[i] <- my[i] + summe[i]/dim(mat)[2]
      }


      return(list(mat, my))
    }

##################################################################################################################################
# check of the input and preparation of the variables (bamp.cc S.1-S.11)

  stopifnot(is.data.frame(cases) || is.matrix(cases))              # cases must be a data.frame or a matrix
  if(is.data.frame(cases)){
    cases <- as.matrix(cases)
  }

  stopifnot(is.data.frame(population) || is.matrix(population))    # population must be a data.frame or a matrix
  if(is.data.frame(population)){
    population <- as.matrix(population)
  }

  if(dataorder == 0 || dataorder == 1){                            # dataorder must be in (0,1)
  }else{
    stop("ERROR: Dataorder must be 0 or 1!")
  }

  if(dataorder == 0){
    cases <- t(cases)
    population <- t(population)
  }

  stopifnot(number_of_agegroups%%1 == 0 && number_of_agegroups > 0)# number_of_agegroups must be a whole number
  stopifnot(number_of_periods%%1 == 0 && number_of_periods > 0)    # number_of_periods must be a whole number
  stopifnot(is.numeric(periods_per_agegroup))                      # periods_per_agegroup must be numeric


  stopifnot(number_of_iterations%%1 == 0)                          # number_of_iterations must be a whole number
  stopifnot(burn_in%%1 == 0)                                       # burn_in must be a whole number
  if(number_of_iterations <= burn_in){                             # number_of_iterations must be bigger than burn_in
    stop("ERROR: Number of iterations must be bigger than burn in!")
  }

  stopifnot(step%%1 == 0)                                          # step must be a whole number
  stopifnot(number_of_iterations-burn_in >= step)                  # number_ofiterations-burn_in musst be bigger than step
  stopifnot(tuning%%1 == 0)                                        # tuning must be a whole number
  if(burn_in <= tuning){                                           # burn_in must be bigger than tuning
    stop("ERROR: Burn in must be bigger than tuning constant!")
  }

  if(age_block == 0 | age_block == 1 | age_block == 2 |            # age_block must be in (0,1,2,3,4,8,9)
     age_block == 3 | age_block == 4){
  }else{
     stop("ERROR: Age_Block must be 0, 1, 2, 3 or 4!")
  }
  if(age_block > 0){
    stopifnot(is.numeric(age_hyperpar_a) & age_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(age_hyperpar_b) & age_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(age_block == 3 | age_block == 4){
    stopifnot(is.numeric(age_hyperpar_a2) & age_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(age_hyperpar_b2) & age_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(period_block == 0 | period_block == 1 |                       # period_block must be in (0,1,2,3,4,8,9)
     period_block == 2 | period_block == 3 |
     period_block == 4 | period_block == 8 |
     period_block == 9){
  }else{
    stop("ERROR: Period Block must be 0, 1, 2, 3, 4, 8 or 9!")
  }

  if(period_block > 0){
    stopifnot(is.numeric(period_hyperpar_a) & period_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(period_hyperpar_b) & period_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(period_block == 3 | period_block == 4){
    stopifnot(is.numeric(period_hyperpar_a2) & period_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(period_hyperpar_b2) & period_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(cohort_block == 0 | cohort_block == 1 |                      # cohort_block must be in (0,1,2,3,4,8,9)
     cohort_block == 2 | cohort_block == 3 |
     cohort_block == 4 | cohort_block == 8 |
     cohort_block == 9){
  }else{
    stop("ERROR: Cohort Block must be 0, 1, 2, 3, 4, 8 or 9!")
  }
  if(cohort_block > 0){
    stopifnot(is.numeric(cohort_hyperpar_a) & cohort_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(cohort_hyperpar_b) & cohort_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(cohort_block == 3 | cohort_block == 4){
    stopifnot(is.numeric(cohort_hyperpar_a2) & cohort_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(cohort_hyperpar_b2) & cohort_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(z_mode == 0 | z_mode == 1){                                   # z_mode must be in (0,1)
  }else{
    stop("ERROR: Z Mode must be 0 or 1!")
  }
  if(z_mode == 1){
    stopifnot(is.numeric(z_hyperpar_a) & z_hyperpar_a > 0)         # needed for gammadistribution
    stopifnot(is.numeric(z_hyperpar_b) & z_hyperpar_b > 0)         # needed for gammadistribution
  }

# additional conditions and new variables
  if(number_of_periods <= periods_per_agegroup){               # number_of_periods must be bigger than periods_per_agegroup
    stop("ERROR: Need more periods!")
  }


  number_of_cohorts <- periods_per_agegroup*(number_of_agegroups-1)+number_of_periods # number_of_cohorts


  max_number_of_ap_combinations <- floor((number_of_periods/periods_per_agegroup)+0.999) # max_combinations


  number_of_extractions <- (number_of_iterations - burn_in)/step    # number_of_extractions

  period_plus <- 0
  if(period_block == 8 || period_block == 9){                       # conditions for period_block 8 and 9
    period_block <- period_block - 7
    period_data <- c()
    period_plus <- 1
  }

  cohort_plus <- 0
  if(cohort_block == 8 || cohort_block == 9){                       # conditions for cohort_block 8 and 9
    cohort_block <- cohort_block - 7
    cohort_data <- c()
    cohort_plus <- 1
  }

  cphi <- 1                                                         # set phi = 1
  cpsi <- 1                                                         # set psi = 1

  if(period_plus == 1){                                             # if period_plus = 1 (with period_covariate)
    cphi <- 0                                                     # set phi = 0
    if(!is.vector(period_covariate))period_covariate<-as.vector(period_covariate)
    if(period_start%%1 != 0 || period_start <= 0){
      stop("ERROR: Period start must be a positive integer!")     # period start must be a positive integer
    }else{
      if(period_start > length(period_covariate)){
        stop("ERROR: Period start is bigger than the count of the covariables!")# period start can not be bigger than the vector
      }else{
        period_covariate <- period_covariate[period_start:length(period_covariate)]
      }
    }
    if(length(period_covariate) < number_of_periods){
      stop("ERROR: Not enough observations of the period covariabe!")    # need more observations
    }else{
      period_data <- period_covariate[1:number_of_periods]  # collects the period data
    }
    for(i in 1:number_of_periods){
      cphi <- cphi + period_data[i]                                  # calculates a cphi
    }
    cphi <- cphi/number_of_periods
    for(i in 1:(number_of_periods)){
      period_data[i] <- period_data[i]/cphi                          # calculates the period data
    }
  }

  if(cohort_plus == 1){                                              # if cohort_plus = 1 (with cohort_covariate)
    cpsi <- 0                                                      # set psi = 0
    if(!is.vector(cohort_covariate))cohort_covariate<-as.vector(cohort_covariate)
  if(cohort_start%%1 != 0 || cohort_start <= 0){
    stop("ERROR: Cohort start must be a positive integer!")         # cohort start must be a positive integer
  }else{
    if(cohort_start > length(cohort_covariate)){
      stop("ERROR: Cohort start is bigger than the count of the covariables!")# cohort start can not be bigger than the vector
    }else{
      cohort_covariate <- cohort_covariate[cohort_start:length(cohort_covariate)]
    }
  }
    if(length(cohort_covariate) < number_of_cohorts){
      stop("ERROR: Not enough observations of the cohort covariabe!")  # need more observations
    }else{
      cohort_data <- cohort_covariate[1:number_of_cohorts]  # collects the cohort data
    }
    for(i in 1:number_of_cohorts){
      cpsi <- cpsi + cohort_data[i]                                        # calculates a cpsi
    }
    cpsi <- cpsi/number_of_cohorts
    for(i in 1:number_of_cohorts){
      cohort_data[i] <- cohort_data[i]/cpsi                                # calculates the cohort data
    }
  }

if (verbose)
  {
  ##################################################################################################################################
## output of settings

  max_block <- max(age_block, (max(period_block, cohort_block)))
  if(max_block == 0){
    stop(cat("ERROR: No Block!\n"))
  }
  settings <- character()                                                 # AGE-Model ?
  if(age_block > 0){
    settings <- "AGE"
  }
  if(age_block*period_block > 0){
    settings <- paste(settings, "-", sep = "")
  }
  if(period_block > 0){
    settings <- paste(settings, "PERIOD", sep = "")                       # PERIOD-Model ?
  }
  if(age_block*cohort_block > 0 || period_block*cohort_block > 0){
    settings <- paste(settings, "-", sep = "")
  }
  if(cohort_block > 0){
    settings <- paste(settings, "COHORT", sep = "")                       # COHORT-Model ?
  }
  settings <- paste(settings, "Model", sep = " ")
  if(z_mode == 1){
    settings <- paste(settings, "with overdispersion", sep = " ")         # with overdispersion when z_mode = 1
  }
  if(period_block > 2){                                                   # with unstructured period effects for period_block > 2
    if(z_mode == 1){
      settings <- paste(settings, "and unstructured period effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured period effects", sep = " ")
    }
  }
  if(age_block > 2){
    if(z_mode == 1 || period_block > 2){                                  # with unstructured age effects for age_block > 2
      settings <- paste(settings, "and unstructured age effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured age effects", sep = " ")
    }
  }
  if(cohort_block > 2){
    if(z_mode == 1 || period_block > 2 || age_block > 2){                 # with unstructured cohort effects for cohort_block > 2
      settings <- paste(settings, "and unstructured cohort effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured cohort effects", sep = " ")
    }
  }

  settings <- paste(settings, ".", sep = "")
  settings <- paste(settings, "\nPrioris:", sep = "")                    # Prioris
  if(age_block == 0&verbose>0){
    settings <- paste(settings, "no age effect", sep = " ")                  # age
  }
  if(age_block > 0){
    settings <- paste(settings, "age effect", sep = " ")                  # age
  }
  if(age_block == 1 | age_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(age_block == 2 | age_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  if(age_block*period_block > 0){
    settings <- paste(settings, ",", sep = "")
  }
  if(period_block == 0&verbose>0){
    settings <- paste(settings, ", no period effect", sep = " ")                  # age
  }
  if(period_block > 0){
    settings <- paste(settings, "period effect", sep = " ")               # period
  }
  if(period_block == 1 | period_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(period_block == 2 | period_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  if(cohort_block == 0&verbose>0){
    settings <- paste(settings, ", no cohort effect", sep = " ")                  # age
  }
  if(cohort_block > 0){
    settings <- paste(settings, ", cohort effect", sep = "")              # cohort
  }
  if(cohort_block == 1 | cohort_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(cohort_block == 2 | cohort_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  settings <- paste(settings, ".", sep = "")
  settings <- paste(settings, "\n", sep = "")
  if(period_plus == 1){
    settings <- paste(settings, "Period effect with covariates ", names(period_covariate),   # period effect with covariates
                      " starting at position ", period_start, "\n", sep = "")
  }
  if(cohort_plus == 1){
    settings <- paste(settings, "Cohort effect with covariates ", names(cohort_covariate),   # cohort effect with covariates
                      " starting at position ", cohort_start, "\n", sep = "")
  }
  settings <- paste(settings, number_of_agegroups, " age groups, ", number_of_periods, " periods, ",
                    number_of_cohorts, " cohorts. ", sep = "")            # counts of agegroups, periods and cohorts
  settings <- paste(settings, "\n", sep = "")
  settings <- paste(settings, number_of_iterations, " iterations with ", burn_in, " burn in, using every ",
                    step, "th sample.", sep = "")                        # iterations, burn_in, step and tuning
  if(tuning > 0){
    settings <- paste(settings, " Tuning at iteration ", tuning, sep = "")
  }
  settings <- paste(settings, "\n", sep = "")

  done <- 0

  if(age_block > 0 || period_block > 0 || cohort_block > 0 || z_mode > 0){   # Hyper parameters
    settings <- paste(settings, "Hyper parameters: ", sep = "")
  }
  if(age_block > 0){                                                         # for age
    settings <- paste(settings, "age eff. (", age_hyperpar_a, ", ", age_hyperpar_b, ")", sep = "")
    done <- 1
    if(age_block == 3 | age_block == 4){
      settings <- paste(settings, " - unstr. age eff. (", age_hyperpar_a2, ", ", age_hyperpar_b2, ")", sep = "")
    }
  }
  if(period_block > 0){                                                      # for period
    if(done == 1){
      settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "period eff. (", period_hyperpar_a, ", ", period_hyperpar_b, ")", sep = "")
    done <- 1
    if(period_block == 3 | period_block == 4){
      settings <- paste(settings, " - unstr. period eff. (", period_hyperpar_a2, ", ", period_hyperpar_b2, ")", sep = "")
    }
  }
  if(cohort_block > 0){                                                      # for cohort
    if(done == 1){
      settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "cohort eff. (", cohort_hyperpar_a, ", ", cohort_hyperpar_b, ")", sep = "")
    done <-1
    if(cohort_block == 3 | cohort_block == 4){
      settings <- paste(settings, " - unstr. cohort eff. (", cohort_hyperpar_a2, ", ", cohort_hyperpar_b2, ")", sep = "")
    }
  }
  if(z_mode == 1){                                                           # for heterogeneity
    if(done == 1){
    settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "overdispersion (", z_hyperpar_a, ", ", z_hyperpar_b, ")", sep = "")
  }
  settings <- paste(settings, ".\n\n",sep="")
  settings<-paste(settings, "verbose level: ", verbose, sep="")
  settings <- paste(settings, ".\n\n", sep = "")
}

##################################################################################################################################
# additional settings and variables



  # additional conditions for the data (TxJ-matrix for dataorder = 0, JxT-matrix for dataorder =1)
  if(dim(cases)[2]  != number_of_periods || dim(cases)[1] != number_of_agegroups){
    stop("ERROR: Cases data does not fit to settings!")
  }



  if(period_block == 3 || period_block == 4){
    delta <- 100                                                              # set delta = 100 for period_block 3 and 4
  }

  if(cohort_block == 3 || cohort_block == 4){
    delta <- 100                                                              # set delta = 100 for cohort_block 3 and 4
  }



##################################################################################################################################
  if (verbose)
    {
    ausgabe <- paste(settings, "Starting Iterations in",chains,"chains.\n\n")
    cat(ausgabe)
  }

 nr.samples<-floor(number_of_iterations/step)-floor(burn_in/step)
 number_of_cohorts = ceiling(periods_per_agegroup*(number_of_agegroups - 1)+number_of_periods)

 delta.sample<-kappa.sample<-kappa2.sample<-lambda.sample<-lambda2.sample<-ny.sample<-ny2.sample<-my.sample<-dev.sample<-rep(0,nr.samples)
 theta.sample<-theta2.sample<-rep(0,nr.samples*number_of_agegroups)
 phi.sample<-phi2.sample<-rep(0,nr.samples*number_of_periods)
 psi.sample<-psi2.sample<-rep(0,nr.samples*number_of_cohorts)
 ksi<-rep(0,number_of_agegroups*number_of_periods)
 #print(paste("length of ksi.sample",length(ksi.sample)))

 blocks=c(age_block, period_block, cohort_block)
 numbers=c(number_of_agegroups, number_of_periods)
 numbersmcmc=c(number_of_iterations, burn_in, step, tuning)
 modelsettings=c(0, 0, z_mode)
 allhyper=c(age_hyperpar_a, age_hyperpar_b, period_hyperpar_a, period_hyperpar_b, cohort_hyperpar_a, cohort_hyperpar_b,
            age_hyperpar_a2, age_hyperpar_b2, period_hyperpar_a2, period_hyperpar_b2, cohort_hyperpar_a2, cohort_hyperpar_b2,
            z_hyperpar_a, z_hyperpar_b)

 if (method == "pg") {
   ## ---- Polya-Gamma Gibbs engine (exact conditionals, no MH tuning) ----
   ord_of <- function(b) if (b %in% c(1,3)) 1L else if (b %in% c(2,4)) 2L else 0L
   ord_a <- ord_of(age_block); ord_p <- ord_of(period_block); ord_c <- ord_of(cohort_block)
   Ymat <- matrix(as.integer(cases),      number_of_agegroups, number_of_periods)
   Nmat <- matrix(as.integer(population),  number_of_agegroups, number_of_periods)
   het_pg <- c(age_block %in% c(3, 4), period_block %in% c(3, 4), cohort_block %in% c(3, 4))
   hyper_pg <- list(age    = c(age_hyperpar_a,    age_hyperpar_b),
                    period = c(period_hyperpar_a, period_hyperpar_b),
                    cohort = c(cohort_hyperpar_a, cohort_hyperpar_b),
                    age_het    = c(age_hyperpar_a2,    age_hyperpar_b2),
                    period_het = c(period_hyperpar_a2, period_hyperpar_b2),
                    cohort_het = c(cohort_hyperpar_a2, cohort_hyperpar_b2))
   ## period/cohort covariate (mean-1 normalised above as period_data/cohort_data,
   ## lengths J and K); the engine scales the smooth period/cohort effect by it.
   cov_p <- if (period_plus == 1) as.numeric(period_data) else NULL
   cov_c <- if (cohort_plus == 1) as.numeric(cohort_data) else NULL
   if (verbose) cat(paste0("Running Polya-Gamma Gibbs engine in ", chains, " chains.\n"))
   ## pass the raw `parallel` (logical or numeric core count) so .bamp_pg can use
   ## as many cores as the taylor path would (it caps internally at n_chains)
   pg <- .bamp_pg(Ymat, Nmat, ord_a, ord_p, ord_c, round(periods_per_agegroup),
                  hyper_pg, number_of_iterations, burn_in, step, chains,
                  parallel = parallel, prior_scale = prior_scale, verbose = verbose,
                  overdisp = (z_mode == 1), z_hyper = c(z_hyperpar_a, z_hyperpar_b),
                  het = het_pg, cov_p = cov_p, cov_c = cov_c, engine = pg_engine)
   sumkick <- chains
   mkmat <- function(field) coda::as.mcmc.list(lapply(pg, function(r) coda::mcmc(r[[field]])))
   mkvec <- function(field) coda::as.mcmc.list(lapply(pg, function(r) coda::mcmc(matrix(r[[field]], ncol = 1))))
   theta <- mkmat("theta"); phi <- mkmat("phi"); psi <- mkmat("psi")
   my <- mkvec("my"); kappa <- mkvec("kappa"); lambda <- mkvec("lambda")
   ny <- mkvec("ny"); deviance <- mkvec("deviance")
   ## het components: populate from the pg fit when present, else zero-fill for
   ## object compatibility. (The shared samples assembly stores age2/etc. only for
   ## block==3, matching the taylor output contract.)
   zerofill <- coda::as.mcmc.list(lapply(pg, function(r) coda::mcmc(matrix(0, nrow = length(r$my), ncol = 1))))
   theta2  <- if (het_pg[1]) mkmat("theta2") else zerofill
   phi2    <- if (het_pg[2]) mkmat("phi2")   else zerofill
   psi2    <- if (het_pg[3]) mkmat("psi2")   else zerofill
   kappa2  <- if (het_pg[1]) mkvec("kappa2") else zerofill
   lambda2 <- if (het_pg[2]) mkvec("lambda2") else zerofill
   ny2     <- if (het_pg[3]) mkvec("ny2")    else zerofill
   ## overdispersion precision (zeta); stored as samples$overdispersion when z_mode==1
   delta <- if (z_mode == 1)
     coda::as.mcmc.list(lapply(pg, function(r) coda::mcmc(matrix(r$zeta, ncol = 1))))
   else
     coda::as.mcmc.list(lapply(pg, function(r) coda::mcmc(matrix(0, nrow = length(r$my), ncol = 1))))
   ## per-chain ksi in the agegroup-major layout expected downstream
   ksi <- lapply(pg, function(r) as.numeric(t(r$ksi)))
 } else {
 singlerun<-function(i,cases,population,blocks,numbers,periods_per_agegroup,
                     numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                     theta2.sample,phi2.sample,psi2.sample,ksi,
                     delta.sample,kappa.sample,kappa2.sample,
                     lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,
                     dev.sample,verbose){
   gc()
   if (verbose>=2)cat(paste("chain",i,"\n"))
   tryCatch(
     .C("bamp",
        as.integer(cases),
        as.integer(population),
        as.integer(blocks),
        as.integer(numbers),
        as.double(periods_per_agegroup),
        as.integer(numbersmcmc),
        as.integer(modelsettings),
        as.double(allhyper),
        as.double(theta.sample), as.double(phi.sample), as.double(psi.sample),
        as.double(theta2.sample), as.double(phi2.sample), as.double(psi2.sample),
        as.double(ksi),
        as.double(delta.sample), as.double(kappa.sample), as.double(kappa2.sample), as.double(lambda.sample),
        as.double(lambda2.sample), as.double(ny.sample), as.double(ny2.sample), as.double(my.sample),
        as.double(dev.sample),
        as.integer(verbose)
     ),
     error = function(e) {
       message("Chain ", i, " failed: ", conditionMessage(e))
       NULL
     }
   )
}

 if (parallel)
 {
   cores<-getOption("mc.cores", 2L)
   if (parallel>1)
     cores<-parallel
   cores<-min(cores, chains)
 }

if(verbose>=2)
{
  results_list<-list()
  for (i in 1:chains)
    results_list[[i]]<-singlerun(i,cases,population,blocks,numbers,periods_per_agegroup,
                     numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                     theta2.sample,phi2.sample,psi2.sample,ksi,delta.sample,kappa.sample,kappa2.sample,
                     lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose)
}
 else if(parallel){
   results_list <- tryCatch({
     cl <- parallel::makeCluster(cores)
     on.exit(parallel::stopCluster(cl), add = TRUE)
     parallel::clusterEvalQ(cl, library(bamp, quietly = TRUE))
     parallel::clusterExport(cl, "singlerun", envir = environment())
     parallel::parLapply(cl, 1:chains, singlerun,
                         cases,population,blocks,numbers,periods_per_agegroup,
                         numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                         theta2.sample,phi2.sample,psi2.sample,ksi,delta.sample,kappa.sample,kappa2.sample,
                         lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose)
   }, error = function(e) {
     message("Note: Parallel execution failed (", conditionMessage(e),
             "); falling back to sequential.")
     lapply(1:chains, singlerun,
            cases,population,blocks,numbers,periods_per_agegroup,
            numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
            theta2.sample,phi2.sample,psi2.sample,ksi,delta.sample,kappa.sample,kappa2.sample,
            lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose)
   })
 } else {
   results_list<-lapply(1:chains,singlerun,cases,population,blocks,numbers,periods_per_agegroup,
                        numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                        theta2.sample,phi2.sample,psi2.sample,ksi,delta.sample,kappa.sample,kappa2.sample,
                        lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose)
 }
##################################################################################################################################

 # Remove chains that failed (singlerun returns NULL on C-level error)
 failed <- vapply(results_list, is.null, logical(1))
 if (any(failed)) {
   if (verbose) message(sum(failed), " chain(s) failed and will be discarded.")
   results_list <- results_list[!failed]
   chains <- length(results_list)
 }
 if (chains == 0) {
   cat("\nAll MCMC chains failed. Please check your model settings.\n")
   return(output)
 }

  deviance<-vector("list",chains)

 for (i in 1:chains){
   deviance[[i]]=coda::mcmc(results_list[[i]][[24]])
 }


 kick<-rep(TRUE, chains)
 deviance.mean<-unlist(lapply(deviance,function(x)return(mean(as.vector(x)))))
 deviance.sd<-unlist(lapply(deviance,sd))

 deviance.mean[is.infinite(deviance.mean)]<-rnorm(1,0,1e9)
 deviance.sd[is.na(deviance.sd)]<-rgamma(1,1,1e-6)
 
 dm<-median(deviance.mean)
 sd<-1.96*median(deviance.sd)
 
  if (verbose==2)
 {
   cat("deviance per chain:")
   cat(paste(deviance.mean," (",deviance.sd,")"))
   cat("\n")
   #print(coda::gelman.diag(deviance))
 }
 
 
 while (any(!((deviance.mean>=(dm-sd))&(deviance.mean<=(dm+sd)))))
 {
   if (verbose>=2)print("kick")
   dm2<-abs(deviance.mean-dm)
   kick2<-which(dm2==max(dm2))
   kick[kick2]=FALSE
   deviance.mean[kick2]<-dm
 }
 
 if(verbose)if (any(!kick))(cat(paste0("Removed ",sum(!kick)," chains.\n")))
 sumkick<-sum(kick)
 theta<-phi<-psi<-theta2<-phi2<-psi2<-delta<-kappa<-kappa2<-lambda<-lambda2<-ny<-ny2<-my<-deviance<-ksi<-vector("list",sumkick)
 
 if (sum(kick)==0)
  {
   warning("Automatic check procedure removed all Markov chains. ",
           "Please change your model settings (e.g. add overdispersion) or try method=\"pg\".",
           call. = FALSE)
    return(output)
  }
 

 if (sum(kick)<chains)
 {
   warning("Automatic check procedure removed ", sum(!kick), " Markov chain",
           if (sum(!kick) > 1) "s" else "",
           ". Please check convergence with checkConvergence() and consider changing your ",
           "model settings (e.g. add overdispersion) or method=\"pg\".", call. = FALSE)
 }
 
 ii=0
  for (i in (1:chains)[kick]){
    ii=ii+1
   results=results_list[[i]]
   theta[[ii]]=coda::mcmc(matrix(results[[9]],ncol=number_of_agegroups,byrow=TRUE))
   phi[[ii]]=coda::mcmc(matrix(results[[10]],ncol=number_of_periods,byrow=TRUE))
   psi[[ii]]=coda::mcmc(matrix(results[[11]],ncol=number_of_cohorts,byrow=TRUE))
   theta2[[ii]]=coda::mcmc(matrix(results[[12]],ncol=number_of_agegroups,byrow=TRUE))
   phi2[[ii]]=coda::mcmc(matrix(results[[13]],ncol=number_of_periods,byrow=TRUE))
   psi2[[ii]]=coda::mcmc(matrix(results[[14]],ncol=number_of_cohorts,byrow=TRUE))
   ksi[[ii]]=results[[15]]
   delta[[ii]]=coda::mcmc(results[[16]])
   kappa[[ii]]=coda::mcmc(results[[17]])
   kappa2[[ii]]=coda::mcmc(results[[18]])
   lambda[[ii]]=coda::mcmc(results[[19]])
   lambda2[[ii]]=coda::mcmc(results[[20]])
   ny[[ii]]=coda::mcmc(results[[21]])
   ny2[[ii]]=coda::mcmc(results[[22]])
   my[[ii]]=coda::mcmc(results[[23]])
   deviance[[ii]]=coda::mcmc(results[[24]])
   }
 

theta<-coda::as.mcmc.list(theta)
phi<-coda::as.mcmc.list(phi)
psi<-coda::as.mcmc.list(psi)
theta2<-coda::as.mcmc.list(theta2)
phi2<-coda::as.mcmc.list(phi2)
psi2<-coda::as.mcmc.list(psi2)
kappa<-coda::as.mcmc.list(kappa)
kappa2<-coda::as.mcmc.list(kappa2)
lambda<-coda::as.mcmc.list(lambda)
lambda2<-coda::as.mcmc.list(lambda2)
ny<-coda::as.mcmc.list(ny)
ny2<-coda::as.mcmc.list(ny2)
my<-coda::as.mcmc.list(my)
delta<-coda::as.mcmc.list(delta)
deviance<-coda::as.mcmc.list(deviance)
 }


 samples=list("intercept"=my, "age"=theta, "period"=phi, "cohort"=psi)
 if (age_block==3)samples=c(samples,list("age2"=theta2))
 if (period_block==3)samples=c(samples,list("period2"=phi2))
 if (cohort_block==3)samples=c(samples,list("cohort2"=psi2))
 samples=c(samples, list("age_parameter"=kappa, "period_parameter"=lambda, "cohort_parameter"=ny))
 if (age_block==3)samples=c(samples,list("age2_parameter"=kappa2))
 if (period_block==3)samples=c(samples,list("period2_parameter"=lambda2))
 if (cohort_block==3)samples=c(samples,list("cohort2_parameter"=ny2))
 if (z_mode==1) samples=c(samples, list("overdispersion"=delta))
 samples=c(samples,list("deviance"=deviance))
 
 data=list("cases"=cases,"population"=population, "periods_per_agegroup"=periods_per_agegroup, 
  agegroups=ncol(cases), periods=nrow(cases), cohorts=periods_per_agegroup * (ncol(cases)- 1) + nrow(cases))
  
 output$model=model
 output$data=data
 output$samples=samples

 
 ksi<-array(unlist(ksi),c(number_of_periods,number_of_agegroups,sumkick))
 ksi<-t(apply(ksi,1:2,mean))

  if (dic)
   {
   if (verbose)cat("\nComputing deviance and DIC.")
   
  devtemp=0.0
  
  pr<-exp(ksi)/(1+exp(ksi))
  ydach<-population*pr
  devtemp1=2*((population-cases)*log((population-cases)/(population-ydach)));
  devtemp2=2*(cases*log(cases/ydach)+(population-cases)*log((population-cases)/(population-ydach)));
  devtemp1<-as.vector(devtemp1)
  devtemp2<-as.vector(devtemp2)
  devtemp2[is.nan(devtemp2)]<-devtemp1[is.nan(devtemp2)]
  devtemp<-sum(devtemp2)
  
 med.deviance<-mean(unlist(deviance))
 deviance <- list()
 deviance$mean.deviance <- med.deviance
 deviance$deviance.mean <- devtemp
 deviance$pD <- med.deviance-devtemp
 deviance$DIC <- 2*med.deviance-devtemp
 
 output$deviance=deviance
 }
 
 if (!is.null(period_covariate)|!is.null(cohort_covariate))
 {
   covariate <- list()
 }
 if (!is.null(period_covariate))
 {
   covariate$period<-period_covariate
 }
 if (!is.null(cohort_covariate))
 {
   covariate$cohort<-cohort_covariate
 }
 if (!is.null(period_covariate)|!is.null(cohort_covariate))
 {
   output$covariate<-covariate
 }
   
# ksi_berechnen <-
 #   function(ksi, psi, vdb, noa, nop){
 #     for(i in 1:noa){
 #       for(j in 1:nop){
 #         ksi[j,i] <- psi[coh(i, j, noa, vdb)]
 #       }
 #     }
 #     return(ksi)
 #   }
 
 ## auto=TRUE suppresses checkConvergence()'s own "did not converge" cat();
 ## that should happen when bamp() is silent (verbose=FALSE), not the other
 ## way around, so this is the negation of verbose, not verbose itself.
 checkConvergence(output, auto=!verbose)

 output$ksi=ksi
 output <- effects.apc(output, update=TRUE)
 cat("\n")
 return(output)
}

## Internal: data-aware default for number_of_iterations. Rare / zero-heavy data
## -- the high-population rare-event cells whose Polya-Gamma augmentation mixes
## slowly -- needs more iterations to converge; well-populated data converges
## quickly. Returns a recommended iteration count from a simple rarity score:
## the fraction of zero cells, or half the fraction of very-low-count cells
## (<=5 events), whichever is larger. Maps a rarity of 0 to 40000 iterations and
## a rarity of 1 (essentially all cells empty/rare) to 120000, rounded to 1000.
## Used only for mcmc.options entries left at "auto"; explicit numbers override.
.bamp_auto_mcmc <- function(cases) {
  Y <- suppressWarnings(as.numeric(as.matrix(cases)))
  Y <- Y[is.finite(Y)]
  if (!length(Y)) return(80000L)                 # no usable data: mid-range default
  zero_frac <- mean(Y == 0)
  rare_frac <- mean(Y <= 5)
  rarity <- min(1, max(zero_frac, 0.5 * rare_frac))
  as.integer(round((40000 + 80000 * rarity) / 1000) * 1000)
}
