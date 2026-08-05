bamp_label_width <- function(label) {
  parts <- as.integer(strsplit(label[1], "-")[[1]])
  parts[2] - parts[1] + 1L
}

bamp_set_to_data <- function(one_set) {
  cases <- round(one_set$counts)
  population <- one_set$exposure

  period_width <- bamp_label_width(one_set$period_labels)
  age_width <- bamp_label_width(one_set$age_labels)

  list(
    cases = cases,
    population = population,
    periods_per_agegroup = age_width / period_width
  )
}

bamp_model_grid <- function(
    age_options = c("rw1", "rw2"),
    period_options = c("rw1", "rw2", " "),
    cohort_options = c("rw1", "rw2", " "),
    overdisp_options = c(FALSE, TRUE)
) {
  grid <- expand.grid(
    age = age_options,
    period = period_options,
    cohort = cohort_options,
    overdisp = overdisp_options,
    stringsAsFactors = FALSE
  )

  label_part <- function(x) ifelse(trimws(x) == "", "none", x)
  grid$model_id <- sprintf(
    "A-%s_P-%s_C-%s_OD-%s",
    label_part(grid$age),
    label_part(grid$period),
    label_part(grid$cohort),
    ifelse(grid$overdisp, "T", "F")
  )

  grid
}

bamp_fit_one <- function(
    dat,
    spec,
    mcmc.options,
    parallel,
    hyperpar
) {
  bamp::bamp(
    cases = dat$cases,
    population = dat$population,
    age = spec$age,
    period = spec$period,
    cohort = spec$cohort,
    overdisp = spec$overdisp,
    periods_per_agegroup = dat$periods_per_agegroup,
    mcmc.options = mcmc.options,
    hyperpar = hyperpar,
    dic = TRUE,
    parallel = parallel,
    verbose = FALSE
  )
}

run_bamp_grid <- function(
    apc_sets_path = "/Users/volkerschmid/projects/Apc-Full-vs-Empirical-Bayes--R.Inla/nord_data_prepare/data/nordcan_processed/nordcan_apc_sets.rds",
    output_dir = "bamp_analysis/data/bamp_results",
    set_ids = NULL,
    model_grid = bamp_model_grid(),
    mcmc.options = list(
      number_of_iterations = "auto",
      burn_in = "auto",
      step = "auto",
      tuning = 500
    ),
    hyperpar = list(
      age = c(1, 0.5),
      period = c(1, 5e-04),
      cohort = c(1, 5e-04),
      overdisp = c(1, 0.05)
    ),
    parallel = TRUE,
    overwrite = FALSE,
    verbose = TRUE
) {
  sets <- data.table::as.data.table(readRDS(apc_sets_path))
  if (!is.null(set_ids)) {
    sets <- sets[sets$set_id %in% set_ids]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_path <- file.path(output_dir, "bamp_dic_summary.csv")

  for (i in seq_len(nrow(sets))) {
    set_id <- sets$set_id[i]
    dat <- bamp_set_to_data(sets$data[[i]])

    for (j in seq_len(nrow(model_grid))) {
      spec <- model_grid[j, ]
      out_file <- file.path(output_dir, sprintf("%s__%s.rds", set_id, spec$model_id))

      if (file.exists(out_file) && !isTRUE(overwrite)) {
        if (isTRUE(verbose)) message(sprintf("skip (done): %s | %s", set_id, spec$model_id))
        next
      }

      if (isTRUE(verbose)) message(sprintf("fitting: %s | %s", set_id, spec$model_id))

      start_time <- Sys.time()
      fit <- tryCatch(
        bamp_fit_one(dat, spec, mcmc.options, parallel, hyperpar),
        error = function(e) {
          message(sprintf("  FAILED: %s", conditionMessage(e)))
          NULL
        }
      )
      runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

      converged <- if (!is.null(fit)) {
        isTRUE(tryCatch(bamp::checkConvergence(fit), error = function(e) NA))
      } else {
        NA
      }

      row <- data.table::data.table(
        set_id = set_id,
        model_id = spec$model_id,
        age = spec$age,
        period = spec$period,
        cohort = spec$cohort,
        overdisp = spec$overdisp,
        dic = if (!is.null(fit)) fit$deviance$DIC else NA_real_,
        pD = if (!is.null(fit)) fit$deviance$pD else NA_real_,
        mean_deviance = if (!is.null(fit)) fit$deviance$mean.deviance else NA_real_,
        converged = converged,
        runtime_sec = runtime_sec,
        fitted_at = as.character(Sys.time())
      )
      data.table::fwrite(row, summary_path, append = file.exists(summary_path))

      if (!is.null(fit)) saveRDS(fit, out_file)
    }
  }

  invisible(data.table::fread(summary_path))
}

if (FALSE) {
  library(bamp)
  source("internaltest/nordstat.R")

  # one set, small grid, short MCMC run - to size the runtime
  test_grid <- bamp_model_grid(
    period_options = c("rw2", " "),
    cohort_options = c("rw2", " "),
    overdisp=TRUE
  )
  summary <- run_bamp_grid(
    set_ids = "Sweden__male__Lip",
    model_grid = test_grid#,
    #mcmc.options = list(number_of_iterations = 2000, burn_in = 1000, step = 2, tuning = 200)
  )
  summary[order(summary$dic), ]

  # full run: all NORDCAN sets x full 36-model grid, default (auto) MCMC length
  summary_full <- run_bamp_grid()
}
