#' Fit A Hidden Markov Movement Model via CmdStan and CmdStanR
#'
#' The Stan model data may be specified by individual arguments or by passing
#' a suitable \code{data} argument. The \code{data} argument takes precedence
#' when it is non-null.
#'
#' When the model data are specified by individual arguments, a number of
#' helper functions with sensible defaults are used to assemble the \code{data}
#' object. This is usually preferable for user convenience.
#'
#' Passing the \code{data} object directly may be preferred for example when the
#' \code{data} object is the output of a simulation model, or when greater
#' flexibility is desired. In either event, care must be taken to ensure the
#' \code{data} object matches the requirements of the underlying Stan model.
#'
#' @param tag_data [data.frame()]
#' @param list_regions [list()]
#' @param list_sizes [list()]
#' @param year_start [integer()] year of initial tag released
#' @param year_end [integer()] year of final tag recovered
#' @param step_interval [character()] one of \code{"month", "quarter", "year"}
#' @param step_duration_max [integer()]
#' @param days_duration_min [integer()] minimum number of days before recovery
#' @param colname_date_released [character()]
#' @param colname_date_recovered [character()]
#' @param colname_region_released [character()]
#' @param colname_region_recovered [character()]
#' @param colname_size_released [character()]
#' @param movement_pattern [integer()]
#' @param movement_allow [integer()][matrix()]
#' @param movement_disallow [integer()][matrix()]
#' @param mu_movement_step_diag [numeric()][vector()]
#' @param sd_movement_step_diag [numeric()][vector()]
#' @param mu_fishing_rate [numeric()][array()]
#' @param cv_fishing_rate [numeric()]
#' @param mu_selectivity [numeric()][array()]
#' @param cv_selectivity [numeric()]
#' @param mu_fishing_weight [numeric()]
#' @param sd_fishing_weight [numeric()]
#' @param mu_natural_mortality_rate [numeric()]
#' @param sd_natural_mortality_rate [numeric()]
#' @param mu_reporting_rate [numeric()]
#' @param sd_reporting_rate [numeric()]
#' @param mu_initial_loss_rate [numeric()]
#' @param sd_initial_loss_rate [numeric()]
#' @param mu_ongoing_loss_rate [numeric()]
#' @param sd_ongoing_loss_rate [numeric()]
#' @param mu_dispersion [numeric()]
#' @param sd_dispersion [numeric()]
#' @param tolerance_expected [numeric()]
#' @param tolerance_fishing [numeric()]
#' @param data [list()] See details
#' @param chains [integer()] number of chains
#' @param step_size [numeric()] initial step size
#' @param adapt_delta [numeric()] the adaptation target acceptance statistic
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param max_treedepth [integer()]
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
#' @param refresh [integer()]
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @details TBD
#'
#' @return [mmmstan::mmmstan()]
#' @export
#'
mmmstan <- function (tag_data,
                     # Tag arguments
                     list_regions,
                     list_sizes,
                     year_start,
                     year_end,
                     step_interval = "quarter",
                     step_duration_max = NULL,
                     days_duration_min = 90,
                     colname_date_released = "date_released",
                     colname_date_recovered = "date_recovered",
                     colname_region_released = "region_released",
                     colname_region_recovered = "region_recovered",
                     colname_size_released = "size_released",
                     # Movement index
                     movement_pattern = 2,
                     movement_allow = NULL,
                     movement_disallow = NULL,
                     # Movement step mean priors
                     mu_movement_step_diag = NULL,
                     sd_movement_step_diag = NULL,
                     # Fishing rate priors
                     mu_fishing_rate = NULL,
                     cv_fishing_rate = NULL,
                     # Selectivity priors
                     mu_selectivity = NULL,
                     cv_selectivity = NULL,
                     # Fishing weight priors
                     mu_fishing_weight = NULL,
                     sd_fishing_weight = NULL,
                     # Natural mortality rate priors
                     mu_natural_mortality_rate = NULL,
                     sd_natural_mortality_rate = NULL,
                     # Fractional (per tag) reporting rate priors
                     mu_reporting_rate = NULL,
                     sd_reporting_rate = NULL,
                     # Fractional (per tag) initial loss rate priors
                     mu_initial_loss_rate = 0.1,
                     sd_initial_loss_rate = 0.01,
                     # Instantaneous ongoing loss rate priors
                     mu_ongoing_loss_rate = 0.02,
                     sd_ongoing_loss_rate = 0.001,
                     # Dispersion priors
                     mu_dispersion = 1.0,
                     sd_dispersion = 0.5,
                     # Tolerance values
                     tolerance_expected = 1e-12,
                     tolerance_fishing =  1e-12,
                     # CmdStanR
                     data = NULL,
                     chains = 1,
                     step_size = 0.01,
                     adapt_delta = 0.95,
                     iter_warmup = 250,
                     iter_sampling = 750,
                     max_treedepth = 10,
                     use_reduce_sum = FALSE,
                     threads_per_chain = parallel::detectCores()/(2*chains),
                     refresh = 100,
                     ...) {

  # Check arguments ------------------------------------------------------------

  # Has data?
  has_data <- ifelse(is.null(data), FALSE, TRUE)
  # Tag data
  checkmate::assert_data_frame(
    tag_data,
    all.missing = FALSE,
    null.ok = has_data
  )
  # Tag arguments
  checkmate::assert_list(
    list_regions,
    any.missing = FALSE,
    min.len = 2,
    null.ok = has_data
  )
  checkmate::assert_list(
    list_sizes,
    any.missing = FALSE,
    min.len = 1,
    null.ok = has_data
  )
  checkmate::assert_integerish(
    year_start,
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_integerish(
    year_end,
    lower = year_start + 1L,
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_choice(
    step_interval,
    choices = c("month", "quarter", "year"),
    null.ok = has_data
  )
  checkmate::assert_integerish(
    step_duration_max,
    lower = 2,
    any.missing = FALSE,
    len = 1,
    null.ok = TRUE
  )
  checkmate::assert_integerish(
    days_duration_min,
    lower = 1,
    any.missing = FALSE,
    len = 1,
    null.ok = FALSE
  )
  checkmate::assert_character(
    colname_date_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_date_recovered,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_region_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_region_recovered,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_size_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  # Movement index
  checkmate::assert_choice(
    movement_pattern,
    choices =  c(1L, 2L),
    null.ok = has_data
  )
  checkmate::assert_matrix(
    movement_allow,
    mode = "integerish",
    any.missing = FALSE,
    ncols = 2,
    null.ok = TRUE
  )
  checkmate::assert_matrix(
    movement_disallow,
    mode = "integerish",
    any.missing = FALSE,
    ncols = 2,
    null.ok = TRUE
  )
  # Movement step mean priors
  checkmate::assert_numeric(
    mu_movement_step_diag,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    len = length(list_regions),
    null.ok = TRUE
  )
  checkmate::assert_numeric(
    sd_movement_step_diag,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions),
    null.ok = TRUE
  )

  # Fishing rate priors
  # Selectivity priors
  # Fishing weight priors
  # Natural mortality rate priors
  # Fractional (per tag) reporting rate priors
  # Fractional (per tag) initial loss rate priors
  # Instantaneous ongoing loss rate priors

  # Dispersion priors
  checkmate::assert_number(
    mu_dispersion,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    sd_dispersion,
    lower = 0,
    null.ok = has_data
  )
  # Tolerance value
  checkmate::assert_number(
    tolerance_expected,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    tolerance_fishing,
    lower = 0,
    null.ok = has_data
  )
  # Data
  checkmate::assert_list(
    data,
    types = c("integer", "double"),
    any.missing = FALSE,
    null.ok = TRUE
  )
  # CmdStanR arguments
  checkmate::assert_integerish(
    chains,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_number(
    step_size,
    lower = 0,
    finite = TRUE
  )
  checkmate::assert_number(
    adapt_delta,
    lower = 0,
    upper = 1,
    finite = TRUE
  )
  checkmate::assert_integerish(
    iter_warmup,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    iter_sampling,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    max_treedepth,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_logical(
    use_reduce_sum,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    threads_per_chain,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    refresh,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )

  # Print messages -------------------------------------------------------------

  if (has_data) {
    cat("\nmmmfit(): data argument present--ignoring most other arguments\n")
  }

  # Assemble data --------------------------------------------------------------

  if (is.null(data)) {

    # Populate null arguments --------------------------------------------------

    # Step duration at large maximum
    if (is.null(step_duration_max)) {
      step_duration_max <- count_model_steps(
        year_start,
        year_end,
        step_interval
      )
    }
    # Movement step prior
    if (is.null(mu_movement_step_diag)) {
      mu_movement_step_diag <- rep(0.98, length(list_regions))
    }
    if (is.null(sd_movement_step_diag)) {
      sd_movement_step_diag <- rep(0.05, length(list_regions))
    }
    # Fishing rate prior
    if (is.null(mu_fishing_rate)) {
      mu_fishing_rate <- matrix(
        0.05,
        nrow = year_end - year_start + 1L,
        ncol = length(list_regions)
      )
    }
    if (is.null(cv_fishing_rate)) {
      cv_fishing_rate <- 0.1
    }
    # Selectivity prior
    if (length(list_sizes) == 1) {
      mu_selectivity_short <- array(1, dim = c(0, length(list_regions))) # [0, X]
      cv_selectivity <- numeric(0) # [0]
    } else {
      if (is.null(mu_selectivity)) {
        mu_selectivity_short <- array(
          0.9,
          dim = c(length(list_sizes) - 1, length(list_regions))
        )
      } else {
        mu_selectivity_short <-
          mu_selectivity[seq_len(length(list_sizes) - 1),, drop = FALSE]
      }
      if (is.null(cv_selectivity)) {
        cv_selectivity <- 0.1
      }
    }
    # Natural mortality rate prior
    if (is.null(mu_natural_mortality_rate)) {
      mu_natural_mortality_rate <- rep(0.1, length(list_regions))
    }
    if (is.null(sd_natural_mortality_rate)) {
      sd_natural_mortality_rate <- rep(0.01, length(list_regions))
    }
    # Reporting rate prior
    if (is.null(mu_reporting_rate)) {
      mu_reporting_rate <- rep(0.5, length(list_regions))
    }
    if (is.null(sd_reporting_rate)) {
      sd_reporting_rate <- rep(0.05, length(list_regions))
    }

    # Assemble tag data --------------------------------------------------------

    tag_array <- create_tag_array(
      tag_data = tag_data,
      list_regions = list_regions,
      list_sizes = list_sizes,
      year_released_start = year_start,
      year_recovered_end = year_end,
      step_duration_max = step_duration_max,
      step_interval = step_interval,
      days_duration_min = days_duration_min,
      colname_date_released = colname_date_released,
      colname_date_recovered = colname_date_recovered,
      colname_region_released = colname_region_released,
      colname_region_recovered = colname_region_recovered,
      colname_size_released = colname_size_released
    )

    # Assemble movement index --------------------------------------------------

    movement_index <- create_movement_index(
      number_of_regions = length(list_regions),
      movement_pattern = movement_pattern,
      movement_allow = movement_allow,
      movement_disallow = movement_disallow
    )

    # Assemble data ------------------------------------------------------------

    data <- list(
      # Model index limits
      N = count_model_steps(year_start, year_end, step_interval),
      D = step_duration_max,
      L = length(list_sizes),
      X = length(list_regions),
      T = year_end - year_start + 1L,
      K = count_intervals(step_interval, "year"),
      # Index arrays
      n_to_t = index_n_to_t(year_start, year_end, step_interval),
      # Tag data
      tags = tag_array, # [N - 1, D, L, X, X]
      # Movement index
      movement_index = movement_index, # [X, X]
      # Movement step priors
      mu_movement_step_diag = mu_movement_step_diag, # [X]
      sd_movement_step_diag = sd_movement_step_diag, # [X]
      # Fishing rate priors
      mu_fishing_rate = mu_fishing_rate, # [T, X]
      cv_fishing_rate = cv_fishing_rate, # [1]
      # Selectivity priors
      mu_selectivity_short = mu_selectivity_short, # [0, X] OR [L - 1, X]
      cv_selectivity = cv_selectivity, # NULL OR [1]
      # Fishing weight priors
      #  mu_fishing_weight = , # [K, X]
      #  cv_fishing_weight = , # [1]
      # Natural mortality rate priors
      mu_natural_mortality_rate = mu_natural_mortality_rate, # [X]
      sd_natural_mortality_rate = sd_natural_mortality_rate, # [X]
      # Fractional (per tag) reporting rate priors
      mu_reporting_rate = mu_reporting_rate, # [X]
      sd_reporting_rate = sd_reporting_rate, # [X]
      # Fractional (per tag) initial loss rate priors
      mu_initial_loss_rate = mu_initial_loss_rate,
      sd_initial_loss_rate = sd_initial_loss_rate,
      # Instantaneous ongoing loss rate priors
      mu_ongoing_loss_rate = mu_ongoing_loss_rate,
      sd_ongoing_loss_rate = sd_ongoing_loss_rate,
      # Dispersion priors
      mu_dispersion = mu_dispersion,
      sd_dispersion = sd_dispersion,
      # Tolerance values
      tolerance_expected = tolerance_expected,
      tolerance_fishing = tolerance_fishing
    )
  }

  # Check data elements --------------------------------------------------------

  # Initialize model -----------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    ifelse(
      use_reduce_sum,
      system.file("stan", "mmm_reduce_sum.stan", package = "mmmstan"),
      system.file("stan", "mmm.stan", package = "mmmstan")
    ),
    include_path = system.file("stan", package = "mmmstan"),
    cpp_options = list(stan_threads = TRUE)
  )

  # Fit model ------------------------------------------------------------------

  fit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    adapt_delta = adapt_delta,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    max_treedepth = max_treedepth,
    threads_per_chain = threads_per_chain,
    refresh = refresh,
    ...
  )

  # Placate R-CMD-check --------------------------------------------------------

  # Values names
  movement_rate <- NULL
  fishing_rate <- NULL
  natural_mortality_rate <- NULL
  reporting_rate <- NULL
  initial_loss_rate <- NULL
  ongoing_loss_rate <- NULL
  selectivity <- NULL
  dispersion <- NULL
  # Dimension names
  x <- NULL
  y <- NULL
  t <- NULL
  k <- NULL
  l <- NULL

  # Compute unconditional fit summaries ----------------------------------------

  # Movement rate
  movement_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(movement_rate[l,x,y]) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Fishing rate
  fishing_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(fishing_rate[t,x]) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Natural mortality rate
  natural_mortality_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(natural_mortality_rate[x]) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Reporting rate
  reporting_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(reporting_rate[x]) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Initial loss rate
  initial_loss_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(initial_loss_rate) %>%
    tidybayes::summarise_draws()
  # Ongoing loss rate
  ongoing_loss_rate_summary <- fit$draws() %>%
    tidybayes::spread_draws(ongoing_loss_rate) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Selectivity
  if (length(list_sizes) > 1) {
    selectivity_summary <- fit$draws() %>%
      tidybayes::spread_draws(selectivity[l, x]) %>%
      tidybayes::summarise_draws() %>%
      dplyr::ungroup()
  } else {
    selectivity_summary <- NULL
  }
  # Dispersion
  dispersion_summary <- fit$draws() %>%
    tidybayes::spread_draws(dispersion) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()

  # Assemble fit summary -------------------------------------------------------

  summary <- list(
    # Movement rate
    movement_rate = movement_rate_summary,
    # Fishing rate
    fishing_rate = fishing_rate_summary,
    # Natural mortality rate
    natural_mortality_rate = natural_mortality_rate_summary,
    # Reporting rate
    reporting_rate = reporting_rate_summary,
    # Tag loss
    initial_loss_rate = initial_loss_rate_summary,
    ongoing_loss_rate = ongoing_loss_rate_summary,
    # Selectivity
    selectivity = selectivity_summary,
    # Dispersion
    dispersion = dispersion_summary
  )

  # Return values --------------------------------------------------------------

  structure(list(
    data = data,
    fit = fit,
    rstanfit = rstan::read_stan_csv(fit$output_files()),
    summary = summary),
    class = "mmmstan")
}
