#' Fit Markov Movement Model via CmdStanR
#'
#' @param data [list()] See details
#' @param chains [integer()] number of chains
#' @param step_size [numeric()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @details TBD
#'
#' @return [mmmstan::fit()]
#' @export
#'
fit <- function (data = NULL,
                 chains = 1,
                 step_size = 0.01,
                 iter_warmup = 250,
                 iter_sampling = 750,
                 use_reduce_sum = FALSE,
                 threads_per_chain = 8,
                 ...) {

  # Check arguments ------------------------------------------------------------

  # Assemble tag data ----------------------------------------------------------

  # if (is.null(data)) {
  #   tags_released <-
  #   tags_recovered <-
  # }

  # Assemble data --------------------------------------------------------------

  # if (is.null(data)) {
  #   data <- list(
  #     # Index limits
  #     I = n_terms, # Number of terms per unit of time
  #     L = n_liberty, # Number of maximum steps at liberty
  #     N = n_steps, # Number of release steps ((T * I) - 1)
  #     P = sum(movement_index), # Number of movement rate mean parameters
  #     S = n_sizes, # Number of size classes
  #     T = n_times, # Number of times (usually years; release only)
  #     X = n_regions, # Number of geographic regions
  #     # Movement index array
  #     movement_index = movement_index, # array[X, X]
  #     # Tag data
  #     tags_released = tags_released, # array[N, S, X]
  #     tags_recovered = tags_recovered, # array[N, S, X, L, X]
  #     # Rates
  #     initial_loss_rate = initial_loss_rate, # Fraction
  #     ongoing_loss_rate = ongoing_loss_rate, # Instantaneous
  #     # Prior means
  #     mu_mortality_rate = mu_mortality_rate, # array[X]
  #     mu_reporting_rate = mu_reporting_rate, # array[X]
  #     mu_fishing_mean_rate = mu_fishing_mean_rate, # array[X]
  #     mu_dispersion = mu_dispersion,
  #     # Prior standard deviations
  #     sd_mortality_rate = sd_mortality_rate, # array[X]
  #     sd_reporting_rate = sd_reporting_rate, # array[X]
  #     sd_fishing_mean_rate = sd_fishing_mean_rate, # array[X]
  #     sd_dispersion = sd_dispersion,
  #     # Autoregressive movement process parameters
  #     cv_movement_deviation_time_rate = cv_movement_deviation_time_rate,
  #     cv_movement_deviation_term_rate = cv_movement_deviation_term_rate,
  #     cv_movement_deviation_size_rate = cv_movement_deviation_size_rate,
  #     # Autoregressive fishing process parameters
  #     cv_fishing_deviation_time_rate = cv_fishing_deviation_time_rate,
  #     cv_fishing_deviation_term_rate = cv_fishing_deviation_term_rate,
  #     cv_fishing_deviation_size_rate = cv_fishing_deviation_size_rate,
  #     # Fudge values
  #     expected_fudge = expected_fudge,
  #   )
  # }

  # Initialize the model -------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    ifelse(
      use_reduce_sum,
      system.file("stan", "mmm_reduce_sum.stan", package = "mmmstan"),
      system.file("stan", "mmm.stan", package = "mmmstan")
    ),
    include_path = system.file("stan", package = "mmmstan"),
    cpp_options = list(stan_threads = TRUE)
  )

  # Fit the model --------------------------------------------------------------

  fit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    ...
  )

  # Return values --------------------------------------------------------------

  structure(list(
    data = data,
    fit = fit),
    class = "mmmstan::fit")
}





#' Fit Markov Movement Model via CmdStanR
#'
#' @param data [list()] See details
#' @param chains [integer()] number of chains
#' @param step_size [numeric()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @details TBD
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data,
                    chains = 1,
                    step_size = 0.01,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    use_reduce_sum = FALSE,
                    threads_per_chain = 8,
                    ...) {

  # Check arguments ------------------------------------------------------------

  # Check T_study
  T_released <- data$T_released
  T_liberty <- data$T_liberty
  T_study <- data$T_study
  if (T_study == T_liberty) {
    if (T_study > T_released) {
      # cat("\nT_study == T_liberty & T_study > T_released\n")
    } else {
      stop("T_study must be greater than T_released\n")
    }
  } else if (T_study == T_liberty + T_released - 1) {
    # cat("T_study == T_liberty + T_released - 1\n")
  } else {
    stop("T_study must equal T_liberty or T_liberty + T_released - 1\n")
  }
  # Check harvest prior parameters
  h_prior_mean <- data$h_prior_mean
  h_prior_sd <- data$h_prior_sd
  if (!all((h_prior_sd * h_prior_sd) < (h_prior_mean * (1 - h_prior_mean)))) {
    stop("all((h_prior_sd^2) < (h_prior_mean*(1-h_prior_mean))) must be true\n")
  }


  # Check that max(rowSums(data$z)) <= 6
  # Check that if rw == 1, P > 2
  # Check dimensions
  # Check lengths of indices


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(
    ifelse(
      use_reduce_sum,
      system.file("stan", "mmm_reduce_sum.stan", package = "mmmstan"),
      system.file("stan", "mmm.stan", package = "mmmstan")
    ),
    include_path = system.file("stan", package = "mmmstan"),
    cpp_options = list(stan_threads = TRUE))
  # Fit the model
  fit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    ...
  )

  # Assemble summaries ---------------------------------------------------------

  summaries <- summarise_posterior_draws(fit, data)

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    fit = fit,
    summaries = summaries),
    class = "mmmfit")
}
