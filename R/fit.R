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
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    ...
  )

  # Assemble parameters --------------------------------------------------------

  parameters <- summary_to_tibbles(
    cmdfit$summary(),
    data$h_prior_mean,
    data$h_prior_sd,
    data$phi_prior_mean,
    data$phi_prior_sd,
    data$sigma_prior_mean,
    data$sigma_prior_sd
  )

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    output = cmdfit$output_files(),
    summary = cmdfit$summary(),
    parameters = parameters),
    class = "mmmfit")
}
