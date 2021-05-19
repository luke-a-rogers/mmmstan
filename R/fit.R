#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data,
                    chains = 1,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    use_reduce_sum = FALSE,
                    threads_per_chain = 8,
                    ...) {

  # Check arguments ------------------------------------------------------------

  # Check that max(rowSums(data$z)) <= 6


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
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    ...
  )
  # Convert to stanfit
  samples <- rstan::read_stan_csv(cmdfit$output_files())

  # Fit the model --------------------------------------------------------------

  # samples <- rstan::sampling(
  #   object = stanmodels$mmm,
  #   data = data,
  #   iter = iter,
  #   chains = chains,
  #   warmup = warmup,
  #   cores = parallel::detectCores(),
  #   ...
  # )

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    samples = samples,
    cmdfit = cmdfit),
    class = "mmmfit")
}
