#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data, iter = 750, chains = 1, warmup = 250, ...) {


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(system.file("stan", "mmm.stan", package = "mmmstan"))
  # Fit the model
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    parallel_chains = parallel::detectCores(),
    threads_per_chain = 1,
    iter_warmup = warmup,
    iter_sampling = iter,
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
    samples = samples),
    class = "mmmfit")
}
