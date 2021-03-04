#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data,
                    chains = 1,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    ...) {


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "mmm.stan", package = "mmmstan"),
    cpp_options = list(stan_threads = TRUE))
  # Fit the model
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
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
