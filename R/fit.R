#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to [$sample()]
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data,
                    chains = 1,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    threads_per_chain = 8,
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
    threads_per_chain = threads_per_chain#,
    #...
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
