#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data, iter = 750, chains = 1, warmup = 250, ...) {

  # Fit the model --------------------------------------------------------------

  samples <- rstan::stan(
    file = here::here("inst", "stan", "mmm.stan"),
    data = data,
    iter = iter,
    chains = chains,
    warmup = warmup,
    cores = parallel::detectCores(),
    ...
  )

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    samples = samples),
    class = "mmmfit")
}
