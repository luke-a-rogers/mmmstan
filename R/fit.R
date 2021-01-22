#' Fit Markov Movement Model via RStan
#'
#' @param data [list()]
#'
#' @return [mmmfit()]
#' @export
#'
mmmfit <- function (data) {

  # Define options -------------------------------------------------------------

  if (is.null(options()$mc.cores)) {
    options(mc.cores = parallel::detectCores())
  }

  # Fit the model --------------------------------------------------------------

  samples <- rstan::stan(
    file = here::here("inst", "stan", "mmm.stan"),
    chains = 1,
    iter = 2000,
    data = data
  )

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    samples = samples),
    class = "mmmfit")
}
