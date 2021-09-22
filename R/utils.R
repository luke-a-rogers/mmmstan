#' Beta Distribution Parameters
#'
#' Compute beta distribution mean, standard deviation, alpha and beta
#' parameters from (array) mean and standard deviation or alpha and beta
#' parameters.
#'
#' @param mu [numeric()] mean values
#' @param sd [numeric()] standard deviations
#' @param alpha [numeric()] alpha parameters
#' @param beta [numeric()] beta parameters
#'
#' @return [list()]
#' @export
#'
#' @examples
#' # Scalar mu
#' mu <- 0.05
#' l <- beta_parameters(mu, 0.01)
#' hist(rbeta(10000, l$alpha, l$beta), breaks = 100)
#'
#' # Array mu
#' mu <- array(0.05, dim = c(2, 10, 3))
#' l <- beta_parameters(mu, 0.001)
#' hist(rbeta(10000, l$alpha[1,1,1], l$beta[1,1,1]), breaks = 100)
#'
beta_parameters <- function(mu = NULL,
                            sd = NULL,
                            alpha = NULL,
                            beta = NULL) {

  # Check arguments
  checkmate::assert_numeric(mu, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(sd, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(alpha, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(beta, lower = 0, finite = TRUE, null.ok = TRUE)

  # Compute parameters
  if ((is.null(mu) | is.null(sd)) & (is.null(alpha) | is.null(beta))) {
    mu <- NULL
    sd <- NULL
    alpha <- NULL
    beta <- NULL
  } else if (!is.null(mu) & !is.null(sd)) {
    # Handle scalar mu
    if (is.null(dim(mu))) dim(mu) <- 1
    # Compute variance
    var <- sd * sd
    # Check parameter condition
    if (!all(var < (mu * (1 - mu)))) {
      stop("all((sd * sd) < (mu * (1 - mu))) must be true")
    }
    # Compute parameters
    nu <- (mu * (1 - mu) / var) - 1
    alpha <- mu * nu
    beta <- (1 - mu) * nu
  } else if (!is.null(alpha) & !is.null(beta)) {
    mu <- alpha / (alpha + beta)
    sd <- sqrt(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1)))
  }
  # Return list
  return(list(mu = mu, sd = sd, alpha = alpha, beta = beta))
}

#' Gamma Distribution Parameters
#'
#' Compute gamma distribution mean, standard deviation, alpha and beta
#' parameters from (array) mean and standard deviation or alpha and beta
#' parameters.
#'
#' @param mu [numeric()] mean values
#' @param sd [numeric()] standard deviations
#' @param alpha [numeric()] alpha parameters
#' @param beta [numeric()] beta parameters
#'
#' @return [list()]
#' @export
#'
#' @examples
#' # Scalar mu
#' mu <- 0.05
#' l <- gamma_parameters(mu, 0.01)
#' hist(rgamma(10000, l$alpha, l$beta), breaks = 100)
#'
#' # Array mu
#' mu <- array(0.05, dim = c(2, 10, 3))
#' l <- gamma_parameters(mu, 0.001)
#' hist(rgamma(10000, l$alpha[1,1,1], l$beta[1,1,1]), breaks = 100)
#'
gamma_parameters <- function (mu = NULL,
                              sd = NULL,
                              alpha = NULL,
                              beta = NULL) {

  # Check arguments
  checkmate::assert_numeric(mu, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(sd, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(alpha, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(beta, lower = 0, finite = TRUE, null.ok = TRUE)

  # Compute parameters
  if ((is.null(mu) | is.null(sd)) & (is.null(alpha) | is.null(beta))) {
    mu <- NULL
    sd <- NULL
    alpha <- NULL
    beta <- NULL
  } else if (!is.null(mu) & !is.null(sd)) {
    # Handle scalar mu
    if (is.null(dim(mu))) dim(mu) <- 1
    # Compute parameters
    alpha <- mu^2 / sd^2
    beta <- mu / sd^2
  } else if (!is.null(alpha) & !is.null(beta)) {
    mu <- alpha / beta
    sd <- sqrt(alpha/beta^2)
  }
  # Return list
  return(list(mu = mu, sd = sd, alpha = alpha, beta = beta))
}
