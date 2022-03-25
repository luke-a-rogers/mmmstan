#' Beta Distribution Parameters
#'
#' Compute beta distribution mean, standard deviation, alpha and beta
#' parameters from (array) mean and standard deviation or alpha and beta
#' parameters.
#'
#' @param mu [numeric()] mean values
#' @param sigma [numeric()] standard deviations
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
                            sigma = NULL,
                            alpha = NULL,
                            beta = NULL) {

  # Check arguments
  checkmate::assert_numeric(mu, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(sigma, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(alpha, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(beta, lower = 0, finite = TRUE, null.ok = TRUE)

  # Compute parameters
  if ((is.null(mu) | is.null(sigma)) & (is.null(alpha) | is.null(beta))) {
    mu <- NULL
    sigma <- NULL
    alpha <- NULL
    beta <- NULL
  } else if (!is.null(mu) & !is.null(sigma)) {
    # Handle scalar mu
    if (is.null(dim(mu))) dim(mu) <- 1
    # Compute variance
    var <- sigma * sigma
    # Check parameter condition
    if (!all(var < (mu * (1 - mu)))) {
      stop("all((sigma * sigma) < (mu * (1 - mu))) must be true")
    }
    # Compute parameters
    nu <- (mu * (1 - mu) / var) - 1
    alpha <- mu * nu
    beta <- (1 - mu) * nu
  } else if (!is.null(alpha) & !is.null(beta)) {
    mu <- alpha / (alpha + beta)
    sigma <- sqrt(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1)))
  }
  # Return list
  return(list(mu = mu, sigma = sigma, alpha = alpha, beta = beta))
}

#' Gamma Distribution Parameters
#'
#' Compute gamma distribution mean, standard deviation, alpha and beta
#' parameters from (array) mean and standard deviation or alpha and beta
#' parameters.
#'
#' @param mu [numeric()] mean values
#' @param sigma [numeric()] standard deviations
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
                              sigma = NULL,
                              alpha = NULL,
                              beta = NULL) {

  # Check arguments
  checkmate::assert_numeric(mu, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(sigma, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(alpha, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(beta, lower = 0, finite = TRUE, null.ok = TRUE)

  # Compute parameters
  if ((is.null(mu) | is.null(sigma)) & (is.null(alpha) | is.null(beta))) {
    mu <- NULL
    sigma <- NULL
    alpha <- NULL
    beta <- NULL
  } else if (!is.null(mu) & !is.null(sigma)) {
    # Handle scalar mu
    if (is.null(dim(mu))) dim(mu) <- 1
    # Compute parameters
    alpha <- mu^2 / sigma^2
    beta <- mu / sigma^2
  } else if (!is.null(alpha) & !is.null(beta)) {
    mu <- alpha / beta
    sigma <- sqrt(alpha/beta^2)
  }
  # Return list
  return(list(mu = mu, sigma = sigma, alpha = alpha, beta = beta))
}

#' Summarise Posterior Draws
#'
#' @param x [cmdstanr::sample()] model fit object
#' @param data [list()] model data
#'
#' @return [list()]
#' @export
#'
summarise_posterior_draws <- function (x, data) {
  # Placate R-CMD-check
  previous_area <- NULL
  current_area <- NULL
  movement_time <- NULL
  released_group <- NULL
  harvest_group <- NULL
  harvest_time <- NULL
  current_area <- NULL
  # lp__
  lp__ <- tidybayes::spread_draws(x, lp__) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Movement rates
  p <- tidybayes::spread_draws(
    x,
    p[previous_area, current_area, movement_time, released_group]
  ) %>%
    tidybayes::summarise_draws() %>%
    dplyr::mutate(
      prior_lower = 0,
      prior_upper = 1
    ) %>%
    dplyr::ungroup()
  # Harvest rates
  h <- tidybayes::spread_draws(
    x,
    h[harvest_group, harvest_time, current_area]
  ) %>%
    tidybayes::summarise_draws() %>%
    dplyr::mutate(
      prior_mean = NA_real_,
      prior_sd = NA_real_
    ) %>%
    dplyr::ungroup()
  for(i in seq_len(nrow(h))) {
    h$prior_mean[i] <- data$h_prior_mean[
      h$harvest_group[i],
      h$harvest_time[i],
      h$current_area[i]]
    h$prior_sd[i] <- data$h_prior_sd[
      h$harvest_group[i],
      h$harvest_time[i],
      h$current_area[i]]
  }
  # Negative binomial dispersion
  phi <- tidybayes::spread_draws(x, phi) %>%
    tidybayes::summarise_draws() %>%
    dplyr::mutate(
      prior_mean = data$phi_prior_mean,
      prior_sd = data$phi_prior_sd
    ) %>%
    dplyr::ungroup()
  # Random walk standard deviation
  n <- nrow(dplyr::filter(x$summary(), startsWith(.data$variable, "sigma")))
  if (n > 0) {
    sigma <- tidybayes::spread_draws(x, sigma[current_area]) %>%
      tidybayes::summarise_draws() %>%
      dplyr::mutate(
        prior_mean = NA_real_,
        prior_sd = NA_real_
      ) %>%
      dplyr::ungroup()
    for(i in seq_len(nrow(sigma))) {
      sigma$prior_mean[i] <- data$sigma_prior_mean[sigma$current_area[i]]
      sigma$prior_sd[i] <- data$sigma_prior_sd[sigma$current_area[i]]
    }
  } else {
    sigma <- NULL
  }
  # Return list
  list(
    lp__ = lp__,
    p = p,
    h = h,
    phi = phi,
    sigma = sigma
  )
}

#' Tag Arrays
#'
#' @param released [data.frame()] released tags. See details.
#' @param recovered [data.frame()] recovered tags. See details.
#' @param released_time_unit [character()] one of \code{year}, \code{quarter},
#'   or \code{month}.
#' @param released_time_max [numeric()] equivalent to \code{T_released}.
#' @param liberty_time_max [numeric()] equivalent to \code{T_liberty}.
#' @param liberty_days_min [numeric()] (not currently implemented).
#' @param colname_released_date [character()] released date column name.
#' @param colname_released_area [character()] released area column name.
#' @param colname_group [character()] released group column name.
#' @param colname_recovered_date [character()] recovered date column name.
#' @param colname_recovered_area [character()] recovered area column name.
#' @param colname_id [character()] tag identifier column name.
#' @param area_list [list()] of named atomic vectors. See details.
#' @param group_list [list()] of named atomic vectors. See details.
#' @param released_date_start [character()] earliest released date.
#' @param released_date_end [character()] latest released date.
#'
#' @details TBD
#'
#' @importFrom rlang .data
#'
#' @return [list()] of released and recovered [array()]s
#' @export
#'
tag_arrays <- function (released,
                        recovered,
                        released_time_unit,
                        released_time_max,
                        liberty_time_max,
                        liberty_days_min,
                        colname_released_date,
                        colname_released_area,
                        colname_group,
                        colname_recovered_date,
                        colname_recovered_area,
                        colname_id,
                        area_list,
                        group_list,
                        released_date_start,
                        released_date_end) {

  # Check arguments ------------------------------------------------------------


  # Released tibble ------------------------------------------------------------

  xt <- released %>%
    dplyr::rename(
      released_date = .data[[colname_released_date]],
      released_area = .data[[colname_released_area]],
      group_raw = .data[[colname_group]],
      id = .data[[colname_id]]
    ) %>%
    # TODO: filter disallowed values
    dplyr::mutate(
      rt = create_tag_release_steps(
        .data$released_date,
        released_date_start,
        released_date_end,
        released_time_unit
      ),
      ra = create_tag_areas(
        .data$released_area,
        area_list
      ),
      rg = create_tag_groups(
        .data$group_raw,
        group_list
      )
    ) %>%
    dplyr::select(.data$rt, .data$ra, .data$rg) %>%
    dplyr::group_by(.data$rt, .data$ra, .data$rg) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(.data$rt, .data$ra, .data$rg) %>%
    tidyr::drop_na()

  # Recovered tibble -----------------------------------------------------------

  yt <- recovered %>%
    dplyr::rename(
      released_date = .data[[colname_released_date]],
      released_area = .data[[colname_released_area]],
      group_raw = .data[[colname_group]],
      recovered_date = .data[[colname_recovered_date]],
      recovered_area = .data[[colname_recovered_area]],
      id = .data[[colname_id]]
    ) %>%
    # TODO: filter disallowed values
    dplyr::mutate(
      rt = create_tag_release_steps(
        .data$released_date,
        released_date_start,
        released_date_end,
        released_time_unit
      ),
      ra = create_tag_areas(
        .data$released_area,
        area_list
      ),
      rg = create_tag_groups(
        .data$group_raw,
        group_list
      ),
      lt = create_tag_liberty_steps(
        .data$released_date,
        .data$recovered_date,
        released_date_start,
        released_date_end,
        released_time_unit,
        liberty_time_max
      ),
      ca = create_tag_areas(
        .data$recovered_area,
        area_list
      )
    ) %>%
    dplyr::select(.data$rt, .data$ra, .data$rg, .data$lt, .data$ca) %>%
    dplyr::group_by(.data$rt, .data$ra, .data$rg, .data$lt, .data$ca) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(.data$rt, .data$ra, .data$rg, .data$lt, .data$ca) %>%
    tidyr::drop_na()

  # Compute array dimensions ---------------------------------------------------

  nt <- released_time_max
  na <- length(area_list)
  ng <- length(group_list)
  nl <- liberty_time_max

  # Populate arrays ------------------------------------------------------------

  # Initialize
  x <- array(0, dim = c(nt, na, ng))
  y <- array(0, dim = c(nt, na, ng, nl, na))

  # Populate from xt
  for (i in seq_len(nrow(xt))) {
    x[xt$rt[i], xt$ra[i], xt$rg[i]] <- xt$count[i]
  }
  # Populate from yt
  for (i in seq_len(nrow(yt))) {
    y[yt$rt[i], yt$ra[i], yt$rg[i], yt$lt[i], yt$ca[i]] <- yt$count[i]
  }

  # Return tag array -----------------------------------------------------------

  return(list(x = x, y = y))
}
