#' Fit A Markov Movement Model via CmdStanR
#'
#' Stan model data may be specified in two ways: via arguments \code{tag_data}
#' through \code{expected_fudge} several of which take advantage of helper
#' functions and many of which have sensible defaults, xor via the \code{data}
#' argument which takes precedent when present and must be compatible with the
#' Stan model directly.
#'
#' The \code{data} argument can offer greater flexibility, for example
#' the times and terms that define the step size used in the model
#' are currently constrained to times = years and terms in
#' \code{c("year", "month", "quarter")} by the helper functions. By specifying
#' \code{tags_released} and \code{tags_recovered} manually in the \code{data}
#' argument and adjusting the corresponding index limits, arbitrary time and
#' step units are possible.
#'
#' @param tag_data [data.frame()]
#' @param list_regions [list()]
#' @param list_sizes [list()]
#' @param date_released_start [character()] amenable to
#'   \code{lubridate::as_date()}
#' @param date_released_end [character()] amenable to
#'   \code{lubridate::as_date()}
#' @param step_liberty_max [integer()]
#' @param term_released [character()] in \code{c("year", "month", "quarter")}
#' @param colname_date_released [character()]
#' @param colname_date_recovered [character()]
#' @param colname_region_released [character()]
#' @param colname_region_recovered [character()]
#' @param colname_size_released [character()]
#' @param movement_pattern [integer()] One of \code{1}: movement between
#'   all pairs of regions at each step, or \code{2}:
#'   movement between numerically sequential regions (the default)
#' @param movement_allow [integer()] [matrix()] Adjusts the movement pattern
#'   to allow additional movement between regions. Each row indicates
#'   directional movement between a pair of regions
#' @param movement_disallow [integer()] [matrix()] As for
#'   \code{movement_allow}, but specified movement is disallowed
#' @param mu_mortality_rate [numeric()]
#' @param mu_reporting_rate [numeric()]
#' @param mu_fishing_mean_rate [numeric()]
#' @param mu_initial_loss_rate [numeric()]
#' @param mu_ongoing_loss_rate [numeric()]
#' @param mu_dispersion [numeric()]
#' @param mu_fishing_deviation_time_rate [numeric()]
#' @param sd_mortality_rate [numeric()]
#' @param sd_reporting_rate [numeric()]
#' @param sd_fishing_mean_rate [numeric()]
#' @param sd_initial_loss_rate [numeric()]
#' @param sd_ongoing_loss_rate [numeric()]
#' @param sd_dispersion [numeric()]
#' @param cv_movement_deviation_time_rate [numeric()]
#' @param cv_movement_deviation_term_rate [numeric()]
#' @param cv_movement_deviation_size_rate [numeric()]
#' @param cv_fishing_deviation_time_rate [numeric()]
#' @param cv_fishing_deviation_term_rate [numeric()]
#' @param cv_fishing_deviation_size_rate [numeric()]
#' @param expected_fudge [numeric()]
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
fit <- function (tag_data,
                 # Tag arguments
                 list_regions,
                 list_sizes,
                 date_released_start,
                 date_released_end,
                 step_liberty_max = NULL,
                 term_released = "quarter",
                 colname_date_released = "date_released",
                 colname_date_recovered = "date_recovered",
                 colname_region_released = "region_released",
                 colname_region_recovered = "region_recovered",
                 colname_size_released = "size_released",
                 # Movement index
                 movement_pattern = 2,
                 movement_allow = NULL,
                 movement_disallow = NULL,
                 # Prior means
                 mu_mortality_rate = rep(0.1, length(list_regions)),
                 mu_reporting_rate = rep(1, length(list_regions)),
                 mu_fishing_mean_rate = rep(0.1, length(list_regions)),
                 mu_initial_loss_rate = 0.1,
                 mu_ongoing_loss_rate = 0.02,
                 mu_dispersion = 1,
                 # Fishing prior mean deviations by time
                 mu_fishing_deviation_time_rate = NULL, # [T, X]
                 # Prior standard deviations
                 sd_mortality_rate = 0.1 * mu_mortality_rate,
                 sd_reporting_rate = 0.1 * mu_reporting_rate,
                 sd_fishing_mean_rate = 0.1 * mu_fishing_mean_rate,
                 sd_initial_loss_rate = 0.1 * mu_initial_loss_rate,
                 sd_ongoing_loss_rate = 0.1 * mu_ongoing_loss_rate,
                 sd_dispersion = 0.25 * mu_dispersion,
                 # Movement prior coefficients of variation
                 cv_movement_deviation_time_rate = 0.1,
                 cv_movement_deviation_term_rate = 0.1,
                 cv_movement_deviation_size_rate = 0.1,
                 # Fishing prior coefficients of variation
                 cv_fishing_deviation_time_rate = 0.1,
                 cv_fishing_deviation_term_rate = 0.1,
                 cv_fishing_deviation_size_rate = 0.1,
                 # Fudge value
                 expected_fudge = 1e-12,
                 # CmdStanR
                 data = NULL,
                 chains = 1,
                 step_size = 0.01,
                 iter_warmup = 250,
                 iter_sampling = 750,
                 use_reduce_sum = FALSE,
                 threads_per_chain = 8,
                 ...) {

  # Check arguments ------------------------------------------------------------

  # Has data?
  has_data <- ifelse(is.null(data), FALSE, TRUE)
  # Tag data
  checkmate::assert_data_frame(
    tag_data,
    all.missing = FALSE,
    null.ok = has_data
  )
  # Tag arguments
  checkmate::assert_list(
    list_regions,
    any.missing = FALSE,
    min.len = 2
  )
  checkmate::assert_list(
    list_sizes,
    any.missing = FALSE,
    min.len = 4
  )
  checkmate::assert_date(
    lubridate::as_date(date_released_start),
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_date(
    lubridate::as_date(date_released_end),
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    step_liberty_max,
    lower = 2,
    any.missing = FALSE,
    len = 1,
    null.ok = TRUE
  )
  checkmate::assert_choice(
    term_released,
    choices = c("year", "quarter", "month")
  )
  checkmate::assert_character(
    colname_date_released,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_character(
    colname_date_recovered,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_character(
    colname_region_released,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_character(
    colname_region_recovered,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_character(
    colname_size_released,
    any.missing = FALSE,
    len = 1
  )
  # Movement index
  checkmate::assert_choice(
    movement_pattern,
    choices =  c(1L, 2L)
  )
  checkmate::assert_matrix(
    movement_allow,
    mode = "integerish",
    any.missing = FALSE,
    ncols = 2,
    null.ok = TRUE
  )
  checkmate::assert_matrix(
    movement_disallow,
    mode = "integerish",
    any.missing = FALSE,
    ncols = 2,
    null.ok = TRUE
  )
  # Prior means
  checkmate::assert_numeric(
    mu_mortality_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_numeric(
    mu_reporting_rate,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_numeric(
    mu_fishing_mean_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_number(
    mu_initial_loss_rate,
    lower = 0
  )
  checkmate::assert_number(
    mu_ongoing_loss_rate,
    lower = 0
  )
  checkmate::assert_number(
    mu_dispersion,
    lower = 0
  )
  # Fishing prior mean deviations by time
  checkmate::assert_array(
    mu_fishing_deviation_time_rate,
    mode = double,
    any.missing = FALSE,
    d = 2,
    null.ok = TRUE
  )
  # Prior standard deviations
  checkmate::assert_numeric(
    sd_mortality_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_numeric(
    sd_reporting_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_numeric(
    sd_fishing_mean_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = length(list_regions)
  )
  checkmate::assert_number(
    sd_initial_loss_rate,
    lower = 0
  )
  checkmate::assert_number(
    sd_ongoing_loss_rate,
    lower = 0
  )
  checkmate::assert_number(
    sd_dispersion,
    lower = 0
  )
  # Movement prior coefficients of variation
  checkmate::assert_number(
    cv_movement_deviation_time_rate,
    lower = 0
  )
  checkmate::assert_number(
    cv_movement_deviation_term_rate,
    lower = 0
  )
  checkmate::assert_number(
    cv_movement_deviation_size_rate,
    lower = 0
  )
  # Fishing prior coefficients of variation
  checkmate::assert_number(
    cv_fishing_deviation_time_rate,
    lower = 0
  )
  checkmate::assert_number(
    cv_fishing_deviation_term_rate,
    lower = 0
  )
  checkmate::assert_number(
    cv_fishing_deviation_size_rate,
    lower = 0
  )
  # Fudge value
  checkmate::assert_number(
    expected_fudge,
    lower = 0
  )
  # Data
  checkmate::assert_list(
    data,
    types = c("integer", "double"),
    any.missing = FALSE,
    null.ok = TRUE
  )
  # CmdStanR arguments
  checkmate::assert_integerish(
    chains,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_number(
    step_size,
    lower = 0,
    finite = TRUE
  )
  checkmate::assert_integerish(
    iter_warmup,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    iter_sampling,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_logical(
    use_reduce_sum,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    threads_per_chain,
    lower = 1,
    any.missing = FALSE,
    len = 1
  )

  # Print messages -------------------------------------------------------------

  if (has_data) {
    cat("\nfit(): data argument present--ignoring most other arguments\n")
  }

  # Assemble tag data ----------------------------------------------------------

  # Released
  if (is.null(data)) {
    tags_released <- create_tags_released(
      tags = tag_data,
      list_regions = list_regions,
      list_sizes = list_sizes,
      date_released_start = date_released_start,
      date_released_end = date_released_end,
      # step_liberty_max = step_liberty_max,
      term_released = term_released,
      colname_date_released = colname_date_released,
      # colname_date_recovered = colname_date_recovered,
      colname_region_released = colname_region_released,
      # colname_region_recovered = colname_region_recovered,
      colname_size_released = colname_size_released
    )
  }
  # Recovered
  if (is.null(data)) {
    tags_recovered <- create_tags_recovered(
      tags = tag_data,
      list_regions = list_regions,
      list_sizes = list_sizes,
      date_released_start = date_released_start,
      date_released_end = date_released_end,
      step_liberty_max = step_liberty_max,
      term_released = term_released,
      colname_date_released = colname_date_released,
      colname_date_recovered = colname_date_recovered,
      colname_region_released = colname_region_released,
      colname_region_recovered = colname_region_recovered,
      colname_size_released = colname_size_released
    )
  }

  # Assemble index limits ------------------------------------------------------

  n_regions <- length(list_regions)
  n_times <- compute_n_times(date_released_start, date_released_end)
  n_terms <- compute_n_terms(term_released)
  n_sizes <- length(list_sizes)
  n_steps <- (n_times * n_terms) - 1
  n_liberty <- ifelse(is.null(step_liberty_max), n_steps, step_liberty_max)

  # Assemble movement index ----------------------------------------------------

  movement_index <- create_movement_index(
    n_regions = n_regions,
    movement_pattern = movement_pattern,
    movement_allow = movement_allow,
    movement_disallow = movement_disallow
  )

  # Assemble prior mean fishing deviation time rate ----------------------------

  if (is.null(data)) {
    if (is.null(mu_fishing_deviation_time_rate)) {
      mu_fishing_deviation_time_rate <- array(0.0, dim = c(n_times, n_regions))
    }
  }

  # Assemble data --------------------------------------------------------------

  if (is.null(data)) {
    data <- list(
      # Index limits
      X = n_regions, # Number of geographic regions
      T = n_times, # Number of times (usually years; release only)
      I = n_terms, # Number of terms per unit of time
      S = n_sizes, # Number of size classes
      N = n_steps, # Number of release steps ((T * I) - 1)
      L = n_liberty, # Number of maximum steps at liberty
      P = sum(movement_index), # Number of movement rate mean parameters
      # Tag data
      tags_released = tags_released, # array[N, S, X]
      tags_recovered = tags_recovered, # array[N, S, X, L, X]
      # Movement index array
      movement_index = movement_index, # array[X, X]
      # Prior means
      mu_mortality_rate = mu_mortality_rate, # array[X]
      mu_reporting_rate = mu_reporting_rate, # array[X]
      mu_fishing_mean_rate = mu_fishing_mean_rate, # array[X]
      mu_initial_loss_rate = mu_initial_loss_rate,
      mu_ongoing_loss_rate = mu_ongoing_loss_rate,
      mu_dispersion = mu_dispersion,
      # Fishing prior mean deviations by time
      mu_fishing_deviation_time_rate = mu_fishing_deviation_time_rate, # [T, X]
      # Prior standard deviations
      sd_mortality_rate = sd_mortality_rate, # array[X]
      sd_reporting_rate = sd_reporting_rate, # array[X]
      sd_fishing_mean_rate = sd_fishing_mean_rate, # array[X]
      sd_initial_loss_rate = sd_initial_loss_rate,
      sd_ongoing_loss_rate = sd_ongoing_loss_rate,
      sd_dispersion = sd_dispersion,
      # Movement prior coefficients of variation
      cv_movement_deviation_time_rate = cv_movement_deviation_time_rate,
      cv_movement_deviation_term_rate = cv_movement_deviation_term_rate,
      cv_movement_deviation_size_rate = cv_movement_deviation_size_rate,
      # Fishing prior coefficients of variation
      cv_fishing_deviation_time_rate = cv_fishing_deviation_time_rate,
      cv_fishing_deviation_term_rate = cv_fishing_deviation_term_rate,
      cv_fishing_deviation_size_rate = cv_fishing_deviation_size_rate,
      # Fudge values
      expected_fudge = expected_fudge
    )
  }

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
