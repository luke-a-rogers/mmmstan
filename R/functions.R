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
#' @param year_released_start [integer()] year that the first tag was released
#' @param year_recovered_end [integer()] year that the last tag was recovered
#' @param step_interval [character()] in \code{c("month", "quarter", "year")}
#' @param step_liberty_max [integer()]
#' @param terms_unique_movement [integer()]
#' @param nest_terms_within_times [logical()] nest terms within times?
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
#' @param mu_fishing_rate [numeric()][array()]
#' @param mu_mortality_rate [numeric()]
#' @param mu_reporting_rate [numeric()]
#' @param mu_selectivity [numeric()]
#' @param mu_initial_loss_rate [numeric()]
#' @param mu_ongoing_loss_rate [numeric()]
#' @param mu_dispersion [numeric()]
#' @param cv_fishing_rate [numeric()]
#' @param cv_mortality_rate [numeric()]
#' @param cv_reporting_rate [numeric()]
#' @param cv_selectivity [numeric()]
#' @param cv_initial_loss_rate [numeric()]
#' @param cv_ongoing_loss_rate [numeric()]
#' @param cv_dispersion [numeric()]
#' @param cv_movement_time_deviation [numeric()]
#' @param cv_movement_term_deviation [numeric()]
#' @param cv_movement_size_deviation [numeric()]
#' @param tolerance_expected [numeric()]
#' @param tolerance_movement [numeric()]
#' @param tolerance_fishing [numeric()]
#' @param data [list()] See details
#' @param chains [integer()] number of chains
#' @param step_size [numeric()] initial step size
#' @param adapt_delta [numeric()] the adaptation target acceptance statistic
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @details TBD
#'
#' @return [mmmstan::mmmstan()]
#' @export
#'
mmmstan <- function (tag_data,
                     # Tag arguments
                     list_regions,
                     list_sizes,
                     year_released_start,
                     year_recovered_end,
                     step_interval = "month",
                     step_liberty_max = NULL,
                     terms_unique_movement = 1L,
                     nest_terms_within_times = FALSE,
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
                     mu_fishing_rate = array(
                       0.1,
                       dim = c(
                         year_recovered_end - year_released_start + 1,
                         length(list_regions)
                       )
                     ),
                     mu_mortality_rate = rep(0.1, length(list_regions)),
                     mu_reporting_rate = rep(1, length(list_regions)),
                     mu_selectivity = rep(1, length(list_sizes) - 1L),
                     mu_initial_loss_rate = 0.1,
                     mu_ongoing_loss_rate = 0.02,
                     mu_dispersion = 1,
                     # Prior coefficients of variation
                     cv_fishing_rate = 0.1,
                     cv_mortality_rate = 0.1,
                     cv_reporting_rate = 0.1,
                     cv_selectivity = 0.1,
                     cv_initial_loss_rate = 0.1,
                     cv_ongoing_loss_rate = 0.1,
                     cv_dispersion = 0.25,
                     # Movement prior coefficients of variation
                     cv_movement_time_deviation = 0.1,
                     cv_movement_term_deviation = 0.1,
                     cv_movement_size_deviation = 0.1,
                     # Tolerance value
                     tolerance_expected = 1e-12,
                     tolerance_movement = 1e-12,
                     tolerance_fishing = 1e-12,
                     # CmdStanR
                     data = NULL,
                     chains = 1,
                     step_size = 0.01,
                     adapt_delta = 0.8,
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
    min.len = 2,
    null.ok = has_data
  )
  checkmate::assert_list(
    list_sizes,
    any.missing = FALSE,
    min.len = 1,
    null.ok = has_data
  )
  checkmate::assert_integerish(
    year_released_start,
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_integerish(
    year_recovered_end,
    lower = year_released_start + 1L,
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_choice(
    step_interval,
    choices = c("year", "quarter", "month"),
    null.ok = has_data
  )
  checkmate::assert_integerish(
    step_liberty_max,
    lower = 2,
    any.missing = FALSE,
    len = 1,
    null.ok = TRUE
  )
  checkmate::assert_integerish(
    terms_unique_movement,
    lower = 1L,
    upper = (year_recovered_end - year_released_start + 1) *
      compute_steps_per_year(step_interval),
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_logical(
    nest_terms_within_times,
    any.missing = FALSE,
    len = 1L,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_date_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_date_recovered,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_region_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_region_recovered,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  checkmate::assert_character(
    colname_size_released,
    any.missing = FALSE,
    len = 1,
    null.ok = has_data
  )
  # Movement index
  checkmate::assert_choice(
    movement_pattern,
    choices =  c(1L, 2L),
    null.ok = has_data
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
  checkmate::assert_array(
    mu_fishing_rate,
    mode = "double",
    any.missing = FALSE,
    d = 2,
    null.ok = has_data
  )
  checkmate::assert_numeric(
    mu_fishing_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    null.ok = has_data
  )
  checkmate::assert_numeric(
    mu_mortality_rate,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    null.ok = has_data
  )
  checkmate::assert_numeric(
    mu_reporting_rate,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    null.ok = has_data
  )
  checkmate::assert_numeric(
    mu_selectivity,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    null.ok = has_data
  )
  checkmate::assert_number(
    mu_initial_loss_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    mu_ongoing_loss_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    mu_dispersion,
    lower = 0,
    null.ok = has_data
  )
  # Fishing prior mean deviations by time
  checkmate::assert_number(
    cv_fishing_rate,
    lower = 0,
    null.ok = has_data
  )
  # Prior standard deviations
  checkmate::assert_number(
    cv_mortality_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_reporting_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_selectivity,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_initial_loss_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_ongoing_loss_rate,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_dispersion,
    lower = 0,
    null.ok = has_data
  )
  # Movement prior coefficients of variation
  checkmate::assert_number(
    cv_movement_time_deviation,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_movement_term_deviation,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    cv_movement_size_deviation,
    lower = 0,
    null.ok = has_data
  )
  # Tolerance value
  checkmate::assert_number(
    tolerance_expected,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    tolerance_movement,
    lower = 0,
    null.ok = has_data
  )
  checkmate::assert_number(
    tolerance_fishing,
    lower = 0,
    null.ok = has_data
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
  checkmate::assert_number(
    adapt_delta,
    lower = 0,
    upper = 1,
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

  # Assemble index limits ------------------------------------------------------

  if (is.null(data)) {
    n_times <- year_recovered_end - year_released_start + 1L
    n_steps <- n_times * compute_steps_per_year(step_interval)
    n_terms <- terms_unique_movement
    n_sizes <- length(list_sizes)
    n_liberty <- ifelse(is.null(step_liberty_max), n_steps-1L, step_liberty_max)
    n_regions <- length(list_regions)
  }

  # Assemble tag data ----------------------------------------------------------

  if (is.null(data)) {
    tag_array <- create_tag_array(
      tag_data = tag_data,
      list_regions = list_regions,
      list_sizes = list_sizes,
      year_released_start = year_released_start,
      year_recovered_end = year_recovered_end,
      step_liberty_max = n_liberty,
      step_interval = step_interval,
      colname_date_released = colname_date_released,
      colname_date_recovered = colname_date_recovered,
      colname_region_released = colname_region_released,
      colname_region_recovered = colname_region_recovered,
      colname_size_released = colname_size_released
    )
  }

  # Assemble movement index ----------------------------------------------------

  if (is.null(data)) {
    movement_index <- create_movement_index(
      n_regions = n_regions,
      movement_pattern = movement_pattern,
      movement_allow = movement_allow,
      movement_disallow = movement_disallow
    )
  }

  # Assemble index vectors -----------------------------------------------------

  if (is.null(data)) {
    step_to_time <- create_step_to_time(n_steps, n_times)
    step_to_term <- create_step_to_term(
      n_steps,
      n_times,
      n_terms,
      nest_terms_within_times
    )
  }

  # Assemble data --------------------------------------------------------------

  if (is.null(data)) {
    data <- list(
      # Index limits
      N = n_steps, # Number of model steps
      T = n_times, # Number of times (years)
      I = n_terms, # Number of unique terms (seasons)
      S = n_sizes, # Number of size classes
      L = n_liberty, # Number of maximum steps at liberty
      X = n_regions, # Number of geographic regions
      P = sum(movement_index) - n_regions, # Number of movement rate mean parameters
      # Tag data
      tags = tag_array, # [N, S, L, X, X]
      # Movement index array
      movement_index = movement_index, # [X, X]
      step_to_time = step_to_time, # [N]
      step_to_term = step_to_term, # [N]
      # Prior means
      mu_fishing_rate = mu_fishing_rate, # [T, X]
      mu_mortality_rate = mu_mortality_rate, # [X]
      mu_reporting_rate = mu_reporting_rate, # [X]
      mu_selectivity = mu_selectivity, # [S - 1]
      mu_initial_loss_rate = mu_initial_loss_rate,
      mu_ongoing_loss_rate = mu_ongoing_loss_rate,
      mu_dispersion = mu_dispersion,
      # Prior coefficients of variation
      cv_fishing_rate = cv_fishing_rate,
      cv_mortality_rate = cv_mortality_rate,
      cv_reporting_rate = cv_reporting_rate,
      cv_selectivity = cv_selectivity,
      cv_initial_loss_rate = cv_initial_loss_rate,
      cv_ongoing_loss_rate = cv_ongoing_loss_rate,
      cv_dispersion = cv_dispersion,
      # Movement prior coefficients of variation
      cv_movement_time_deviation = cv_movement_time_deviation,
      cv_movement_term_deviation = cv_movement_term_deviation,
      cv_movement_size_deviation = cv_movement_size_deviation,
      # Tolerance values
      tolerance_expected = tolerance_expected,
      tolerance_movement = tolerance_movement,
      tolerance_fishing = tolerance_fishing
    )
  }

  # Check dimensions data elements ---------------------------------------------

  checkmate::assert_true(data$P == sum(data$movement_index) - data$X)

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
    adapt_delta = adapt_delta,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    # refresh = 100 # ,
    # ...
  )

  # Return values --------------------------------------------------------------

  structure(list(
    data = data,
    fit = fit),
    class = "mmmstan")
}
