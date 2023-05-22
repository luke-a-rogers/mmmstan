#' Count Model Steps
#'
#' @param year_start [integer()] year of first tag release
#' @param year_end [integer()] year of final tag recovery
#' @param step_interval [character()] one of \code{"month"}, \code{"quarter"}
#'   or \code{"year"}
#'
#' @return [integer()]
#' @export
#'
#' @examples
#'
#' count_model_steps(2011, 2018, "month")
#' count_model_steps(2011, 2018, "quarter")
#' count_model_steps(2011, 2018, "year")
#'
count_model_steps <- function (year_start,
                               year_end,
                               step_interval) {

  # Check arguments ------------------------------------------------------------

  # Compute model steps --------------------------------------------------------

  # Count steps per year
  steps_per_year <- count_intervals(step_interval, "year")
  # Count model years
  years <- year_end - year_start + 1L
  # Count model steps
  model_steps <- as.integer(years * steps_per_year)

  # Return model steps ---------------------------------------------------------

  return(model_steps)
}

#' Count Intervals
#'
#' @param a [character()] one of \code{"month"}, \code{"quarter"}
#' @param b [character()] one of \code{"month"}, \code{"quarter"}
#'
#' @return [integer()]
#' @export
#'
#' @examples
#'
#' count_intervals("month", "month")
#' count_intervals("month", "quarter")
#' count_intervals("month", "year")
#' count_intervals("quarter", "quarter")
#' count_intervals("quarter", "year")
#' count_intervals("year", "year")
#'
count_intervals <- function (a, b) {

  # Define reference options ---------------------------------------------------

  options <- c("month", "quarter", "year")
  offsets <- c(1L, 3L, 12L)

  # Check arguments ------------------------------------------------------------

  checkmate::assert_choice(a, options)
  checkmate::assert_choice(b, options)

  # Check nesting of intervals -------------------------------------------------

  if (match(a, options) > match(b, options)) stop("a must be <= b")

  # Compute value --------------------------------------------------------------

  value <- offsets[match(b, options)] %/% offsets[match(a, options)]

  # Return value ---------------------------------------------------------------

  return(value)
}

#' Index N to T
#'
#' @param year_start [integer()] year of first tag release
#' @param year_end [integer()] year of final tag recovery
#' @param step_interval [character()] one of \code{"month"}, \code{"quarter"}
#'   or \code{"year"}
#'
#' @return [integer()][vector()]
#' @export
#'
#' @examples
#'
#' index_n_to_t(2011, 2018, "month")
#' index_n_to_t(2011, 2018, "quarter")
#' index_n_to_t(2011, 2018, "year")
#'
index_n_to_t <- function (year_start, year_end, step_interval) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(
    year_start,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    year_end,
    lower = year_start + 1L,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_choice(step_interval, c("month", "quarter", "year"))

  # Compute values -------------------------------------------------------------

  model_steps <- count_model_steps(year_start, year_end, step_interval)
  model_years <- year_end - year_start + 1
  steps_per_year <- count_intervals(step_interval, "year")

  # Populate index -------------------------------------------------------------

  index <- rep(seq_len(model_years), each=steps_per_year)[seq_len(model_steps)]

  # Return index ---------------------------------------------------------------

  return(index)
}

#' Index N to K
#'
#' @param year_start [integer()] year of first tag release
#' @param year_end [integer()] year of final tag recovery
#' @param step_interval [character()] one of \code{"month"}, \code{"quarter"}
#'   or \code{"year"}
#' @param term_interval [character()] one of \code{"month"}, \code{"quarter"}
#'   or \code{"year"}
#'
#' @return [integer()][vector()]
#' @export
#'
#' @examples
#'
#' index_n_to_k(2011, 2018, "month", "month")
#' index_n_to_k(2011, 2018, "month", "quarter")
#' index_n_to_k(2011, 2018, "month", "year")
#' index_n_to_k(2011, 2018, "quarter", "quarter")
#' index_n_to_k(2011, 2018, "quarter", "year")
#' index_n_to_k(2011, 2018, "year", "year")
#'
index_n_to_k <- function (year_start, year_end, step_interval, term_interval) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(
    year_start,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_integerish(
    year_end,
    lower = year_start + 1L,
    any.missing = FALSE,
    len = 1
  )
  checkmate::assert_choice(step_interval, c("month", "quarter", "year"))
  checkmate::assert_choice(term_interval, c("month", "quarter", "year"))

  # Compute values -------------------------------------------------------------

  model_steps <- count_model_steps(year_start, year_end, step_interval)
  model_terms <- count_intervals(term_interval, "year")
  model_years <- year_end - year_start + 1L
  steps_per_term <- count_intervals(step_interval, term_interval)

  # Populate index -------------------------------------------------------------

  index_within_year <- rep(seq_len(model_terms), each = steps_per_term)
  index <- rep(index_within_year, model_years)[seq_len(model_steps)]

  # Return index ---------------------------------------------------------------

  return(index)
}


#' Compute The Number of Years
#'
#' @param a [character()] date as \code{"\%Y-\%m-\%d"}
#' @param b [character()] date as \code{"\%Y-\%m-\%d"}
#'
#' @return [integer()]
#' @export
#'
#' @examples
#'
#' compute_n_years("2011-01-01", "2014-12-31")
#'
compute_n_years <- function (a, b) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_date(lubridate::date(a))
  checkmate::assert_date(lubridate::date(b))

  # Return the number of years -------------------------------------------------

  return(as.integer(abs(lubridate::year(b) - lubridate::year(a)) + 1))
}

#' Create A Movement Index Matrix
#'
#' @param number_of_regions [integer()] Number of regions
#' @param movement_pattern [integer()] One of \code{1}: movement between
#'   all pairs of regions at each step, or \code{2}:
#'   movement between numerically sequential regions (the default)
#' @param movement_allow [integer()] [matrix()] Adjusts the movement pattern
#'   to allow additional movement between regions. Each row indicates
#'   directional movement between a pair of regions
#' @param movement_disallow [integer()] [matrix()] As for
#'   \code{movement_allow}, but specified movement is disallowed
#'
#' @return A square matrix of zeros and ones
#' @export
#'
#' @examples
#'
#' # Neighbors (default)
#' create_movement_index(6)
#'
#' # All pairwise movement
#' create_movement_index(6, 1)
#'
#' # Neighbors plus 1-6, 2-5, and 5-2, but not 6-1
#' movement_allow <- matrix(c(1,6,2,5,5,2), ncol = 2, byrow = TRUE)
#' create_movement_index(6, 2, movement_allow)
#'
#' # All pairwise movement but not 1-6
#' movement_disallow <- matrix(c(1,6), ncol = 2, byrow = TRUE)
#' create_movement_index(6, 1, movement_disallow = movement_disallow)
#'
create_movement_index <- function (number_of_regions,
                                   movement_pattern = 2,
                                   movement_allow = NULL,
                                   movement_disallow = NULL) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(
    number_of_regions,
    lower = 2,
    len = 1,
    any.missing = FALSE
  )
  checkmate::assert_integerish(
    movement_pattern,
    lower = 0,
    len = 1,
    any.missing = FALSE
  )
  checkmate::assert_integerish(
    movement_allow,
    lower = 1,
    upper = number_of_regions,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_matrix(
    movement_allow,
    min.rows = 1,
    ncols = 2,
    null.ok = TRUE
  )
  checkmate::assert_integerish(
    movement_disallow,
    lower = 1,
    upper = number_of_regions,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_matrix(
    movement_disallow,
    min.rows = 1,
    ncols = 2,
    null.ok = TRUE
  )

  # Initialize matrix ----------------------------------------------------------

  movement_index <- diag(
    1L,
    nrow = number_of_regions,
    ncol = number_of_regions,
    names = FALSE
  )

  # Add movement pattern -------------------------------------------------------

  # All pairwise
  if (movement_pattern == 1) {
    movement_index <- matrix(
      1L,
      nrow = number_of_regions,
      ncol = number_of_regions
    )
  }
  # Neighbors only
  if (movement_pattern == 2) {
    for (i in seq_len(number_of_regions - 1L)) {
      movement_index[i, i + 1L] <- 1L
      movement_index[i + 1L, i] <- 1L
    }
  }

  # Allow additional pairwise movement -----------------------------------------

  if (!is.null(movement_allow)) {
    movement_index[movement_allow] <- 1L
  }

  # Disallow specified pairwise movement ---------------------------------------

  if (!is.null(movement_disallow)) {
    movement_index[movement_disallow] <- 0L
  }

  # Return movement index ------------------------------------------------------

  return(movement_index)
}

#' Create Groups from Index List
#'
#' @param x [atomic()] [vector()] of values
#' @param list_x [list()] of groups of x
#'
#' @return [integer()] [vector()] of numeric index of groups
#' @export
#'
#' @examples
#'
#' # Numeric
#' x <- c(1, 4, 2, 3, NA, 2, 7)
#' g <- list(a = 1:2, b = 3:4)
#' create_group(x, g)
#'
#' # Character
#' x <- c("M", "F", "F", NA, "M", "N", "F")
#' g <- list(m = "M", f = "F")
#' create_group(x, g)
#'
#' # Factor
#' x <- factor(c("L", "S", "M", NA, "M", "H"))
#' g <- list(s = "S", m = "M", l = "L")
#' create_group(x, g)
#'
create_group <- function (x, list_x) {

  # Check arguments ------------------------------------------------------------

  # Define size released -------------------------------------------------------

  x <- data.frame(value = x)
  list_sequence <- seq_along(list_x)
  x_key <- data.frame(
    value = unlist(list_x),
    index = rep(list_sequence, lengths(list_x))
  )
  x_groups <- as.integer(
    dplyr::left_join(x, x_key, by = "value")$index
  )

  # Return size classes --------------------------------------------------------

  return(x_groups)
}

#' Create Step Duration At Large Max
#'
#' @param date_released [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_recovered [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_recovered_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param step_interval [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#' @param step_duration_max [integer()] maximum steps duration at large
#'
#' @return [integer()] [vector()]
#' @export
#'
#' @examples
#'
#' # Yearly
#' u <- rep("2011-01-01", 5)
#' v <- c("2011-01-01", "2012-01-01", NA, "2020-01-01", "2021-01-01")
#' w <- c("2011-01-01")
#' x <- c("2020-01-01")
#' y  <- "year"
#' z <- 10
#' create_step_duration_max(u, v, w, x, y, z)
#'
#' # Quarterly
#' u <- rep("2011-01-01", 5)
#' v <- c("2011-01-01", "2012-12-31", NA, "2020-12-31", "2021-01-01")
#' w <- c("2011-01-01")
#' x <- c("2020-01-01")
#' y <- "quarter"
#' z <- 40
#' create_step_duration_max(u, v, w, x, y, z)
#'
#' # Monthly
#' u <- rep("2011-01-01", 5)
#' v <- c("2011-01-01", "2011-02-01", NA, "2020-12-31", "2021-01-01")
#' w <- c("2011-01-01")
#' x <- c("2020-01-01")
#' y <- "month"
#' z <- 120
#' create_step_duration_max(u, v, w, x, y, z)
#'
create_step_duration_max <- function (date_released,
                                      date_recovered,
                                      date_released_start,
                                      date_recovered_end,
                                      step_interval,
                                      step_duration_max) {

  # Check arguments ------------------------------------------------------------

  # Convert to dates -----------------------------------------------------------

  date_released <- lubridate::as_date(date_released)
  date_recovered <- lubridate::as_date(date_recovered)
  date_released_start <- lubridate::as_date(date_released_start)
  date_recovered_end <- lubridate::as_date(date_recovered_end)

  # Replace disallowed dates by NAs --------------------------------------------

  date_released[which(date_released < date_released_start)] <- NA_real_
  date_released[which(date_released > date_recovered_end)] <- NA_real_

  # Create date released parts -------------------------------------------------

  # Released
  year_released <- lubridate::year(date_released)
  quarter_released <- lubridate::quarter(date_released)
  month_released <- lubridate::month(date_released)
  # Recovered
  year_recovered <- lubridate::year(date_recovered)
  quarter_recovered <- lubridate::quarter(date_recovered)
  month_recovered <- lubridate::month(date_recovered)

  # Create step liberty --------------------------------------------------------

  if (step_interval == "year") {
    step_duration <- year_recovered - year_released + 1
  } else if (step_interval == "quarter") {
    step_duration <- 4 * (year_recovered - year_released) +
      (quarter_recovered - quarter_released) + 1
  } else if (step_interval == "month") {
    step_duration <- 12 * (year_recovered - year_released) +
      (month_recovered - month_released) + 1
  } else {
    stop("step_interval must be 'year', 'quarter', or 'month'")
  }

  # Replace disallowed dates by NAs --------------------------------------------

  step_duration[which(step_duration < 2)] <- NA_real_
  step_duration[which(step_duration > step_duration_max)] <- NA_real_

  # Return step liberty --------------------------------------------------------

  return(as.integer(step_duration))
}

#' Create Step Released
#'
#' @param date_released [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_recovered_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param step_interval [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#'
#' @return [integer()] [vector()] step released
#' @export
#'
#' @examples
#'
#' x <- c("2010-01-01", "2011-01-01", NA, "2020-12-31", "2021-01-01")
#' d0 <- c("2011-01-01")
#' d1 <- c("2020-12-31")
#' create_step_released(x, d0, d1, "year")
#' create_step_released(x, d0, d1, "quarter")
#' create_step_released(x, d0, d1, "month")
#'
create_step_released <- function (date_released,
                                  date_released_start,
                                  date_recovered_end,
                                  step_interval) {

  # Check arguments ------------------------------------------------------------

  # Convert to dates -----------------------------------------------------------

  date_released <- lubridate::as_date(date_released)
  date_released_start <- lubridate::as_date(date_released_start)
  date_recovered_end <- lubridate::as_date(date_recovered_end)

  # Replace disallowed dates by NAs --------------------------------------------

  date_released[which(date_released < date_released_start)] <- NA_real_
  date_released[which(date_released > date_recovered_end)] <- NA_real_

  # Create date released parts -------------------------------------------------

  # Released
  year_released <- lubridate::year(date_released)
  quarter_released <- lubridate::quarter(date_released)
  month_released <- lubridate::month(date_released)
  # Released start
  year_start <- lubridate::year(date_released_start)
  quarter_start <- lubridate::quarter(date_released_start)
  month_start <- lubridate::month(date_released_start)

  # Create step released -------------------------------------------------------

  if (step_interval == "year") {
    step_released <- year_released - year_start + 1
  } else if (step_interval == "quarter") {
    step_released <- 4 * (year_released - year_start) +
      (quarter_released - quarter_start) + 1
  } else if (step_interval == "month") {
    step_released <- 12 * (year_released - year_start) +
      (month_released - month_start) + 1
  } else {
    stop("step_interval must be 'year', 'quarter', or 'month'")
  }

  # Return step released -------------------------------------------------------

  return(as.integer(step_released))
}

#' Create Tag Array
#'
#' @param tag_data [data.frame()]
#' @param list_regions [list()]
#' @param list_sizes [list()]
#' @param year_released_start [integer()] year that the first tag was released
#' @param year_recovered_end [integer()] year the last tag was recovered
#' @param step_interval [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#' @param step_duration_max [integer()] optional constraint on the maximum
#'   steps duration at large
#' @param colname_date_released [character()]
#' @param colname_date_recovered [character()]
#' @param colname_region_released [character()]
#' @param colname_region_recovered [character()]
#' @param colname_size_released [character()]
#'
#' @importFrom rlang .data
#'
#' @return [array()]
#' @export
#'
create_tag_array <- function (tag_data,
                              list_regions,
                              list_sizes,
                              year_released_start,
                              year_recovered_end,
                              step_interval = "month",
                              step_duration_max = NULL,
                              colname_date_released = "date_released",
                              colname_date_recovered = "date_recovered",
                              colname_region_released = "region_released",
                              colname_region_recovered = "region_recovered",
                              colname_size_released = "size_released") {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_data_frame(tag_data)
  checkmate::assert_list(list_sizes)
  checkmate::assert_list(list_regions)
  checkmate::assert_integerish(
    year_released_start,
    any.missing = FALSE,
    len = 1L
  )
  checkmate::assert_integerish(
    year_recovered_end,
    lower = year_released_start + 1L,
    any.missing = FALSE,
    len = 1L
  )
  checkmate::assert_choice(
    x = step_interval,
    choices = c("year", "quarter", "month")
  )
  checkmate::assert_integerish(
    step_duration_max,
    lower = 2,
    len = 1,
    any.missing = FALSE,
    null.ok = TRUE
  )
  checkmate::assert_choice(
    x = colname_date_released,
    choices = colnames(tag_data)
  )
  checkmate::assert_choice(
    x = colname_date_recovered,
    choices = colnames(tag_data)
  )
  checkmate::assert_choice(
    x = colname_region_released,
    choices = colnames(tag_data)
  )
  checkmate::assert_choice(
    x = colname_region_recovered,
    choices = colnames(tag_data)
  )
  checkmate::assert_choice(
    x = colname_size_released,
    choices = colnames(tag_data)
  )

  # Placate R-CMD-check --------------------------------------------------------

  date_released <- NULL
  date_recovered <- NULL
  region_released <- NULL
  region_recovered <- NULL
  size_released <- NULL
  count <- NULL

  # Compute dates --------------------------------------------------------------

  # Initial date of the first released step
  date_released_start <- lubridate::as_date(
    paste0(year_released_start, "-01-01")
  )
  # Final date of the last released step
  date_recovered_end <- lubridate::as_date(
    paste0(year_recovered_end, "-12-31")
  )

  # Compute index limits -------------------------------------------------------

  n_times <- compute_n_years(date_released_start, date_recovered_end)
  n_steps <- n_times * count_intervals(step_interval, "year")
  n_sizes <- length(list_sizes)
  n_regions <- length(list_regions)
  n_duration <- ifelse(is.null(step_duration_max), n_steps-1, step_duration_max)

  # Assemble released tibble ---------------------------------------------------

  tags_released_tibble <- tag_data %>%
    dplyr::select(
      date_released = colname_date_released,
      region_released = colname_region_released,
      size_released = colname_size_released
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      n = create_step_released(
        date_released = .data$date_released,
        date_released_start = date_released_start,
        date_recovered_end = date_recovered_end,
        step_interval = step_interval
      ),
      l = create_group(x = .data$size_released, list_x = list_sizes),
      x = create_group(x = .data$region_released, list_x = list_regions)
    ) %>%
    dplyr::filter(.data$n < n_steps) %>% # n_steps = N
    dplyr::select(
      .data$n,
      .data$l,
      .data$x
    ) %>%
    dplyr::group_by(
      .data$n,
      .data$l,
      .data$x
    ) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(
      .data$n,
      .data$l,
      .data$x
    ) %>%
    tidyr::drop_na()

  # Assemble recovered tibble --------------------------------------------------

  tags_recovered_tibble <- tag_data %>%
    dplyr::select(
      date_released = colname_date_released,
      date_recovered = colname_date_recovered,
      region_released = colname_region_released,
      region_recovered = colname_region_recovered,
      size_released = colname_size_released
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      n = create_step_released(
        date_released = .data$date_released,
        date_released_start = date_released_start,
        date_recovered_end = date_recovered_end,
        step_interval = step_interval
      ),
      d = create_step_duration_max(
        date_released = .data$date_released,
        date_recovered = .data$date_recovered,
        date_released_start = date_released_start,
        date_recovered_end = date_recovered_end,
        step_interval = step_interval,
        step_duration_max = step_duration_max
      ),
      l = create_group(x = .data$size_released, list_x = list_sizes),
      x = create_group(x = .data$region_released, list_x = list_regions),
      y = create_group(x = .data$region_recovered, list_x = list_regions)
    ) %>%
    tidyr::drop_na() %>%
    dplyr::filter(.data$n < n_steps) %>% # n_steps = N
    dplyr::filter(.data$d > 1L) %>%
    dplyr::filter(.data$d <= n_duration) %>%
    dplyr::filter(.data$n + .data$d - 2L < n_steps) %>% # One beyond released
    dplyr::select(
      .data$n,
      .data$d,
      .data$l,
      .data$x,
      .data$y
    ) %>%
    dplyr::group_by(
      .data$n,
      .data$d,
      .data$l,
      .data$x,
      .data$y
    ) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(
      .data$n,
      .data$d,
      .data$l,
      .data$x,
      .data$y
    ) %>%
    tidyr::drop_na()

  # Initialize array -----------------------------------------------------------

  tag_array <- array(
    0L,
    dim = c(n_steps - 1L, n_duration, n_sizes, n_regions, n_regions)
  )

  # Populate released ----------------------------------------------------------

  for (i in seq_len(nrow(tags_released_tibble))) {
    tag_array[
      tags_released_tibble$n[i],
      1L,
      tags_released_tibble$l[i],
      tags_released_tibble$x[i],
      tags_released_tibble$x[i]
    ] <- as.integer(tags_released_tibble$count[i])
  }

  # Populate recovered ---------------------------------------------------------

  for (i in seq_len(nrow(tags_recovered_tibble))) {
    tag_array[
      tags_recovered_tibble$n[i],
      tags_recovered_tibble$d[i],
      tags_recovered_tibble$l[i],
      tags_recovered_tibble$x[i],
      tags_recovered_tibble$y[i]
    ] <- as.integer(tags_recovered_tibble$count[i])
  }

  # Return array ---------------------------------------------------------------

  return(tag_array)
}


# create_fishing_weight <- function (tag_array,
#                                    year_start,
#                                    year_end,
#                                    step_interval,
#                                    term_interval) {
#
#   # Check arguments ------------------------------------------------------------
#
#   checkmate::assert_array(
#     tag_array,
#     mode = "integer",
#     any.missing = FALSE,
#     d = 5L
#   )
#
#   # Get dimensions -------------------------------------------------------------
#
#   N <- dim(tag_array)[1L] # Stand-in for [N - 1]
#   D <- dim(tag_array)[2L]
#   L <- dim(tag_array)[3L]
#   X <- dim(tag_array)[4L]
#
#   # Compute constants ----------------------------------------------------------
#
#   terms_per_year <- count_intervals(term_interval, "year")
#
#   # Assemble term index --------------------------------------------------------
#
#   n_to_k <- index_n_to_k(year_start, year_end, step_interval, term_interval)
#
#   # Assemble fishing weight ----------------------------------------------------
#
#   fishing_weight <- array(0, dim = c(terms))
#
#   # Return value ---------------------------------------------------------------
#
#   return(fishing_weight)
# }

#' Read File From Path
#'
#' @param path [character()] file path
#'
#' @return the object at the file path
#' @export
#'
read_from_path <- function (path) {
  envir <- environment()
  data_name <- load(path, envir = envir)
  get(data_name)
}
