#' Compute The Number of Terms Per Year
#'
#' @param x [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#'
#' @return [integer()]
#' @export
#'
#' @examples
#'
#' compute_n_terms("year")
#' compute_n_terms("quarter")
#' compute_n_terms("month")
#'
compute_n_terms <- function (x) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_choice(x = x, choices = c("year", "quarter", "month"))

  # Return the number of terms -------------------------------------------------

  ifelse(x == "year", 1L, ifelse(x == "quarter", 4L, 12L))
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
#' compute_n_times("2011-01-01", "2014-12-31")
#'
compute_n_times <- function (a, b) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_date(lubridate::date(a))
  checkmate::assert_date(lubridate::date(b))

  # Return the number of years -------------------------------------------------

  return(as.integer(abs(lubridate::year(b) - lubridate::year(a)) + 1))
}

#' Create A Movement Index Matrix
#'
#' @param n_regions [integer()] Number of regions
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
create_movement_index <- function (n_regions,
                                   movement_pattern = 2,
                                   movement_allow = NULL,
                                   movement_disallow = NULL) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(
    n_regions,
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
    upper = n_regions,
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
    upper = n_regions,
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

  movement_index <- diag(1L, nrow = n_regions, ncol = n_regions, names = FALSE)

  # Add movement pattern -------------------------------------------------------

  # All pairwise
  if (movement_pattern == 1) {
    movement_index <- matrix(1L, nrow = n_regions, ncol = n_regions)
  }
  # Neighbors only
  if (movement_pattern == 2) {
    for (i in seq_len(n_regions - 1L)) {
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

#' Create Step Liberty
#'
#' @param date_released [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_recovered [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_released_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param term_released [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#' @param step_liberty_max [integer()] maximum steps at liberty
#'
#' @return [integer()] [vector()] step liberty
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
#' create_step_liberty(u, v, w, x, y, z)
#'
#' # Quarterly
#' u <- rep("2011-01-01", 5)
#' v <- c("2011-01-01", "2012-12-31", NA, "2020-12-31", "2021-01-01")
#' w <- c("2011-01-01")
#' x <- c("2020-01-01")
#' y <- "quarter"
#' z <- 40
#' create_step_liberty(u, v, w, x, y, z)
#'
#' # Monthly
#' u <- rep("2011-01-01", 5)
#' v <- c("2011-01-01", "2011-02-01", NA, "2020-12-31", "2021-01-01")
#' w <- c("2011-01-01")
#' x <- c("2020-01-01")
#' y <- "month"
#' z <- 120
#' create_step_liberty(u, v, w, x, y, z)
#'
create_step_liberty <- function (date_released,
                                 date_recovered,
                                 date_released_start,
                                 date_released_end,
                                 term_released,
                                 step_liberty_max) {

  # Check arguments ------------------------------------------------------------

  # Convert to dates -----------------------------------------------------------

  date_released <- lubridate::as_date(date_released)
  date_recovered <- lubridate::as_date(date_recovered)
  date_released_start <- lubridate::as_date(date_released_start)
  date_released_end <- lubridate::as_date(date_released_end)

  # Replace disallowed dates by NAs --------------------------------------------

  date_released[which(date_released < date_released_start)] <- NA_real_
  date_released[which(date_released > date_released_end)] <- NA_real_

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

  if (term_released == "year") {
    step_liberty <- year_recovered - year_released + 1
  } else if (term_released == "quarter") {
    step_liberty <- 4 * (year_recovered - year_released) +
      (quarter_recovered - quarter_released) + 1
  } else if (term_released == "month") {
    step_liberty <- 12 * (year_recovered - year_released) +
      (month_recovered - month_released) + 1
  } else {
    stop("term_released must be 'year', 'quarter', or 'month'")
  }

  # Replace disallowed dates by NAs --------------------------------------------

  step_liberty[which(step_liberty < 2)] <- NA_real_
  step_liberty[which(step_liberty > step_liberty_max)] <- NA_real_

  # Return step liberty --------------------------------------------------------

  return(as.integer(step_liberty))
}

#' Create Step Released
#'
#' @param date_released [character()] [vector()] dates as \code{"\%Y-\%m-\%d"}
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_released_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param term_released [character()] one of \code{"year"}, \code{"quarter"}
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
                                  date_released_end,
                                  term_released) {

  # Check arguments ------------------------------------------------------------

  # Convert to dates -----------------------------------------------------------

  date_released <- lubridate::as_date(date_released)
  date_released_start <- lubridate::as_date(date_released_start)
  date_released_end <- lubridate::as_date(date_released_end)

  # Replace disallowed dates by NAs --------------------------------------------

  date_released[which(date_released < date_released_start)] <- NA_real_
  date_released[which(date_released > date_released_end)] <- NA_real_

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

  if (term_released == "year") {
    step_released <- year_released - year_start + 1
  } else if (term_released == "quarter") {
    step_released <- 4 * (year_released - year_start) +
      (quarter_released - quarter_start) + 1
  } else if (term_released == "month") {
    step_released <- 12 * (year_released - year_start) +
      (month_released - month_start) + 1
  } else {
    stop("term_released must be 'year', 'quarter', or 'month'")
  }

  # Return step released -------------------------------------------------------

  return(as.integer(step_released))
}

#' Create Tags Recovered
#'
#' @param tags [data.frame()]
#' @param list_regions [list()]
#' @param list_sizes [list()]
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_released_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param step_liberty_max [integer()]
#' @param term_released [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#' @param colname_date_released [character()]
#' @param colname_date_recovered [character()]
#' @param colname_region_released [character()]
#' @param colname_region_recovered [character()]
#' @param colname_size_released [character()]
#'
#' @importFrom rlang .data
#'
#' @return [tibble::tibble()]
#' @export
#'
create_tags_recovered <- function (tags,
                                   list_regions,
                                   list_sizes,
                                   date_released_start,
                                   date_released_end,
                                   step_liberty_max,
                                   term_released = "quarter",
                                   colname_date_released = "date_released",
                                   colname_date_recovered = "date_recovered",
                                   colname_region_released = "region_released",
                                   colname_region_recovered = "region_recovered",
                                   colname_size_released = "size_released") {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_data_frame(tags)
  checkmate::assert_list(list_sizes)
  checkmate::assert_list(list_regions)
  checkmate::assert_date(lubridate::date(date_released_start))
  checkmate::assert_date(lubridate::date(date_released_end))
  checkmate::assert_integerish(
    step_liberty_max,
    lower = 2,
    len = 1,
    any.missing = FALSE
  )
  checkmate::assert_choice(
    x = term_released,
    choices = c("year", "quarter", "month")
  )
  checkmate::assert_choice(
    x = colname_date_released,
    choices = colnames(tags)
  )
  checkmate::assert_choice(
    x = colname_date_recovered,
    choices = colnames(tags)
  )
  checkmate::assert_choice(
    x = colname_region_released,
    choices = colnames(tags)
  )
  checkmate::assert_choice(
    x = colname_region_recovered,
    choices = colnames(tags)
  )
  checkmate::assert_choice(
    x = colname_size_released,
    choices = colnames(tags)
  )

  # Placate R-CMD-check --------------------------------------------------------

  date_released <- NULL
  date_recovered <- NULL
  region_released <- NULL
  region_recovered <- NULL
  size_released <- NULL

  # Compute index limits -------------------------------------------------------

  n_times <- compute_n_times(date_released_start, date_released_end)
  n_terms <- compute_n_terms(term_released)
  n_steps <- n_times * n_terms - 1
  n_sizes <- length(list_sizes)
  n_regions <- length(list_regions)
  n_liberty <- step_liberty_max

  # Assemble tibble ------------------------------------------------------------

  tags_recovered_tibble <- tags %>%
    dplyr::rename(
      date_released = .data[[colname_date_released]],
      date_recovered = .data[[colname_date_recovered]],
      region_released = .data[[colname_region_released]],
      region_recovered = .data[[colname_region_recovered]],
      size_released = .data[[colname_size_released]]
    ) %>%
    dplyr::select(
      date_released,
      date_recovered,
      region_released,
      region_recovered,
      size_released
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      n = create_step_released(
        date_released = .data$date_released,
        date_released_start = date_released_start,
        date_released_end = date_released_end,
        term_released = term_released
      ),
      s = create_group(x = .data$size_released, list_x = list_sizes),
      x = create_group(x = .data$region_released, list_x = list_regions),
      l = create_step_liberty(
        date_released = .data$date_released,
        date_recovered = .data$date_recovered,
        date_released_start = date_released_start,
        date_released_end = date_released_end,
        term_released = term_released,
        step_liberty_max = step_liberty_max
      ),
      y = create_group(x = .data$region_recovered, list_x = list_regions)
    ) %>%
    dplyr::filter(.data$n + .data$l - 1L <= n_times * n_terms) %>%
    dplyr::select(
      .data$n,
      .data$s,
      .data$x,
      .data$l,
      .data$y
    ) %>%
    dplyr::group_by(
      .data$n,
      .data$s,
      .data$x,
      .data$l,
      .data$y
    ) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(
      .data$n,
      .data$s,
      .data$x,
      .data$l,
      .data$y
    ) %>%
    tidyr::drop_na()

  # Initialize array -----------------------------------------------------------

  tags_recovered <- array(
    0L,
    dim = c(n_steps, n_sizes, n_regions, n_liberty, n_regions)
  )

  # Populate array -------------------------------------------------------------

  for (i in seq_len(nrow(tags_recovered_tibble))) {
    tags_recovered[
      tags_recovered_tibble$n[i],
      tags_recovered_tibble$s[i],
      tags_recovered_tibble$x[i],
      tags_recovered_tibble$l[i],
      tags_recovered_tibble$y[i]
    ] <- as.integer(tags_recovered_tibble$count[i])
  }

  # Return array ---------------------------------------------------------------

  return(tags_recovered)
}

#' Create Tags Released
#'
#' @param tags [data.frame()]
#' @param list_regions [list()]
#' @param list_sizes [list()]
#' @param date_released_start [character()] date as \code{"\%Y-\%m-\%d"}
#' @param date_released_end [character()] date as \code{"\%Y-\%m-\%d"}
#' @param term_released [character()] one of \code{"year"}, \code{"quarter"}
#'   or \code{"month"}
#' @param colname_date_released [character()]
#' @param colname_region_released [character()]
#' @param colname_size_released [character()]
#'
#' @importFrom rlang .data
#'
#' @return [tibble::tibble()]
#' @export
#'
create_tags_released <- function (tags,
                                  list_regions,
                                  list_sizes,
                                  date_released_start,
                                  date_released_end,
                                  # step_liberty_max,
                                  term_released = "quarter",
                                  colname_date_released = "date_released",
                                  # colname_date_recovered = "date_recovered",
                                  colname_region_released = "region_released",
                                  # colname_region_recovered = "region_recovered",
                                  colname_size_released = "size_released") {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_data_frame(tags)
  checkmate::assert_list(list_sizes)
  checkmate::assert_list(list_regions)
  checkmate::assert_date(lubridate::date(date_released_start))
  checkmate::assert_date(lubridate::date(date_released_end))
  # checkmate::assert_integerish(
  #   step_liberty_max,
  #   lower = 2,
  #   len = 1,
  #   any.missing = FALSE
  # )
  checkmate::assert_choice(
    x = term_released,
    choices = c("year", "quarter", "month")
  )
  checkmate::assert_choice(
    x = colname_date_released,
    choices = colnames(tags)
  )
  # checkmate::assert_choice(
  #   x = colname_date_recovered,
  #   choices = colnames(tags)
  # )
  checkmate::assert_choice(
    x = colname_region_released,
    choices = colnames(tags)
  )
  # checkmate::assert_choice(
  #   x = colname_region_recovered,
  #   choices = colnames(tags)
  # )
  checkmate::assert_choice(
    x = colname_size_released,
    choices = colnames(tags)
  )

  # Placate R-CMD-check --------------------------------------------------------

  date_released <- NULL
  # date_recovered <- NULL
  region_released <- NULL
  # region_recovered <- NULL
  size_released <- NULL

  # Compute index limits -------------------------------------------------------

  n_times <- compute_n_times(date_released_start, date_released_end)
  n_terms <- compute_n_terms(term_released)
  n_steps <- n_times * n_terms - 1
  n_sizes <- length(list_sizes)
  n_regions <- length(list_regions)
  # n_liberty <- step_liberty_max

  # Assemble tibble ------------------------------------------------------------

  tags_released_tibble <- tags %>%
    dplyr::rename(
      date_released = .data[[colname_date_released]],
      # date_recovered = .data[[colname_date_recovered]],
      region_released = .data[[colname_region_released]],
      # region_recovered = .data[[colname_region_recovered]],
      size_released = .data[[colname_size_released]]
    ) %>%
    dplyr::select(
      date_released,
      # date_recovered,
      region_released,
      # region_recovered,
      size_released
    ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(
      n = create_step_released(
        date_released = .data$date_released,
        date_released_start = date_released_start,
        date_released_end = date_released_end,
        term_released = term_released
      ),
      s = create_group(x = .data$size_released, list_x = list_sizes),
      x = create_group(x = .data$region_released, list_x = list_regions)
      # l = create_step_liberty(
      #   date_released = .data$date_released,
      #   date_recovered = .data$date_recovered,
      #   date_released_start = date_released_start,
      #   date_released_end = date_released_end,
      #   term_released = term_released,
      #   step_liberty_max = step_liberty_max
      # ),
      # y = create_group(x = .data$region_recovered, list_x = list_regions)
    ) %>%
    # dplyr::filter(.data$n + .data$l - 1L <= n_times * n_terms) %>%
    dplyr::select(
      .data$n,
      .data$s,
      .data$x
      # .data$l,
      # .data$y
    ) %>%
    dplyr::group_by(
      .data$n,
      .data$s,
      .data$x
      # .data$l,
      # .data$y
    ) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(
      .data$n,
      .data$s,
      .data$x
      # .data$l,
      # .data$y
    ) %>%
    tidyr::drop_na()

  # Initialize array -----------------------------------------------------------

  tags_released <- array(
    0L,
    dim = c(n_steps, n_sizes, n_regions)
  )

  # Populate array -------------------------------------------------------------

  for (i in seq_len(nrow(tags_released_tibble))) {
    tags_released[
      tags_released_tibble$n[i],
      tags_released_tibble$s[i],
      tags_released_tibble$x[i]
    ] <- as.integer(tags_released_tibble$count[i])
  }

  # Return array ---------------------------------------------------------------

  return(tags_released)
}

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
