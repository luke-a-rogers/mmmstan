#' Create Harvest Rate Annual
#'
#' @param h_step [numeric()] step harvest rate
#' @param y [integer()] number of time steps in one year
#'
#' @return [numeric()]
#' @export
#'
#' @examples
#' create_h_annual(0.05, 1)
#' create_h_annual(0.01274146, 4)
#'
create_h_annual <- function (h_step, y) {
  return(1 - (1 - h_step)^(y))
}

#' Create Harvest Rate Step
#'
#' @param h_annual [numeric()] annual harvest rate
#' @param y [integer()] number of time steps in one year
#'
#' @return [numeric()]
#' @export
#'
#' @examples
#' create_h_step(0.05, 1)
#' create_h_step(0.05, 4)
#'
create_h_step <- function (h_annual, y) {
  return(1 - (1 - h_annual)^(1/y))
}

#' Movement Index Matrix
#'
#' @param n [integer()] Number of areas.
#' @param pattern [integer()] One of \code{1}: full movement, or \code{0}:
#'   direct movement between numerically sequential neighbours only.
#' @param allow [integer()] [matrix()] Each row indicates directional
#'   movement between a pair of areas indexed from one.
#' @param disallow [integer()] [matrix()] As for \code{allow}, but
#'   specified movement is disallowed.
#'
#' @return A square binary matrix
#' @export
#'
#' @examples
#'
#' # Neighbours (default)
#' movement_index(6)
#'
#' # Full movement
#' movement_index(6, 1)
#'
#' # Neighbours plus 1-6, 2-5, and 5-2, but not 6-1
#' allow <- matrix(c(1,6,2,5,5,2), ncol = 2, byrow = TRUE)
#' movement_index(6, 0, allow)
#'
#' # Full minus 1-6
#' disallow <- matrix(c(1,6), ncol = 2, byrow = TRUE)
#' movement_index(6, 1, disallow = disallow)
#'
movement_index <- function (n, pattern = NULL, allow = NULL, disallow = NULL) {

  # Check arguments
  # checkmate::assert_integerish(n, lower = 2, len = 1, any.missing = FALSE)
  # checkmate::assert_integerish(pattern, len = 1, null.ok = TRUE)
  # checkmate::assert_integerish(pattern, lower = 0, upper = 1, null.ok = TRUE)
  # checkmate::assert_integerish(pattern, any.missing = FALSE, null.ok = TRUE)
  # checkmate::assert_matrix(allow, mode = "integerish", null.ok = TRUE)
  # checkmate::assert_matrix(allow, any.missing = FALSE, null.ok = TRUE)
  # checkmate::assert_matrix(allow, ncols = 2, null.ok = TRUE)
  # Initialize matrix
  z <- matrix(0L, nrow = n, ncol = n)
  diag(z) <- 1L
  # Default
  if (is.null(pattern) & is.null(allow)) {
    # Immediate neighbours only
    for (i in seq_len(n - 1L)) {
      z[i, i + 1L] <- 1L
      z[i + 1L, i] <- 1L
    }
  }
  # Add pattern
  if (!is.null(pattern)) {
    if (pattern) {
      z <- matrix(1L, nrow = n, ncol = n)
    } else {
      for (i in seq_len(n - 1L)) {
        z[i, i + 1L] <- 1L
        z[i + 1L, i] <- 1L
      }
    }
  }
  # Allow indexes
  if (!is.null(allow)) {
    z[allow] <- 1L
  }
  # Disallow indexes
  if (!is.null(disallow)) {
    z[disallow] <- 0L
  }
  # Return
  return(z)
}

#' Create Tag Areas
#'
#' @param x [atomic()] vector of area identifiers
#' @param area_list [list()] of named vectors of area identifiers
#'
#' @return [integer()] vector
#' @export
#'
#' @examples
#'
#' # Numeric area identifiers
#' x <- c(1, 4, 2, 3, NA, 2, 7)
#' a <- list(area_1 = 1:2, area_2 = 3:4)
#' create_tag_areas(x, a)
#'
#' # Character area identifiers
#' x <- c("A", "B", "B", NA, "A", "C", "B")
#' a <- list(area_1 = "A", area_2 = "B")
#' create_tag_areas(x, a)
#'
#' # Factor
#' x <- factor(c("C", "A", "B", NA, "B", "D"))
#' a <- list(area_1 = "A", area_2 = "B", area_3 = "C")
#' create_tag_areas(x, a)
#'
create_tag_areas <- function (x, area_list) {
  # Check arguments

  # Create areas
  s <- seq_along(area_list)
  x <- data.frame(val = x)
  y <- data.frame(val = unlist(area_list), ind = rep(s, lengths(area_list)))
  # Return group vector
  return(as.integer(dplyr::left_join(x, y, by = "val")$ind))
}

#' Create Tag Groups
#'
#' @param x [atomic()] vector of values to be grouped
#' @param group_list [list()] of named vectors of group elements
#'
#' @return [integer()] vector
#' @export
#'
#' @examples
#'
#' # Numeric
#' x <- c(1, 4, 2, 3, NA, 2, 7)
#' g <- list(a = 1:2, b = 3:4)
#' create_tag_groups(x, g)
#'
#' # Character
#' x <- c("M", "F", "F", NA, "M", "N", "F")
#' g <- list(m = "M", f = "F")
#' create_tag_groups(x, g)
#'
#' # Factor
#' x <- factor(c("L", "S", "M", NA, "M", "H"))
#' g <- list(s = "S", m = "M", l = "L")
#' create_tag_groups(x, g)
#'
create_tag_groups <- function (x, group_list) {
  # Check arguments

  # Create groups
  s <- seq_along(group_list)
  x <- data.frame(val = x)
  y <- data.frame(val = unlist(group_list), ind = rep(s, lengths(group_list)))
  # Return group vector
  return(as.integer(dplyr::left_join(x, y, by = "val")$ind))
}

#' Create Tag Liberty Steps
#'
#' @param release_date [character()] vector of release dates
#' @param recovery_date [character()] vector of recovery dates
#' @param release_start [character()] release date start
#' @param release_end [character()] release date end
#' @param time_step [character()] one of \code{c("year", "quarter", "month")}
#' @param max_steps_liberty [integer()] maximum number of time steps at liberty
#'
#' @return [integer()] vector
#' @export
#'
#' @examples
#'
#' # Yearly
#' rel_date <- rep("2011-01-01", 5)
#' rec_date <- c("2011-01-01", "2012-01-01", NA, "2020-01-01", "2021-01-01")
#' rel_start <- c("2011-01-01")
#' rel_end <- c("2020-01-01")
#' time_step <- "year"
#' max_steps_liberty <- 10
#' create_tag_liberty_steps(rel_date,
#'                          rec_date,
#'                          rel_start,
#'                          rel_end,
#'                          time_step,
#'                          max_steps_liberty)
#'
#' # Quarterly
#' rel_date <- rep("2011-01-01", 5)
#' rec_date <- c("2011-01-01", "2012-12-31", NA, "2020-12-31", "2021-01-01")
#' rel_start <- c("2011-01-01")
#' rel_end <- c("2020-01-01")
#' time_step <- "quarter"
#' max_steps_liberty <- 40
#' create_tag_liberty_steps(rel_date,
#'                          rec_date,
#'                          rel_start,
#'                          rel_end,
#'                          time_step,
#'                          max_steps_liberty)
#'
#' # Monthly
#' rel_date <- rep("2011-01-01", 5)
#' rec_date <- c("2011-01-01", "2011-02-01", NA, "2020-12-31", "2021-01-01")
#' rel_start <- c("2011-01-01")
#' rel_end <- c("2020-01-01")
#' time_step <- "month"
#' max_steps_liberty <- 120
#' create_tag_liberty_steps(rel_date,
#'                          rec_date,
#'                          rel_start,
#'                          rel_end,
#'                          time_step,
#'                          max_steps_liberty)
#'
create_tag_liberty_steps <- function (release_date,
                                      recovery_date,
                                      release_start,
                                      release_end,
                                      time_step,
                                      max_steps_liberty) {
  # Check arguments

  # Convert to dates
  release_date <- lubridate::as_date(release_date)
  recovery_date <- lubridate::as_date(recovery_date)
  release_start <- lubridate::as_date(release_start)
  release_end <- lubridate::as_date(release_end)
  # Replace disallowed release dates by NAs
  release_date[which(release_date < release_start)] <- NA
  release_date[which(release_date > release_end)] <- NA
  # Define year and months
  release_year <- lubridate::year(release_date)
  release_quart <- lubridate::quarter(release_date)
  release_month <- lubridate::month(release_date)
  recovery_year <- lubridate::year(recovery_date)
  recovery_quart <- lubridate::quarter(recovery_date)
  recovery_month <- lubridate::month(recovery_date)
  # Create liberty steps
  if (time_step == "year") {
    liberty_steps <- recovery_year - release_year + 1
  } else if (time_step == "quarter") {
    liberty_years <- recovery_year - release_year
    liberty_quarters <- recovery_quart - release_quart
    liberty_steps <- 4 * liberty_years + liberty_quarters + 1
  } else if (time_step == "month") {
    liberty_years <- recovery_year - release_year
    liberty_months <- recovery_month - release_month
    liberty_steps <- 12 * liberty_years + liberty_months + 1
  } else {
    stop("time_step must be 'year' or 'month'")
  }
  # Replace disallowed liberty steps by NAs
  liberty_steps[which(liberty_steps < 2)] <- NA
  liberty_steps[which(liberty_steps > max_steps_liberty)] <- NA
  # Return release step
  return(as.integer(liberty_steps))
}

#' Create Tag Release Step
#'
#' @param release_date [character()] vector of release dates
#' @param release_start [character()] release date start
#' @param release_end [character()] release date end
#' @param time_step [character()] one of \code{c("year", "quarter", "month")}
#'
#' @return [integer()] vector of tag release time steps
#' @export
#'
#' @examples
#'
#' x <- c("2010-01-01", "2011-01-01", NA, "2020-12-31", "2021-01-01")
#' d0 <- c("2011-01-01")
#' d1 <- c("2020-12-31")
#' create_tag_release_steps(x, d0, d1, "year")
#' create_tag_release_steps(x, d0, d1, "quarter")
#' create_tag_release_steps(x, d0, d1, "month")
#'
create_tag_release_steps <- function (release_date,
                                      release_start,
                                      release_end,
                                      time_step) {
  # Check arguments

  # Convert to dates
  release_date <- lubridate::as_date(release_date)
  release_start <- lubridate::as_date(release_start)
  release_end <- lubridate::as_date(release_end)
  # Replace disallowed dates by NAs
  release_date[which(release_date < release_start)] <- NA
  release_date[which(release_date > release_end)] <- NA
  # Define year and months
  release_year <- lubridate::year(release_date)
  release_month <- lubridate::month(release_date)
  release_quart <- lubridate::quarter(release_date)
  start_year <- lubridate::year(release_start)
  start_quart <- lubridate::quarter(release_start)
  start_month <- lubridate::month(release_start)
  # Create time step
  if (time_step == "year") {
    s <- release_year - start_year + 1
  } else if (time_step == "quarter") {
    s <- 4 * (release_year - start_year) + (release_quart - start_quart) + 1
  } else if (time_step == "month") {
    s <- 12 * (release_year - start_year) + (release_month - start_month) + 1
  } else {
    stop("time_step must be 'year', 'quarter', or 'month'")
  }
  # Return release step
  return(as.integer(s))
}

#' Extract Indices From A Character String
#'
#' @param x [character()]
#' @param index [numeric()]
#'
#' @return [numeric()] [vector()]
#' @export
#'
#' @examples
#' x <- c("h[1,1,1]", "h[1,1,2]")
#' strings_to_indices(x)
#' strings_to_indices(x, index = 3)
#'
#' x <- list("h[1,1,1]", "h[1,1,2]")
#' strings_to_indices(x)
#' strings_to_indices(x, index = 3)
#'
strings_to_indices <- function (x, index = NULL) {

  # Check arguments
  # index must be NULL or numeric length 1

  if (is.null(index)) {
    l <- list()
    for (i in seq_along(x)) {
      l[[i]] <- as.numeric(
        strsplit(
          stringr::str_extract_all(
            x[[i]],
            "(?<=\\[).+?(?=\\])"
          )[[1]],
          ","
        )[[1]]
      )
    }
  } else {
    l <- list()
    for (i in seq_along(x)) {
      l[[i]] <- as.numeric(
        strsplit(
          stringr::str_extract_all(
            x[[i]],
            "(?<=\\[).+?(?=\\])"
          )[[1]],
          ","
        )[[1]]
      )[[index]]
    }
  }
  return(l)
}

#' Convert Movement Model Summary To List Of Tibbles
#'
#' @param x [cmdstanr::summary()] from [mmmfit()] object
#' @param h_mean [numeric()] [array()] of harvest prior mean values
#' @param h_sd [numeric()] [array()] of harvest prior standard deviations
#' @param phi_mean [numeric()] dispersion prior mean value
#' @param phi_sd [numeric()] dispersion prior standard deviation
#' @param sigma_mean [numeric()] [vector()] of random walk standard deviation
#'   prior mean values
#' @param sigma_sd [numeric()] [vector()] of random walk standard deviation
#'   prior standard deviations
#'
#' @return [list()] of [tibble::tibble()]
#' @export
#'
summary_to_tibbles <- function (x,
                                h_mean,
                                h_sd,
                                phi_mean,
                                phi_sd,
                                sigma_mean,
                                sigma_sd) {

  # Movement probabilities
  p <- x %>%
    dplyr::filter(startsWith(.data$variable, "p[")) %>%
    dplyr::mutate(
      name = "p",
      previous_area = as.numeric(strings_to_indices(.data$variable, 1L)),
      current_area  = as.numeric(strings_to_indices(.data$variable, 2L)),
      movement_time  = as.numeric(strings_to_indices(.data$variable, 3L)),
      released_group = as.numeric(strings_to_indices(.data$variable, 4L)),
      prior_lower = 0,
      prior_upper = 1
    ) %>%
    dplyr::select(
      .data$name,
      .data$variable,
      .data$previous_area,
      .data$current_area,
      .data$movement_time,
      .data$released_group,
      .data$prior_lower,
      .data$prior_upper,
      .data$mean:.data$ess_tail
    )
  # Harvest rates
  h <- x %>%
    dplyr::filter(startsWith(.data$variable, "h[")) %>%
    dplyr::mutate(
      name = "h",
      harvest_group = as.numeric(strings_to_indices(.data$variable, 1L)),
      harvest_time = as.numeric(strings_to_indices(.data$variable, 2L)),
      current_area = as.numeric(strings_to_indices(.data$variable, 3L)),
      prior_mean = NA_real_,
      prior_sd = NA_real_
    ) %>%
    dplyr::select(
      .data$name,
      .data$variable,
      .data$harvest_group,
      .data$harvest_time,
      .data$current_area,
      .data$prior_mean,
      .data$prior_sd,
      .data$mean:.data$ess_tail
    )
  for(i in seq_len(nrow(h))) {
    h$prior_mean[i] <- h_mean[
      h$harvest_group[i],
      h$harvest_time[i],
      h$current_area[i]]
    h$prior_sd[i] <- h_sd[
      h$harvest_group[i],
      h$harvest_time[i],
      h$current_area[i]]
  }

  # Negative binomial dispersion parameter
  phi <- x %>%
    dplyr::filter(startsWith(.data$variable, "phi")) %>%
    dplyr::mutate(
      name = "phi",
      prior_mean = phi_mean,
      prior_sd = phi_sd
    ) %>%
    dplyr::select(
      .data$name,
      .data$variable,
      .data$prior_mean,
      .data$prior_sd,
      .data$mean:.data$ess_tail
    )
  # Random walk standard deviation
  sigma <- x %>%
    dplyr::filter(startsWith(.data$variable, "sigma")) %>%
    dplyr::mutate(
      name = "sigma",
      current_area = as.numeric(strings_to_indices(.data$variable, 1L)),
      prior_mean = NA_real_,
      prior_sd = NA_real_
    ) %>%
    dplyr::select(
      .data$name,
      .data$variable,
      .data$current_area,
      .data$prior_mean,
      .data$prior_sd,
      .data$mean:.data$ess_tail
    )
  for(i in seq_len(nrow(sigma))) {
    sigma$prior_mean[i] <- sigma_mean[sigma$current_area[i]]
    sigma$prior_sd[i] <- sigma_sd[sigma$current_area[i]]
  }

  # Return list
  return(list(
    p = p,
    h = h,
    phi = phi,
    sigma = sigma)
  )
}
