#' Create Tag Array
#'
#' @param x [data.frame()] tag release observation data. See details
#' @param y [data.frame()] tag recovery observation data. See details
#' @param time_step [character()] desired time step for tag array. One of
#'    \code{c("year", "quarter", "month")}
#' @param max_steps_liberty [integer()] maximum time steps at liberty
#' @param release_date_colname [character()] column name of the release date
#'   in \code{x} and \code{y}
#' @param release_area_colname [character()] column name of the release area
#'   in \code{x} and \code{y}
#' @param group_colname [character()] column name of the group indicator
#'   in \code{x} and \code{y}
#' @param recovery_date_colname [character()] column name of the recovery date
#'   in \code{y}
#' @param recovery_area_colname [character()] column name of the recovery area
#'   in \code{y}
#' @param id_colname [character()] column name of the tag identifier
#'   in \code{x} and \code{y}
#' @param area_list [list()] named atomic vectors. See details
#' @param group_list [list()] named atomic vectors. See details
#' @param release_date_start [character()] release start date
#' @param release_date_end [character()] release end date
#'
#' @details TBD
#'
#' @importFrom magrittr `%>%`
#'
#' @return [array()] aggregated tag observation data.
#' @export
#'
create_tag_array <- function (x,
                              y,
                              time_step,
                              max_steps_liberty,
                              release_date_colname,
                              release_area_colname,
                              group_colname,
                              recovery_date_colname,
                              recovery_area_colname,
                              id_colname,
                              area_list,
                              group_list,
                              release_date_start,
                              release_date_end) {

  # Check arguments ------------------------------------------------------------


  # Rename release columns -----------------------------------------------------

  colnames(x)[which(colnames(x) == release_date_colname)] <- "rel_date"
  colnames(x)[which(colnames(x) == release_area_colname)] <- "rel_area"
  colnames(x)[which(colnames(x) == group_colname)] <- "raw_group"
  colnames(x)[which(colnames(x) == id_colname)] <- "id"

  # Rename recovery columns ----------------------------------------------------

  colnames(y)[which(colnames(y) == release_date_colname)] <- "rel_date"
  colnames(y)[which(colnames(y) == release_area_colname)] <- "rel_area"
  colnames(y)[which(colnames(y) == group_colname)] <- "raw_group"
  colnames(y)[which(colnames(y) == recovery_date_colname)] <- "rec_date"
  colnames(y)[which(colnames(y) == recovery_area_colname)] <- "rec_area"
  colnames(y)[which(colnames(y) == id_colname)] <- "id"

  # Filter disallowed release values -------------------------------------------
  # TODO: filter disallowed values
  x0 <- x

  # Filter disallowed recovery values ------------------------------------------
  # TODO: filter disallowed values

  y0 <- y

  # Aggregate release data -----------------------------------------------------

  x1 <- x0 %>%
    dplyr::mutate(
      mt = create_tag_release_steps(
        rel_date,
        release_date_start,
        release_date_end,
        time_step),
      ma = create_tag_areas(rel_area, area_list),
      mg = create_tag_groups(raw_group, group_list)) %>%
    dplyr::select(mt, ma, mg) %>%
    dplyr::group_by(mt, ma, mg) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(mt, ma, mg) %>%
    tidyr::drop_na()

  # Aggregate recovery data ----------------------------------------------------

  y1 <- y0 %>%
    dplyr::mutate(
      mt = create_tag_release_steps(
        rel_date,
        release_date_start,
        release_date_end,
        time_step),
      ma = create_tag_areas(rel_area, area_list),
      mg = create_tag_groups(raw_group, group_list),
      cl = create_tag_liberty_steps(
        rel_date,
        rec_date,
        release_date_start,
        release_date_end,
        time_step,
        max_steps_liberty),
      ca = create_tag_areas(rec_area, area_list)) %>%
    dplyr::select(mt, ma, mg, cl, ca) %>%
    dplyr::group_by(mt, ma, mg, cl, ca) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(mt, ma, mg, cl, ca) %>%
    tidyr::drop_na()

  # Compute array dimensions ---------------------------------------------------

  nt <- max(x1$mt, y1$mt, na.rm = TRUE)
  na <- length(area_list)
  ng <- length(group_list)
  nl <- max_steps_liberty

  # Populate array -------------------------------------------------------------

  # Initialize
  a <- array(0, dim = c(nt, na, ng, nl, na))

  # Populate from y1
  for (i in seq_len(nrow(y1))) {
    a[y1$mt[i], y1$ma[i], y1$mg[i], y1$cl[i], y1$ca[i]] <- y1$count[i]
  }

  # Populate from x1
  for (i in seq_len(nrow(x1))) {
    a[x1$mt[i], x1$ma[i], x1$mg[i], 1, x1$ma[i]] <- x1$count[i]
  }

  # Return tag array -----------------------------------------------------------

  return(a)
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

#' Create Beta Parameters
#'
#' @param mu [numeric()] Mean
#' @param sd [numeric()] Standard deviation
#' @param n [integer()] Number of samples for histogram
#' @param b [numeric()] Breaks passed to [hist()]
#'
#' @return [numeric()] vector of length two.
#' @export
#'
#' @examples
#' create_beta_parameters(0.01, 0.001)
#'
create_beta_parameters <- function (mu, sd, n = 10000, b = 100) {
  # Compute var
  var <- sd * sd
  # Check parameter condition
  stopifnot(var < mu * (1 - mu))
  # Compute parameters
  nu <- (mu * (1 - mu) / var) - 1
  alpha <- mu * nu
  beta <- (1 - mu) * nu
  # Plot distribution
  x <- stats::rbeta(n, alpha, beta)
  graphics::hist(x, breaks = b)
  # Return vector
  return(c(alpha = alpha, beta = beta))
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

