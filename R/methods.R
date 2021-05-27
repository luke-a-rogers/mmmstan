#' Plot Posterior Summaries
#'
#' @param x [mmmfit()] model fit
#' @param pars [character()] [vector()]
#' @param ci [numeric()] [vector()]
#' @param area_names [character()] [vector()]
#' @param group_names [character()] [vector()]
#' @param ... not currently implemented
#'
#' @return [ggplot2::ggplot()] [ggpubr::ggarrange()]
#' @export
#'
plot.mmmfit <- function (x = NULL,
                         pars = c("p", "h", "phi"),
                         ci = c(0.8, 0.95),
                         area_names = NULL,
                         group_names = NULL,
                         ...) {

  # Check arguments ------------------------------------------------------------


  # Create sample summary ------------------------------------------------------

  tbl <- create_sample_tibble(x, pars = pars)
  s <- create_sample_summary(tbl, ci = ci)

  # Create p plot --------------------------------------------------------------

  if (is.element("p", pars)) {
    # Filter sample summary
    s_p <- dplyr::filter(s, .data$parameter == "p")
    # Plot parameter summary
    p_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step, group = .data$step),
        dplyr::filter(s_p, .data$type == "outer")) +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step, group = .data$step),
        dplyr::filter(s_p, .data$type == "inner"),
        color = "red",
        size = 1.5) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$value, y = .data$step),
        dplyr::filter(s_p, .data$type == "mean")) +
      ggplot2::facet_grid(
        cols = ggplot2::vars(.data$area),
        rows = ggplot2::vars(.data$group))
  } else {
    p_plot <- NULL
  }

  # Create h plot --------------------------------------------------------------

  if (is.element("h", pars)) {
    # Filter sample summary
    s_h <- dplyr::filter(s, .data$parameter == "h")
    # Plot parameter summary
    h_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step, group = .data$step),
        dplyr::filter(s_h, .data$type == "outer")) +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step, group = .data$step),
        dplyr::filter(s_h, .data$type == "inner"),
        color = "red",
        size = 1.5) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$value, y = .data$step),
        dplyr::filter(s_h, .data$type == "mean")) +
      ggplot2::facet_grid(
        cols = ggplot2::vars(.data$area),
        rows = ggplot2::vars(.data$group))
  } else {
    h_plot <- NULL
  }

  # Create phi plot ------------------------------------------------------------

  if (is.element("phi", pars)) {
    # Filter sample summary
    s_phi <- dplyr::filter(s, .data$parameter == "phi")
    # Plot parameter summary
    phi_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step),
        dplyr::filter(s_phi, .data$type == "outer")) +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$value, y = .data$step),
        dplyr::filter(s_phi, .data$type == "inner"),
        color = "red",
        size = 1.5) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$value, y = .data$step),
        dplyr::filter(s_phi, .data$type == "mean"))
  } else {
    phi_plot <- NULL
  }

  # Create composite plot ------------------------------------------------------

  mmm_plot <- ggpubr::ggarrange(
    p_plot, h_plot, phi_plot,
    labels = c("p", "h", "phi"),
    ncol = 1)

  # Return composite plot ------------------------------------------------------

  return(mmm_plot)
}
