# #' @export
# #' @import methods
# plot.mmmstan <- function (x) {
#
#   # Check arguments ------------------------------------------------------------
#
#   checkmate::assert_class(x, "mmmstan")
#
#   # Define values --------------------------------------------------------------
#
#   movement_rate <- x$summary$movement_rate
#   x_max <- max(movement_rate$x)
#   i_max <- max(movement_rate$i)
#
#   # Assemble plots -------------------------------------------------------------
#
#   par(
#     mfrow = c(x_max, x_max),
#     mar = c(3, 3, 0.5, 0.5)
#   )
#   for (x in seq_len(x_max)) {
#     for (y in seq_len(x_max)) {
#       data <- movement_rate %>%
#         dplyr::filter(.data$x == .env$x, .data$y == .env$y)
#       plot(
#         x = data$i,
#         y = data$mean,
#         pch = 20,
#         ylim = c(0, 1),
#         xlab = "",
#         ylab = "",
#         xaxt = "n",
#         yaxt = "n"
#       )
#       for (i in seq_len(i_max)) {
#         lines(
#           x = rep(data$i[i], 2),
#           y = c(data$q5[i], data$q95[i])
#         )
#       }
#       lines(
#         x = data$i,
#         y = data$mean
#       )
#       title(
#         xlab = ifelse(x == x_max, "Movement term (i)", ""),
#         ylab = ifelse(y == 1, "Movement rate", ""),
#         line = 2
#       )
#       if (x == x_max) {
#         axis(
#           side = 1,
#           padj = -0.5
#         )
#       }
#       if (y == 1) {
#         axis(
#           side = 2,
#           at = seq(0, 1, 0.5),
#           hadj = 0.75,
#           las = 1
#         )
#       }
#     }
#   }
#   par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
# }
