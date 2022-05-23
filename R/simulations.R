

# sim <- function (n_steps,
#                  n_liberty,
#                  n_regions,
#                  movement_step,
#                  mean_released = 500,
#                  dispersion_released = Inf, # Poisson
#                  dispersion = Inf, # Poisson
#                  data = NULL) {
#
#   # Check arguments ------------------------------------------------------------
#
#   # Has data?
#   has_data <- ifelse(is.null(data), FALSE, TRUE)
#
#   # Construct movement step ----------------------------------------------------
#
#   # Construct survival step ----------------------------------------------------
#
#   survival_step <- diag(0.85, nrow = n_regions)
#
#   # Construct observed step ----------------------------------------------------
#
#   observed_step <- diag(0.05, nrow = n_regions)
#
#   # Populate abundance ---------------------------------------------------------
#
#   # Declare
#   abundance <- array(
#     data = 0L,
#     dim = c(n_steps - 1L, n_liberty, n_regions, n_regions)
#   )
#   tags <- abundance
#   # Released
#   released <- array(
#     stats::rnbinom(
#       n = (n_steps - 1L) * n_regions,
#       size = dispersion_released, # var = mu + mu^2/dispersion
#       mu = mean_released
#     ),
#     dim = c(n_steps - 1L, n_regions)
#   )
#   # Populate released
#   for (n in seq_len(n_steps - 1)) {
#     for (x in n_regions) {
#       abundance[n, 1L, x, x] <- released[n, x] * (1 - initial_loss_step)
#       tags[n, 1L, x, x] <- released[n, x]
#     }
#   }
#
#   # Perform simulation ---------------------------------------------------------
#
#   for (n in seq_len(n_steps - 1L)) {
#     for (l in 2:min(n_steps - n + 1, n_liberty)) {
#       # Abundance
#       abundance[n, l,,] <- as.matrix(abundance[n, l - 1,,]) %*%
#         survival_step %*%
#         movement_step
#       # Observed
#       tags[n, l,,] <- array(
#         stats::rnbinom(
#           n = n_regions * n_regions,
#           size = dispersion,
#           mu = as.vector(
#             as.matrix(abundance[n, l,,]) %*%
#               observed_step
#           ),
#         ),
#         dim = c(n_regions, n_regions)
#       )
#     }
#   }
#
#   # Assemble data --------------------------------------------------------------
#
#   # Return values --------------------------------------------------------------
#
#   return(sim_data)
# }
