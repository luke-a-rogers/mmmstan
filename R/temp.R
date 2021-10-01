#
# # data
# harvest_rates <- sablefishdata::harvest_rates
# tags_released <- sablefishdata::tags_released
# tags_recovered <- sablefishdata::tags_recovered
#
# # lists
# area_list <- list(ak = 1, bc = 2, cc = 3)
# group_list <- list(small = 400:549, large = 550:800)
#
# # data_list arguments
# n_areas <- length(area_list)
# released_group_max <- length(group_list)
# harvest_group_max <- length(group_list)
# released_time_max <- 152
# liberty_time_max <- 160
# study_time_max <- liberty_time_max
# movement_time_max <- 40
# harvest_time_max <- 40
# reporting_time_max <- 1
# year_time_max <- 4
# natural_mortality_rate <- 0.1
# tag_loss_rate_initial <- 0.1
# tag_loss_rate_ongoing <- 0.02
# reporting_rates <- array(c(.4, .5, .3), c(1, 3))
# z <- mmmstan::movement_index(n = 3, pattern = 0)
# random_walk <- 1
# harvest_prior_mean <- array(0.05, dim = c(harvest_group_max, harvest_time_max, n_areas))
# harvest_prior_sd <- array(0.01, dim = c(harvest_group_max, harvest_time_max, n_areas))
# phi_prior_mean <- 0.6
# phi_prior_sd <- 0.2
# sigma_prior_mean <- rep(0.01, n_areas)
# sigma_prior_sd <- rep(0.005, n_areas)
# movement_rate_fudge <- 1e-12
# predicted_tags_fudge <- 1e-12
#
# # indexes
# movement_time_index <- rep(
#   x = seq_len(ceiling(study_time_max / year_time_max)),
#   each = year_time_max
# )[seq_len(study_time_max)]
#
# harvest_time_index <- rep(
#   x = seq_len(ceiling(study_time_max / year_time_max)),
#   each = year_time_max
# )[seq_len(study_time_max)]
#
# reporting_time_index <- rep(1, study_time_max)
# harvest_group_index <- seq_len(length(group_list))
#
# # tags
# tags <- mmmstan::tag_arrays(
#   released = tags_released,
#   recovered = tags_recovered,
#   released_time_unit = "quarter",
#   released_time_max = released_time_max,
#   liberty_time_max = liberty_time_max,
#   liberty_days_min = 30,
#   colname_released_date = "released_date",
#   colname_released_area = "released_region",
#   colname_group = "released_length",
#   colname_recovered_date = "recovered_date",
#   colname_recovered_area = "recovered_region",
#   colname_id = "tag_id",
#   area_list = area_list,
#   group_list = group_list,
#   released_date_start = "1979-01-01",
#   released_date_end = "2016-12-31"
# )
#
# # data_list
# data_list <- list(
#   A = n_areas,
#   G_released = released_group_max,
#   G_harvest = harvest_group_max,
#   T_released = released_time_max,
#   T_liberty = liberty_time_max,
#   T_study = study_time_max,
#   T_movement = movement_time_max,
#   T_harvest = harvest_time_max,
#   T_reporting = reporting_time_max,
#   T_year = year_time_max,
#   m = natural_mortality_rate,
#   u = (1 - tag_loss_rate_initial),
#   v = tag_loss_rate_ongoing,
#   w = reporting_rates,
#   x = tags$x,
#   y = tags$y,
#   z = z,
#   p_time_index = movement_time_index,
#   h_time_index = harvest_time_index,
#   w_time_index = reporting_time_index,
#   h_group_index = harvest_group_index,
#   random_walk = random_walk,
#   h_prior_mean = harvest_prior_mean,
#   h_prior_sd = harvest_prior_sd,
#   phi_prior_mean = phi_prior_mean,
#   phi_prior_sd = phi_prior_sd,
#   sigma_prior_mean = sigma_prior_mean,
#   sigma_prior_sd = sigma_prior_sd,
#   p_fudge = movement_rate_fudge,
#   y_fudge = predicted_tags_fudge
# )
#
# # fit arguments
# data <- data_list
# chains <- 1
# step_size <- 0.01
# iter_warmup <- 150
# iter_sampling <- 200
# use_reduce_sum <- TRUE
# threads_per_chain <- 16
#
# # fit model
# mmmstan::mmmfit(
#   data = data,
#   chains = chains,
#   step_size = step_size,
#   iter_warmup = iter_warmup,
#   iter_sampling = iter_sampling,
#   use_reduce_sum = use_reduce_sum,
#   threads_per_chain = threads_per_chain
# )
# # Error in sigma_mean[sigma$harvest_group[i], sigma$harvest_time[i], sigma$current_area[i]] :
# #   incorrect number of dimensions
# # In addition: Warning messages:
# #   1: Unknown or uninitialised column: `harvest_group`.
# # 2: Unknown or uninitialised column: `harvest_time`.
#
# # Create model object
# mod <- cmdstanr::cmdstan_model(
#   ifelse(
#     use_reduce_sum,
#     system.file("stan", "mmm_reduce_sum.stan", package = "mmmstan"),
#     system.file("stan", "mmm.stan", package = "mmmstan")
#   ),
#   include_path = system.file("stan", package = "mmmstan"),
#   cpp_options = list(stan_threads = TRUE))
# # Fit the model
# cmdfit <- mod$sample(
#   data = data,
#   chains = chains,
#   step_size = step_size,
#   iter_warmup = iter_warmup,
#   iter_sampling = iter_sampling,
#   threads_per_chain = threads_per_chain
# )
# # Summary
# tibble::view(cmdfit$summary())
# # Assemble parameters
# parameters <- mmmstan::summary_to_tibbles(
#   cmdfit$summary(),
#   data$h_prior_mean,
#   data$h_prior_sd,
#   data$phi_prior_mean,
#   data$phi_prior_sd,
#   data$sigma_prior_mean,
#   data$sigma_prior_sd
# )
# # Error in sigma_mean[sigma$harvest_group[i], sigma$harvest_time[i], sigma$current_area[i]] :
# #   incorrect number of dimensions
# # In addition: Warning messages:
# #   1: Unknown or uninitialised column: `harvest_group`.
# # 2: Unknown or uninitialised column: `harvest_time`.
#
# # summary_to_tibbles arguments
# x <- cmdfit$summary()
# h_mean <- data$h_prior_mean
# h_sd <- data$h_prior_sd
# phi_mean <- data$phi_prior_mean
# phi_sd <- data$phi_prior_sd
# sigma_mean <- data$sigma_prior_mean
# sigma_sd <- data$sigma_prior_sd
#
# # harvest rates
# library(rlang)
# library(magrittr)
# h <- x %>%
#   dplyr::filter(startsWith(.data$variable, "h[")) %>%
#   dplyr::mutate(
#     name = "h",
#     harvest_group = as.numeric(mmmstan::strings_to_indices(.data$variable, 1L)),
#     harvest_time = as.numeric(mmmstan::strings_to_indices(.data$variable, 2L)),
#     current_area = as.numeric(mmmstan::strings_to_indices(.data$variable, 3L)),
#     prior_mean = NA_real_,
#     prior_sd = NA_real_
#   ) %>%
#   dplyr::select(
#     .data$name,
#     .data$variable,
#     .data$harvest_group,
#     .data$harvest_time,
#     .data$current_area,
#     .data$prior_mean,
#     .data$prior_sd,
#     .data$mean:.data$ess_tail
#   )
# for(i in seq_len(nrow(h))) {
#   h$prior_mean[i] <- h_mean[
#     h$harvest_group[i],
#     h$harvest_time[i],
#     h$current_area[i]]
#   h$prior_sd[i] <- h_sd[
#     h$harvest_group[i],
#     h$harvest_time[i],
#     h$current_area[i]]
# }
# # random walk standard deviation
# sigma <- x %>%
#   dplyr::filter(startsWith(.data$variable, "sigma")) %>%
#   dplyr::mutate(
#     name = "sigma",
#     current_area = as.numeric(mmmstan::strings_to_indices(.data$variable, 1L)),
#     prior_mean = NA_real_,
#     prior_sd = NA_real_
#   ) %>%
#   dplyr::select(
#     .data$name,
#     .data$variable,
#     .data$current_area,
#     .data$prior_mean,
#     .data$prior_sd,
#     .data$mean:.data$ess_tail
#   )
# for(i in seq_len(nrow(sigma))) {
#   sigma$prior_mean[i] <- sigma_mean[
#     sigma$harvest_group[i],
#     sigma$harvest_time[i],
#     sigma$current_area[i]]
#   sigma$prior_sd[i] <- sigma_sd[
#     sigma$harvest_group[i],
#     sigma$harvest_time[i],
#     sigma$current_area[i]]
# }
# # Error in sigma_mean[sigma$harvest_group[i], sigma$harvest_time[i], sigma$current_area[i]] :
# #   incorrect number of dimensions
# # In addition: Warning messages:
# #   1: Unknown or uninitialised column: `harvest_group`.
# # 2: Unknown or uninitialised column: `harvest_time`.
#
# for(i in seq_len(nrow(sigma))) {
#   sigma$prior_mean[i] <- sigma_mean[sigma$current_area[i]]
#   sigma$prior_sd[i] <- sigma_sd[sigma$current_area[i]]
# }
