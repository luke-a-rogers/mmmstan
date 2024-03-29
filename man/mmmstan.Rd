% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mmmstan}
\alias{mmmstan}
\title{Fit A Hidden Markov Movement Model via CmdStan and CmdStanR}
\usage{
mmmstan(
  tag_data,
  list_regions,
  list_sizes,
  year_start,
  year_end,
  step_interval = "quarter",
  step_duration_max = NULL,
  days_duration_min = 90,
  colname_date_released = "date_released",
  colname_date_recovered = "date_recovered",
  colname_region_released = "region_released",
  colname_region_recovered = "region_recovered",
  colname_size_released = "size_released",
  movement_pattern = 2,
  movement_allow = NULL,
  movement_disallow = NULL,
  mu_movement_step_diag = NULL,
  sd_movement_step_diag = NULL,
  mu_fishing_rate = NULL,
  cv_fishing_rate = NULL,
  mu_selectivity = NULL,
  cv_selectivity = NULL,
  mu_fishing_weight = NULL,
  sd_fishing_weight = NULL,
  mu_natural_mortality_rate = NULL,
  sd_natural_mortality_rate = NULL,
  mu_reporting_rate = NULL,
  sd_reporting_rate = NULL,
  mu_initial_loss_rate = 0.1,
  sd_initial_loss_rate = 0.01,
  mu_ongoing_loss_rate = 0.02,
  sd_ongoing_loss_rate = 0.001,
  mu_dispersion = 1,
  sd_dispersion = 0.5,
  tolerance_expected = 1e-12,
  tolerance_fishing = 1e-12,
  data = NULL,
  chains = 1,
  step_size = 0.01,
  adapt_delta = 0.95,
  iter_warmup = 250,
  iter_sampling = 750,
  max_treedepth = 10,
  use_reduce_sum = FALSE,
  threads_per_chain = parallel::detectCores()/(2 * chains),
  refresh = 100,
  ...
)
}
\arguments{
\item{tag_data}{\code{\link[=data.frame]{data.frame()}}}

\item{list_regions}{\code{\link[=list]{list()}}}

\item{list_sizes}{\code{\link[=list]{list()}}}

\item{year_start}{\code{\link[=integer]{integer()}} year of initial tag released}

\item{year_end}{\code{\link[=integer]{integer()}} year of final tag recovered}

\item{step_interval}{\code{\link[=character]{character()}} one of \code{"month", "quarter", "year"}}

\item{step_duration_max}{\code{\link[=integer]{integer()}}}

\item{days_duration_min}{\code{\link[=integer]{integer()}} minimum number of days before recovery}

\item{colname_date_released}{\code{\link[=character]{character()}}}

\item{colname_date_recovered}{\code{\link[=character]{character()}}}

\item{colname_region_released}{\code{\link[=character]{character()}}}

\item{colname_region_recovered}{\code{\link[=character]{character()}}}

\item{colname_size_released}{\code{\link[=character]{character()}}}

\item{movement_pattern}{\code{\link[=integer]{integer()}}}

\item{movement_allow}{\link[=matrix]{integer()}}

\item{movement_disallow}{\link[=matrix]{integer()}}

\item{mu_movement_step_diag}{\link[=vector]{numeric()}}

\item{sd_movement_step_diag}{\link[=vector]{numeric()}}

\item{mu_fishing_rate}{\link[=array]{numeric()}}

\item{cv_fishing_rate}{\code{\link[=numeric]{numeric()}}}

\item{mu_selectivity}{\link[=array]{numeric()}}

\item{cv_selectivity}{\code{\link[=numeric]{numeric()}}}

\item{mu_fishing_weight}{\code{\link[=numeric]{numeric()}}}

\item{sd_fishing_weight}{\code{\link[=numeric]{numeric()}}}

\item{mu_natural_mortality_rate}{\code{\link[=numeric]{numeric()}}}

\item{sd_natural_mortality_rate}{\code{\link[=numeric]{numeric()}}}

\item{mu_reporting_rate}{\code{\link[=numeric]{numeric()}}}

\item{sd_reporting_rate}{\code{\link[=numeric]{numeric()}}}

\item{mu_initial_loss_rate}{\code{\link[=numeric]{numeric()}}}

\item{sd_initial_loss_rate}{\code{\link[=numeric]{numeric()}}}

\item{mu_ongoing_loss_rate}{\code{\link[=numeric]{numeric()}}}

\item{sd_ongoing_loss_rate}{\code{\link[=numeric]{numeric()}}}

\item{mu_dispersion}{\code{\link[=numeric]{numeric()}}}

\item{sd_dispersion}{\code{\link[=numeric]{numeric()}}}

\item{tolerance_expected}{\code{\link[=numeric]{numeric()}}}

\item{tolerance_fishing}{\code{\link[=numeric]{numeric()}}}

\item{data}{\code{\link[=list]{list()}} See details}

\item{chains}{\code{\link[=integer]{integer()}} number of chains}

\item{step_size}{\code{\link[=numeric]{numeric()}} initial step size}

\item{adapt_delta}{\code{\link[=numeric]{numeric()}} the adaptation target acceptance statistic}

\item{iter_warmup}{\code{\link[=integer]{integer()}} number of warmup iterations}

\item{iter_sampling}{\code{\link[=integer]{integer()}} number of sampling iterations}

\item{max_treedepth}{\code{\link[=integer]{integer()}}}

\item{use_reduce_sum}{\code{\link[=logical]{logical()}} use within chain parallel threading}

\item{threads_per_chain}{\code{\link[=integer]{integer()}} number of threads per chain}

\item{refresh}{\code{\link[=integer]{integer()}}}

\item{...}{additional arguments to pass to \code{$sample()} method}
}
\value{
\code{\link[=mmmstan]{mmmstan()}}
}
\description{
The Stan model data may be specified by individual arguments or by passing
a suitable \code{data} argument. The \code{data} argument takes precedence
when it is non-null.
}
\details{
When the model data are specified by individual arguments, a number of
helper functions with sensible defaults are used to assemble the \code{data}
object. This is usually preferable for user convenience.

Passing the \code{data} object directly may be preferred for example when the
\code{data} object is the output of a simulation model, or when greater
flexibility is desired. In either event, care must be taken to ensure the
\code{data} object matches the requirements of the underlying Stan model.

TBD
}
