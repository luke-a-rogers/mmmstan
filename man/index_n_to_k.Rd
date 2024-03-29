% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{index_n_to_k}
\alias{index_n_to_k}
\title{Index N to K}
\usage{
index_n_to_k(year_start, year_end, step_interval, term_interval)
}
\arguments{
\item{year_start}{\code{\link[=integer]{integer()}} year of first tag release}

\item{year_end}{\code{\link[=integer]{integer()}} year of final tag recovery}

\item{step_interval}{\code{\link[=character]{character()}} one of \code{"month"}, \code{"quarter"}
or \code{"year"}}

\item{term_interval}{\code{\link[=character]{character()}} one of \code{"month"}, \code{"quarter"}
or \code{"year"}}
}
\value{
\link[=vector]{integer()}
}
\description{
Index N to K
}
\examples{

index_n_to_k(2011, 2018, "month", "month")
index_n_to_k(2011, 2018, "month", "quarter")
index_n_to_k(2011, 2018, "month", "year")
index_n_to_k(2011, 2018, "quarter", "quarter")
index_n_to_k(2011, 2018, "quarter", "year")
index_n_to_k(2011, 2018, "year", "year")

}
