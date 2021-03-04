# mmmstan
An R package to fit Markov movement models (MMM) via CmdStan

<!-- badges: start -->
[![R-CMD-check](https://github.com/luke-a-rogers/mmmstan/workflows/R-CMD-check/badge.svg)](https://github.com/luke-a-rogers/mmmstan/actions)
<!-- badges: end -->

## R-CMD-check

The mmmstan package currently depends on the cmdstanr package to take advantage
of within-chain parallel threading proviced by the `reduce_sum()` function. The 
cmdstanr package is not currently available on CRAN, so R-CMD-check will fail 
until (I imagine) I specify the appropriate repository in the DESCRIPTION and 
convince the system that performs the check to install CmdStan.

## Installation

1. Install the R package `cmdstanr` (see <https://mc-stan.org/cmdstanr/index.html>).

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

2. Install CmdStan (see <https://mc-stan.org/cmdstanr/articles/cmdstanr.html>).

``` r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(cores = parallel::detectCores())
```

3. Install `mmmstan`.

``` r
devtools::install_github("luke-a-rogers/mmmstan")
```
