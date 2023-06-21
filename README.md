# mmmstan
An R package to fit Markov movement models (MMM) via CmdStan and CmdStanR

<!-- badges: start -->
[![R-CMD-check](https://github.com/luke-a-rogers/mmmstan/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luke-a-rogers/mmmstan/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## CmdStan / CmdStanR

The mmmstan package currently depends on the cmdstanr package to take advantage
of within-chain parallel threading provided by the `reduce_sum()` function in
CmdStan >= 2.23.

## Installation

1. Install the R package `cmdstanr` (see <https://mc-stan.org/cmdstanr/index.html>).

``` r
remotes::install_github("stan-dev/cmdstanr")
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
