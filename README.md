# mmmstan
An R package to fit Markov movement models (MMM) via RStan

<!-- badges: start -->
[![R-CMD-check](https://github.com/luke-a-rogers/mmmstan/workflows/R-CMD-check/badge.svg)](https://github.com/luke-a-rogers/mmmstan/actions)
<!-- badges: end -->

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
