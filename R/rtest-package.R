#' The 'rtest' package.
#'
#' @docType package
#' @name    rtest-package
#' @aliases rtest
#' @useDynLib rtest, .registration = TRUE
#'
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @import stats

#' @importFrom rstan sampling extract stanc rstan_options traceplot stan_rhat
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines par plot points text arrows grid rect
#' @importFrom parallel detectCores
#' @importFrom utils as.roman
#'
#' @description A template for R package that uses rstan and Rcpp
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#'
NULL
