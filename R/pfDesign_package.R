#' The 'pfDesign' package.
#'
#' @docType package
#' @name    pfDesign-package
#' @aliases pfDesign
#' @useDynLib pfDesign, .registration = TRUE
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
#' @importFrom mvtnorm   rmvnorm
#'
#' @description An Analysis Toolbox for Platform Design
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan.
#' R package version 2.18.2. http://mc-stan.org
#'
NULL

#' Parameters
#'
#' @name parameters
#'
#' @param epsilon mixing proportion from the posterior distribution
#' @param nsmps number of samples to draw
#' @param mix_ab shape parameters for Beta
#'
#'
NULL
