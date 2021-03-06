% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pfDesign_filter.R
\name{pdTran}
\alias{pdTran}
\title{HMM transition function}
\usage{
pdTran(smps, epsilon = 1, f_pi0 = pdBinPi0, ...)
}
\arguments{
\item{epsilon}{mixing proportion from the posterior distribution}
}
\value{
column 1: theta; column 2: mixture indicator
}
\description{
HMM transition function
}
