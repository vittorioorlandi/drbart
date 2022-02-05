\name{drbart-package}
\alias{drbart-package}
\docType{package}
\title{
  An implemention of Density Regression with Bayesian Additive Regression Trees (DR-BART)
}
\description{
  This package implements the DR-BART algorithm for flexible and accurate density regression using Bayesian Additive Regression Trees.
}
\details{
  \code{\link{drbart}} is the workhorse of the package; it fits a density regression model to supplied data and returns an object of class \code{drbart}. To estimate and visualize (functionals of) conditional densities, \code{predict} and \code{plot} methods are provided for this class. 
}
\author{
Vittorio Orlandi
Jared Murray

Maintainer: Vittorio Orlandi <vdo@duke.edu>
}
\references{
  Orlandi, Murray, Linero, and Volfovsky (2021) \emph{Density Regression with Bayesian Additive Regression Trees}.
  \url{https://arxiv.org/abs/2112.12259}.
}
\keyword{ package }
\seealso{
  \code{\link{drbart}}
}
