#' Density Regression with Bayesian Additive Regression Trees
#'
#' Fits a density regression model using DR-BART.
#'
#' \code{drbart} fits a density regression model using DR-BART as described in
#' Orlandi, Murray, Linero, and Volfovsky (2021)
#' \href{https://arxiv.org/abs/2112.12259}{here}. DR-BART uses a continuous
#' latent variable model that generalizes covariate-dependent discrete mixture
#' models to flexibly and accurately model conditional densities and perform
#' density regression. The most flexible version of DR-BART models both the mean
#' and variance functions using BART priors that are functions of observed
#' covariates and a latent \eqn{U}. Two restricted models are also available and
#' can be specified by the \code{variance} argument. DR-BART-LH only models the
#' variance as a flexible function of covariates (\code{variance = 'x'}) and
#' DR-BART-L specifies a constant variance (\code{variance = 'const'}). See the
#' paper for more details.
#'
#' Given a continuous covariate \eqn{x}, there are infinitely many conditional
#' densities \eqn{p(y | x)} that may be of interest. Because of this, the notion
#' of fitted values or an in-sample fit becomes less important in a density
#' regression setting. For this reason, \code{drbart} stores information about
#' the fits (namely the trees) in files specified by \code{mean_file} and
#' \code{prec_file} (created if not already existing). This information is then
#' used to estimate certain conditional densities requested by the user via
#' \code{\link{predict.drbart}} and so should not be deleted until all desired
#' predictions have been made. This behavior may change in the future.
#'
#' Hyperparameters for the BART prior
#' can be modified and are described briefly above. For full details, see CGM
#' 2010 \href{https://arxiv.org/pdf/0806.3286.pdf}{here}.
#'
#' @param y A vector of observed responses.
#' @param x A matrix of observed covariates. Rows correspond to observations and
#'   columns to different covariates.
#' @param nburn Number of MCMC burn-in iterations
#' @param nsim Number of MCMC iterations to be returned.
#' @param nthin Thinning parameter for the MCMC -- the number of iterations run
#'   before one is stored.
#' @param printevery How often should the number of MCMC iterations be printed.
#' @param m_mean,m_var Number of trees used to model the mean, variance
#'   functions.
#' @param alpha,beta Hyperparameters for the tree splitting prior. The
#'   probability that a node at depth d splits is given by \eqn{\alpha / (1 + d)
#'   ^ \beta}. Default values are \eqn{\alpha = 0.95} and \eqn{\beta = 2}.
#' @param lambda,nu Hyperparameters in the prior for sigma.
#' @param kfac Hyperparameter in the variance term of the prior for the mu.
#'   Larger values imply larger shrinkage of the fit towards 0.
#' @param phi0 Baseline precision if \code{variance = 'x'} or \code{variance =
#'   'ux'}.
#' @param variance One of \code{'ux'}, \code{'x'}, or \code{'const'}. If
#'   \code{'ux'}, the variance term is treated as a function of both the latent
#'   u and the observed x, which can be helpful for flexible learning in finite
#'   sample settings. If \code{'x'}, the variance is solely a function of x,
#'   allowing for heteroscedasticity in the observed covariates. If
#'   \code{'const'}, an assumption of constant variance is made.
#' @param censor A logical vector denoting which values of \code{y}, if any, are
#'   censored from below. Such observations will be imputed within the sampler.
#'   Defaults to assuming no observations are censored.
#' @param mean_file,prec_file File location for information about the mean,
#'   variance fit, respectively (trees, number of covariates, etc.). Primarily
#'   used for prediction after model fitting.
#' @param mean_cuts,prec_cuts Optional. Cut points to consider when building
#'   trees to model the mean, variance function, respectively. If supplied, a
#'   list of length \code{ncol(x)} containing the split points associated with
#'   each covariate. Otherwise, an intelligent choice based off observed
#'   covariate values will be made.
#'
#' @return An object of class `drbart`, containing:
#'
#' @importFrom utils flush.console
#' @importFrom graphics legend lines points arrows
#' @importFrom stats approxfun integrate quantile runif
#' @importFrom methods new
#' @export
#'
#' @seealso \code{\link{predict.drbart}}, \code{\link{plot.drbart}}.
#'
drbart <- function(y, x,
                   nburn = 5000, nsim = 5000,  nthin = 1,
                   printevery = round((nburn + nsim * nthin) / 20),
                   m_mean = 200, m_var = 100, alpha = 0.95, beta = 2,
                   lambda = 1,
                   nu = 2, kfac = 2, phi0 = 1,
                   variance = c('ux', 'x', 'const'),
                   censor = logical(length(y)),
                   mean_file = 'dr_bart_mean.txt',
                   prec_file = 'dr_bart_prec.txt',
                   mean_cuts, prec_cuts) {

  x <-
    check_args(x, y, nburn, nsim, nthin, m_mean,
             m_var, alpha, beta, lambda, nu, kfac, censor,
             mean_file, prec_file)

  # No actual way of preventing people from passing in (u, x)
  variance <- match.arg(variance)

  n <- dim(x)[1]
  p <- dim(x)[2]

  ux <- cbind(runif(n), x)

  if (missing(mean_cuts)) {
    mean_cuts <- lapply(data.frame(x), .cp_quantile) # check just apply
  }
  mean_cuts <- c(list((1:9999) / 10000), mean_cuts)

  if (missing(prec_cuts)) {
    if (variance == 'ux') {
      prec_cuts <- mean_cuts
    }
    else {
      prec_cuts <- mean_cuts[-1]
    }
  }
  # the below can be misleading if they pass in u -- but works so fine for now
  stopifnot(length(mean_cuts) == p + 1)
  if (variance == 'ux') {
    stopifnot(length(prec_cuts) == p + 1)
  }
  else {
    stopifnot(length(prec_cuts) == p)
  }

  if (variance == 'ux') {
    out <- drbartRcppHeteroClean(y, t(ux), t(ux),
                                 mean_cuts, prec_cuts,
                                 nburn, nsim, nthin, printevery,
                                 m_mean, m_var, alpha, beta,
                                 lambda, nu, kfac, phi0,
                                 TRUE,
                                 censor,
                                 mean_file, prec_file)
  }
  else if (variance == 'x') {
    out <- drbartRcppHeteroClean(y, t(ux), t(x),
                                 mean_cuts, prec_cuts,
                                 nburn, nsim, nthin, printevery,
                                 m_mean, m_var, alpha, beta,
                                 lambda, nu, kfac, phi0,
                                 FALSE,
                                 censor,
                                 mean_file, prec_file)
  }
  else {
    out <- drbartRcppClean(y, t(ux), t(ux[1, ]),
                           mean_cuts,
                           nburn, nsim, nthin, printevery,
                           m_mean, alpha, beta,
                           lambda, nu, kfac,
                           censor, mean_file)
  }
  out <- list(fit = out,
              variance = variance,
              mean_file = mean_file)

  if (variance != 'const') {
    out <- c(out, list(prec_file = prec_file))
  }

  class(out) <- 'drbart'
  return(out)
}
