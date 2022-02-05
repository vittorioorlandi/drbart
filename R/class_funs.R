Rcpp::loadModule("TreeSamples", TRUE)

get_q_from_cdf <- function(p, grid, cdf) approxfun(cdf, grid, ties = 'ordered')(p)

get_mean_from_pdf <- function(grid, pdf) {
  pdf_fun <- approxfun(grid, pdf, ties = 'ordered') 
  return(integrate(function(z) pdf_fun(z) * z, 
                   lower = min(grid), upper = max(grid))$value)
}

#' Posterior Predictive Quantities from DR-BART
#' 
#' Compute and plot conditional density functions, distribution functions,
#' quantiles, and means.
#' 
#' \code{predict.drbart} returns posterior predictive draws of functionals of the
#' conditional densities and is called by \code{plot.drbart} to generate appropriate
#' plots. The functional is determined by the value of \code{type}; either entire
#' density or distribution functions, selected quantiles, or means may be
#' estimated. Because quantiles and means are derived from the relevant
#' densities themselves, in all cases, prediction is a relatively time intensive
#' procedure, dependent on the number of posterior samples and on the size of
#' \code{ygrid}. Note also that estimated quantities will be inaccurate should
#' \code{ygrid} not fully capture the high density regions of the conditional
#' densities.
#' 
#' 
#' @param object,x An object of class \code{drbart}, returned by a call to
#'   \code{\link{drbart}}.
#' @param xpred A matrix of points describing which conditional densities /
#'   distributions should be estimated. Rows correspond to different conditional
#'   densities and columns to different covariates.
#' @param ygrid A numeric vector of y values at which the conditional density /
#'   distribution should be evaluated.
#' @param type Type of predictions to be returned. If \code{'density'}, returns an
#'   estimate of the conditional density functions (pdfs) evaluated at points
#'   in \code{ygrid}. If \code{'distribution'}, returns an estimate of the conditional
#'   distribution functions (cdfs), evaluated at points in \code{'ygrid'}. If
#'   \code{'quantiles'}, returns the quantiles of the conditional density specified
#'   by \code{quantiles}. If \code{mean}, returns the conditional mean.
#' @param quantiles If \code{type = 'quantiles'}, the quantiles of the conditional
#'   densities that should be estimated.
#' @param ... Ignored.
#' @name Methods
#' @return If predicting, an array containing posterior draws of the conditional
#'   densities, conditional distributions, conditional quantiles, or conditional
#'   means, depending on the value of `type`. If plotting, the original `drbart`
#'   object is invisibly returned.
#' @export
predict.drbart <- function(object, xpred, ygrid, 
                           type = c('density', 'distribution', 
                                    'quantiles', 'mean'),
                           quantiles = c(0.025, 0.5, 0.975), ...) {
  
  
  type <- match.arg(type)
  
  mean_file <- object$mean_file
  stopifnot(file.exists(mean_file))
  
  variance <- object$variance
  if (variance != 'const') {
    prec_file <- object$prec_file
    stopifnot(file.exists(prec_file))
  }
  
  ts_mean <- TreeSamples$new()
  ts_mean$load(mean_file)
  
  if (variance != 'const') {
    ts_prec <- TreeSamples$new()
    ts_prec$load(prec_file)
  }
  
  fit <- object$fit
  
  logprobs <- lapply(fit$ucuts, function(u) log(diff(c(0, u, 1))))
  mids <- lapply(fit$ucuts, function(u) c(0, u) + diff(c(0, u, 1)) / 2)
  
  nsim <- length(fit$ucuts)
  
  if (type == 'mean') {
    preds <- array(dim = c(nrow(xpred), 1, nsim))
  }
  else if (type == 'quantiles') {
    preds <- array(dim = c(nrow(xpred), length(quantiles), nsim))
  }
  else {
    preds <- array(dim = c(nrow(xpred), length(ygrid), nsim))
  }
  
  for (j in seq_len(nrow(xpred))) {
    message(paste0('Predicting conditional density ', j, ' of ', nrow(xpred), 
                   ' (', round(100 * (j - 1) / nrow(xpred)), '%)', '\r'),
            appendLF = FALSE)
    flush.console()
    
    
    des <- lapply(mids, function(m) do.call(data.frame, c(list(m), as.list(xpred[j, ]))))
    mu <- lapply(seq_along(fit$ucuts), function(i) {
      c(ts_mean$predict_i(t(des[[i]]), i - 1))
    })
    
    if (variance == 'const' | variance == 'x') {
      if (type == 'quantiles') {
        post_fun <- pmixnorm0_post
      }
      else {
        post_fun <- dmixnorm0_post
      }
    }
    else {
      if (type == 'quantiles') {
        post_fun <- pmixnorm_post
      }
      else {
        post_fun <- dmixnorm_post
      }
    }
    
    if (variance == 'const') {
      sigma <- fit$sigma
    }
    else if (variance == 'x') {
      sigma <- 1 / sqrt(fit$phistar * ts_prec$predict_prec(matrix(xpred[j, ])))
    }
    else {
      phi <- lapply(seq_along(fit$ucuts), function(i) {
        fit$phistar[i] * c(ts_prec$predict_prec_i(t(des[[i]]), i - 1))
      })
      sigma <- lapply(phi, function(x) 1 / sqrt(x))
    }
    post <- exp(post_fun(ygrid, mu, sigma, logprobs))
    
    if (type == 'mean') {
      preds[j, 1, ] <- apply(post, 2, function(pdf) get_mean_from_pdf(ygrid, pdf))
    }
    else if (type == 'quantiles') {
      preds[j, , ] <- apply(post, 2, function(cdf) get_q_from_cdf(quantiles, ygrid, cdf))
    }
    else {
      preds[j, , ] <- post
    }
  }
  
  return(preds)
}

#' Plot Conditional Densities
#' @param CI Whether credible intervals should be plotted.
#' @param alpha If \code{CI = TRUE}, plots \code{100(1 - alpha)\%} credible intervals.
#' @param legend_position A legend position keyword as accepted by
#'   \code{\link{legend}} or \code{'none'} if the legend should be omitted.
#' @rdname Methods
#' @export
plot.drbart <- function(x, xpred, ygrid, 
                        type = c('density', 'distribution', 'quantiles', 'mean'), 
                        quantiles = c(0.025, 0.5, 0.975),
                        CI = FALSE, alpha = 0.05, legend_position) {
  
  type <- match.arg(type)
  
  if (!is.matrix(xpred)) {
    if (length(xpred) != 1) {
      stop(paste0('`xpred` must be a matrix, unless there is one covariate ',
                  'and one conditional density or distribution to be plotted, in which case a scalar is acceptable.'))
    }
    xpred <- matrix(xpred)
  }
  
  n_unique <- apply(xpred, 2, function(col) length(unique(col)))
  
  if (type == 'density' || type == 'distribution') {
    n_colors <- nrow(xpred)
  }
  else if (type == 'quantiles') {
    n_colors <- length(quantiles)
  }
  else {
    n_colors <- 1
  }
  
  if (nrow(xpred) > 1) {
    non_const <- which(n_unique > 1)
    if (length(non_const) != 1) {
      stop('Only one non-constant column may be provided.')
    }
    if (n_colors > 9) {
      warning(paste0('The used color palette only allows for 9 distinct ', 
                     'colors, while you are asking to plot ', nrow(xpred), ' ',
                     'different conditional densities. Coloring and an ',
                     'associated legend to identify the conditional densities ',
                     'will be suppressed.'))
    }
  }
  else {
    non_const <- 1
  }
  
  vals <- sort(xpred[, non_const], index.return = TRUE)
  xpred <- xpred[vals$ix, , drop = FALSE]
  vals <- vals$x
  
  colors <- 
    tryCatch(
      RColorBrewer::brewer.pal(n_colors, 'Blues'),
      warning = function(cnd) {
        if (n_colors > 9) {
          rep('black', n_colors)
        }
        else { # <= 2
          RColorBrewer::brewer.pal(3, 'Blues')[1 + seq_len(n_colors)]
        }
      }
    )
  
  if (type == 'density' | type == 'distribution') {
    estimand_size <- length(ygrid)
  }
  else if (type == 'quantiles') {
    estimand_size <- length(quantiles)
  }
  else {
    estimand_size <- 1
  }
  
  all_preds <- predict(x, xpred, ygrid, type = type, quantiles = quantiles)
  
  if (CI) {
    summary_preds <- array(dim = c(nrow(xpred), estimand_size, 3))
  }
  else {
    summary_preds <- array(dim = c(nrow(xpred), estimand_size, 1))
  }
  
  summary_preds[, , 1] <- apply(all_preds, 1:2, mean)
  if (CI) {
    summary_preds[, , 2] <- apply(all_preds, 1:2, quantile, alpha / 2)
    summary_preds[, , 3] <- apply(all_preds, 1:2, quantile, 1 - alpha / 2)
  }
  
  limits <- range(summary_preds) + c(0, 0.05)
  
  if (type == 'density' | type == 'distribution') {
    if (type == 'density') {
      ylab <- 'p(y|x)'
    }
    else {
      ylab <- 'P(y|x)'
    }
    for (i in seq_len(nrow(xpred))) {
      if (i == 1) {
        plot(ygrid, summary_preds[i, , 1], type = 'l', col = colors[i], 
             xlab = 'y', ylab = ylab, ylim = limits)
      }
      else {
        lines(ygrid, summary_preds[i, , 1], col = colors[i])
      }
      if (CI) {
        lines(ygrid, summary_preds[i, , 2], lty = 'dashed', col = colors[i])
        lines(ygrid, summary_preds[i, , 3], lty = 'dashed', col = colors[i])
      }
    }
  }
  else if (type == 'quantiles') {
    for (i in seq_along(quantiles)) {
      if (i == 1) {
        plot(vals, summary_preds[, i, 1], type = 'l', col = colors[i],
             xlab = 'x', ylab = 'Q(y|x)', ylim = limits)
      }
      else {
        lines(vals, summary_preds[, i, 1], type = 'l', col = colors[i])
      }
      if (CI) {
        lines(vals, summary_preds[, i, 2], lty = 'dashed', col = colors[i])
        lines(vals, summary_preds[, i, 3], lty = 'dashed', col = colors[i])
      }
    }
  }
  else {
    plot(vals, summary_preds[, 1, 1], type = 'l', col = colors[1], 
         xlab = 'x', ylab = 'E(y|x)')
    if (CI) {
      lines(vals, summary_preds[, 1, 2], lty = 'dashed', col = colors[1])
      lines(vals, summary_preds[, 1, 3], lty = 'dashed', col = colors[1])
    }
  }
  
  if (missing(legend_position)) {
    legend_position <- 'topleft'
  }
  
  if (type == 'mean' | legend_position == 'none' | n_colors > 9) {
    return(invisible(x))
  }
  
  if (type == 'quantiles') {
    legend(legend_position, legend = quantiles, fill = colors)  
  }
  else {
    legend(legend_position, legend = vals, fill = colors)  
  }

  return(invisible(x))
}
# 
# if (object$variance == 'const') {
#   hist(object$fit$sigma ^ 2, xlab = expression(sigma ^ 2), main = '')
# }
