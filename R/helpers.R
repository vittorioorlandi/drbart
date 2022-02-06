check_args <- function(x, y,
                       nburn, nsim, nthin,
                       m_mean, m_var, alpha, beta, lambda, nu, kfac, censor,
                       mean_file, prec_file) {

  stopifnot(is.vector(y) && is.atomic(y))
  if (is.vector(x) && is.atomic(x)) {
    x <- matrix(x, ncol = 1)
  }
  stopifnot(is.matrix(x))
  stopifnot(dim(x)[1] == length(y))

  stopifnot(0 <= nburn)
  stopifnot(0 < nsim)
  stopifnot(nthin <= nsim)
  stopifnot(0 < m_mean)
  stopifnot(0 < m_var)
  stopifnot(0 < alpha & alpha < 1)
  stopifnot(0 <= beta)
  stopifnot(0 < lambda)
  stopifnot(0 < nu)
  stopifnot(0 < kfac)
  stopifnot(length(censor) == length(y))
  return(x)
}


.cp_quantile <- function(x, num = 10000, cat_levels = 8) {
  # BCF function for supplying BART split points

  nobs <- length(x)
  nuniq <- length(unique(x))

  if (nuniq == 1) {
    ret <- x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if (nuniq < cat_levels) {
    xx <- sort(unique(x))
    ret <- xx[-length(xx)] + diff(xx) / 2
  } else {
    q <- approxfun(sort(x), quantile(x, p = 0:(nobs - 1) / nobs))
    ind <- seq(min(x), max(x), length.out = num)
    ret <- q(ind)
  }

  return(ret)
}

get_q_from_cdf <- function(p, grid, cdf) {
  approxfun(cdf, grid, ties = 'ordered')(p)
}

get_mean_from_pdf <- function(grid, pdf) {
  pdf_fun <- approxfun(grid, pdf, ties = 'ordered')
  return(integrate(function(z) pdf_fun(z) * z,
                   lower = min(grid), upper = max(grid),
                   stop.on.error = FALSE)$value)
}

preprocess_plot_args <- function(xpred, ygrid, type, quantiles, CI, alpha,
                                 legend_position) {

  match.arg(legend_position,
            c("bottomright", "bottom", "bottomleft",
              "left", "topleft", "top", "topright", "right", "center", "none"))

  stopifnot(is.logical(CI) & length(CI) == 1)
  stopifnot(is.numeric(alpha) & is.null(dim(alpha))
            & length(alpha) == 1 & 0 < alpha & alpha < 1)
  stopifnot(is.numeric(quantiles) & is.null(dim(quantiles))
            & all(0 < quantiles) & all(quantiles < 1))

  if (!is.matrix(xpred)) {
    if (length(xpred) != 1) {
      stop(paste0('`xpred` must be a matrix, unless there is one covariate ',
                  'and one conditional density or distribution to be plotted, ',
                  'in which case a scalar is acceptable.'))
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
                     'colors, while you are asking to plot information from',
                     nrow(xpred), ' different conditional densiies. Coloring ',
                     'and an associated legend to identify the conditional ',
                     'densities will be suppressed.'))
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

  return(list(colors = colors, estimand_size = estimand_size,
              vals = vals, xpred = xpred))
}
