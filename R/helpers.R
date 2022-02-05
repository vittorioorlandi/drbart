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
  
  mean_file_name <- basename(mean_file)
  mean_file_dir <- dirname(mean_file)
  
  prec_file_name <- basename(prec_file)
  prec_file_dir <- dirname(prec_file)
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

get_q_from_cdf <- function(p, grid, cdf) approxfun(cdf, grid, ties = 'ordered')(p)

get_mean_from_pdf <- function(grid, pdf) {
  pdf_fun <- approxfun(grid, pdf, ties = 'ordered') 
  return(integrate(function(z) pdf_fun(z) * z, 
                   lower = min(grid), upper = max(grid),
                   stop.on.error = FALSE)$value)
}
