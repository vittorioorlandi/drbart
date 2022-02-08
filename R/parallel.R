predict_serial <- function(xpred, mids, fit, ts_mean, ts_prec, type, variance, quantiles, ygrid, preds, logprobs, post_fun) {
  for (j in seq_len(nrow(xpred))) {
    message(paste0('Predicting conditional density ', j, ' of ', nrow(xpred),
                   ' (', round(100 * (j - 1) / nrow(xpred)), '%)', '\r'),
            appendLF = FALSE)
    flush.console()
    
    
    des <- lapply(mids, function(m) {
      do.call(data.frame, c(list(m), as.list(xpred[j, ])))
    })
    mu <- lapply(seq_along(fit$ucuts), function(i) {
      c(ts_mean$predict_i(t(des[[i]]), i - 1))
    })
    
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
    
    preds[j, , ] <- post
  }
  return(preds)
}

predict_parallel <- function(xpred, mids, fit, ts_mean, ts_prec, type, variance, quantiles, ygrid, logprobs, post_fun) {
  
  if (!requireNamespace("abind", quietly = TRUE)) {
    warning(paste0('We recommend that you install package `abind` if ',
                   'parallelizing predictions.'), call. = FALSE)
    combine_fun <- function(x, y) {
      if (length(dim(x)) == 2 & length(dim(y)) == 2) {
        arr <- array(dim = c(2, dim(x)))
        arr[1, , ] <- x
        arr[2, , ] <- y
      }
      else if (length(dim(x)) == 2) {
        arr <- array(dim = c(dim(y)[1] + 1, dim(y)[2:3]))
        arr[1, , ] <- x
        arr[2:(dim(y)[1] + 1), , ] <- y
      }
      else if (length(dim(y)) == 2) {
        arr <- array(dim = c(dim(x)[1] + 1, dim(x)[2:3]))
        arr[1:dim(x)[1], , ] <- x
        arr[dim(x)[1] + 1, , ] <- y
      }
      else {
        arr <- array(dim = c(dim(x)[1] + dim(y)[1], dim(x)[2:3]))
        arr[1:dim(x)[1], , ] <- x
        arr[(dim(x)[1] + 1):(dim(x)[1] + dim(y)[1]), , ] <- y
      }
      return(arr)
    }
    multicombine <- FALSE
    packages <- NULL
  }
  else {
    combine_fun <- function(...) {
      return(abind::abind(..., along = 0))
    }
    multicombine <- TRUE
    packages <- 'abind'
  }
  
  j <- NULL # To appease R CMD CHECK 
  
  preds <- 
    foreach::`%dopar%`(foreach::foreach(j=1:nrow(xpred), 
                                        .combine = combine_fun, 
                                        .multicombine = multicombine,
                                        .packages = packages),{
                       des <- lapply(mids, function(m) {
                         do.call(data.frame, c(list(m), as.list(xpred[j, ])))
                       })
                       mu <- lapply(seq_along(fit$ucuts), function(i) {
                         c(ts_mean$predict_i(t(des[[i]]), i - 1))
                       })
                       
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
                       post <- exp(post_fun(ygrid, mu, sigma, logprobs))}
    )
  return(preds)
}
