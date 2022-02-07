predict_serial <- function(xpred, mids, fit, ts_mean, ts_prec, type, variance, quantiles, ygrid, preds) {
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
    
    if (variance == 'const' | variance == 'x') {
      if (type == 'quantiles' | type == 'distribution') {
        post_fun <- pmixnorm0_post
      }
      else {
        post_fun <- dmixnorm0_post
      }
    }
    else {
      if (type == 'quantiles' | type == 'distribution') {
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
      preds[j, 1, ] <-
        apply(post, 2, function(pdf) get_mean_from_pdf(ygrid, pdf))
    }
    else if (type == 'quantiles') {
      preds[j, , ] <-
        apply(post, 2, function(cdf) get_q_from_cdf(quantiles, ygrid, cdf))
    }
    else {
      preds[j, , ] <- post
    }
  }
  return(preds)
}

predict_parallel <- function(xpred, mids, fit, ts_mean, ts_prec, type, variance, quantiles, ygrid, preds) {
  foreach(j=1:nrow(xpred)) %dopar% {
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
    
    if (variance == 'const' | variance == 'x') {
      if (type == 'quantiles' | type == 'distribution') {
        post_fun <- pmixnorm0_post
      }
      else {
        post_fun <- dmixnorm0_post
      }
    }
    else {
      if (type == 'quantiles' | type == 'distribution') {
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
      preds[j, 1, ] <-
        apply(post, 2, function(pdf) get_mean_from_pdf(ygrid, pdf))
    }
    else if (type == 'quantiles') {
      preds[j, , ] <-
        apply(post, 2, function(cdf) get_q_from_cdf(quantiles, ygrid, cdf))
    }
    else {
      preds[j, , ] <- post
    }
  }
  return(preds)
}
