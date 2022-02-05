#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

template <class T>
double logsumexp(T &x) {
  double m = *std::max_element(x.begin(), x.end());
  double s = 0.0;
  typename T::iterator it;
  for (it = x.begin(); it != x.end(); ++it) {
    s += exp(*it - m);
  }
  return(m + log(s));
}

NumericVector dmixnorm0(NumericVector& x, NumericVector& logprob, 
                        NumericVector& mu, double& sd) {
  NumericVector out(x.size());
  std::vector<double> tmp(logprob.size());
  for (int i = 0; i < x.size(); ++i) {
    for (int h = 0; h < logprob.size(); ++h) {
      tmp[h] = logprob(h) + R::dnorm(x(i), mu(h), sd, 1);
    }
    out(i) = logsumexp(tmp);
  }
  return(out);
}

NumericVector pmixnorm0(NumericVector& x, NumericVector& logprob, 
                        NumericVector& mu, double& sd) {
  NumericVector out(x.size());
  std::vector<double> tmp(logprob.size());
  for (int i = 0; i < x.size(); ++i) {
    for (int h = 0; h < logprob.size(); ++h) {
      tmp[h] = logprob(h) + R::pnorm(x(i), mu(h), sd, 1, 1);
    }
    out(i) = logsumexp(tmp);
  }
  return(out);
}

NumericVector dmixnorm(NumericVector& x, NumericVector& logprob, 
                       NumericVector& mu, NumericVector& sd) {
  NumericVector out(x.size());
  std::vector<double> tmp(logprob.size());
  for(int i=0; i<x.size(); ++i) {
    for(int h=0; h<logprob.size(); ++h) {
      tmp[h] = logprob(h) + R::dnorm(x(i), mu(h), sd(h), 1);
    }
    out(i) = logsumexp(tmp);
  }
  return(out);
}

NumericVector pmixnorm(NumericVector& x, NumericVector& logprob, 
                       NumericVector& mu, NumericVector& sd) {
  NumericVector out(x.size());
  std::vector<double> tmp(logprob.size());
  for(int i=0; i<x.size(); ++i) {
    for(int h=0; h<logprob.size(); ++h) {
      tmp[h] = logprob(h) + R::pnorm(x(i), mu(h), sd(h), 1, 1);
    }
    out(i) = logsumexp(tmp);
  }
  return(out);
}

//[[Rcpp::export]]
NumericMatrix dmixnorm0_post(NumericVector x, List mus, NumericVector sd, List logprobs) {
  NumericVector mu, logprob;
  NumericMatrix out(x.size(), mus.size());
  for (int j = 0; j < mus.size(); ++j) {
    mu = as<NumericVector>(mus[j]);
    //sd = as<NumericVector>(sds[j]);
    logprob = as<NumericVector>(logprobs[j]);
    out(_, j) = dmixnorm0(x, logprob, mu, sd(j));
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix pmixnorm0_post(NumericVector x, List mus, NumericVector sd, List logprobs) {
  NumericVector mu, logprob;
  NumericMatrix out(x.size(), mus.size());
  for (int j = 0; j < mus.size(); ++j) {
    mu = as<NumericVector>(mus[j]); // checked: mu and logprobs are the same length
    logprob = as<NumericVector>(logprobs[j]);
    out(_,j) = pmixnorm0(x, logprob, mu, sd(j));
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix dmixnorm_post(NumericVector x, List mus, List sds, List logprobs) {
  NumericVector mu, sd, logprob;
  NumericMatrix out(x.size(), mus.size());
  for(int j=0; j<mus.size(); ++j) {
    mu = as<NumericVector>(mus[j]);
    sd = as<NumericVector>(sds[j]);
    logprob = as<NumericVector>(logprobs[j]);
    out(_,j) = dmixnorm(x, logprob, mu, sd);
  }
  return out;
}


//[[Rcpp::export]]
NumericMatrix pmixnorm_post(NumericVector x, List mus, List sds, List logprobs) {
  NumericVector mu, sd, logprob;
  NumericMatrix out(x.size(), mus.size());
  for(int j=0; j<mus.size(); ++j) {
    mu = as<NumericVector>(mus[j]);
    sd = as<NumericVector>(sds[j]);
    logprob = as<NumericVector>(logprobs[j]);
    out(_,j) = pmixnorm(x, logprob, mu, sd);
  }
  return out;
}
