// // [[Rcpp::depends(RcppParallel)]]
// #include <Rcpp.h>
// #include <RcppParallel.h>
// 
// using namespace std;
// using namespace Rcpp;
// using namespace RcppParallel;
// 
// template <class T>
// double logsumexp(T &x) {
//   double m = *std::max_element(x.begin(), x.end());
//   double s = 0.0;
//   typename T::iterator it;
//   for (it = x.begin(); it != x.end(); ++it) {
//     s += exp(*it - m);
//   }
//   return(m + log(s));
// }
// 
// NumericVector dmixnorm0(NumericVector& x, NumericVector& logprob, 
//                         NumericVector& mu, double& sd) {
//   NumericVector out(x.size());
//   std::vector<double> tmp(logprob.size());
//   for (int i = 0; i < x.size(); ++i) {
//     for (int h = 0; h < logprob.size(); ++h) {
//       tmp[h] = logprob(h) + R::dnorm(x(i), mu(h), sd, 1);
//     }
//     out(i) = logsumexp(tmp);
//   }
//   return(out);
// }
// 
// 
// struct mydmix : public Worker
// {
//   // source matrix
//   const RMatrix<double> input;
//   
//   // destination matrix
//   RVector<double> output;
//   
//   // initialize with source and destination
//   mydmix(const NumericMatrix input, NumericMatrix output) 
//     : input(input), output(output) {}
//   
//   // take the square root of the range of elements requested
//   void operator()(std::size_t begin, std::size_t end) {
//     std::transform(input.begin() + begin, 
//                    input.begin() + end, 
//                    output.begin() + begin, 
//                    ::sqrt);
//   }
// };
// 
// //[[Rcpp::export]]
// NumericMatrix dmixnorm0_post(NumericVector x, List mus, NumericVector sd, List logprobs) {
//   NumericVector mu, logprob;
//   NumericMatrix out(x.size(), mus.size());
//   for (int j = 0; j < mus.size(); ++j) {
//     mu = as<NumericVector>(mus[j]);
//     logprob = as<NumericVector>(logprobs[j]);
//     out(_, j) = dmixnorm0(x, logprob, mu, sd(j));
//   }
//   return out;
// }