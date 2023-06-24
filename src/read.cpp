#include <Rcpp.h>
#include <vector>
#include <ctime>

#include "funs.h"

using namespace Rcpp;

xinfo load_cutpoints(List xinfo_list, int p) {
  xinfo xi;
  xi.resize(p);
  for (int i = 0; i < p; ++i) {
    NumericVector tmp = xinfo_list[i];
    std::vector<double> tmp2;
    for (size_t j = 0; j < tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi[i] = tmp2;
  }
  return xi; 
}

std::vector<xinfo> load_cutpoints(List xinfo_lists, IntegerVector p) {
  int n_groups = p.size();
  std::vector<xinfo> all_xi (n_groups);
  for (size_t i = 0; i < n_groups; ++i) {
    all_xi[i] = load_cutpoints(xinfo_lists[i], p[i]);
  }
  return all_xi; 
}

std::vector<double> load_x(NumericVector x_) {
  std::vector<double> x; 
  for (NumericVector::iterator it = x_.begin(); it != x_.end(); ++it) {
    x.push_back(*it);
  }
  return x;
}

std::vector<std::vector<double>> load_x(List x_, int n_groups) {
  std::vector<std::vector<double>> x (n_groups);  
  for (size_t i = 0; i < n_groups; ++i) {
    x[i] = load_x(x_[i]);
  }
  return x;
}