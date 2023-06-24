#include <Rcpp.h>
#include <vector>
#include <ctime>
#include "tree.h"

using namespace Rcpp;

xinfo load_cutpoints(List xinfo_list, int p);
std::vector<xinfo> load_cutpoints(List xinfo_lists, IntegerVector p);

std::vector<double> load_x(NumericVector x_);
std::vector<std::vector<double>> load_x(List x_, int n_groups);

