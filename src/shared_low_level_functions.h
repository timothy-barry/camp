#include <Rcpp.h>
using namespace Rcpp;

#ifndef COMPUTE_EMPIRICAL_P_VALUE
#define COMPUTE_EMPIRICAL_P_VALUE
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side);
#endif
