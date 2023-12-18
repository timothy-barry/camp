#include <Rcpp.h>
using namespace Rcpp;

#ifndef COMPUTE_EMPIRICAL_P_VALUE
#define COMPUTE_EMPIRICAL_P_VALUE
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side);
#endif

#ifndef FIT_SKEW_NORMAL_FUNCT
#define FIT_SKEW_NORMAL_FUNCT
std::vector<double> fit_skew_normal_funct(const std::vector<double>& y);
#endif

#ifndef FIT_AND_EVALUATE_SKEW_NORMAL
#define FIT_AND_EVALUATE_SKEW_NORMAL
std::vector<double> fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics, int side_code);
#endif
