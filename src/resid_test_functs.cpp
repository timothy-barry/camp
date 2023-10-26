#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
#include "shared_low_level_functions.h"


double compute_observed_resid_statistic(const NumericVector& resids, int s, IntegerVector trt_idxs) {
  double sum = 0, stat;
  for (int i = 0; i < s; i ++) {
    sum += resids[trt_idxs[i]];
  }
  stat = 1/(sqrt(s)) * sum;
  return(stat);
}


std::vector<double> compute_null_resid_statistics(const NumericVector& resids, int start_pos, int n_stats_to_compute, int s, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);
  
  // initialize variables
  int* curr_vect;
  std::vector<double> out(n_stats_to_compute);
  double sum, stat;
  
  for (int k = start_pos; k < n_stats_to_compute + start_pos; k ++) {
    sum = 0;
    curr_vect = &(*synth_idx_list)[k][0];
    for (int j = 0; j < s; j ++) {
      sum += resids[curr_vect[j]];
    }
    stat = 1/(sqrt(s)) * sum;
    out[k - start_pos] = stat;
  }
  return(out);
}


// [[Rcpp::export]]
SEXP run_test_resid_stat_binary_cpp(IntegerVector trt_idxs,
                                    SEXP synthetic_idxs,
                                    NumericVector resids,
                                    int s,
                                    int B,
                                    int side,
                                    bool adaptive,
                                    int B_0,
                                    double p_thresh,
                                    bool return_null_distribution) {
  // initialize variables
  double p, p_0;
  std::vector<double> null_statistics;
  List out;
  
  // compute the original statistic
  double z_orig = compute_observed_resid_statistic(resids, s, trt_idxs);
  
  // compute the null statistics
  if (!adaptive) {
    null_statistics = compute_null_resid_statistics(resids, 0, B, s, synthetic_idxs);
    p = compute_empirical_p_value(null_statistics, z_orig, side);
  } else {
    null_statistics = compute_null_resid_statistics(resids, 0, B_0, s, synthetic_idxs);
    p_0 = compute_empirical_p_value(null_statistics, z_orig, side);
    if (p_0 < p_thresh) {
      null_statistics = compute_null_resid_statistics(resids, B_0, B - B_0, s, synthetic_idxs);
      p = compute_empirical_p_value(null_statistics, z_orig, side);
    } else {
      p = p_0;
    }
  }
  
  if (return_null_distribution) {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("null_statistics") = null_statistics);
  } else {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig);
  }
  return(out);
}