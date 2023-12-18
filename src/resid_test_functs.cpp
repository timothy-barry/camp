#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
#include "shared_low_level_functions.h"


double compute_observed_resid_statistic(const NumericVector& resids, int s, IntegerVector trt_idxs) {
  double tot_sum = 0, trt_sum = 0, control_sum = 0, stat = 0;
  int n_control = resids.size() - s;
  for (int i = 0; i < resids.size(); i ++) tot_sum += resids[i];
  for (int i = 0; i < s; i ++) trt_sum += resids[trt_idxs[i]];
  control_sum = tot_sum - trt_sum;
  stat = trt_sum/s - control_sum/n_control;
  return(stat);
}


std::vector<double> compute_null_resid_statistics(const NumericVector& resids, int start_pos, int n_stats_to_compute, int s, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);

  // initialize variables
  int* curr_vect;
  std::vector<double> out(n_stats_to_compute);
  double tot_sum = 0, trt_sum = 0, control_sum = 0, stat = 0;

  // compute control_sum and n_control
  int n_control = resids.size() - s;
  for (int i = 0; i < resids.size(); i ++) tot_sum += resids[i];

  for (int k = start_pos; k < n_stats_to_compute + start_pos; k ++) {
    trt_sum = 0;
    curr_vect = &(*synth_idx_list)[k][0];
    for (int j = 0; j < s; j ++) {
      trt_sum += resids[curr_vect[j]];
    }
    control_sum = tot_sum - trt_sum;
    stat = trt_sum/s - control_sum/n_control;
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
                                    int B_0,
                                    double p_thresh,
                                    bool return_null_distribution,
                                    bool fit_sn) {
  // initialize variables
  double p, p_0;
  std::vector<double> null_statistics, fit_sn_out;
  List out;

  // compute the original statistic
  double z_orig = compute_observed_resid_statistic(resids, s, trt_idxs);

  // compute the null statistics
  null_statistics = compute_null_resid_statistics(resids, 0, B_0, s, synthetic_idxs);
  p_0 = compute_empirical_p_value(null_statistics, z_orig, side);
  if (p_0 < p_thresh) {
    null_statistics = compute_null_resid_statistics(resids, B_0, B - B_0, s, synthetic_idxs);
    if (fit_sn) {
      fit_sn_out = fit_and_evaluate_skew_normal(z_orig, null_statistics, side);
      p = fit_sn_out[3];
    } else {
      p = compute_empirical_p_value(null_statistics, z_orig, side);
    }
  } else {
    p = p_0;
  }

  if (return_null_distribution) {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("null_statistics") = null_statistics);
  } else {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig);
  }
  return(out);
}
