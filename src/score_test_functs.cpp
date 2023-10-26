#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>
#include "shared_low_level_functions.h"


double compute_observed_full_statistic(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int s, const IntegerVector& trt_idxs) {
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), D_ncol = D.ncol();

  // iterate over the rows of D
  for (int i = 0; i < D_nrow; i ++) {
    inner_sum = 0;
    for (int j = 0; j < s; j ++) {
      inner_sum += D(i, trt_idxs[j]);
    }
    lower_right += inner_sum * inner_sum;
  }

  // second, compute the lower-left hand of the denominator; also, compute the top
  for (int j = 0; j < s; j ++) {
    top += a[trt_idxs[j]];
    lower_left += w[trt_idxs[j]];
  }

  // finally, compute the z-score
  return(top/sqrt(lower_left - lower_right));
}


std::vector<double> compute_null_full_statistics(const NumericVector& a, const NumericVector& w, const NumericMatrix& D,
                                                 int start_pos, int n_stats_to_compute, int s, SEXP synthetic_idxs) {
  // dereference the synthetic treatment idxs
  Rcpp::XPtr<std::vector<std::vector<int>>> synth_idx_list(synthetic_idxs);

  // initialize variables
  int* curr_vect;
  std::vector<double> out(n_stats_to_compute);
  double lower_right = 0, lower_left = 0, top = 0, inner_sum;
  int D_nrow = D.nrow(), D_ncol = D.ncol(), idx;

  // iterate over the n_stats_to_compute tests
  for (int k = start_pos; k < n_stats_to_compute + start_pos; k ++) {
    curr_vect = &(*synth_idx_list)[k][0];
    lower_right = 0;
    inner_sum = 0;

    // iterate over the rows of D
    for (int i = 0; i < D_nrow; i ++) {
      inner_sum = 0;
      for (int j = 0; j < s; j ++) {
        inner_sum += D(i, curr_vect[j]);
      }
      lower_right += inner_sum * inner_sum;
    }

    // second, compute the lower-left hand of the denominator; also, compute the top
    lower_left = 0;
    top = 0;
    for (int j = 0; j < s; j ++) {
      idx = curr_vect[j];
      top += a[idx];
      lower_left += w[idx];
    }

    // compute the z-score
    out[k - start_pos] = top/sqrt(lower_left - lower_right);
  }

  return(out);
}


// [[Rcpp::export]]
SEXP run_test_score_stat_binary_cpp(IntegerVector trt_idxs,
                                    SEXP synthetic_idxs,
                                    NumericVector a,
                                    NumericVector w,
                                    NumericMatrix D,
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
  double z_orig = compute_observed_full_statistic(a, w, D, s, trt_idxs);

  // compute the null statistics
  if (!adaptive) {
    null_statistics = compute_null_full_statistics(a, w, D, 0, B, s, synthetic_idxs);
    p = compute_empirical_p_value(null_statistics, z_orig, side);
  } else {
    null_statistics = compute_null_full_statistics(a, w, D, 0, B_0, s, synthetic_idxs);
    p_0 = compute_empirical_p_value(null_statistics, z_orig, side);
    if (p_0 < p_thresh) {
      null_statistics = compute_null_full_statistics(a, w, D, B_0, B - B_0, s, synthetic_idxs);
      p = compute_empirical_p_value(null_statistics, z_orig, side);
    } else {
      p = p_0;
    }
  }

  // construct output
  if (return_null_distribution) {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig, Named("null_statistics") = null_statistics);
  } else {
    out = List::create(Named("p") = p, Named("z_orig") = z_orig);
  }
  return(out);
}
