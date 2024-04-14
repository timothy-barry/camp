// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fisher_yates_samlper
SEXP fisher_yates_samlper(int n, int s, int B);
RcppExport SEXP _camp_fisher_yates_samlper(SEXP nSEXP, SEXP sSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fisher_yates_samlper(n, s, B));
    return rcpp_result_gen;
END_RCPP
}
// run_test_resid_stat_binary_cpp
SEXP run_test_resid_stat_binary_cpp(IntegerVector trt_idxs, SEXP synthetic_idxs, NumericVector resids, int s, int B, int side, int B_0, double p_thresh, bool return_null_distribution, bool fit_sn);
RcppExport SEXP _camp_run_test_resid_stat_binary_cpp(SEXP trt_idxsSEXP, SEXP synthetic_idxsSEXP, SEXP residsSEXP, SEXP sSEXP, SEXP BSEXP, SEXP sideSEXP, SEXP B_0SEXP, SEXP p_threshSEXP, SEXP return_null_distributionSEXP, SEXP fit_snSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type trt_idxs(trt_idxsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type synthetic_idxs(synthetic_idxsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    Rcpp::traits::input_parameter< int >::type B_0(B_0SEXP);
    Rcpp::traits::input_parameter< double >::type p_thresh(p_threshSEXP);
    Rcpp::traits::input_parameter< bool >::type return_null_distribution(return_null_distributionSEXP);
    Rcpp::traits::input_parameter< bool >::type fit_sn(fit_snSEXP);
    rcpp_result_gen = Rcpp::wrap(run_test_resid_stat_binary_cpp(trt_idxs, synthetic_idxs, resids, s, B, side, B_0, p_thresh, return_null_distribution, fit_sn));
    return rcpp_result_gen;
END_RCPP
}
// compute_observed_full_statistic
double compute_observed_full_statistic(const NumericVector& a, const NumericVector& w, const NumericMatrix& D, int s, const IntegerVector& trt_idxs);
RcppExport SEXP _camp_compute_observed_full_statistic(SEXP aSEXP, SEXP wSEXP, SEXP DSEXP, SEXP sSEXP, SEXP trt_idxsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type trt_idxs(trt_idxsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_observed_full_statistic(a, w, D, s, trt_idxs));
    return rcpp_result_gen;
END_RCPP
}
// run_test_score_stat_binary_cpp
SEXP run_test_score_stat_binary_cpp(IntegerVector trt_idxs, SEXP synthetic_idxs, NumericVector a, NumericVector w, NumericMatrix D, int s, int B, int side, int B_0, double p_thresh, bool return_null_distribution, bool fit_sn);
RcppExport SEXP _camp_run_test_score_stat_binary_cpp(SEXP trt_idxsSEXP, SEXP synthetic_idxsSEXP, SEXP aSEXP, SEXP wSEXP, SEXP DSEXP, SEXP sSEXP, SEXP BSEXP, SEXP sideSEXP, SEXP B_0SEXP, SEXP p_threshSEXP, SEXP return_null_distributionSEXP, SEXP fit_snSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type trt_idxs(trt_idxsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type synthetic_idxs(synthetic_idxsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    Rcpp::traits::input_parameter< int >::type B_0(B_0SEXP);
    Rcpp::traits::input_parameter< double >::type p_thresh(p_threshSEXP);
    Rcpp::traits::input_parameter< bool >::type return_null_distribution(return_null_distributionSEXP);
    Rcpp::traits::input_parameter< bool >::type fit_sn(fit_snSEXP);
    rcpp_result_gen = Rcpp::wrap(run_test_score_stat_binary_cpp(trt_idxs, synthetic_idxs, a, w, D, s, B, side, B_0, p_thresh, return_null_distribution, fit_sn));
    return rcpp_result_gen;
END_RCPP
}
// compute_empirical_p_value
double compute_empirical_p_value(const std::vector<double>& null_statistics, double z_orig, int side);
RcppExport SEXP _camp_compute_empirical_p_value(SEXP null_statisticsSEXP, SEXP z_origSEXP, SEXP sideSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type null_statistics(null_statisticsSEXP);
    Rcpp::traits::input_parameter< double >::type z_orig(z_origSEXP);
    Rcpp::traits::input_parameter< int >::type side(sideSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_empirical_p_value(null_statistics, z_orig, side));
    return rcpp_result_gen;
END_RCPP
}
// fit_skew_normal_funct
std::vector<double> fit_skew_normal_funct(const std::vector<double>& y);
RcppExport SEXP _camp_fit_skew_normal_funct(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(fit_skew_normal_funct(y));
    return rcpp_result_gen;
END_RCPP
}
// fit_and_evaluate_skew_normal
std::vector<double> fit_and_evaluate_skew_normal(double z_orig, std::vector<double>& null_statistics, int side_code);
RcppExport SEXP _camp_fit_and_evaluate_skew_normal(SEXP z_origSEXP, SEXP null_statisticsSEXP, SEXP side_codeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z_orig(z_origSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type null_statistics(null_statisticsSEXP);
    Rcpp::traits::input_parameter< int >::type side_code(side_codeSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_and_evaluate_skew_normal(z_orig, null_statistics, side_code));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_camp_fisher_yates_samlper", (DL_FUNC) &_camp_fisher_yates_samlper, 3},
    {"_camp_run_test_resid_stat_binary_cpp", (DL_FUNC) &_camp_run_test_resid_stat_binary_cpp, 10},
    {"_camp_compute_observed_full_statistic", (DL_FUNC) &_camp_compute_observed_full_statistic, 5},
    {"_camp_run_test_score_stat_binary_cpp", (DL_FUNC) &_camp_run_test_score_stat_binary_cpp, 12},
    {"_camp_compute_empirical_p_value", (DL_FUNC) &_camp_compute_empirical_p_value, 3},
    {"_camp_fit_skew_normal_funct", (DL_FUNC) &_camp_fit_skew_normal_funct, 1},
    {"_camp_fit_and_evaluate_skew_normal", (DL_FUNC) &_camp_fit_and_evaluate_skew_normal, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_camp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
