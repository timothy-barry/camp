#' @export
run_perm_test_score_stat_binary_trt <- function(permutations, precomputation, side = 0L,
                                                return_null_distribution = FALSE, fit_sn = TRUE,
                                                p_thresh = 0.1, B_0 = 499L) {
  a <- precomputation$a
  w <- precomputation$w
  D <- precomputation$D
  synthetic_idxs <- permutations$synthetic_idxs
  trt_idxs <- permutations$trt_idxs
  B <- permutations$B
  s <- permutations$s

  out <- run_test_score_stat_binary_cpp(trt_idxs = trt_idxs,
                                        synthetic_idxs = synthetic_idxs,
                                        a = a,
                                        w = w,
                                        D = D,
                                        s = s,
                                        B = B,
                                        B_0 = B_0,
                                        fit_sn = fit_sn,
                                        side = side,
                                        p_thresh = p_thresh,
                                        return_null_distribution = return_null_distribution)
  return(out)
}

#' @export
run_perm_test_resid_stat_binary_trt <- function(permutations, precomputation, side = 0L,
                                                return_null_distribution = FALSE, fit_sn = TRUE,
                                                p_thresh = 0.1, B_0 = 499L) {
  resids <- precomputation
  synthetic_idxs <- permutations$synthetic_idxs
  trt_idxs <- permutations$trt_idxs
  B <- permutations$B
  s <- permutations$s
  out <- run_test_resid_stat_binary_cpp(trt_idxs = trt_idxs,
                                        synthetic_idxs = synthetic_idxs,
                                        resids = resids,
                                        s = s,
                                        B = B,
                                        side = side,
                                        B_0 = B_0,
                                        p_thresh = p_thresh,
                                        fit_sn = fit_sn,
                                        return_null_distribution = return_null_distribution)
  return(out)
}
