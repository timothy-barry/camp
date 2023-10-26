run_perm_test_score_stat_binary_trt <- function(x, precomputation, permutations, return_null_distribution = FALSE,
                                                adaptive = TRUE, p_0 = 0.1, B_0 = 500L) {
  a <- precomputation$a
  w <- precomputation$w
  D <- precomputation$D
  out <- run_perm_test_score_stat_binary_cpp(x = x,
                                             permutations = permutations,
                                             a = a,
                                             w = w,
                                             D = D,
                                             adaptive = adaptive,
                                             p_0 = p_0,
                                             return_null_distribution = return_null_distribution)
  return(out)
}
