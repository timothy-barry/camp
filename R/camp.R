#' camp
#'
#' `camp` is an R package for powerful and fast permutation testing and conditional randomization testing.
#'
#' @useDynLib camp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import BH
#' @docType package
#' @name camp
#' @examples
#' # simulate data
#' n <- 5000
#' x <- rbinom(n = n, size = 1, prob = 0.1)
#' z <- MASS::mvrnorm(n = n, mu = c(-0.5, 0.5), Sigma = toeplitz(c(1, 0.5)))
#' family_object <- poisson()
#' design_matrix <- cbind(x, z)
#' y <- generate_glm_data(
#'   design_matrix = design_matrix,
#'   coefficients = c(0.6, 0.1, 0.1, 0.3),
#'   family_object = family_object,
#'   add_intercept = TRUE
#' )
#'
#' # 1. fit a GLM excluding the treatment
#' fit <- glm(
#'   y ~ design_matrix[,-x],
#'   family = family_object,
#' )
#'
#' # 2. obtain the permuted treatment vectors
#' permutations <- permute_bernoulli_treatment_vector(x)
#'
#' # 3. perform the score statistic precomputation
#' precomputation_score <- run_score_stat_precomputation(fit)
#'
#' # 4. run the permutation test
#' res_score <- run_perm_test_score_stat_binary_trt(permutations, precomputation_score)
#'
#' # 5. perform deviance resid precomputation
#' precomputation_resid <- run_resid_precomputation(fit)
#'
#' # 6. run the permutation test
#' res_resid <- run_perm_test_resid_stat_binary_trt(permutations, precomputation_resid)
#'
#' # 6. compare to standard GLM p-value
#' fit_full <- glm(
#'   y ~ design_matrix,
#'   family = family_object,
#' )
#' summary(fit_full)$coefficients[2,"Pr(>|z|)"]
NULL
