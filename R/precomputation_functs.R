#' Perform score stat precomputation
run_score_stat_precomputation <- function(fitted_glm) {
  # extract pieces from GLM: family object, mu, y, and covariate matrix
  family_object <- fitted_glm$family
  fam_str <- gsub(x = family_object$family, replacement = "", pattern = "\\([0-9]+\\)")
  mu <- fitted_glm$fitted.values
  y <- fitted_glm$y
  covariate_matrix <- stats::model.matrix(fitted_glm$formula, data = fitted_glm$model)

  # special case for NB reg and Poisson reg with log link (for our genomics peeps!)
  if (fam_str == "Negative Binomial" && family_object$link == "log") {
    theta <- get(x = ".Theta", envir = environment(family_object$variance))
    denom <- 1 + mu/theta
    w <- mu/denom
    a <- (y - mu)/denom
  } else if (fam_str == "poisson" && family_object$link == "log") {
    w <- mu
    a <- y - mu
  } else {
    w <- fitted_glm$weights
    eta <- fitted_glm$linear.predictors
    M <- 1/family_object$mu.eta(eta)
    a <- M * w * (y - mu)
  }

  # compute the D matrix
  wZ <- w * covariate_matrix
  Zt_wZ <- t(covariate_matrix) %*% wZ
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  precomp_list <- list(a = a, w = w, D = D)
  return(precomp_list)
}


run_resid_precomputation <- function(fit, type = "deviance") {
 resids <- stats::residuals(fit, type = type)
 return(resids)
}


permute_bernoulli_treatment_vector <- function(x, B = 5498L) {
  if (!all(x %in% c(0, 1))) stop("x is not a bernoulli vector.")
  n <- length(x)
  trt_idxs <- which(x == 1) - 1L
  s <- length(trt_idxs)
  synthetic_idxs <- fisher_yates_samlper(n, s, B)
  return(list(synthetic_idxs = synthetic_idxs, trt_idxs = trt_idxs, B = B, s = s))
}
