#' Extract a Variance-Covariance Matrix for SPVIM Estimates
#' 
#' Extract a variance-covariance matrix based on the efficient influence function
#' for each of the estimated SPVIMs.
#' 
#' @param spvim_ests estimated SPVIMs
#' 
#' @return a variance-covariance matrix
#' @export
spvim_vcov <- function(spvim_ests) {
  if (is.null(names(spvim_ests$ic))) {
    # cross-fitted SEs were used; create a vector where each column is when
    # that "observation" was in the validation fold
    ic_v <- do.call(cbind, lapply(spvim_ests$ic, function(x) x$contrib_v[-1, ]))
    ic_s <- do.call(cbind, lapply(spvim_ests$ic, function(x) x$contrib_s[-1, ]))
    var_contrib_v <- ic_v %*% t(ic_v)
    var_contrib_s <- (1 / spvim_ests$gamma) * ic_s %*% t(ic_s)
  } else {
    var_contrib_v <- spvim_ests$ic$contrib_v[-1, ] %*%
      t(spvim_ests$ic$contrib_v[-1, ])
    var_contrib_s <- (1 / spvim_ests$gamma) * spvim_ests$ic$contrib_s[-1, ] %*%
      t(spvim_ests$ic$contrib_s[-1, ]) 
  }
  cov_mat <- var_contrib_v + var_contrib_s
  return(cov_mat)
}