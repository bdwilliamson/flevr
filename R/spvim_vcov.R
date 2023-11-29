#' Extract a Variance-Covariance Matrix for SPVIM Estimates
#' 
#' Extract a variance-covariance matrix based on the efficient influence function
#' for each of the estimated SPVIMs.
#' 
#' @param spvim_ests estimated SPVIMs
#' 
#' @return a variance-covariance matrix
#' @examples
#' \donttest{
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))]
#' feature_nms <- names(x)
#' # estimate SPVIMs (using simple library and V = 2 for illustration only)
#' set.seed(20231129)
#' library("SuperLearner")
#' est <- vimp::sp_vim(Y = y, X = x, V = 2, type = "auc", SL.library = "SL.glm", 
#'                     cvControl = list(V = 2))
#' # get variance-covariance matrix
#' vcov <- spvim_vcov(spvim_ests = est)
#' }
#' @export
spvim_vcov <- function(spvim_ests = NULL) {
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