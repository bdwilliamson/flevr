#' Extract the learner-specific importance from a fitted SuperLearner algorithm
#'
#' Extract the individual-algorithm extrinsic importance from one fitted
#' algorithm within the Super Learner, along with the importance rank.
#'
#' @param fit the specific learner (e.g., from the Super Learner's
#'     \code{fitLibrary} list).
#' @param coef the Super Learner coefficient associated with the learner.
#' @param feature_names the feature names.
#' @param ... other arguments to pass to algorithm-specific importance extractors.
#'
#' @return a tibble, with columns \code{algorithm} (the fitted algorithm),
#'   \code{feature} (the feature), \code{importance} (the algorithm-specific
#'   extrinsic importance of the feature), \code{rank} (the feature importance
#'   rank, with 1 indicating the most important feature), and \code{weight}
#'   (the algorithm's weight in the Super Learner)
#' @export
extract_importance_SL_learner <- function(fit, coef, feature_names, ...) {
  if (inherits(fit, "xgboost") | inherits(fit, "xgb.Booster")) {
    imp_dt <- extract_importance_xgboost(fit = fit, feature_names = feature_names,
                                         coef = coef)
  } else if (inherits(fit, "ranger")) {
    imp_dt <- extract_importance_ranger(fit = fit, feature_names = feature_names,
                                        coef = coef)
  } else if (inherits(fit, "glmnet") | inherits(fit, "cv.glmnet")) {
    imp_dt <- extract_importance_glmnet(fit = fit, feature_names = feature_names,
                                        coef = coef)
  } else if (inherits(fit, "glm")) {
    imp_dt <- extract_importance_glm(fit = fit, feature_names = feature_names,
                                     coef = coef)
  } else if (inherits(fit, "numeric")) {
    imp_dt <- extract_importance_mean(fit = fit, feature_names = feature_names,
                                      coef = coef)
  } else if (inherits(fit, "ksvm")) {
    L <- list(...)
    imp_dt <- extract_importance_svm(fit = fit, feature_names = feature_names,
                                     coef = coef, x = L$x, y = L$y)
  } else if (inherits(fit, "polymars") | inherits(fit, "polyclass")) {
    imp_dt <- extract_importance_polymars(fit = fit, feature_names = feature_names,
                                          coef = coef)
  } else {
    stop("One of the algorithms in fitLibrary is unsupported at this time.")
  }
  imp_dt
}
