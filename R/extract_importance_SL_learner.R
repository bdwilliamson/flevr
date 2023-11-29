#' Extract the learner-specific importance from a fitted SuperLearner algorithm
#'
#' Extract the individual-algorithm extrinsic importance from one fitted
#' algorithm within the Super Learner, along with the importance rank.
#'
#' @param fit the specific learner (e.g., from the Super Learner's
#'     \code{fitLibrary} list).
#' @inheritParams extract_importance_glm
#' @param ... other arguments to pass to algorithm-specific importance extractors.
#'
#' @inherit extract_importance_glm return
#' 
#' @examples
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))]
#' feature_nms <- names(x)
#' # get the fit (using a simple library and 2 folds for illustration only)
#' library("SuperLearner")
#' set.seed(20231129)
#' fit <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = c("SL.glm", "SL.mean"), 
#'                                   cvControl = list(V = 2))
#' # extract importance
#' importance <- extract_importance_SL_learner(fit = fit$fitLibrary[[1]]$object, 
#'                                             feature_names = feature_nms, coef = fit$coef[1])
#' importance
#' 
#' @export
extract_importance_SL_learner <- function(fit = NULL, coef = 0, feature_names = "", ...) {
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
