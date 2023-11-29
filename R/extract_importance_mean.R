#' Extract the learner-specific importance from a mean object
#' 
#' Extract the individual-algorithm extrinsic importance from a mean object, 
#' along with the importance rank.
#' 
#' @param fit the \code{mean} object.
#' @inheritParams extract_importance_glm
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
#' # get the mean outcome
#' fit <- mean(y)
#' # extract importance
#' importance <- extract_importance_mean(fit = fit, feature_names = feature_nms)
#' importance
#' 
#' @export 
extract_importance_mean <- function(fit = NULL, feature_names = "", coef = 0) {
  imp_dt <- tibble::tibble(algorithm = "mean", feature = feature_names, 
                           importance = NA, 
                           rank = rep(mean(1:length(feature_names)), 
                                      length(feature_names)), 
                           weight = coef)
  imp_dt
}
