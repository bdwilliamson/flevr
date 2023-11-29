#' Extract the learner-specific importance from an svm object
#'
#' Extract the individual-algorithm extrinsic importance from a glm object,
#' along with the importance rank.
#'
#' @param fit the \code{svm} object.
#' @param x the features
#' @param y the outcome
#' @inheritParams extract_importance_glm
#' 
#' @examples
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- as.data.frame(dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))])
#' x_mat <- as.matrix(x)
#' feature_nms <- names(x)
#' # get the fit 
#' set.seed(20231129)
#' fit <- kernlab::ksvm(x_mat, y)
#' # extract importance
#' importance <- extract_importance_svm(fit = fit, feature_names = feature_nms, x = x, y = y)
#' importance
#'
#' @inherit extract_importance_glm return
#' @importFrom kernlab kpar kernelf param
#' @export
extract_importance_svm <- function(fit = NULL, feature_names = "", coef = 0, x = NULL, y = NULL) {
  if (!inherits(fit, "ksvm")) {
    stop("This is not an svm object. Please use a different importance extraction function.")
  } else {
    param_df <- data.frame(sigma = kpar(kernelf(fit))$sigma, C = param(fit)$C)
    if (inherits(kernelf(fit), "rbfkernel")) {
      caret_mod <- "svmRadial"
    } else {
      stop("The entered kernel is not currently supported.")
    }
    svm_train <- caret::train(x = x, y = y,
                              method = caret_mod, tuneGrid = param_df)
    svm_imp <- caret::varImp(svm_train, useModel = FALSE, value = "gcv")
    imp_dt <- tibble::tibble(algorithm = "svm", feature = rownames(svm_imp$importance),
                             importance = svm_imp$importance$Overall, rank = rank(-abs(svm_imp$importance)),
                             weight = coef)
    imp_dt[order(imp_dt$rank), ]
  }
}
