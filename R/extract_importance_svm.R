#' Extract the learner-specific importance from an svm object
#'
#' Extract the individual-algorithm extrinsic importance from a glm object,
#' along with the importance rank.
#'
#' @param fit the \code{svm} object.
#' @param feature_names the feature names
#' @param coef the Super Learner coefficient associated with the learner.
#' @param x the features
#' @param y the outcome
#'
#' @return a tibble, with columns \code{algorithm} (the fitted algorithm),
#'   \code{feature} (the feature), \code{importance} (the algorithm-specific
#'   extrinsic importance of the feature), \code{rank} (the feature importance
#'   rank, with 1 indicating the most important feature), and \code{weight}
#'   (the algorithm's weight in the Super Learner)
#' @export
extract_importance_svm <- function(fit, feature_names, coef = 0, x = NULL, y = NULL) {
  if (!any(grepl("svm", class(fit)))) {
    stop("This is not a glm object. Please use a different importance extraction function.")
  } else {
    param_df <- data.frame(sigma = kpar(kernelf(obj))$sigma, C = param(obj)$C)
    if (class(kernelf(obj)) == "rbfkernel") {
      caret_mod <- "svmRadial"
    } else {
      stop("The entered kernel is not currently supported.")
    }
    svm_train <- caret::train(x = x, y = y,
                              method = caret_mod, tuneGrid = param_df)
    svm_imp <- caret::varImp(svm_train, useModel = FALSE, value = "gcv")
    imp_dt <- tibble::tibble(algo = "svm", feature = rownames(svm_imp$importance),
                             importance = svm_imp$importance$Overall, rank = rank(-importance),
                             weight = coef)
  }
}
