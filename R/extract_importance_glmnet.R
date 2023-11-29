#' Extract the learner-specific importance from a glmnet object
#'
#' Extract the individual-algorithm extrinsic importance from a glmnet object,
#' along with the importance rank.
#'
#' @inheritParams extract_importance_glm
#' @param fit the \code{glmnet} or \code{cv.glmnet} object
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
#' # get the fit (using only 3 CV folds for illustration only)
#' set.seed(20231129)
#' fit <- glmnet::cv.glmnet(x = as.matrix(x), y = y, 
#'                          family = "binomial", nfolds = 3)
#' # extract importance
#' importance <- extract_importance_glmnet(fit = fit, feature_names = feature_nms)
#' importance
#' 
#' @export
extract_importance_glmnet <- function(fit = NULL, feature_names = "", coef = 0) {
  if (!inherits(fit, "cv.glmnet") & !inherits(fit, "glmnet")) {
    stop("This is not a cv.glmnet object. Please use a different importance extraction function.")
  } else {
    p <- length(feature_names)
    imp_dt_init <- fit$glmnet.fit$beta[, which(fit$lambda == fit$lambda.min)]
    ranks <- rank(-abs(imp_dt_init))
    if (length(ranks) < p) {
      current_length <- length(imp_dt_init)
      current_nms <- names(imp_dt_init)
      avg_remaining_rank <- mean((length(imp_dt_init) + 1):p)
      remaining_features <- feature_names[!(feature_names %in% current_nms)]
      imp_dt_init <- c(imp_dt_init, rep(NA, p - current_length))
      names(imp_dt_init) <- c(current_nms, remaining_features)
      ranks <- c(ranks, rep(avg_remaining_rank, p - current_length))
    }
    imp_dt <- tibble::tibble(algorithm = "lasso", feature = names(imp_dt_init),
                             importance = abs(imp_dt_init), rank = ranks,
                             weight = coef)
    row.names(imp_dt) <- NULL
    imp_dt[order(imp_dt$rank), ]
  }
}
