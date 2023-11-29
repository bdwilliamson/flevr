#' Extract the learner-specific importance from a ranger object
#'
#' Extract the individual-algorithm extrinsic importance from a ranger object,
#' along with the importance rank.
#'
#' @param fit the \code{ranger} object.
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
#' # get the fit
#' set.seed(20231129)
#' fit <- ranger::ranger(y ~ ., data = data.frame(y = y, x), importance = "impurity")
#' # extract importance
#' importance <- extract_importance_ranger(fit = fit, feature_names = feature_nms)
#' importance
#' 
#' @export
extract_importance_ranger <- function(fit = NULL, feature_names = "", coef = 0) {
  if (!inherits(fit, "ranger")) {
    stop("This is not a ranger object. Please use a different importance extraction function.")
  } else {
    p <- length(feature_names)
    imp_dt_init <- ranger::importance(fit)
    ranks <- rank(-abs(imp_dt_init))
    if (length(ranks) < p) {
      avg_remaining_rank <- mean((length(imp_dt_init) + 1):p)
      current_nms <- names(imp_dt_init)
      remaining_features <- feature_names[!(feature_names %in% current_nms)]
      current_length <- length(imp_dt_init)
      imp_dt_init <- c(imp_dt_init, rep(NA, p - current_length))
      names(imp_dt_init) <- c(current_nms, remaining_features)
      ranks <- c(ranks, rep(avg_remaining_rank, p - current_length))
    }
    imp_dt <- tibble::tibble(algorithm = "rf", feature = names(imp_dt_init),
                             importance = imp_dt_init, rank = ranks,
                             weight = coef)
    row.names(imp_dt) <- NULL
    imp_dt[order(imp_dt$rank), ]
  }
}
