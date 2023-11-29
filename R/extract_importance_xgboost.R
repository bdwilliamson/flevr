#' Extract the learner-specific importance from an xgboost object
#'
#' Extract the individual-algorithm extrinsic importance from an xgboost object,
#' along with the importance rank.
#'
#' @param fit the \code{xgboost} object.
#' @inheritParams extract_importance_glm
#' 
#' @examples
#' \donttest{
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- as.matrix(dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))])
#' feature_nms <- names(x)
#' set.seed(20231129)
#' xgbmat <- xgboost::xgb.DMatrix(data = x, label = y)
#' # get the fit, using a small number of rounds for illustration only 
#' fit <- xgboost::xgboost(data = xgbmat, objective = "binary:logistic", nthread = 1, nrounds = 10)
#' # extract importance
#' importance <- extract_importance_xgboost(fit = fit, feature_names = feature_nms)
#' importance
#' }
#'
#' @inherit extract_importance_glm return
#' @export
extract_importance_xgboost <- function(fit = NULL, feature_names = "", coef = 0) {
  if (!inherits(fit, "xgboost") & !inherits(fit, "xgb.Booster")) {
    stop("This is not an xgboost object. Please use a different importance extraction function.")
  } else {
    p <- length(feature_names)
    imp_dt_init <- xgboost::xgb.importance(model = fit)
    colnames(imp_dt_init)[1] <- "Feature"
    imp_dt_init$rank <- seq_len(nrow(imp_dt_init))
    if (nrow(imp_dt_init) < p) {
      avg_remaining_rank <- mean((nrow(imp_dt_init) + 1):p)
      na_mat <- matrix(NA, nrow = p - nrow(imp_dt_init), ncol = ncol(imp_dt_init))
      na_df <- as.data.frame(na_mat)
      names(na_df) <- names(imp_dt_init)
      na_df$Feature <- feature_names[!(feature_names %in% imp_dt_init$Feature)]
      na_df$rank <- avg_remaining_rank
      imp_dt_init <- dplyr::bind_rows(imp_dt_init, na_df)
    }
    imp_dt <- tibble::tibble(algorithm = "xgboost", feature = imp_dt_init$Feature,
                             importance = imp_dt_init$Gain, rank = imp_dt_init$rank,
                             weight = coef)
    row.names(imp_dt) <- NULL
    imp_dt[order(imp_dt$rank), ]
  }
}
