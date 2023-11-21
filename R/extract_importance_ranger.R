#' Extract the learner-specific importance from a ranger object
#'
#' Extract the individual-algorithm extrinsic importance from a ranger object,
#' along with the importance rank.
#'
#' @param fit the \code{ranger} object.
#' @param feature_names the feature names
#' @param coef the Super Learner coefficient associated with the learner.
#'
#' @return a tibble, with columns \code{algorithm} (the fitted algorithm),
#'   \code{feature} (the feature), \code{importance} (the algorithm-specific
#'   extrinsic importance of the feature), \code{rank} (the feature importance
#'   rank, with 1 indicating the most important feature), and \code{weight}
#'   (the algorithm's weight in the Super Learner)
#' @export
extract_importance_ranger <- function(fit, feature_names, coef = 0) {
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
