#' Extract the learner-specific importance from a polymars object
#'
#' Extract the individual-algorithm extrinsic importance from a polymars object,
#' along with the importance rank.
#'
#' @param fit the \code{polymars} object.
#' @param feature_names the feature names
#' @param coef the Super Learner coefficient associated with the learner.
#'
#' @return a tibble, with columns \code{algorithm} (the fitted algorithm),
#'   \code{feature} (the feature), \code{importance} (the algorithm-specific
#'   extrinsic importance of the feature), \code{rank} (the feature importance
#'   rank, with 1 indicating the most important feature), and \code{weight}
#'   (the algorithm's weight in the Super Learner)
#' @export
extract_importance_polymars <- function(fit, feature_names, coef = 0) {
  if (!inherits(fit, "polymars")) {
    stop("This is not a polymars object. Please use a different importance extraction function.")
  } else {
    p <- length(feature_names)
    mod <- fit$model
    all_preds <- c(mod$pred1, mod$pred2)
    nonzero_preds <- all_preds[all_preds > 0]
    unique_preds <- unique(nonzero_preds)
    coeffs <- matrix(nrow = length(unique_preds), ncol = 2)
    for (i in seq_len(nrow(coeffs))) {
      these_coefs <- mod[mod$pred1 == unique_preds[i] | mod$pred2 == unique_preds[i], ]
      coeffs[i, ] <- c(unique_preds[i], sum(abs(these_coefs$coefs / these_coefs$SE)))
    }
    summ2 <- as.data.frame(coeffs[rank(-abs(coeffs[, 2]), ties.method = "last"), ])
    names(summ2) <- c("feature", "importance")
    summ2$rank <- rank(-abs(summ2$importance), ties.method = "last")
    if (nrow(summ2) < p) {
      current_length <- nrow(summ2)
      current_nms <- row.names(summ2)
      avg_remaining_rank <- mean((current_length + 1):p)
      remaining_features <- feature_names[!(feature_names %in% current_nms)]
      na_mat <- matrix(NA, nrow = p - nrow(summ2), ncol = ncol(summ2))
      na_df <- as.data.frame(na_mat)
      names(na_df) <- names(summ2)
      na_df$feature <- remaining_features
      na_df$rank <- avg_remaining_rank
      summ2 <- dplyr::bind_rows(summ2, na_df)
    }
    imp_dt  <- tibble::tibble(algorithm = "polymars", feature = feature_names[summ2$feature],
                              importance = summ2$importance,
                              rank = summ2$rank,
                              weight = coef)
    imp_dt[order(imp_dt$rank), ]
  }
}
