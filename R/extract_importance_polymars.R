#' Extract the learner-specific importance from a polymars object
#'
#' Extract the individual-algorithm extrinsic importance from a polymars object,
#' along with the importance rank.
#'
#' @param fit the \code{polymars} object.
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
#' x_mat <- as.matrix(x)
#' # get the fit
#' set.seed(20231129)
#' fit <- polspline::polyclass(y, x_mat)
#' # extract importance
#' importance <- extract_importance_polymars(fit = fit, feature_names = feature_nms)
#' importance
#' 
#' @export
extract_importance_polymars <- function(fit = NULL, feature_names = "", coef = 0) {
  if (!inherits(fit, "polymars") & !inherits(fit, "polyclass")) {
    stop("This is not a polymars object. Please use a different importance extraction function.")
  } else {
    p <- length(feature_names)
    if (inherits(fit, "polymars")) {
      mod <- fit$model
      all_preds <- c(mod$pred1, mod$pred2)
      nonzero_preds <- all_preds[all_preds > 0]
      unique_preds <- unique(nonzero_preds)
      coeffs <- matrix(nrow = length(unique_preds), ncol = 2)
      for (i in seq_len(nrow(coeffs))) {
        these_coefs <- mod[mod$pred1 == unique_preds[i] | mod$pred2 == unique_preds[i], ]
        coeffs[i, ] <- c(unique_preds[i], sum(abs(these_coefs$coefs / these_coefs$SE)))
      }
    } else {
      mod <- fit$fcts
      all_preds <- as.numeric(c(mod[, 1], mod[, 3]))
      nonzero_preds <- all_preds[all_preds > 0 & !is.na(all_preds)]
      unique_preds <- sort(unique(nonzero_preds))
      coeffs <- matrix(nrow = length(unique_preds), ncol = 2)
      for (i in seq_len(nrow(coeffs))) {
        these_coefs <- mod[(mod[, 1] == unique_preds[i] & !is.na(mod[, 1])) |
                             (mod[, 3] == unique_preds[i] & !is.na(mod[, 3])), , drop = FALSE]
        coeffs[i, ] <- c(unique_preds[i], ifelse(
          all(these_coefs[, 6] == 0), sum(abs(these_coefs[, 5])),
          sum(abs(these_coefs[, 5] / these_coefs[, 6]))
        ))
      }
    }
    summ2 <- as.data.frame(coeffs[rank(-abs(coeffs[, 2]), ties.method = "last"), ])
    names(summ2) <- c("feature", "importance")
    summ2$rank <- rank(-abs(summ2$importance), ties.method = "last")
    if (nrow(summ2) < p) {
      current_length <- nrow(summ2)
      current_nms <- summ2$feature
      avg_remaining_rank <- mean((current_length + 1):p)
      remaining_features <- (1:p)[!(1:p %in% current_nms)]
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
