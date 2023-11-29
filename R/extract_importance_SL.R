#' Extract extrinsic importance from a Super Learner object
#'
#' Extract the individual-algorithm extrinsic importance from each fitted
#' algorithm within the Super Learner; compute the average weighted rank of the
#' importance scores, with weights specified by each algorithm's weight in the
#' Super Learner.
#'
#' @param fit the fitted Super Learner ensemble
#' @param feature_names the names of the features
#' @param import_type the level of granularity for importance: \code{"all"} is the
#'     importance based on the weighted average of ranks across algorithmrithms
#'     (weights are SL coefs); \code{"best"} is the importance based on the algorithmrithm
#'     with highest weight. Defaults to \code{"all"}.
#' @param ... other arguments to pass to individual-algorithm extractors.
#'
#' @return a tibble, with columns \code{feature} (the feature) and
#'   \code{rank} (the weighted feature importance rank, with 1 indicating the
#'   most important feature).
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
#' set.seed(20231129)
#' library("SuperLearner")
#' fit <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = c("SL.glm", "SL.mean"), 
#'                                   cvControl = list(V = 2))
#' # extract importance using all learners
#' importance <- extract_importance_SL(fit = fit, feature_names = feature_nms)
#' importance
#' # extract importance of best learner
#' best_importance <- extract_importance_SL(fit = fit, feature_names = feature_nms, 
#'                                          import_type = "best")
#' best_importance
#'
#' @importFrom dplyr case_when filter group_by mutate summarize
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @export
extract_importance_SL <- function(fit, feature_names, import_type = "all", ...) {
  if (import_type == "all") {
    nonzero_weights <- (1:length(fit$coef))[fit$coef > 0]
    importances <- sapply(
      nonzero_weights,
      function(i) extract_importance_SL_learner(
        fit = fit$fitLibrary[[i]]$object, coef = fit$coef[i],
        feature_names = feature_names, ...
      ),
      simplify = FALSE
    )
    imp_dt <- do.call(rbind, importances) %>%
      dplyr::filter(.data$weight != 0) %>%
      dplyr::group_by(.data$feature, .data$algorithm) %>%
      dplyr::mutate(wgt_rank = .data$weight * .data$rank, .groups = "drop") %>%
      dplyr::group_by(.data$feature) %>%
      dplyr::summarize(rank = sum(.data$wgt_rank), .groups = "drop")
  } else {
    biggest_weight <- which.max(fit$coef)
    best_obj <- fit$fitLibrary[[biggest_weight]]$object
    imp_dt <- extract_importance_SL_learner(
      fit = best_obj, coef = fit$coef[biggest_weight],
      feature_names = feature_names, ...
    ) %>%
      dplyr::select("feature", "rank")
  }
  imp_dt[order(imp_dt$rank), ]
}
