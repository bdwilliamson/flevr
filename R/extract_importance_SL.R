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
#' 
#' @return a tibble, with columns \code{feature} (the feature) and 
#'   \code{rank} (the weighted feature importance rank, with 1 indicating the 
#'   most important feature).
#' 
#' @importFrom dplyr case_when filter group_by mutate summarize
#' @importFrom data.table rbindlist
extract_importance_SL <- function(fit, feature_names, import_type = "all") {
  if (import_type == "all") {
    nonzero_weights <- (1:length(fit$coef))[fit$coef > 0]
    importances <- sapply(
      nonzero_weights,
      function(i) extract_importance_SL_learner(
        fit = fit$fitLibrary[[i]]$object, coef = fit$coef[i], 
        feature_names = feature_names
      ),
      simplify = FALSE
    )
    imp_dt <- importances %>%
      data.table::rbindlist() %>%
      dplyr::filter(weight != 0) %>%
      dplyr::group_by(feature, algorithm) %>%
      dplyr::mutate(wgt_rank = weight * rank, .groups = "drop") %>%
      dplyr::group_by(feature) %>%
      dplyr::summarize(rank = sum(wgt_rank), .groups = "drop")
  } else {
    biggest_weight <- which.max(fit$coef)
    best_obj <- fit$fitLibrary[[biggest_weight]]$object
    imp_dt <- extract_importance_SL_learner(
      fit = best_obj, coef = fit$coef[biggest_weight], 
      feature_names = feature_names
    ) %>% 
      dplyr::select(feature, rank)
  } 
  imp_dt[order(imp_dt$rank), ]
}
