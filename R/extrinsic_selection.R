#' Perform extrinsic, ensemble-based variable selection
#'
#' Based on a fitted Super Learner ensemble, extract extrinsic
#' variable importance estimates, rank them, and do variable
#' selection using the specified rank threshold.
#'
#' @param fit the fitted Super Learner ensemble.
#' @param feature_names the names of the features (a character vector of 
#'    length \code{p} (the total number of features)); only used if the 
#'    fitted Super Learner ensemble was fit on a \code{matrix} rather than on a 
#'    \code{data.frame}, \code{tibble}, etc.
#' @param threshold the threshold for selection based on ranked
#'    variable importance; rank 1 is the most important. Defaults
#'    to 20 (though this is arbitrary, and really should be
#'    specified for the task at hand).
#' @param import_type the type of extrinsic importance (either \code{"all"}, 
#'   the default, for a weighted combination of the individual-algorithm importance;
#'   or \code{"best"}, for the importance from the algorithm with the highest
#'   weight in the Super Learner).
#'
#' @return a tibble with the estimated extrinsic variable importance,
#'    the corresponding variable importance ranks, and the selected
#'    variables.
#' @seealso \code{\link[SuperLearner]{SuperLearner}} for specific usage of
#'    the \code{SuperLearner} function and package.
#' @importFrom magrittr `%>%`
#' @export
extrinsic_selection <- function(fit, feature_names, threshold = 20,
                                import_type = "all") {
  if (!any(grepl("SuperLearner", class(fit)))) {
    stop("The entered fitted ensemble must be a `SuperLearner` object.")
  } 

  var_nms <- fit$varNames
  if (is.null(var_nms)) {
    var_nms <- feature_names
  }
  importance_df <- extract_importance_SL(fit = fit, 
                                         feature_names = feature_names,
                                         import_type = import_type) %>%
    mutate(selected = rank < threshold)
  importance_df
}
