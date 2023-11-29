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
#' @param ... other arguments to pass to algorithm-specific importance extractors.
#'
#' @return a tibble with the estimated extrinsic variable importance,
#'    the corresponding variable importance ranks, and the selected
#'    variables.
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
#' library("SuperLearner")
#' set.seed(20231129)
#' fit <- SuperLearner::SuperLearner(Y = y, X = x, SL.library = c("SL.glm", "SL.mean"), 
#'                                   cvControl = list(V = 2))
#' # extract importance
#' importance <- extrinsic_selection(fit = fit, feature_names = feature_nms, threshold = 1.5, 
#'                                   import_type = "all")
#' importance
#' 
#' @seealso \code{\link[SuperLearner]{SuperLearner}} for specific usage of
#'    the \code{SuperLearner} function and package.
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate
#' @export
extrinsic_selection <- function(fit = NULL, feature_names = "", threshold = 20,
                                import_type = "all", ...) {
  if (!any(grepl("SuperLearner", class(fit)))) {
    stop("The entered fitted ensemble must be a `SuperLearner` object.")
  }

  var_nms <- fit$varNames
  if (is.null(var_nms)) {
    var_nms <- feature_names
  }
  importance_df <- extract_importance_SL(fit = fit,
                                         feature_names = feature_names,
                                         import_type = import_type, ...) %>%
    dplyr::mutate(selected = rank < threshold)
  importance_df
}
