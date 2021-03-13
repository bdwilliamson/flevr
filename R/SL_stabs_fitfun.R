#' Wrapper for using Super Learner-based extrinsic selection within stability selection
#' 
#' A wrapper function for Super Learner-based extrinsic variable selection within
#' stability selection, using the \code{stabs} package.
#' 
#' @param x the features.
#' @param y the outcome of interest.
#' @param q the number of features to select on average.
#' @param ... other arguments to pass to \code{SuperLearner}.
#' 
#' @return a named list, with elements: \code{selected} (a logical vector 
#'   indicating whether or not each variable was selected); and \code{path} (
#'   a logical matrix indicating which variable was selected at each step).
#'   
#' @seealso \code{\link[stabs]{stabsel}} for general usage of stability selection.
#' @export
SL_stabs_fitfun <- function(x, y, q, ...) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package ", sQuote("SuperLearner"), " needed but not available")
  }
  x_df <- as.data.frame(x)
  if (is.null(names(x))) {
    names(x_df) <- paste0("V", 1:ncol(x_df))  
  }
  fit <- SuperLearner::SuperLearner(Y = y, X = x_df, verbose = FALSE, ...)
  selected <- extrinsic_selection(fit = fit, feature_names = names(x_df), 
                                  threshold = q + 1, import_type = "all")
  list(selected = selected$selected, path = NULL)
}