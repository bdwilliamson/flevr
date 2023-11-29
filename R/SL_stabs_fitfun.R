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
#' @examples
#' \donttest{
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))]
#' feature_nms <- names(x)
#' # use stability selection with SL (using small number of folds for CV, 
#' # small SL library and small number of bootstrap replicates for illustration only)
#' set.seed(20231129)
#' library("SuperLearner")
#' sl_stabs <- stabs::stabsel(x = x, y = y,
#'                            fitfun = SL_stabs_fitfun,
#'                            args.fitfun = list(SL.library = "SL.glm", cvControl = list(V = 2)),
#'                            q = 2, B = 5, PFER = 5)
#' sl_stabs
#' }
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
                                  threshold = q + 1, import_type = "all",
                                  x = x, y = y)
  list(selected = selected$selected, path = NULL)
}
