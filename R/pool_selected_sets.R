#' Pool selected sets from multiply-imputed data
#'
#' Pool the selected sets from multiply-imputed or bootstrap + imputed data. Uses
#' the "stability" of the variables over the multiple selected sets to select
#' variables that are stable across the sets, where stability is determined by
#' presence in a certain fraction of the selected sets (and the fraction must be
#' above the specified threshold to be "stable").
#'
#' @param sets a list of sets of selected variables from the multiply-imputed datasets.
#'     Expects each set of selected variables to be a binary vector, where 1 denotes
#'     that the variable was selected.
#' @param threshold a numeric threshold between 0 and 1 detemining the "stability"
#'     of a feature; only features with stability above the threshold after pooling
#'     will be in the final selected set of variables.
#'
#' @return a vector denoting the final set of selected variables (1 denotes
#'     selected, 0 denotes not selected)
#' @export
pool_selected_sets <- function(sets = list(), threshold = 0.8) {
  set_tib <- do.call(rbind, sets)
  fractions <- colMeans(set_tib)
  fractions >= threshold
}
