#' Control parameters for intrinsic variable selection
#' 
#' Control parameters for SPVIM-based intrinsic variable selection.
#' 
#' @param quantity the desired quantity for error-rate control: possible values 
#'    are \code{"gFWER"} (the generalized family-wise error rate), 
#'    \code{"PFP"} (the proportion of false positives), and \code{"FDR"} (the 
#'    false discovery rate).
#' @param base_method the family-wise error rate controlling method to use for 
#'    obtaining the initial set of selected variables. Possible values are
#'    \code{"maxT"} and \code{"minP"} (for step-down procedures based on the 
#'    test statistics ranked from largest to smallest or the p-values ranked from 
#'    smallest to largest, respectively) or \code{"Holm"} for a procedure based
#'    on Holm-adjusted p-values.
#' @param fdr_method the method for controlling the FDR (if 
#'    \code{quantity = "FDR"}); possible values are \code{"BY"} (for 
#'    Benjamini-Yekutieli) or one of the \code{base_method}s.
#' @param q the desired proportion of false positives (only used if 
#'    \code{quantity = "PFP"} or \code{"FDR"}; a fraction between 0 and 1).
#' @param k the desired number of family-wise errors (an integer, greater than 
#'    or equal to zero.) 
#'
#' @return a list with the control parameters.
#' @examples
#' control <- intrinsic_control(quantity = "gFWER", base_method = "Holm", fdr_method = "Holm", 
#'                              k = 1)
#' control
#' @export
intrinsic_control <- function(quantity = "gFWER", base_method = "Holm",
                              fdr_method = "Holm", q = 0.2, k = 5) {
  list(quantity = quantity, base_method = base_method, fdr_method = fdr_method,
       q = q, k = k)
}