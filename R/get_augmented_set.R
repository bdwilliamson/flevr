#' Get an augmented set based on the next-most significant variables
#'
#' Based on the adjusted p-values from a FWER-controlling procedure and a
#' more general error rate for which control is desired (e.g., generalized
#' FWER, proportion of false positives, or FDR), augment the set based on FWER
#' control with the next-most significant variables.
#'
#' @param p_values the adjusted p-values.
#' @param num_rejected the number of rejected null hypotheses from the base
#'     FWER-controlling procedure.
#' @param alpha the significance level.
#' @param quantity the quantity to control (i.e., \code{"gFWER"}, \code{"PFP"},
#'     or \code{"FDR"}).
#' @param q the proportion for FDR or PFP control.
#' @param k the number of false positives for gFWER control.
#'
#' @return a list of the variables selected into the augmentation set. Contains the following values:
#' \itemize{
#'   \item \code{set}, a numeric vector where 1 denotes that the variable was selected and 0 otherwise
#'   \item \code{k}, the value of k used
#'   \item \code{q_star}, the value of q-star used
#' }
#' 
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
#' # estimate SPVIMs (using simple library and V = 2 for illustration only)
#' set.seed(20231129)
#' library("SuperLearner")
#' est <- vimp::sp_vim(Y = y, X = x, V = 2, type = "auc", SL.library = "SL.glm", 
#'                     cvControl = list(V = 2))
#' # get base set
#' base_set <- get_base_set(test_statistics = est$test_statistic, p_values = est$p_value, 
#'                          alpha = 0.2, method = "Holm")
#' # get augmented set
#' augmented_set <- get_augmented_set(p_values = base_set$p_values, 
#'                                    num_rejected = sum(base_set$decision), alpha = 0.2, 
#'                                    quantity = "gFWER", k = 1)
#' augmented_set$set
#' }
#' @export
get_augmented_set <- function(p_values = NULL, num_rejected = 0, alpha = 0.05,
                              quantity = "gFWER", q = 0.05, k = 1) {
  # rank the test statistics and p-values
  ranks_p <- rank(p_values, ties.method = "first")
  ranked_p <- p_values[order(ranks_p)]
  # get the reversed order
  reverse_ord <- rank(ranks_p)
  p <- length(ranked_p)
  if (quantity == "PFP" || quantity == "FDR") {
    if (num_rejected == 0) {
      k <- 0
      q_star <- 0
    } else {
      k <- max(which(sapply(0:(p - num_rejected),
                      function(j) j / (j + num_rejected) <= q)), 
               na.rm = TRUE)
      q_star <- k / (k + num_rejected)
    }
  } else {
    q_star <- q
  }
  augmented_set <- rep(0, length(ranked_p))
  if (k > 0) {
    augmented_set[num_rejected + 1:k] <- 1
  }
  list(set = augmented_set[reverse_ord], k = k, q_star = q_star)
}
