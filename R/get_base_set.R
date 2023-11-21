#' Get an initial selected set based on intrinsic importance and a base method
#' 
#' Using the estimated intrinsic importance and a base method 
#' designed to control the family-wise error rate (e.g., Holm), 
#' obtain an initial selected set.
#' 
#' @param test_statistics the test statistics (used with "maxT")
#' @param p_values (used with "minP" or "Holm")
#' @param alpha the alpha level
#' @param method the method (one of "maxT", "minP", or "Holm")
#' @param B the number of resamples (for minP or maxT)
#' @param Sigma the estimated covariance matrix for the test statistics
#' @param q the false discovery rate (for method = "BY")
#' 
#' @importFrom mvtnorm rmvnorm
#' 
#' @return the initial selected set
#' @export
get_base_set <- function(test_statistics, p_values, alpha = 0.05, 
                         method = "maxT", B = 1e4, 
                         Sigma = diag(1, nrow = length(test_statistics)),
                         q = NULL) {
  # rank the test statistics and p-values
  ranks_t <- rank(test_statistics, ties.method = "first")
  ranks_p <- rank(p_values, ties.method = "first")
  ranked_t <- test_statistics[order(-ranks_t)]
  ranked_p <- p_values[order(ranks_p)]
  # get the reversed order
  reverse_ord <- switch((method == "maxT") + 1, rank(ranks_p), rank(-ranks_t))
  p <- length(ranked_t)
  if (method == "BY") {
    # sort p-values in increasing order
    ranked_p <- sort(p_values, decreasing = FALSE)
    denom <- sum(sapply(1:p, function(i) 1 / i))
    # find k
    all_k <- sapply(1:length(ranked_p), function(i) ranked_p[i] <= (i * q) / (p * denom))
    if (!any(all_k)) {
      k <- 0
    } else {
      k <- max(which(all_k))
    }
    decision <- rep(0, length(p_values))
    if (k > 0) {
      decision <- as.numeric(p_values <= sort(p_values, decreasing = FALSE)[k])
    }
    ret_lst <- list(decision = decision, p_values = p_values)
  } else {
    # compute the adjusted p-values according to the base method
    if (method == "Holm") {
      # triangular array of p-values, scaled by (p - k + 1), where k is the column
      pval_array <- matrix(rep(ranked_p * (p - (1:p) + 1), p), nrow = p, ncol = p, 
                           byrow = TRUE)
      pval_array[upper.tri(pval_array)] <- 0
      adj_p <- apply(pmin(pval_array, 1), 1, max)
    } else if (method == "maxT" || method == "minP") {
      # generate some random numbers from the distribution of test statistics
      z <- mvtnorm::rmvnorm(n = B, mean = rep(0, length(test_statistics)), 
                            sigma = Sigma)
      stat <- switch((method == "maxT") + 1, ranked_p, ranked_t)
      if (method == "minP") {
        z <- 1 - stats::pnorm(z)
      }
      comparisons <- sapply(1:B, function(i) {
        z_tri <- matrix(rep(z[i, ], p), nrow = p, ncol = p, byrow = TRUE)
        z_tri[lower.tri(z_tri)] <- NA
        z_summary <- switch((method == "maxT") + 1, 
                            apply(z_tri, 1, min, na.rm = TRUE),
                            apply(z_tri, 1, max, na.rm = TRUE))
        comparison <- switch((method == "maxT") + 1,
                             z_summary <= stat,
                             z_summary >= stat)
      })
      mean_comparison <- rowMeans(comparisons)
      pval_array <- matrix(rep(mean_comparison, p), nrow = p, ncol = p, byrow = TRUE)
      pval_array[upper.tri(pval_array)] <- 0
      adj_p <- apply(pval_array, 1, max)
    } else {
      stop("The requested method has not yet been implemented. Please enter one of 'Holm', 'minP', or 'maxT'.")
    }
    # re-order to match the original indices
    reordered_adj_p <- adj_p[reverse_ord]
    # return
    ret_lst <- list(decision = as.numeric(reordered_adj_p <= alpha),
                    p_values = reordered_adj_p)
  }
  ret_lst
}