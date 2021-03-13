#' Perform intrinsic, ensemble-based variable selection
#'
#' Based on estimated SPVIM values, do variable selection using the 
#' specified error-controlling method.
#'
#' @param spvim_ests the estimated SPVIM values (a tibble, resulting from a 
#'   call to \code{vimp::sp_vim}).
#' @param sample_size the number of independent observations used to estimate
#'   the SPVIM values.
#' @param alpha the nominal generalized family-wise error rate, proportion of 
#'    false positives, or false discovery rate level to control at (e.g., 0.05).
#' @param control a list of parameters to control the variable selection process.
#'    Parameters include \code{quantity}, \code{base_method}, \code{q}, and 
#'    \code{k}. See \code{\link[flevr]{intrinsic_control}} for details.
#'
#' @return a tibble with the estimated intrinsic variable importance,
#'    the corresponding variable importance ranks, and the selected
#'    variables.
#' @seealso \code{\link[vimp]{sp_vim}} for specific usage of
#'    the \code{sp_vim} function and the \code{vimp} package for estimating 
#'    intrinsic variable importance.
#' @importFrom magrittr `%>%`
#' @export
intrinsic_selection <- function(spvim_ests, sample_size, 
                                alpha = 0.05, control = list()) {
  # make sure that control is in the correct format; if nothing was entered,
  # get default values
  control <- do.call(intrinsic_control, control)
  
  if (control$fdr_method == "BY") {
    # get k
    k <- benjamini_yekutieli(p_values = spvim_ests$p_value, q = control$q)
    selected_set <- rep(0, length(spvim_ests$p_value))
    if (k > 0) {
      # return indices with with p-value less than or equal to p[k]
      selected_set <- as.numeric(
        spvim_ests$p_value <= sort(spvim_ests$p_value, decreasing = FALSE)[k]
      )
    }
  } else {
    var_contrib_v <- spvim_ests$ic$contrib_v[-1, ] %*% 
      t(spvim_ests$ic$contrib_v[-1, ])
    var_contrib_s <- (1 / spvim_ests$gamma) * spvim_ests$ic$contrib_s[-1, ] %*% 
      t(spvim_ests$ic$contrib_s[-1, ])
    cov_mat <- var_contrib_v + var_contrib_s
    test_statistics <- spvim_ests$test_statistic
    p_values <- spvim_ests$p_values
    initial_set_lst <- get_base_set(test_statistics = test_statistics, 
                                    p_values = p_values, alpha = alpha,
                                    method = control$method, B = control$B, 
                                    Sigma = cov_mat / sample_size)
  }
}
