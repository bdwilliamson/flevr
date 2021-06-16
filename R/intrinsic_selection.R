#' Perform intrinsic, ensemble-based variable selection
#'
#' Based on estimated SPVIM values, do variable selection using the
#' specified error-controlling method.
#'
#' @param spvim_ests the estimated SPVIM values (an object of class \code{vim},
#'   resulting from a call to \code{vimp::sp_vim}).
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
#' @importFrom dplyr mutate select
#' @export
intrinsic_selection <- function(spvim_ests, sample_size, feature_names,
                                alpha = 0.05,
                                control = list(quantity = "gFWER",
                                               base_method = "Holm",
                                               fdr_method = NULL, q = NULL,
                                               k = NULL)) {
  # make sure that control is in the correct format; if nothing was entered,
  # get default values
  control <- do.call(intrinsic_control, control)
  test_statistics <- spvim_ests$test_statistic
  p_values <- spvim_ests$p_value
  if (control$fdr_method == "BY") {
    selected_set_lst <- get_base_set(test_statistics = test_statistics,
                                     p_values = p_values, alpha = alpha,
                                     method = control$fdr_method, B = control$B,
                                     Sigma = NULL, q = control$q)
    selected_set <- selected_set_lst$decision
    adj_p <- selected_set$p_values
  } else {
    if (is.null(names(spvim_ests$ic))) {
      # cross-fitted SEs were used; create a vector where each column is when
      # that "observation" was in the validation fold
      ic_v <- do.call(cbind, lapply(spvim_ests$ic, function(x) x$contrib_v[-1, ]))
      ic_s <- do.call(cbind, lapply(spvim_ests$ic, function(x) x$contrib_s[-1, ]))
      var_contrib_v <- ic_v %*% t(ic_v)
      var_contrib_s <- (1 / spvim_ests$gamma) * ic_s %*% t(ic_s)
    } else {
      var_contrib_v <- spvim_ests$ic$contrib_v[-1, ] %*%
        t(spvim_ests$ic$contrib_v[-1, ])
      var_contrib_s <- (1 / spvim_ests$gamma) * spvim_ests$ic$contrib_s[-1, ] %*%
        t(spvim_ests$ic$contrib_s[-1, ]) 
    }
    cov_mat <- var_contrib_v + var_contrib_s
    initial_set_lst <- get_base_set(test_statistics = test_statistics,
                                    p_values = p_values, alpha = alpha,
                                    method = control$base_method, B = control$B,
                                    Sigma = cov_mat / sample_size)
    initial_selected_set <- initial_set_lst$decision
    adj_p <- initial_set_lst$p_values
    if (control$quantity == "FDR") {
      q_1 <- alpha
      alpha <- q <- 1 - sqrt(1 - q_1)
    } else {
      q <- control$q
    }
    augmented_set <- get_augmented_set(p_values, sum(initial_selected_set),
                                       alpha, control$quantity, q, control$k)
    selected_set <- ifelse((initial_selected_set == 1) | (augmented_set$set == 1),
                           1, 0)
  }
  importance_df <- spvim_ests$mat %>%
    dplyr::mutate(feature = feature_names, adjusted_p_value = adj_p,
           rank = rank(-abs(est)),
           selected = (selected_set == 1)) %>%
    dplyr::select(feature, est, p_value, adjusted_p_value, rank, selected)
  importance_df
}
