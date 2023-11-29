#' Perform intrinsic, ensemble-based variable selection
#'
#' Based on estimated SPVIM values, do variable selection using the
#' specified error-controlling method.
#'
#' @param spvim_ests the estimated SPVIM values (an object of class \code{vim},
#'   resulting from a call to \code{vimp::sp_vim}). Can also be a list of
#'   estimated SPVIMs, if multiple imputation was used to handle missing data; in
#'   this case, Rubin's rules will be used to combine the estimated SPVIMs, and
#'   then selection will be based on the combined SPVIMs.
#' @param sample_size the number of independent observations used to estimate
#'   the SPVIM values.
#' @param feature_names the names of the features (a character vector of
#'    length \code{p} (the total number of features)); only used if the
#'    fitted Super Learner ensemble was fit on a \code{matrix} rather than on a
#'    \code{data.frame}, \code{tibble}, etc.
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
#' # do intrinsic selection
#' intrinsic_set <- intrinsic_selection(spvim_ests = est, sample_size = nrow(dat_cc), alpha = 0.2, 
#'                                      feature_names = feature_nms, 
#'                                      control = list(quantity = "gFWER", base_method = "Holm", 
#'                                                     k = 1))
#' intrinsic_set
#' }
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @importFrom dplyr mutate select
#' @export
intrinsic_selection <- function(spvim_ests = NULL, sample_size = NULL, feature_names = "",
                                alpha = 0.05,
                                control = list(quantity = "gFWER",
                                               base_method = "Holm",
                                               fdr_method = NULL, q = NULL,
                                               k = NULL)) {
  # make sure that control is in the correct format; if nothing was entered,
  # get default values
  control <- do.call(intrinsic_control, control)
  # if a list of SPVIMs was entered, use Rubin's rules to combine the estimates
  # and get variance estimators
  if (!any(grepl("vim", class(spvim_ests))) & any(grepl("list", class(spvim_ests)))) {
    # pool them
    pooled_results <- pool_spvims(spvim_ests)
    test_statistics <- pooled_results$test_statistics
    p_values <- pooled_results$p_values
    tau <- pooled_results$tau_n
    cov_mat <- pooled_results$vcov
    importance_df_init <- tibble::tibble(
      s = seq_len(length(p_values)), est = pooled_results$est,
      se = pooled_results$se, tau = pooled_results$tau_n,
      p_value = p_values, feature = feature_names) %>%
      mutate(rank = rank(-abs(.data$est)))
  } else {
    test_statistics <- spvim_ests$test_statistic
    p_values <- spvim_ests$p_value
    tau <- 0
    cov_mat <- spvim_vcov(spvim_ests)
    importance_df_init <- spvim_ests$mat %>%
      dplyr::mutate(feature = feature_names, rank = rank(-abs(.data$est)))
  }
  if (control$fdr_method == "BY") {
    selected_set_lst <- get_base_set(test_statistics = test_statistics,
                                     p_values = p_values, alpha = alpha,
                                     method = control$fdr_method, B = control$B,
                                     Sigma = NULL, q = control$q)
    selected_set <- selected_set_lst$decision
    adj_p <- selected_set$p_values
  } else {
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
  importance_df <- importance_df_init %>%
    dplyr::mutate(adjusted_p_value = adj_p, selected = (selected_set == 1)) %>%
    dplyr::select("feature", "est", "p_value", "adjusted_p_value", "rank", "selected")
  importance_df
}
