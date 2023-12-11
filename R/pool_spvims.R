#' Pool SPVIM Estimates Using Rubin's Rules
#'
#' If multiple imputation was used due to the presence of missing data,
#' pool SPVIM estimates from individual imputed datasets using Rubin's rules.
#' Results in point estimates averaged over the imputations, along with
#' within-imputation variance estimates and across-imputation variance estimates;
#' and test statistics and p-values for hypothesis testing.
#'
#' @param spvim_ests a list of estimated SPVIMs (of class \code{vim})
#' @return a list of results containing the following:
#' \itemize{
#'   \item \code{est}, the average SPVIM estimate over the multiply-imputed datasets
#'   \item \code{se}, the average of the within-imputation SPVIM variance estimates
#'   \item \code{test_statistics}, the test statistics for hypothesis tests of zero importance, using the Rubin's rules standard error estimator and average SPVIM estimate
#'   \item \code{p_values}, p-values computed using the above test statistics
#'   \item \code{tau_n}, the across-imputation variance estimates
#'   \item \code{vcov}, the overall variance-covariance matrix
#' }
#' @examples
#' \donttest{
#' data("biomarkers")
#' library("dplyr")
#' # do multiple imputation (with a small number for illustration only)
#' library("mice")
#' n_imp <- 2
#' set.seed(20231129)
#' mi_biomarkers <- mice::mice(data = biomarkers, m = n_imp, printFlag = FALSE)
#' imputed_biomarkers <- mice::complete(mi_biomarkers, action = "long") %>%
#'   rename(imp = .imp, id = .id)
#' # estimate SPVIMs for each imputed dataset, using simple library for illustration only
#' library("SuperLearner")
#' est_lst <- lapply(as.list(1:n_imp), function(l) {
#'   this_x <- imputed_biomarkers %>%
#'     filter(imp == l) %>%
#'     select(starts_with("lab"), starts_with("cea"))
#'   this_y <- biomarkers$mucinous
#'   suppressWarnings(
#'     vimp::sp_vim(Y = this_y, X = this_x, V = 2, type = "auc", 
#'                  SL.library = "SL.glm", gamma = 0.1, alpha = 0.05, delta = 0,
#'                  cvControl = list(V = 2), env = environment())
#'   )
#' })
#' # pool the SPVIMs using Rubin's rules
#' pooled_spvims <- pool_spvims(spvim_ests = est_lst)
#' pooled_spvims
#' }
#' @export
 pool_spvims <- function(spvim_ests = NULL) {
   M <- length(spvim_ests)
   delta <- spvim_ests[[1]]$delta
   # compute the mean SPVIM estimates over the imputations
   all_ests <- do.call(cbind, lapply(spvim_ests, function(l) l$est))
   mean_ests <- rowMeans(all_ests)
   all_zero_ests <- do.call(cbind, lapply(spvim_ests, function(l) l$test_statistic_computation$est_0))
   mean_zero_est <- rowMeans(all_zero_ests)
   # compute across-imputation variance
   across_impute_vars <- rowSums(sweep(all_ests, MARGIN = 1, STATS = mean_ests, FUN = "-") ^ 2) / (M - 1)
   across_impute_zero_vars <- sum((all_zero_ests - mean_zero_est) ^ 2) / (M - 1)
   # add the zero-variable value on, for hypothesis testing
   ests_plus <- mean_ests + mean_ests[1]
   # compute within-imputation variance
   within_impute_vars <- colMeans(do.call(rbind, lapply(spvim_ests, function(l) {
     l$test_statistic_computation$ses ^ 2
   })))
   within_impute_vars_0 <- colMeans(do.call(rbind, lapply(spvim_ests, function(l) {
     l$test_statistic_computation$se_0 ^ 2
   })))
   # compute overall variance-covariance matrix
   vcov_arry <- sapply(seq_len(length(spvim_ests)), function(i) spvim_vcov(spvim_ests[[i]]),
                       simplify = "array")
   vcov_mat <- sweep(apply(vcov_arry, c(1, 2), mean), MARGIN = 1,
                     STATS = ((M + 1) / M) * across_impute_vars[-1], FUN = "+")
   # compute test statistics and p-values based on within-, across-imputation variance
   test_statistics <- unlist(lapply(2:length(mean_ests), function(j) {
     var_n_j <- within_impute_vars[j] + within_impute_vars_0
     tau_n_j <- across_impute_vars[j] + across_impute_zero_vars
     (ests_plus[j] - mean_zero_est - delta) / sqrt(var_n_j + ((M + 1) / M) * tau_n_j)
   }))
   p_values <- 1 - stats::pnorm(test_statistics)
   return(list(est = mean_ests[-1], se = sqrt(within_impute_vars[-1]),
               test_statistics = test_statistics, p_values = p_values,
               tau_n = across_impute_vars[-1], vcov = vcov_mat))
 }
