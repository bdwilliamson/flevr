#' Pool SPVIM Estimates Using Rubin's Rules
#' 
#' If multiple imputation was used due to the presence of missing data,
#' pool SPVIM estimates from individual imputed datasets using Rubin's rules.
#' Results in point estimates averaged over the imputations, along with
#' within-imputation variance estimates and across-imputation variance estimates;
#' and test statistics and p-values for hypothesis testing.
#' 
#' @param spvim_ests a list of estimated SPVIMs (of class \code{vim})
#' @return a list of results
#' 
#' @export
 pool_spvims <- function(spvim_ests) {
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
   p_values <- 1 - pnorm(test_statistics)
   return(list(est = mean_ests[-1], se = sqrt(within_impute_vars[-1]),
               test_statistics = test_statistics, p_values = p_values, 
               tau_n = across_impute_vars[-1], vcov = vcov_mat))
 }