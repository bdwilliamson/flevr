#' Extract the learner-specific importance from a mean object
#' 
#' Extract the individual-algorithm extrinsic importance from a mean object, 
#' along with the importance rank.
#' 
#' @param fit the \code{mean} object.
#' @param feature_names the feature names
#' @param coef the Super Learner coefficient associated with the learner.
#' 
#' @return a tibble, with columns \code{algorithm} (the fitted algorithm),
#'   \code{feature} (the feature), \code{importance} (the algorithm-specific 
#'   extrinsic importance of the feature), \code{rank} (the feature importance 
#'   rank, with 1 indicating the most important feature), and \code{weight} 
#'   (the algorithm's weight in the Super Learner)
#' @export 
extract_importance_mean <- function(fit, feature_names, coef = 0) {
  imp_dt <- tibble::tibble(algorithm = "mean", feature = feature_names, 
                           importance = NA, 
                           rank = rep(mean(1:length(feature_names)), 
                                      length(feature_names)), 
                           weight = coef)
  imp_dt
}
