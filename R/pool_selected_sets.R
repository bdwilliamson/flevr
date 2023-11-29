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
#' @examples
#' \donttest{
#' data("biomarkers")
#' x <- biomarkers[, !(names(biomarkers) %in% c("mucinous", "high_malignancy"))]
#' feature_nms <- names(x)
#' library("dplyr")
#' library("SuperLearner")
#' # do multiple imputation (with a small number for illustration only)
#' library("mice")
#' n_imp <- 2
#' set.seed(20231129)
#' mi_biomarkers <- mice::mice(data = biomarkers, m = n_imp, printFlag = FALSE)
#' imputed_biomarkers <- mice::complete(mi_biomarkers, action = "long") %>%
#'   rename(imp = .imp, id = .id)
#' # set up a list to collect selected sets
#' all_selected_vars <- vector("list", length = 5)
# for each imputed dataset, do extrinsic selection
#' for (i in 1:n_imp) {
#'   # fit a Super Learner using simple library for illustration only
#'   these_data <- imputed_biomarkers %>%
#'     filter(imp == i)
#'   this_y <- these_data$mucinous
#'   this_x <- these_data %>%
#'     select(starts_with("lab"), starts_with("cea"))
#'   this_x_df <- as.data.frame(this_x)
#'   fit <- SuperLearner::SuperLearner(Y = this_y, X = this_x_df,
#'                                   SL.library = "SL.glm",
#'                                   cvControl = list(V = 2),
#'                                   family = "binomial")
#'   # do extrinsic selection
#'   all_selected_vars[[i]] <- extrinsic_selection(
#'     fit = fit, feature_names = feature_nms, threshold = 5, import_type = "all"
#'   )$selected
#' }
#' # perform extrinsic variable selection
#' selected_vars <- pool_selected_sets(sets = all_selected_vars, threshold = 1 / n_imp)
#' feature_nms[selected_vars]
#' }
#' @export
pool_selected_sets <- function(sets = list(), threshold = 0.8) {
  set_tib <- do.call(rbind, sets)
  fractions <- colMeans(set_tib)
  fractions >= threshold
}
