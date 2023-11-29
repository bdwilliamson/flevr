#' Super Learner wrapper for a ranger object with variable importance
#'
#' @inheritParams SuperLearner::SL.ranger
#' @inheritParams ranger::ranger
#' 
#' @return a named list with elements \code{pred} (predictions on \code{newX}) and \code{fit} (the fitted \code{ranger} object).
#' @inherit SuperLearner::SL.ranger seealso
#' @inherit SuperLearner::SL.ranger references
#' @examples
#' data("biomarkers")
#' # subset to complete cases for illustration
#' cc <- complete.cases(biomarkers)
#' dat_cc <- biomarkers[cc, ]
#' # use only the mucinous outcome, not the high-malignancy outcome
#' y <- dat_cc$mucinous
#' x <- dat_cc[, !(names(dat_cc) %in% c("mucinous", "high_malignancy"))]
#' feature_nms <- names(x)
#' # get the fit
#' set.seed(20231129)
#' fit <- SL.ranger.imp(Y = y, X = x, newX = x, family = binomial())
#' fit
#' @importFrom stats predict
#' @export
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)),
                           num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = FALSE, importance = "impurity", ...) {
  if (!requireNamespace("ranger", quietly = FALSE)) {
    stop("loading required package ranger failed", call. = FALSE)
  }
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose, importance = importance)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}
