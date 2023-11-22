#' flevr: Flexible, Ensemble-Based Variable Selection with Potentially Missing Data
#'
#' A framework for flexible, ensemble-based variable selection using either
#' extrinsic or intrinsic variable importance. You provide
#' the data and a library of candidate algorithms for estimating the
#' conditional mean outcome given covariates; `flevr` handles the rest.
#'
#' @section Author(s):
#' \bold{Maintainer}: Brian Williamson \url{https://bdwilliamson.github.io/}
#'
#' Methodology authors:
#' \itemize{
#'   \item{Brian D. Williamson}
#'   \item{Ying Huang}
#' }
#'
#' @section See Also:
#' Papers:
#' \itemize{
#'   \item{\url{https://arxiv.org/abs/2202.12989}}
#' }
#'
#' Other useful links:
#' \itemize{
#'   \item{\url{https://bdwilliamson.github.io/flevr/}}
#'   \item{\url{https://github.com/bdwilliamson/flevr}}
#'   \item{Report bugs at \url{https://github.com/bdwilliamson/flevr/issues}}
#' }
#'
#' @section Imports:
#' The packages that we import either make the internal code nice
#' (dplyr, magrittr, tibble) or are directly relevant for estimating
#' variable importance (SuperLearner, caret).
#'
#' We suggest several other packages: xgboost, ranger, glmnet, kernlab, polspline
#' and quadprog allow a flexible library of candidate learners in the Super
#' Learner; stabs allows importance to be embedded within stability selection;
#' testthat and covr help with unit tests; and
#' knitr, rmarkdown,and RCurl help with the vignettes and examples.
#'
#' @docType package
#' @name flevr
NULL
