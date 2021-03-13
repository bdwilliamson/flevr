#' flevr: Flexible, ensemble-based variable selection in the presence of missing data
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
#' Preprints:
#' \itemize{
#'   \item{ (to come)}
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
#' (dplyr, magrittr, tibble, rlang, MASS), are directly relevant to estimating
#' the conditional mean (SuperLearner),
#' or are necessary for hypothesis testing (stats).
#'
#' We suggest several other packages: xgboost, ranger, gam, glmnet, polspline,
#' and quadprog allow a flexible library of candidate learners in the Super
#' Learner; testthat and covr help with unit tests; and
#' knitr, rmarkdown,and RCurl help with the vignettes and examples.
#'
#' @docType package
#' @name vimp
NULL
