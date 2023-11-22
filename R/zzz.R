.onAttach <- function(...) {
  pkg_desc <- utils::packageDescription("flevr")
  packageStartupMessage(paste0(
    "flevr version ", pkg_desc$Version,
    ": ", pkg_desc$Title
  ))
}

.onLoad <- function(...) {
  # set thread limit for xgboost to pass CRAN check
  Sys.setenv("OMP_THREAD_LIMIT" = 1)
}
