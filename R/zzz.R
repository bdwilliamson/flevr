.onAttach <- function(...) {
  pkg_desc <- utils::packageDescription("flevr")
  packageStartupMessage(paste0(
    "flevr version ", pkg_desc$Version,
    ": ", pkg_desc$Title
  ))
  packageStartupMessage(paste0(
    "Package created on ", utils::packageDate("flevr")
  ))
}
