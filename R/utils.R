# -------------------------------------
# release questions
# -------------------------------------
# @keywords internal
release_questions <- function() {
  c(
    "Have you run cran_prep <- rhub::check_for_cran(env_vars = c(R_COMPILE_AND_INSTALL_PACKAGES = 'always'), show_status = FALSE)?",
    "Have you run devtools::check_win_devel() and devtools::check_win_release()?"
  )
}