# test that intrinsic selection works

# load required functions and packages
library("testthat")
library("SuperLearner")
library("vimp")

# generate the data -- note that this is a simple setting, for speed
set.seed(4747)
p <- 4
n <- 5e3
x <- replicate(p, stats::rnorm(n, 0, 1))
x_names <- paste0("V", 1:p)
# apply the function to the x's
y <- 1 + 0.5 * x[, 1] + 0.75 * x[, 2] + stats::rnorm(n, 0, 1)

true_set <- c(1, 1, 0, 0)

# set up a library for SuperLearner
learners <- "SL.glm"
univariate_learners <- "SL.glm"
V <- 2
B <- 1e2

# estimate the SPVIMs
set.seed(1234)
est <- suppressWarnings(
  sp_vim(Y = y, X = x, V = V, type = "r_squared",
              SL.library = learners, gamma = .1, alpha = 0.05, delta = 0,
              cvControl = list(V = V), env = environment())
)

# get a base set, based on FWER-controlling procedures
test_that("obtaining an FWER-controlling base set works", {
  holm_set <- get_base_set(test_statistics = est$test_statistic,
                           p_values = est$p_value, alpha = 0.2,
                           method = "Holm")$decision
  maxT_set <- get_base_set(test_statistics = est$test_statistic,
                           p_values = est$p_value, alpha = 0.2,
                           method = "maxT", B = B)$decision
  minP_set <- get_base_set(test_statistics = est$test_statistic,
                           p_values = est$p_value, alpha = 0.2,
                           method = "minP", B = B)$decision
  by_set <- get_base_set(test_statistics = est$test_statistic,
                         p_values = est$p_value, alpha = 0.2,
                         method = "BY", q = 0.2)$decision

  expect_equal(holm_set, true_set)
  expect_equal(maxT_set, true_set)
  expect_equal(minP_set, true_set)
  expect_equal(by_set, true_set)
})

# get an augmented set
holm_set <- get_base_set(test_statistics = est$test_statistic,
                         p_values = est$p_value, alpha = 0.2,
                         method = "Holm")
test_that("obtaining an augmented set works", {
  gfwer_set <- get_augmented_set(p_values = holm_set$p_values,
                                 num_rejected = sum(holm_set$decision),
                                 alpha = 0.2, quantity = "gFWER", k = 1)$set
  expect_equal(gfwer_set, c(0, 0, 1, 0))

  pfp_set <- get_augmented_set(p_values = holm_set$p_values,
                               num_rejected = sum(holm_set$decision),
                               alpha = 0.2, quantity = "PFP", k = 1,
                               q = 0.05)$set
  expect_equal(pfp_set, c(0, 0, 1, 0))

  fdr_set <- get_augmented_set(p_values = holm_set$p_values,
                               num_rejected = sum(holm_set$decision),
                               alpha = 0.2, quantity = "FDR", k = 1,
                               q = 0.05)$set
  expect_equal(fdr_set, c(0, 0, 1, 0))
})

# do the whole procedure
test_that("doing intrinsic selection works", {
  intrinsic_set <- intrinsic_selection(spvim_ests = est, sample_size = n,
                                       alpha = 0.2, feature_names = x_names,
                                       control = list(
                                         quantity = "gFWER", base_method = "Holm",
                                         k = 1
                                       ))
  expect_equal(intrinsic_set$selected, c(TRUE, TRUE, FALSE, TRUE))
})
