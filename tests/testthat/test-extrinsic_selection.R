# test that extrinsic selection works

# load required functions and packages
library("testthat")
library("SuperLearner")
library("ranger")
library("xgboost")
library("glmnet")
library("kernlab")
library("caret")
library("polspline")

# generate the data -- note that this is a simple setting, for speed
set.seed(4747)
p <- 2
n <- 500
x <- replicate(p, stats::rnorm(n, 0, 1))
x_df <- as.data.frame(x)
x_names <- names(x_df)
# apply the function to the x's
y <- 1 + 0.5 * x[, 1] + 0.75 * x[, 2] + stats::rnorm(n, 0, 1)
y_bin <- as.numeric(0.5 * x[, 1] + 0.75 * x[, 2] + stats::rnorm(n, 0, 1) > 0)

# fit a Super Learner ensemble
set.seed(1234)
learners <- c("SL.xgboost", "SL.ranger.imp", "SL.glm", "SL.glmnet", "SL.ksvm",
              "SL.polymars", "SL.mean")
V <- 2
fit <- SuperLearner::SuperLearner(Y = y, X = x_df,
                                  SL.library = learners,
                                  cvControl = list(V = V))
coef_nms <- names(fit$coef)
fit_bin <- SuperLearner::SuperLearner(Y = y_bin, X = x_df,
                                      family = "binomial",
                                      SL.library = learners,
                                      cvControl = list(V = V))
coef_nms_bin <- names(fit_bin$coef)

# extract algorithm-specific importance
test_that("algorithm-specific importance extraction works", {
  mean_importance <- extract_importance_mean(
    fit = fit$fitLibrary$SL.mean_All$object, coef = fit$coef[grepl("mean", coef_nms)],
    feature_names = x_names
  )
  expect_equal(mean_importance$rank, c(1.5, 1.5))
  glm_importance <- extract_importance_glm(
    fit = fit$fitLibrary$SL.glm_All$object, coef = fit$coef[grepl("glm_All", coef_nms)],
    feature_names = x_names
  )
  expect_equal(glm_importance$feature, c("V2", "V1"))
  glmnet_importance <- extract_importance_glmnet(
    fit = fit$fitLibrary$SL.glmnet_All$object, coef = fit$coef[grepl("glmnet", coef_nms)],
    feature_names = x_names
  )
  expect_equal(glmnet_importance$feature, c("V2", "V1"))
  ranger_importance <- extract_importance_ranger(
    fit = fit$fitLibrary$SL.ranger.imp_All$object, coef = fit$coef[grepl("ranger", coef_nms)],
    feature_names = x_names
  )
  expect_equal(ranger_importance$feature, c("V2", "V1"))
  xgboost_importance <- extract_importance_xgboost(
    fit = fit$fitLibrary$SL.xgboost_All$object, coef = fit$coef[grepl("xgboost", coef_nms)],
    feature_names = x_names
  )
  expect_equal(xgboost_importance$feature, c("V2", "V1"))
  svm_importance <- extract_importance_svm(
    fit = fit$fitLibrary$SL.ksvm_All$object, coef = fit$coef[grepl("svm", coef_nms)],
    feature_names = x_names, x = x_df, y = y
  )
  expect_equal(svm_importance$feature, c("V2", "V1"))
  polymars_importance <- extract_importance_polymars(
    fit = fit$fitLibrary$SL.polymars_All$object, coef = fit$coef[grepl("polymars", coef_nms)],
    feature_names = x_names
  )
  expect_equal(polymars_importance$feature, c("V1", "V2"))
  # note different name, issue raised at SL
  polyclass_importance <- extract_importance_polymars(
    fit = fit_bin$fitLibrary$SL.polymars_All$fit, coef = fit$coef[grepl("polymars", coef_nms)],
    feature_names = x_names
  )
  expect_equal(polyclass_importance$feature, c("V2", "V1"))
})

# extract importance for the whole Super Learner ensemble
test_that("SL importance extraction works", {
  sl_importance_all <- extract_importance_SL(
    fit = fit, feature_names = x_names, import_type = "all", x = x_df, y = y
  )
  expect_equal(sl_importance_all$feature, c("V2", "V1"))
  sl_importance_best <- extract_importance_SL(
    fit = fit, feature_names = x_names, import_type = "best"
  )
  expect_equal(sl_importance_best$feature, c("V2", "V1"))
})

# do extrinsic variable selection
test_that("Extrinsic variable selection works", {
  extrinsic_selected <- extrinsic_selection(
    fit = fit, feature_names = x_names, threshold = 1.5, import_type = "all"
  )
  expect_equal(extrinsic_selected$selected, c(TRUE, FALSE))
})
