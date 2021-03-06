# v0.0.2

## Major changes

* Updated PFP and FDR control for intrinsic selection
* Added functions `spvim_vcov` (compute variance-covariance matrix for SPVIMs for intrinsic selection); `pool_spvims` (pool point estimates and obtain Rubin's rules-based variance estimates for multiple imputation-based procedures).
* Added a dataset, `biomarkers`, for use in vignettes

# v0.0.1

## Major changes

* Added function `extract_importance_svm`, which computes extrinsic importance from a `ksvm` object
* Added a vignette on basic usage of the package

## Minor changes

None

# v0.0.0.9000

## Major changes

The first version of the `flevr` package! This package provides functions for variable selection using flexible ensembles. These functions can be used on both complete-case and multiply-imputed data. The two procedures considered are based on extrinsic variable importance (estimated using the Super Learner) and intrinsic variable importance (estimated using the Shapley Population Variable Importance Measure). Both can also be embedded within stability selection using `stabs`.

## Minor changes

None
