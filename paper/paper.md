---
# Example from https://joss.readthedocs.io/en/latest/submitting.html
title: 'flevr: Flexible, Ensemble-based Variable Selection in R'
tags:
  - R
  - Variable selection
  - Machine learning
  - Variable importance
authors:
  - name: Brian D. Williamson
    orcid: 0000-0002-7024-548X
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Ying Huang
    orcid: 0000-0002-9655-7502
    affiliation: "2, 3"
affiliations:
 - name: Biostatistics Division, Kaiser Permanente Washington Health Research Institute, Seattle, USA
   index: 1
 - name: Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, USA
   index: 2
 - name: Department of Biostatistics, University of Washington, Seattle, USA
   index: 3
citation_author: Williamson et al.
date: 21 November 2023
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

Variable selection is a common goal in biomedical research, among other fields. Traditional tools for variable selection are generally not robust to departures from possibly restrictive modelling assumptions [@leng2006]. Additionally, in settings with missing data, procedures to combine variable selection results across imputed datasets become more difficult to interpret when machine learning tools are used [@heymans2007]. 

The `flevr` package for `R` [@r_software] provides researchers with the tools necessary to perform variable selection in settings with missing data based on flexible prediction approaches. One approach is to use the Super Learner [@vanderlaan2007super] to perform variable selection based on an ensemble of individual prediction algorithms. The Super Learner has asymptotic and finite-sample guarantees on its prediction performance [@vanderlaan2007super]. The second approach is to base variable selection on _intrinsic variable importance_ [@williamson2021vim2;@williamson2020spvim], which quantifies the population prediction potential of features. The package facilitates both handling missing data using multiple imputation [@mice2011] and incorporating state-of-the-art machine learning algorithms to perform variable selection.

# Statement of need

Variable selection in missing-data contexts is complicated by two issues: (i) a desire for the results to be robust to departures from simplifying modelling assumptions, and (ii) the need to incorporate a missing-data approach. Often, researchers turn to machine learning to increase robustness, but the options for missing-data variable selection when machine learning tools have been used are not easily interpretable. The `flevr` `R` package provides an open-source tool for performing variable selection using flexible, ensemble machine learning algorithms. The _intrinsic selection_ approach (`intrinsic_selection`) is based on population variable importance [@williamson2023flevr]; this approach outputs a single set of selected variables, in contrast to other procedures that require post-hoc harmonization of multiple selected sets [@heymans2007]. We also provide functions to perform variable selection based on the Super Learner [@vanderlaan2007super], improving robustness to model misspecification.

# Examples

This section should serve as a quick guide to using the `flevr` package --- we will cover the main functions for doing extrinsic and intrinsic variable selection using a simulated data example. 

First, we create some data:
```r
# generate the data -- note that this is a simple setting, for speed
set.seed(4747)
p <- 2
n <- 500
# generate features
x <- replicate(p, stats::rnorm(n, 0, 1))
x_df <- as.data.frame(x)
x_names <- names(x_df)
# generate outcomes
y <- 1 + 0.5 * x[, 1] + 0.75 * x[, 2] + stats::rnorm(n, 0, 1)
```

This creates a matrix of covariates `x` with 2 columns and a vector `y` of normally-distributed outcome values for a sample of `n = 500` study participants.

There are two main types of variable selection available in `flevr`: extrinsic and intrinsic. Extrinsic selection is the most common type of variable selection: in this approach, a given algorithm (and perhaps its associated algorithm-specific variable importance) is used for variable selection. The lasso is a widely-used example of extrinsic selection. Intrinsic selection, on the other hand, uses estimated intrinsic variable importance (a population quantity) to perform variable selection. This intrinsic importance is both defined and estimated in a model-agnostic manner.

## Extrinsic variable selection

We recommend using the Super Learner [@vanderlaan2007super] to do extrinsic variable selection to protect against model misspecification. This requires specifying a _library_ of _candidate learners_ (e.g., lasso, random forests). We can do this in `flevr` using the following code:
```r
library("flevr")
set.seed(1234)
# fit a Super Learner ensemble; note its simplicity, for speed
library("SuperLearner")
learners <- c("SL.glm", "SL.mean")
V <- 2
fit <- SuperLearner::SuperLearner(Y = y, X = x_df,
                                  SL.library = learners,
                                  cvControl = list(V = V))
# extract importance based on the whole Super Learner
sl_importance_all <- extract_importance_SL(
  fit = fit, feature_names = x_names, import_type = "all"
)
sl_importance_all

# A tibble: 2 × 2
#   feature  rank
#   <chr>   <dbl>
# 1 V2       1.01
# 2 V1       1.99
```
These results suggest that feature 2 is more important than feature 1 within the Super Learner ensemble (since a lower rank is better). If we want to scrutinize the importance of features within the best-fitting algorithm in the Super Learner ensemble, we can do the following:
```r
# check the best-fitting learner
fit
# Call:  SuperLearner::SuperLearner(Y = y, X = x_df, SL.library = learners, cvControl = list(V = V)) 
 
 
 
#                  Risk       Coef
# SL.glm_All  0.9614026 0.98842732
# SL.mean_All 1.8183525 0.01157268

# importance within the best-fitting learner
sl_importance_best <- extract_importance_SL(
  fit = fit, feature_names = x_names, import_type = "best"
)
sl_importance_best
# # A tibble: 2 × 2
#   feature  rank
#   <chr>   <int>
# 1 V2          1
# 2 V1          2
```
Finally, to do variable selection, we need to select a threshold for the variable importance rank (ideally before looking at the data). In this case, since there are only two variables, we choose a threshold of 1.5, which means we will select only one variable:
```r
extrinsic_selected <- extrinsic_selection(
  fit = fit, feature_names = x_names, threshold = 1.5, import_type = "all"
)
extrinsic_selected
# A tibble: 2 × 3
#   feature  rank selected
#   <chr>   <dbl> <lgl>   
# 1 V2       1.01 TRUE    
# 2 V1       1.99 FALSE  
```
In this case, we select only variable 2.

## Intrinsic variable selection

Intrinsic variable selection is based on population variable importance [@williamson2021vim2;@williamson2020spvim]. Intrinsic selection [@williamson2023flevr] also uses the Super Learner under the hood, and requires specifying a useful _measure of predictiveness_ (e.g., R-squared or classification accuracy). The first step in doing intrinsic selection is estimating the variable importance:
```r
set.seed(1234)

# set up a library for SuperLearner
learners <- "SL.glm"
univariate_learners <- "SL.glm"
V <- 2

# estimate the SPVIMs
library("vimp")
est <- suppressWarnings(
  sp_vim(Y = y, X = x, V = V, type = "r_squared",
              SL.library = learners, gamma = .1, alpha = 0.05, delta = 0,
              cvControl = list(V = V), env = environment())
)
est
# Variable importance estimates:
#       Estimate  SE         95% CI                  VIMP > 0 p-value     
# s = 1 0.1515809 0.06090463 [0.03221005, 0.2709518] TRUE     1.330062e-03
# s = 2 0.2990449 0.06565597 [0.17036157, 0.4277282] TRUE     6.863052e-09
```
This procedure again shows (correctly) that variable 2 is more important than variable 1 in this population, since the point estimate of variable importance is larger. 

The next step is to choose an error rate to control and a method for controlling the family-wise error rate. Here, we choose to control the generalized family-wise error rate, setting `k = 1` and `alpha = 0.2` (so controlling the probability of making more than 1 type I error at 20%), and choose Holm-adjusted p-values to control the initial family-wise error rate: 
```r
intrinsic_set <- intrinsic_selection(
  spvim_ests = est, sample_size = n, alpha = 0.2, feature_names = x_names,
  control = list( quantity = "gFWER", base_method = "Holm", k = 1)
)
intrinsic_set
# A tibble: 2 × 6
#   feature   est       p_value adjusted_p_value  rank selected
#   <chr>   <dbl>         <dbl>            <dbl> <dbl> <lgl>   
# 1 V1      0.152 0.00133           0.00133          2 TRUE    
# 2 V2      0.299 0.00000000686     0.0000000137     1 TRUE
```
In this case, we select both variables.

## Variable selection with missing data

Settings with missing data render variable selection more challenging. In this section, we consider a simulated dataset inspired by data collected by the Early Detection Research Network. Biomarkers developed at six "labs" are validated at at least one of four "validation sites" on 306 cysts. The data also include two binary outcome variables: whether or not the cyst was classified as mucinous, and whether or not the cyst was determined to have high malignant potential. The data contain some missing information, which complicates variable selection; only 212 cysts have complete information. We will use AUC to measure intrinsic importance. We begin by loading the data and performing multiple imputation.

```r
# load the dataset
data("biomarkers")
library("dplyr")
# set up vector "y" of outcomes and matrix "x" of features
y <- biomarkers$mucinous
x <- biomarkers %>%
  na.omit() %>%
  select(starts_with("lab"), starts_with("cea"))
x_names <- names(x)
library("mice")
set.seed(20231121)
mi_biomarkers <- mice::mice(data = biomarkers, m = 5, printFlag = FALSE)
imputed_biomarkers <- mice::complete(mi_biomarkers, action = "long") %>%
  rename(imp = .imp, id = .id)
```

### Extrinsic selection with missing data

We can perform extrinsic variable selection using the imputed data. First, we fit a Super Learner and perform extrinsic variable selection for each imputed dataset. Then, we select a final set of variables based on those that are selected in a pre-specified number of imputed datasets (e.g., 3 of 5) [@heymans2007]. Again, we use a rank of 5 for each imputed dataset to select variables.

```r
extrinsic_learners <- c("SL.glm", "SL.ranger.imp", "SL.xgboost")
set.seed(20231121)
# set up a list to collect selected sets
all_selected_vars <- vector("list", length = 5)
# for each imputed dataset, do extrinsic selection
for (i in 1:5) {
  # fit a Super Learner
  these_data <- imputed_biomarkers %>%
    filter(imp == i)
  this_y <- these_data$mucinous
  this_x <- these_data %>%
    select(starts_with("lab"), starts_with("cea"))
  this_x_df <- as.data.frame(this_x)
  fit <- SuperLearner::SuperLearner(Y = this_y, X = this_x_df,
                                  SL.library = extrinsic_learners,
                                  cvControl = list(V = V),
                                  family = "binomial")
  # do extrinsic selection
  all_selected_vars[[i]] <- extrinsic_selection(
    fit = fit, feature_names = x_names, threshold = 5, import_type = "all"
  )$selected
}
# perform extrinsic variable selection
selected_vars <- pool_selected_sets(sets = all_selected_vars, threshold = 3 / 5)
x_names[selected_vars]
# [1] "lab1_actb"             "lab1_molecules_score"  "lab1_telomerase_score"
```

### Intrinsic selection with missing data

To perform intrinsic variable selection using the imputed data, we first estimate variable importance for each imputed dataset. 

```r
set.seed(20231121)
est_lst <- lapply(as.list(1:5), function(l) {
  this_x <- imputed_biomarkers %>%
    filter(imp == l) %>%
    select(starts_with("lab"), starts_with("cea"))
  this_y <- biomarkers$mucinous
  suppressWarnings(
    sp_vim(Y = this_y, X = this_x, V = V, type = "auc", 
    SL.library = learners, gamma = 0.1, alpha = 0.05, delta = 0,
    cvControl = list(V = V), env = environment())
  )
})
```
Next, we use Rubin's rules [@rubin2018] to combine the variable importance estimates, and use this to perform variable selection. Here, we control the probability of making more than 5 errors (the generalized family-wise error rate) at 5%, using Holm-adjusted p-values to control the initial family-wise error rate.
```r
intrinsic_set <- intrinsic_selection(
  spvim_ests = est_lst, sample_size = nrow(biomarkers),
  feature_names = x_names, alpha = 0.05, 
  control = list(quantity = "gFWER", base_method = "Holm", k = 5)
)
library("dplyr")
intrinsic_set %>%
  filter(selected)
# A tibble: 5 × 6
#   feature                            est p_value adjusted_p_value  rank selected
#   <chr>                            <dbl>   <dbl>            <dbl> <dbl> <lgl>   
# 1 lab1_actb                       0.0334   0.409                1     2 TRUE    
# 2 lab1_telomerase_score           0.0216   0.439                1     5 TRUE    
# 3 lab3_muc3ac_score               0.0238   0.433                1     4 TRUE    
# 4 lab6_ab_score                   0.0241   0.433                1     3 TRUE    
# 5 lab2_fluorescence_mucinous_call 0.0354   0.396                1     1 TRUE
```
We select five variables, here those with the top-5 estimated variable importance. The point estimates and p-values have been computed using Rubin's rules. Two of the variables are the same as those selected using extrinsic selection.

# Availability

The `flevr` package is publicly available on [GitHub](https://github.com/bdwilliamson/flevr); stable releases are publicly available on the [Comprehensive R Archive Network](https://cran.r-project.org/package=flevr). Documentation and examples may be found in the package manual pages and vignettes, and on the `pkgdown` website [@wickham_pkgdown] at [https://bdwilliamson.github.io/flevr](https://bdwilliamson.github.io/flevr).

<!-- # Citations

Citations to entries in paper.bib should be in
[rMarkdown](https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"
-->

# Acknowledgements

The authors gratefully acknowledge bug reports and feature requests submitted by Bhavesh Borate.

# References

