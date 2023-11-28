<!-- badges: start -->
  [![R-CMD-check](https://github.com/bdwilliamson/flevr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bdwilliamson/flevr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# R/`flevr`: flexible, ensemble-based variable selection

Software author: [Brian Williamson](https://bdwilliamson.github.io)

Methodology authors: [Brian Williamson](https://bdwilliamson.github.io), [Ying Huang](https://www.fredhutch.org/en/faculty-lab-directory/huang-ying.html)

## Introduction

`flevr` is an `R` package for doing variable selection based on flexible ensembles. The package provides functions for extrinsic variable selection using the [Super Learner](https://github.com/ecpolley/SuperLearner) and for intrinsic variable selection using the Shapley Population Variable Importance Measure ([SPVIM](https://github.com/bdwilliamson/vimp)).

The author and maintainer of the `flevr` package is  [Brian Williamson](https://bdwilliamson.github.io). For details on the method, check out our [preprint](https://arxiv.org/abs/2202.12989).

## Installation

You can install a development release of `flevr` from GitHub via `devtools` by running the following code:
```r
# install devtools if you haven't already
# install.packages("devtools", repos = "https://cloud.r-project.org")
devtools::install_github(repo = "bdwilliamson/flevr")
```

## Example

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/bdwilliamson/flevr/issues).

## Citation

After using the `flevr` package, please cite the following:

```
@article{williamson2023,
    author={Williamson, BD and Huang, Y},
    title={Flexible variable selection in the presence of missing data},
    journal={International Journal of Biostatistics},
    year={2023},
    url={https://arxiv.org/abs/2202.12989}
}
```

## License

The contents of this repository are distributed under the MIT license. See below for details:
```
MIT License

Copyright (c) [2021--present] [Brian D. Williamson]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
