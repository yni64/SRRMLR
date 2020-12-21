
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SRRMLR

<!-- badges: start -->
<!-- badges: end -->

The goal of SRRMLR is to provides an tractable algorithm for building a
multinomial logistic regression model with the model proposed in \[Wen
et al. (2020)\]. The proposed method combines the merits of rank
reduction and variable selection, constructing a parameter matrix with
row sparsity and low-rank. Four other baseline methods are also included
for comparison (grouped/ungrouped multinomial lasso \[Simon et
al. (2013)\]; NNET \[Ripley et al. (2016)\]; RRVGLM \[Yee and Hastie
(2003)\]). As a representative example, the HAM10000 dataset is included
as well. See Simultaneous Dimension Reduction and Variable Selection for
Multinomial Logistic Regression’ for further information of this
algorithm and related backgrounds.

“Function.R” file includes the prerequisite functions which perform
SRRMLR with specific sparsity, rank and start points. “Simulation.R”
file generates simulation data and do the simulation analysis.
“cv\_srrmlr.R” chooses the best sparsity and rank for our target dataset
X and Y, providing the logistic model and corresponding prediction
error. Meanwhile, the results of other 4 methods are listed in output.
“cmp.R” performs an averaged result for “cv\_srrmlr.R” and presents the
summary of results as output.

## Installation

You can install the released version of SRRMLR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SRRMLR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SRRMLR)
## basic example code

fit = cv_srrmlr(HAM$hamimgs, HAM$hamlab)
fit$tab
#>           Pred r.hat s.hat
#> [1,]  1.937945     1    27
#> [2,]  1.916948    NA    49
#> [3,]  1.949405    NA    34
#> [4,] 13.968729    NA   192
#> [5,]        NA    NA    NA
```
