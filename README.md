# phylomd

[![Travis-CI Build Status](https://travis-ci.org/dunleavy005/phylomd.svg?branch=master)](https://travis-ci.org/dunleavy005/phylomd)

This package implements a post-order tree traversal meta-algorithm that computes phylogenetic stochastic mapping moments and likelihood derivatives of any order and scales linearly in the number of phylogeny tips.
Utility functions related to simulation-based stochastic mapping and higher-order continuous-time Markov chain (CTMC) moment and derivative calculations are also provided.

## Installation

First, install the `devtools` package (if it isn't already installed) using `install.packages("devtools")`.
Check to make sure you have a working development environment by running `devtools::has_devel()`.
If it returns `TRUE`, then your development environment is correctly set up; otherwise, you'll need to install additional tools.
More information on this subject can be found at https://github.com/hadley/devtools.
In addition, a C++11 compliant compiler is needed to build this package.
If you're using Mac OS X, please be aware of [Rcpp-FAQ](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf) question 2.16.

The `phylomd` package can then be installed from GitHub using the `install_github` function.

```r
library(devtools)
install_github("dunleavy005/phylomd")
```

To install the `phylomd` package with vignettes, run `install_github("dunleavy005/phylomd", build_vignettes = TRUE)`.
The vignettes can be accessed by executing `browseVignettes("phylomd")`.

## References

1. Nielsen R (2002) "Mapping Mutations on Phylogenies", *Systematic Biology*, 51(5):729-739.

2. Minin VN and Suchard MA (2008) "Counting labeled transitions in continuous-time Markov models of evolution", *Journal of Mathematical Biology*, 56(3):391-412.

3. Minin VN and Suchard MA (2008) "Fast, Accurate and Simulation-Free Stochastic Mapping", *Philosophical Transactions of the Royal Society B: Biological Sciences*, 363(1512):3985-3995.

4. Kenney T and Gu H (2012) "Hessian Calculation for Phylogenetic Likelihood based on the Pruning Algorithm and its Applications", *Statistical Applications in Genetics and Molecular Biology*, 11(4).

5. Dhar A and Minin VN (2017) "Calculating Higher-Order Moments of Phylogenetic Stochastic Mapping Summaries in Linear Time", *Journal of Computational Biology*, 24(5):377-399.

6. OUR PAPER!!