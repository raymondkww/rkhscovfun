# rkhscovfun: Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data

This package implements the covariance function estimation proposed in Wong and Zhang (2017).

## Installation
This package can be installed via function `install_github` in R package `devtools`:

``` r
install.packages("devtools")
devtools::install_github("raymondkww/rkhscovfun")

```

## Example
```r
library("rkhscovfun")

#### generate demo data ####
set.seed(1234)
dat <- generate.demo.data(n=50, m=10, sigma=0.1, nk=2)


#### fitting ####
# compute a trace-norm regularized estimator (with pre-chosen tuning parameter;
# see rkhscov.cv for tuning parameter selection via k-fold cross-validation)
res <- rkhscov(time=dat$times, x=dat$y, subject=dat$ids, lam=2e-05, gam=1,
               centered=FALSE, pos=T)

# evaluate the fitted covariance function
tgrid <- seq(0,1,len=50)
C <- predict.rkhscov(tgrid, res)
image(C)

#### FPCA ####
fres <- fpca.rkhscov(res)
fres$values # eigenvalues
efun <- compute.fpca(tgrid, fres) # eigenfunctions
```

## Installation issue on MacOS
If you are using a precompiled binary distribution of R for Mac OS X, you may experience a "-lgfortran" error during the installation of the package. This error is due to incorrect version of gfortran binary. You can find a solution [\[here\]](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/).


## References
* R. K. W. Wong and X. Zhang. (2017) "Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data". [\[arXiv\]](http://arxiv.org/abs/1701.06263)
