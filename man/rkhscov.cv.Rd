\name{rkhscov.cv}
\alias{rkhscov.cv}
\title{
  Nonparametric operator-regularized covariance function estimation
  with k-fold cross-validation for tuning parameter selection
}
\description{
  This function compute  nonparametric operator-regularized covariance function estimatior
  with k-fold cross-validation for tuning parameter selection.
}
\usage{
rkhscov.cv(time, x, subject, nfold=5, weight=NULL, centered=FALSE, ncpu=nfold,
           pos=TRUE, control.alg=list(), control.cv=list())
}
\arguments{
  \item{time}{
    a vector of time points scaled to [0,1]
  }
  \item{x}{
    a vector of the corresponding observed function values
  }
  \item{subject}{
    a vector of the corresponding subject indices where the index takes value from {1, 2, ..., n}
    with n being the number of subjects
  }
  \item{nfold}{
    number of folds in cross-validation
  }
  \item{weight}{
    internal parameters. Just set to default value: NULL.
  }
  \item{centered}{
    whether x is de-meaned or not.
  }
  \item{ncpu}{
  number of CPU cores
}
  \item{pos}{
    whether semi-postivity of the covariance function estimator is enforced.
  }
  \item{control.alg}{
    a list of parameters for the proximal gradient descent algorithm.
  }
  \item{control.cv}{
    a list of parameters for the cross-validation.
  }
}
\references{R. K. W. Wong and X. Zhang. (2017) "Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data".}
\seealso{
\code{\link{predict.rkhscov, fpca.rkhscov}}
}

\examples{
library(rkhscovfun)

#### generate data ####
set.seed(1234)
dat <- generate.demo.data(n=50, m=10, sigma=0.1, nk=2)


#### fitting ####
ncpu <- 5 # number of cpu; for efficiency, it should be the number of folds in CV.

# compute a trace-norm regularized estimator
res <- rkhscov.cv(time=dat$times, x=dat$y, subject=dat$ids, nfold=5, ncpu=ncpu,
                  centered=FALSE, pos=T, control.alg=list(maxit=1000), control.cv=list(gams=1))

# evaluate the fitted covariance function
tgrid <- seq(0,1,len=50)
C <- predict.rkhscov(tgrid, res$sobj)
image(C)

# FPCA
fres <- fpca.rkhscov(res$sobj)
fres$values # eigenvalues

efun <- compute.fpca(tgrid, fres) # eigenfunctions

}
