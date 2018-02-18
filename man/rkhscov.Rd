\name{rkhscov}
\alias{rkhscov}
\title{
  Nonparametric operator-regularized covariance function estimation
}
\description{
  This function computes a nonparametric operator-regularized covariance function estimatior.
}
\usage{
rkhscov.cv(time, x, subject, lam, gam=1, weight=NULL, centered=FALSE, B0v=NULL,
           pos=TRUE, control=list())
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
  \item{gam}{
    weight between trace norm (gam=1) and (squared) Hilbert-Schmidt norm (gam=0).
    For instance, the penalty takes the form: lam * (gam * sum_j w_j |eig_j| + (1-gam)/2 * sum_j eig_j^2) if pos=false
  }
  \item{weight}{
    internal parameter. Just set to default value: NULL.
  }
  \item{centered}{
    whether x is de-meaned or not.
  }
  \item{B0v}{
    internal parameter. Just set to default value: NULL.
  }
  \item{pos}{
    whether semi-postivity of the covariance function estimator is enforced.
  }
  \item{control}{
    a list of parameters for the proximal gradient descent algorithm.
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

# compute a trace-norm regularized estimator
res <- rkhscov(time=dat$times, x=dat$y, subject=dat$ids, lam=2e-05, gam=1,
               centered=FALSE, pos=T)

# evaluate the fitted covariance function
tgrid <- seq(0,1,len=50)
C <- predict.rkhscov(tgrid, res)
image(C)

# FPCA
fres <- fpca.rkhscov(res)
fres$values # eigenvalues

efun <- compute.fpca(tgrid, fres) # eigenfunctions

}
