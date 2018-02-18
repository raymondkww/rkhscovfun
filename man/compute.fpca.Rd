\name{compute.fpca}
\alias{compute.fpca}
\title{
  Evaluate the eigenfunctions based on the output of \link{fpca.rkhscov}
}
\description{
  This function evaluate the eigenfunction of (L_2) functional principal component analysis.
}
\usage{
compute.fpca(tt, fpca.obj)
}
\arguments{
  \item{tt}{
    a vector of time points at which the eigenfunctions are evaluated
  }
  \item{sobj}{
    the return of \link{fpca.rkhscov}
  }
}
\value{
  \item{values}{
    L_2 eigenvalues
  }
}
\references{R. K. W. Wong and X. Zhang. (2017) "Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data".}
