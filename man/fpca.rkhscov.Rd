\name{fpca.rkhscov}
\alias{fpca.rkhscov}
\title{
  Functional principal component analysis
}
\description{
  This function performs (L_2) functional principal component analysis.
}
\usage{
fpca.rkhscov(sobj, Q=NULL)
}
\arguments{
  \item{sobj}{
    the return of \link{rkhscov} or the sobj of the return \link{rkhscov.cv}.
  }
  \item{Q}{
    internal parameters. Just set to default value: NULL.
  }
}
\value{
  \item{values}{
    L_2 eigenvalues
  }
}
\references{R. K. W. Wong and X. Zhang. (2017) "Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data".}

\seealso{
\code{\link{compute.fpca}}
}
