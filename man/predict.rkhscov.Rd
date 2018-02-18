\name{predict.rkhscov}
\alias{predict.rkhscov}
\title{
  Fitted covariance function
}
\description{
  This function evaluates the fitted covariance function.
}
\usage{
predict.rkhscov(newtt, sobj)
}
\arguments{
  \item{newtt}{
    time points for evaluation.
  }
  \item{sobj}{
    the return of \link{rkhscov} or the sobj of the return \link{rkhscov.cv}.
  }
}
\value{
  fitted value of the covariance function over outer(newtt, newtt).
}
\references{R. K. W. Wong and X. Zhang. (2017) "Nonparametric Operator-Regularized Covariance Function Estimation for Functional Data".}
