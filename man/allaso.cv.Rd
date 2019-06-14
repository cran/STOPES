\name{alasso.cv}
\alias{alasso.cv}

\title{ALASSO variable selection via cross-validation regularization
parameter selection}

\description{ \code{alasso.cv} computes the ALASSO estimator.}

\usage{
alasso.cv(x, y)
}

\arguments{

  \item{x}{n x p covariate matrix}

  \item{y}{n x 1 response vector}

}

\value{

\code{alasso.cv} returns the ALASSO estimate

\item{alasso}{the ALASSO estimator}

}

\references{Hui Zou, (2006). "The adaptive LASSO and its oracle properties", JASA, 101 (476), 1418-1429}

\examples{

p <- 5
n <- 100
beta <- c(2, 1, 0.5, rep(0, p - 3))
x <- matrix(nrow = n, ncol = p, rnorm(n * p))
y <- rnorm(n) + crossprod(t(x), beta)
alasso.cv(x, y)

}
