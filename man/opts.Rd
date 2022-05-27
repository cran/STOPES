\name{opts}
\alias{opts}

\title{Optimization via Subsampling (OPTS)}

\description{ \code{opts} computes the OPTS MLE in low dimensional
case.}

\usage{
opts(X, Y, m, crit = "aic", prop_split = 0.5, cutoff = 0.75, ...)
}

\arguments{

  \item{X}{n x p covariate matrix (without intercept)}

  \item{Y}{n x 1 binary response vector}

  \item{m}{number of subsamples}

  \item{crit}{information criterion to select the variables: (a) aic =
  minimum AIC and (b) bic = minimum BIC}

  \item{prop_split}{proportion of subsample size and sample size,
  default value = 0.5}

  \item{cutoff}{cutoff used to select the variables using the
  stability selection criterion, default value = 0.75}

  \item{...}{other arguments passed to the glm function, e.g., family
  = "binomial"}

}

\value{

\code{opts} returns a list:

\item{betahat}{OPTS MLE of regression parameter vector}

\item{Jhat}{estimated set of active predictors (TRUE/FALSE)
corresponding to the OPTS MLE}

\item{SE}{standard error of OPTS MLE}

\item{freqs}{relative frequency of selection for all variables}

}

\examples{
require(MASS)
P = 15
N = 100
M = 20
BETA_vector = c(0.5, rep(0.5, 2), rep(0.5, 2), rep(0, P - 5))
MU_vector = numeric(P)
SIGMA_mat = diag(P)

X <- mvrnorm(N, MU_vector, Sigma = SIGMA_mat)
linearPred <- cbind(rep(1, N), X) %*% c(1, BETA_vector)
Y <- rbinom(N, 1, plogis(linearPred))

# OPTS-AIC MLE
opts(X, Y, 10, family = "binomial")

}
