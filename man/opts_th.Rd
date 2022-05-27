\name{opts_th}
\alias{opts_th}

\title{Threshold OPTimization via Subsampling (OPTS_TH)}

\description{ \code{opts_th} computes the threshold OPTS MLE in low
dimensional case.}

\usage{
opts_th(X, Y, m, crit = "aic", type = "binseg", prop_split = 0.5,
  prop_trim = 0.2, q_tail = 0.5, ...)
}

\arguments{

  \item{X}{n x p covariate matrix (without intercept)}

  \item{Y}{n x 1 binary response vector}

  \item{m}{number of subsamples}

  \item{crit}{information criterion to select the variables: (a) aic =
  minimum AIC and (b) bic = minimum BIC}

  \item{type}{method used to minimize the trimmed and averaged
  information criterion: (a) min = observed minimum subsampling
  trimmed average information, (b) sd = observed minimum using the
  0.25sd rule (corresponding to OPTS-min in the paper), (c) pelt =
  PELT changepoint algorithm (corresponding to OPTS-PELT in the
  paper), (d) binseg = binary segmentation changepoint algorithm
  (corresponding to OPTS-BinSeg in the paper), (e) amoc = AMOC
  method.}

  \item{prop_split}{proportion of subsample size of the sample size;
  default value is 0.5}

  \item{prop_trim}{proportion that defines the trimmed mean; default
  value = 0.2}

  \item{q_tail}{quantiles for the minimum and maximum p-values across
  the subsample cutpoints used to define the range of cutpoints}

  \item{...}{other arguments passed to the glm function, e.g., family
  = "binomial"}

}

\value{

\code{opts_th} returns a list:

\item{betahat}{STOPES MLE of regression parameters}

\item{SE}{SE of STOPES MLE}

\item{Jhat}{set of active predictors (TRUE/FALSE) corresponding to
STOPES MLE}

\item{cuthat}{estimated cutpoint for variable selection}

\item{pval}{marginal p-values from univariate fit}

\item{cutpoits}{subsample cutpoints}

\item{aic_mean}{mean subsample AIC}

\item{bic_mean}{mean subsample BIC}

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

# Threshold OPTS-BinSeg MLE
opts_th(X, Y, M, family = "binomial")

}
