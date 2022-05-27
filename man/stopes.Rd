\name{stopes}
\alias{stopes}

\title{Selection of Threshold OPtimized Empirically via Splitting (STOPES)}

\description{ \code{stopes} computes the STOPES estimator.}

\usage{
stopes(x, y, m = 20, prop_split = 0.50, prop_trim = 0.20, q_tail = 0.90)
}

\arguments{

  \item{x}{n x p covariate matrix}

  \item{y}{n x 1 response vector}

	\item{m}{number of split samples, with default value = 20}

	\item{prop_split}{proportion of data used for training samples, default value = 0.50} 

	\item{prop_trim}{proportion of trimming, default prop_trim = 0.20}

	\item{q_tail}{proportion of truncation samples across the split samples, default values = 0.90}

}

\value{

\code{stopes} returns a list with the STOPE estimates via data splitting using 0.25 method and the PELT method:

\item{beta_stopes}{the STOPE estimate via data splitting}

\item{J_stopes}{the set of active predictors corresponding to STOPES via data splitting}

\item{final_cutpoints}{the final cutpoint for STOPES}

\item{beta_pelt}{the STOPE estimate via PELT}

\item{J_pelt}{the set of active predictors corresponding to STOPES via PELT}

\item{final_cutpoints_PELT}{the final cutpoint for PELT}

\item{quan_NA}{test if the vector of trimmed cutpoints has length 0, with 1 if TRUE and 0 otherwise}

}

\author{
	Marinela Capanu, Mihai Giurcanu, Colin Begg, and Mithat Gonen
}

\examples{

p <- 5
n <- 100
beta <- c(2, 1, 0.5, rep(0, p - 3))
x <- matrix(nrow = n, ncol = p, rnorm(n * p))
y <- rnorm(n) + crossprod(t(x), beta)
stopes(x, y)

}
