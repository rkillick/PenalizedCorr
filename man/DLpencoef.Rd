% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DLpencoef.R
\name{DLpencoef}
\alias{DLpencoef}
\title{DLpencoef}
\usage{
DLpencoef(
  x,
  lag.max = NULL,
  na.action = na.fail,
  penalized = TRUE,
  lh = NULL,
  return.mat = FALSE,
  ...
)
}
\arguments{
\item{x}{a univariate or multivariate numeric time series.}

\item{lag.max}{maximum lag at which to calculate the coefficients. Defaults to the smaller of \eqn{N-1} and \eqn{10log_{10}(N/nser)} where N is the number of non-missing observations and nser is the number of series.}

\item{na.action}{function to be called to handle missing values. 'na.pass' can be used.}

\item{penalized}{'logical'. If 'TRUE' (the default) the penalized PACF is computed; if 'FALSE' the sample PACF is computed.}

\item{lh}{sequence of threshold values across h. Could be a single value (repeated for all h), a single vector of length h (repeated for all nser), or a h x nser matrix. Default is data driven.}

\item{return.mat}{'logical'. If 'TRUE' the return form is a matrix.  If 'FALSE' the return form is an acf object.}

\item{...}{additional arguments for penalized PACF estimation.}
}
\value{
A coefficients vector to estimate a penalized autocorrelation function which is positive definite.
}
\description{
Compute coefficients via the Durbin Levinson algorithm using the penalized partial autocorrelation function (PACF) estimation to obtain an penalized autocorrelation function (ACF) estimation which is positive definite.
}
\examples{
\dontrun{
data <- arima.sim(n=100, model=list(ar=0.5))

DLpencoef(data)
DLpencoef(data, lag.max=10)
}
}
\references{
Gallagher, C., Killick, R., Tan, X. (2024+) Penalized M-estimation 
for Autocorrelation. \emph{Submitted.}
}
\keyword{internal}
