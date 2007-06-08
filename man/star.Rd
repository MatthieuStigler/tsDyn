\name{STAR}
\alias{STAR}
\alias{star}
\alias{isLinear.star}
\title{STAR model}
\description{
STAR model fitting
}
\usage{
star(x, m=3, d = 1, steps = d, series, rob = FALSE,
	mTh, thDelay=1, thVar, 
	sig=0.05, trace=TRUE, control=list(), ...)
\method{isLinear}{star}(object, rob=FALSE, sig=0.05, trace=TRUE,...)
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{rob}{perform robust test?}
\item{thDelay}{'time delay' for the threshold variable (as multiple of
  embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the
  threshold variable}
\item{thVar}{external threshold variable}
\item{sig}{significance level for no. of regimes testing}
\item{control}{further arguments to be passed as \code{control} list to
  \code{\link{optim}}}
\item{trace}{should additional infos be printed out?}
\item{object}{fitted \code{star} object}
\item{...}{currently unused}
}
\value{
  \code{star} returns an object of class \code{nlar}, subclass \code{star}, i.e. a list
  with fitted model informations.

  The \code{isLinear} method performs a 3rd order Taylor expansion LM test of linearity
  against 2-regimes star model.
}
\details{
TODO
}
\seealso{\code{\link{addRegime}}}
\author{ J. L. Aznarte }
\examples{
##TODO
}
\keyword{ ts }
