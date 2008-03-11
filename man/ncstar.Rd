\name{NCSTAR}
\alias{NCSTAR}
\alias{ncstar}
\title{NCSTAR model}
\description{
TODO
}
\usage{
ncstar(x, m=2, noRegimes, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay, thVar, sig=0.05, trace=TRUE, control=list(), ...)
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{noRegimes}{ max number of regimes }
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
\item{...}{currently unused}
}
\value{
  A fitted NCSTAR model.
}
\details{
TODO
}
\seealso{\code{\link{addRegime}}, \code{\link{star}}}
\author{ J. L. Aznarte }
\examples{
#TODO
}
\keyword{ ts }
