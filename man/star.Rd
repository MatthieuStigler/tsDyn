\name{STAR}
\alias{STAR}
\alias{star}
\title{STAR model}
\description{
STAR model fitting with automatic selection of the number of regimes
based on LM tests.
}
\usage{
star(x, m=2, noRegimes, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay, thVar, sig=0.05, trace=TRUE, control=list(), ...)
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{noRegimes}{ max number of regimes }
\item{series}{ time series name (optional) }
\item{rob}{perform robust test (not implemented)}
\item{thDelay}{'time delay' for the threshold variable (as multiple of
  embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the
  threshold variable}
\item{thVar}{external threshold variable}
\item{sig}{significance level for the tests to select the number of regimes.}
\item{control}{further arguments to be passed as \code{control} list to
  \code{\link{optim}}}
\item{trace}{should additional infos be printed out?}
\item{...}{currently unused}
}
\value{
  \code{star} returns an object of class \code{nlar}, subclass
  \code{star}, i.e. a list with informations about the fitted model.
}
\details{
The function star implements the iterative building strategy described
in [1] to identify and estimate Smooth Transition AutoRegressive
models.

[1] T. Terasvirta, "Specification, estimation and evaluation of smooth
transition autoregresive models", J. Am. Stat. Assoc. 89 (1994):
208-218. 
}
\seealso{\code{\link{addRegime}}}
\author{ J. L. Aznarte M. }
\examples{
mod.star <- star(log10(lynx), mTh=c(0,1), control=list(maxit=3000))
mod.star

addRegime(mod.star)
}
\keyword{ ts }

