\name{selectSETAR}
\alias{selectSETAR}
\title{Automatic selection of SETAR hyper-parameters}
\description{
Automatic selection of SETAR hyper-parameters
}
\usage{
selectSETAR(x, m, d=1, steps=d, thSteps=7, 
	mL = 1:m, mH = 1:m, 
	th=quantile(x, prob=seq(0.15, 0.85, length=thSteps) ), 
	thDelay=0:(m-1), criterion=c("pooled-AIC","AIC"))
}

\arguments{
  \item{x}{ time series }
  \item{m, d, steps}{ embedding parameters. For their meanings, see help about \code{\link{nlar}} }
  \item{thSteps}{ Number of steps along different values of threshold
    (if \code{th} omitted) }
  \item{th}{ Vector of threshold values }
  \item{mL,mH}{ Vector of \sQuote{low} and \sQuote{high} regimes autoregressive orders}
  \item{thDelay}{ Vector of \sQuote{threshold delay} values}
  \item{criterion}{ Model selection criterion  }
}
\details{
  Routine for automatic selection of SETAR models hyper parameters.
  
  An exhaustive search over all possible combinations of values of
  specified hyper-parameters is performed.
  
  Embedding parameters \code{m,d,steps} are kept fixed.

  Possible criteria are the usual AIC and a pooled AIC formula:
  \eqn{AIC(low regime model) + AIC(high regime model)}. The default
  criterion is the pooled AIC formula.
}
\value{
  A data-frame, with columns giving hyper-parameter values and the
  computed AIC for each row (only the best 10s are returned)
}
\seealso{
  \code{\link{selectLSTAR}}, \code{\link{selectNNET}}
}
\author{ Antonio, Fabio Di Narzo }
\examples{
llynx <- log10(lynx)
selectSETAR(llynx, m=2)
#Suggested model is the following:
setar(llynx, m=2, thDelay=1, th=3.4)
}
\keyword{ ts }
