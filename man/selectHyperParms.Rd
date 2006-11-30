\name{selectHyperParms}
\alias{selectLSTAR}
\alias{selectNNET}
\title{Automatic selection of model hyper-parameters}
\description{
Automatic selection of model hyper-parameters
}
\usage{
selectLSTAR(x, m, d=1, steps=d,	mL = 1:m, mH = 1:m, thDelay=0:(m-1))
selectNNET(x, m, d=1, steps=d, size=1:(m+1), maxit=1e3)
}

\arguments{
  \item{x}{ time series }
  \item{m, d, steps}{ embedding parameters. For their meanings, see help about \code{\link{nlar}} }
  \item{mL,mH}{ Vector of \sQuote{low} and \sQuote{high} regimes autoregressive orders}
  \item{thDelay}{ Vector of \sQuote{threshold delay} values}
  \item{size}{ Vector of numbers of hidden units in the nnet model }
  \item{maxit}{ Max. number of iterations for each model estimation }
}
\details{
Functions for automatic selection of LSTAR and NNET models hyper parameters. An exhaustive search over all possible combinations of values of specified hyper-parameters is performed.
Embedding parameters \code{m,d,steps} are kept fixed.

Selection criterion is the usual AIC.
}
\value{
  A data-frame, with columns giving hyper-parameter values and the
  computed AIC for each row (only the best 10s are returned)
}
\author{ Antonio, Fabio Di Narzo }
\examples{
llynx <- log10(lynx)
selectLSTAR(llynx, m=2)
selectNNET(llynx, m=3, size=1:5)
}
\keyword{ ts }
