\name{TVAR}
\alias{TVAR}
\alias{OlsTVAR}
\title{Multivariate Treshold Autoregressive model}
\description{
Estimate a multivariate Threshold VAR
}
\usage{
OlsTVAR(data, lag, type=c("level", "difference"), trend=TRUE, nthresh=1, trim=0.1, ngrid, gamma=NULL, thDelay=1, mTh=1, thVar, around, plot=TRUE)
}
\value{
Fitted model data
}
\arguments{
\item{data}{ time series }
\item{lag}{TODO}
\item{type}{TODO}
\item{trend}{TODO}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d)}
\item{thVar}{external threshold variable}
\item{mTh}{combinations of simulatneous variables for the trehhold. either a single value (indicating which variable to take) or a combination}
\item{gamma}{prespecified threshold values}
\item{nthresh}{number of threshold (see details)}
\item{ngrid}{number of elements of the grid, especially for \code{nthresh=3}}
\item{around}{especially for \code{nthresh=3}. The grid search is made around the values given (with length \var{ngrid})}
\item{Trim}{trimming parameter indicating the minimal percentage of observations in each regime}
}
\details{
nthresh=1: estimation upon a grid of thresholds and possible delays

nthresh=2: conditional estimation: conditional on the threshold found in model where nthresh=1, the second thresghold is searced. When both are foun, a second grid search is made with 30 values around each threshold

nthresh=3: DOES NOT estimate a 3 regimes, but a 2 regimes with a whole grid over the parameters, is there rather to check the consistency of the method nthresh=2
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{
data(zeroyld)
data <- zeroyld[,c(36,19)]
colnames(data) <- c("short", "long")

OlsTVAR(data, lag=2, nthresh=3, thDelay=1,trim=0.1, mTh=1, ngrid, plot=FALSE)
#lag2, 2 thresh, trim00.05: 561.46
}
\keyword{ ts }
\references{
TODO
}

