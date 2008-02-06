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
\item{lag}{Number of lags to include in each regime}
\item{type}{If the variables are to be first differenced (difference) or not (level)}
\item{trend}{If a trend should be included in the model}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d) PLEASE NOTE that the notation is currently different to univariate models in tsDyn. The left side variable is taken at time t, and not t+1 as in univariate cases.  }
\item{thVar}{external transition variable}
\item{mTh}{combination of variables with same lag order for the transition variable. Either a single value (indicating which variable to take) or a combination}
\item{gamma}{prespecified threshold values}
\item{nthresh}{number of threshold (see details)}
\item{ngrid}{number of elements of the grid, especially for \code{nthresh=3}}
\item{around}{The grid search is restricted to \var{ngrid} values around this point. Especially useful for \code{nthresh=3}. }
\item{Trim}{trimming parameter indicating the minimal percentage of observations in each regime}
}
\details{
For fixed \code{th} and threshold variable, the model is linear, so estimation can be done directly by CLS (Conditional Least Squares). The search of the parameters values is made upon a grid of potential values. So it is pretty slow.

nthresh=1: estimation of one threshold model (two regimes) upon a grid of \var{ngrid} values (default to ALL) possible thresholds and delays values. 

nthresh=2: conditional estimation: conditional on the threshold found in model where nthresh=1, the second thresghold is searced. When both are foun, a second grid search is made with 30 values around each threshold

nthresh=3: DOES NOT estimate a 3 thresholds model, but a 2 thresholds model with a whole grid over the parameters, is there rather to check the consistency of the method nthresh=2
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{
data(zeroyld)
data <- zeroyld[,c(36,19)]
colnames(data) <- c("short", "long")

OlsTVAR(data, lag=2, nthresh=2, thDelay=1, trim=0.1, mTh=1, plot=TRUE)

##The one threshold (two regimes) gives a value of 10.698 for the threshold and 1 for the delay. Conditional on this values, the search for a second threshold (three regimes) gives 8.129. Starting from this values, a full grid search finds the same values and confims the first step estimation. 
}
\keyword{ ts }
\references{
TODO
}

