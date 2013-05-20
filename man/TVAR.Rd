\name{TVAR}
\alias{OlsTVAR}
\alias{TVAR}
\title{Multivariate Threshold Autoregressive model}
\usage{
  TVAR(data, lag,
    include = c("const", "trend", "none", "both"),
    model = c("TAR", "MTAR"), commonInter = FALSE,
    nthresh = 1, thDelay = 1, mTh = 1, thVar, trim = 0.1,
    ngrid, gamma = NULL, around, plot = FALSE,
    dummyToBothRegimes = TRUE, trace = TRUE, trick = "for",
    max.iter = 2)
}
\arguments{
  \item{data}{time series}

  \item{lag}{Number of lags to include in each regime}

  \item{include}{Type of deterministic regressors to
  include}

  \item{model}{Whether the transition variable is taken in
  levels (TAR) or difference (MTAR)}

  \item{commonInter}{Whether the deterministic regressors
  are regime specific (commonInter=FALSE) or not.}

  \item{nthresh}{Number of thresholds}

  \item{thDelay}{'time delay' for the threshold variable
  (as multiple of embedding time delay d) PLEASE NOTE that
  the notation is currently different to univariate models
  in tsDyn. The left side variable is taken at time t, and
  not t+1 as in univariate cases.}

  \item{mTh}{combination of variables with same lag order
  for the transition variable. Either a single value
  (indicating which variable to take) or a combination}

  \item{thVar}{external transition variable}

  \item{trim}{trimming parameter indicating the minimal
  percentage of observations in each regime}

  \item{ngrid}{number of elements of the grid, especially
  for \code{nthresh=3}}

  \item{gamma}{prespecified threshold values}

  \item{around}{The grid search is restricted to
  \var{ngrid} values around this point. Especially useful
  for \code{nthresh=3}.}

  \item{plot}{Whether a plot showing the results of the
  grid search should be printed}

  \item{dummyToBothRegimes}{Whether the dummy in the one
  threshold model is applied to each regime or not.}

  \item{trace}{should additional infos be printed out?}

  \item{trick}{type of R function called: \code{for} or
  \code{mapply}}

  \item{max.iter}{Number of iterations for the algorithm}
}
\value{
  An object of class TVAR, with standard methods.
}
\description{
  Estimate a multivariate Threshold VAR
}
\details{
  For fixed \code{th} and threshold variable, the model is
  linear, so estimation can be done directly by CLS
  (Conditional Least Squares). The search of the parameters
  values is made upon a grid of potential values. So it is
  pretty slow.

  nthresh=1: estimation of one threshold model (two
  regimes) upon a grid of \var{ngrid} values (default to
  ALL) possible thresholds and delays values.

  nthresh=2: estimation of two thresholds model (three
  regimes) Conditional on the threshold found in model
  where nthresh=1, the second threshold is searched. When
  both are found, a second grid search is made with 30
  values around each threshold.

  nthresh=3: DOES NOT estimate a 3 thresholds model, but a
  2 thresholds model with a whole grid over the thresholds
  parameters (so is really slow) with a given delay, is
  there rather to check the consistency of the method
  nthresh=2
}
\examples{
data(zeroyld)

data<-zeroyld

TVAR(data, lag=2, nthresh=2, thDelay=1, trim=0.1, mTh=1, plot=TRUE)

##The one threshold (two regimes) gives a value of 10.698 for the threshold and 1 for the delay.
Conditional on this values, the search for a second threshold (three regimes) gives 8.129.
Starting from this values, a full grid search finds the same values and confims the first step estimation.
}
\author{
  Matthieu Stigler
}
\references{
  Lo and Zivot (2001) "Threshold Cointegration and
  Nonlinear Adjustment to the Law of One Price,"
  Macroeconomic Dynamics, Cambridge University Press, vol.
  5(4), pages 533-76, September.
}
\seealso{
  \code{\link{lineVar}} for the linear VAR/VECM,
  \code{\link{TVAR.LRtest}} to test for TVAR,
  \code{\link{TVAR.sim}} to simulate/bootstrap a TVAR.
}
\keyword{ts}

