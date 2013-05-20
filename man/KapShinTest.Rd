\name{KapShinTest}
\alias{KapShinTest}
\alias{KapShinTest}
\title{Test of unit root against SETAR alternative with}
\description{Test of unit root against a stationnary 3 regime SETAR alternative with random walk in the inner regime
}
\usage{
 KapShinTest(x, m=1, series, include = c("none","const", "trend", "both"), 
  c=3, delta=0.5, points=NULL,minObsMid=10, 
  trick=c("for", "apply", "mapply"), trace=FALSE)
}


\arguments{
\item{x}{ time series }
\item{m}{Number of lags under the alternative}
\item{series}{ time series name (optional) }
\item{include}{Whether data should be raw, de-meaned or de-meaned and de-trended}
\item{c}{Argument for the grid search. See details}
\item{delta}{Argument for the grid search. See details}
\item{points}{Points for the grid search. See details}
\item{minObsMid}{Minimal number of observations in the inner regime}
\item{trick}{type of internal function used}
\item{trace}{should additional infos be printed? (logical)}
}

\value{A object of class \code{KapShin2006Test} containing:
\item{statistic}{The three (SupW, AvgW, ExpW) test statistics computed}
\item{case}{Whether the data was transformed, corresponds to input argument include}
\item{series}{The name of the series}

}

\details{
This function is currently spurious.  

}
\seealso{
\code{\link{BBCTest}} for a similar test. \code{\link{setarTest}} for a test with stationarity as a null.
}

\author{ Matthieu Stigler}
\examples{
KapShinTest(lynx, m=1, trace=FALSE, include="none", points=10)
}
