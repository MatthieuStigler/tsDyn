\name{BBCTest}
\alias{BBCTest}
\title{Test of unit root against SETAR alternative}
\usage{
BBCTest(x, m, series, testStat = c("LR", "Wald", "LM"), trim = 0.1,
  grid = c("minPerc", "minObs"))
}
\arguments{
  \item{x}{time series}

  \item{m}{Number of lags under the alternative}

  \item{series}{time series name (optional)}

  \item{testStat}{Type of test statistic to use}

  \item{trim}{trimming parameter indicating the minimal
  percentage of observations in each regime}

  \item{grid}{Whether a minimal number of percentage or
  observations should be imposed. See details}
}
\value{
A object of class "BBC2004Test" containing:

-The value of the sup Test

-The version of test used (either Wald, LM or LR).
}
\description{
Test of unit root against a stationnary three regime SETAR
alternative
}
\details{
TODO
}
\examples{
BBCTest(lynx, m=3, test="Wald", grid="minPerc")
}
\author{
Matthieu Stigler
}
\seealso{
\code{\link{setarTest}} for a test with stationarity as a
null.
}

