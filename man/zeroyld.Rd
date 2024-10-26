\name{zeroyld}
\alias{zeroyld}
\alias{zeroyldMeta}
\title{zeroyld time series}
\description{
U.S. Term Structure Data, 1951-1991. Dataset used by Hansen and Seo (2002). The data contains the 12 month short rate and 120 month long rate. 
}

\usage{
zeroyld
zeroyldMeta
}

\source{
Hansen, B. and Seo, B. (2002), Testing for two-regime threshold cointegration in vector error-correction models, Journal of Econometrics, 110, pages 293 - 318

The data can be downloaded from: \url{https://www.ssc.wisc.edu/~bhansen/progs/joe_02r.zip}.

The authors themselves took the data from the webpage of Huston McCulloch:
\url{https://www.asc.ohio-state.edu/mcculloch.2/ts/mcckwon/mccull.htm}
% old: \url{https://www.econ.ohio-state.edu/jhm/ts/mcckwon/mccull.htm}
}

\format{
 \code{zeroyld} contains two variables, while \code{zeryldMeta} contains also \code{Year} and \code{Month} columns. 
  

  \code{zeroyld} is a data frame with 482 observations and 2 variables:
 \tabular{lll}{
    \code{short.run}  \tab numeric  \tab Short term, 12 month\cr
    \code{long.run}   \tab numeric  \tab Long term, 120 month
  }
}

\seealso{
\code{\link{TVECM.HStest}}: Hansen and Seo test. 

\code{\link{TVECM}} for estimating a TVECM.
}

\examples{
data(zeroyldMeta)
plot(zeroyldMeta$Date, zeroyldMeta$short.run, type = "l", xlab = "Date", ylab ="Rate")
lines(zeroyldMeta$Date, zeroyldMeta$long.run, lty = 2)
legend("topleft", lty = c(1, 2), legend = c("Short rate: 12 months", "Long rate: 120 months"))
}

\keyword{datasets}
