\name{zeroyld}
\alias{zeroyld}
\title{zeroyld time series}
\description{
U.S. Term Structure Data, 1947-1991. Dataset used by Hansen and Seo (2002). The data contains the 12 month short rate and 120 month long rate. 
}
\source{
Hansen, B. and Seo, B. (2002), Testing for two-regime threshold cointegration in vector error-correction models, Journal of Econometrics, 110, pages 293 - 318

The data can be downloaded from: \url{http://www.ssc.wisc.edu/~bhansen/progs/joe_02r.zip}.

The authors themselves took the data from the wepage of Huston McCulloch:
\url{http://www.econ.ohio-state.edu/jhm/ts/mcckwon/mccull.htm}
}

\format{
  A data frame with 482 observations on 2 variables.
 \tabular{rlll}{
    \code{[,1]}  \tab \code{short.run}  \tab numeric  \tab Short term, 12 month\cr
    \code{[,2]}  \tab \code{long.run}   \tab numeric  \tab Long term, 120 month
  }
}

\seealso{
\code{\link{TVECM.HStest}}: Hansen and Seo test. 

\code{\link{TVECM}} for estimating a TVECM.
}


\keyword{datasets}
