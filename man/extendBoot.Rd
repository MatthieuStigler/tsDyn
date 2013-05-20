\name{extendBoot}
\alias{extendBoot}
\title{extension of the bootstrap replications}
\usage{
  extendBoot(x, nboot)
}
\arguments{
  \item{x}{A object from setarTest (hence of class
  Hansen99Test)}

  \item{nboot}{The number of new bootstrap replications}
}
\value{
  Returns an object of the same class with same objects but
  updated values.
}
\description{
  This function updates an existing bootstrap test with new
  bootstrap replications.
}
\details{
  The plot function will draw the old and new distribution,
  hence allowing to test the sensitivity of the results
  obtained.
}
\examples{
\dontrun{
# test with 10 bootstrap replications:
a<-setarTest(sun[1:100], m=1, nboot=10)
plot(a)

#use old results and compue 20 new replications
b<-extendBoot(a, n=20)
#see the different distributions:
plot(b)
}
}
\author{
  Matthieu Stigler
}
\seealso{
  \code{\link{BBCTest}} for a similar test.
  \code{\link{setarTest}} for a test with stationarity as a
  null.
}

