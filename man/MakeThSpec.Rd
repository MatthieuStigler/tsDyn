\name{MakeThSpec}
\alias{MakeThSpec}
\alias{makeThSpec}
\title{Specification of the threshold search}
\usage{
MakeThSpec(ngrid = c("All", "Half", "Third", "Quarter"), exact = NULL,
  int = c("from", "to"), around = "val", ...)
}
\arguments{
  \item{exact}{The user give an exact threshold value}

  \item{int}{The user gives an interval to search inside}

  \item{around}{The user gives an point to search around}

  \item{ngrid}{The number of values to search for}

  \item{...}{currently unused}
}
\value{
The input values are given as output after checking for
consistency (only one of exact/int/around should be given).
}
\description{
This optional function allows the user to set different
restrictions for the threshold grid search in function
\code{\link{selectSETAR}}.
}
\details{
This function is just to check the inputs for the
specification of the grid search. If not provided, the
search will be in the biggest interval (\verb{ngrid
="All"}) between the minimum and maximum values. The user
can reduce it by giving setting "Half" (only every two
points is taken) and so on, or setting a number.

The search can also be made around a point, or between two
points. When between a point, the argument \code{ngrid} is
still used, whereas for around, a value of 30 is taken as
default value if \code{ngrid} is not specified by user.
}
\examples{
sun<-(sqrt(sunspot.year+1)-1)*2
selectSETAR(sun, m=3, th=MakeThSpec(exact=10.40967),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=1)
#when pre-sepcified value does not correspond, function will search nearest value
selectSETAR(sun, m=3, th=MakeThSpec(exact=10.4),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=1)
#search around:
selectSETAR(sun, m=3, th=MakeThSpec(around=10.40967, ngrid=20),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=1)
#search in an interval
selectSETAR(sun, m=3, th=MakeThSpec(int=c(10, 11), ngrid=20),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=1)
#reduce size of the grid:
selectSETAR(sun, m=3, th=MakeThSpec(ngrid="Half"),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=1)


# 2 thresholds:
selectSETAR(sun, m=3, th=MakeThSpec(ngrid="Half"),criterion="SSR", d=1, thDelay=0:2,
           plot=FALSE, nthresh=2)
}
\author{
Matthieu Stigler
}
\seealso{
\code{\link{selectSETAR}}
}

