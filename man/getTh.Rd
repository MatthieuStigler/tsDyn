\name{getTh}
\alias{getTh}
\alias{getTh.default}
\title{ Extract threshold(s) coefficient }
\description{
Extract threshold coefficient(s)
}
\usage{
getTh(object, ...)

## Default S3 method:
getTh(object, ...)
}
\arguments{
  \item{object}{ object of class \code{setar}, \code{summary.setar}, \code{nlVar} }
  \item{\dots}{ additional arguments to \code{getTh} }
}
\value{
Threshold value. 
}
\examples{
set<-setar(lynx, m=3)
getTh(set)
getTh(summary(set))
}
\author{Matthieu Stigler}
\keyword{ ts }
