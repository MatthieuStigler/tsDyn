\name{getTh}
\alias{getTh}
\alias{getTh.default}
\title{Extract threshold(s) coefficient}
\usage{
getTh(object, ...)

\method{getTh}{default}(object, ...)
}
\arguments{
  \item{object}{object of class \code{setar},
  \code{summary.setar}, \code{nlVar}}

  \item{\dots}{additional arguments to \code{getTh}}
}
\value{
Threshold value.
}
\description{
Extract threshold coefficient(s)
}
\examples{
set<-setar(lynx, m=3)
getTh(set)
getTh(summary(set))
}
\author{
Matthieu Stigler
}
\keyword{ts}

