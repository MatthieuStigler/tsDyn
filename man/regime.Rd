\name{regime}
\alias{regime}
\alias{regime.default}
\title{ Extract variable showing regime }
\description{
This function allows to extract the indicator variable specifying the regime in which the process is at time t. 
}
\usage{
regime(object, ...)

## Default S3 method:
regime(object, ...)
}
\arguments{
  \item{object}{ object of class \code{setar} }
}
\value{
Time series of same attributes as input to setar. 
}
\examples{
set<-setar(lynx, m=3)
regime(set)
}
\author{Matthieu Stigler}
\keyword{ ts }
