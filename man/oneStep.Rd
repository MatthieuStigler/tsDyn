\name{oneStep}
\alias{oneStep}
\title{oneStep}
\usage{
oneStep(object, newdata, ...)
}
\arguments{
  \item{object}{fitted \sQuote{nlar} object}

  \item{newdata}{data from which to step forward}

  \item{\dots}{further arguments to be passed to and from
  other methods}
}
\value{
Computed value(s)
}
\description{
Doing one step forward within a \code{NLAR} model
}
\details{
If in \code{object} is encapsulated the \code{NLAR} map,
say, \code{F(X[t], X[t-d], ..., X[t-(m-1)d])}, this
function should return the value of \code{F} (with already
fitted parameters) applied to given new data, which can be
a single vector of length \code{m} or a matrix with
\code{m} columns.
}
\note{
This is an internal function, and should not be called by
the user
}
\examples{
tsDyn:::oneStep.linear
}
\author{
Antonio, Fabio Di Narzo
}
\keyword{internal}
\keyword{ts}

