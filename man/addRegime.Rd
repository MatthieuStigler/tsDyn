\name{addRegime}
\alias{addRegime}
\alias{addRegime.star}
\title{addRegime test}
\description{
addRegime test
}
\usage{
addRegime(object, ...)
\method{addRegime}{star}(object, G, rob=FALSE, sig=0.05, trace = TRUE, ...)
}
\arguments{
\item{object}{ fitted model object with at least 2 regimes }
\item{G}{ gradient matrix as computed by \code{\link{computeGradient}} }
\item{rob}{boolean indicating if robust tests should be used or not}
\item{sig}{significance level of the tests}
\item{trace}{should additional infos be printed?}
\item{...}{arguments to and from other methods}
}
\value{
A list containing the p-value of the F statistic and a boolean,
true if there is some remaining nonlinearity and false otherwise.
}
\seealso{\code{\link{star}}}
\author{ J. L. Aznarte }
\examples{
##TODO
}
\keyword{ ts }
\references{
TODO
}
