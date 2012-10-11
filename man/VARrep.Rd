\name{VARrep}
\alias{VARrep}
\title{VAR representation}
\description{Show the VAR representation of a VECM}
\usage{
VARrep(vecm)
}


\arguments{
\item{vecm}{An object of class \sQuote{VECM} created by \code{\link{VECM}}}
}

\value{A matrix containing the parameters of the VECM under their VAR representation. 
}


\references{
Hamilton (1994) \emph{Time Series Analysis}, Princeton University Press
}

\author{Matthieu Stigler}

\examples{

data(barry)

# VECM model:
mod_vecm <- VECM(barry, lag=2, estim="ML")
VARrep(mod_vecm)

}

\keyword{ ts }
\keyword{VECM}
\keyword{VAR}
\keyword{cointegration}
