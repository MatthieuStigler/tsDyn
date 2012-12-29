\name{VARrep}
\alias{VARrep}
\alias{VARrep.VECM}
\alias{VARrep.VAR}

\title{VAR representation}

\description{Show the VAR representation of a VECM}

\usage{
\method{VARrep}{VECM}(object, ...)
\method{VARrep}{VAR}(object, ...)
}


\arguments{
\item{object}{An object of class \sQuote{VECM} created by \code{\link{VECM}}, or of class \sQuote{VAR} created by \code{\link{lineVar}}}
\item{\dots}{Currently not used}
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

# VAR model:
mod_var <- lineVar(barry, lag=2, I="diff")
VARrep(mod_var)

}

\keyword{ ts }
\keyword{VECM}
\keyword{VAR}
\keyword{cointegration}
