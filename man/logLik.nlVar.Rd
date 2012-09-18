\name{logLik}
\alias{logLik}
\alias{logLik.default}

\title{ Extract Log-Likelihood}
\description{
Log-Likelihood method for VAR models. 
}
\usage{
logLik(object, ...)

## S3 method for class 'VAR'
logLik(object, ...)

}
\arguments{
  \item{object}{ object of class \code{VAR} computed by \code{lineVar}.}
  \item{\dots}{ additional arguments to \code{logLik}.}
}
\details{
  The Log-Likelihood is computed as in Luetkepohl (2006) equ. 3.4.5 (p. 89) and Juselius (2006) p. 56:

  \deqn{
    LL = -(TN/2) \log(2\pi) - (T/2) \log|\Sigma| - (1/2)  \sum^{T} \left [ (y_t - A^{'}x_t)^{'} \Sigma^{-1} (y_t - A^{'}x_t) \right ]
  }
Where \eqn{\Sigma} is the Variance matrix of residuals, and \eqn{x_t} is the matrix stacking the regressors (lags and deterministic). 

(Note that Hamilton (1996) 11.1.10, p. 293 gives \eqn{+ (T/2) \log|\Sigma^{-1}|}, which is the same as \eqn{-(T/2) \log|\Sigma|)}. 
}

\value{
Log-Likelihood value. 
}

\references{
Hamilton (1996) \emph{Time Series Analysis}, Princeton University Press

Juselius (2006) \emph{The Cointegrated VAR model: methodology and Applications}, Oxford Univesity Press

Luetkepohl (2006) \emph{New Introduction  to Multiple Time Series Analysis}, Springer

}

\examples{
data(zeroyld)
data<-zeroyld

#Fit a VAR
VAR<-lineVar(data, lag=1)
logLik(VAR)
}
\author{Matthieu Stigler}
\keyword{ ts }
