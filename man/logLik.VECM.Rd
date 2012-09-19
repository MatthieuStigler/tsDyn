\name{logLik.VECM}

\alias{logLik.VECM}

\title{ Extract Log-Likelihood}
\description{
Log-Likelihood method for VECM models. 
}
\usage{
\method{logLik}{VECM}(object,r=1,...)
}
\arguments{
  \item{object}{ object of class \code{VECM} computed by \code{VECM}.}
  \item{r}{The cointegraing rank.}
  \item{\dots}{ additional arguments to \code{logLik}.}
}
\details{
  The Log-Likelihood is computed in two dfferent ways, depending on whether the \code{VECM} was estimated with ML (Johansen) or 2OLS (Engle and Granger). 

When the model is estimated with ML, the LL is computed as in Hamilton (1994) 20.2.10 (p. 637):

  \deqn{
    LL = -(TN/2) \log(2\pi) - (TN/2) -(T/2) \log|\hat\Sigma_{UU}| - (T/2)  \sum_{i=1}^{r} \log (1-\hat\lambda_{i})
  }
Where \eqn{\Sigma_{UU}} is the variance matrix of residuals from the first auxiliary regression, i.e. regressing \eqn{\Delta y_t} on a constant and lags, \eqn{\Delta y_{t-1}, \ldots, \Delta y_{t-p}}. \eqn{\lambda_{i}} are the eigenvalues from the \eqn{\Sigma_{VV}^{-1}\Sigma_{VU}\Sigma_{UU}^{-1}\Sigma_{UV}}, see 20.2.9 in Hamilton (1994).

When the model is estimated with 2OLS, the LL is computed as:
  \deqn{
    LL = \log|\Sigma|
  }

Where \eqn{\Sigma} is the variance matrix of residuals from the the VECM model. There is hence no correspondance between the LL from the VECM computed with 2OLS or  ML. 
}

\value{
Log-Likelihood value. 
}

\references{
Hamilton (1994) \emph{Time Series Analysis}, Princeton University Press

}

\examples{
data(zeroyld)
data<-zeroyld

#Fit a VAR
vecm<-VECM(data, lag=1,r=1, estim="ML")
logLik(vecm)
}
\author{Matthieu Stigler}
\keyword{ ts }
