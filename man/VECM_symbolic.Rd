\name{VECM_symbolic}
\alias{VECM_symbolic}
\title{Virtua VECM model}
\description{Pedagogical tool to create a symbolic VECM model, i.e. just for representation purpose.}
\usage{
VECM_symbolic(alpha, beta, lags, inc, include = c("none", "const", "trend", "both"))
}


\arguments{
\item{alpha}{Matrix of alpha speed adjustment coefficients.}
\item{beta}{Matrix of alpha,  cointegrating coefficients.}
\item{lags}{Matrix containg the lags coefficients.}
\item{inc}{Matrix containg the include (see following arg) coefficients.}
\item{include}{Character indicating the type of deterministic term included, if any. }
}

\value{An object of class \sQuote{VECM}, without however any data. 
}



\examples{

  a<-matrix(c(-0.4, 0.1), ncol=1)
  b<-matrix(c(1, -2), ncol=2)

  # VECM_symb(alpha=a, beta=t(b))
  d<- VECM_symbolic(alpha=a, beta=t(b))
  VARrep(d)
  d<- VECM_symbolic(alpha=a, beta=t(b), lags=matrix(0, ncol=2, nrow=2))
  VARrep(d)
  d3<- VECM_symbolic(alpha=a, beta=t(b), lags=matrix(c(0.1, 0.3, 0.1, 0.2), ncol=2, nrow=2), inc=matrix(c(0.5, 0.1), ncol=1), include="const")
  VARrep(d3)

}

\keyword{ ts }
\keyword{VECM}
\keyword{cointegration}
