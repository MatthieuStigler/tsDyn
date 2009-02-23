\name{resVar}
\alias{resVar}
\title{Residual variance}
\description{Extracts the global and regime-dependant variance of the residuals
}
\usage{
resVar(x, adj=c("OLS", "ML")
}


\arguments{
\item{x}{setar object}
\item{adj}{Degrees of freedom adjustment for the variance}
}

\value{A vector containing:
\item{Total}{The residual variance of the full sample}
\item{L, (M), H}{The residual variance of the lower (L), middle (if two thresholds) (M) and higher (H) regimes}
}

\details{
The degree of freedom adjustment in the formula for the variance is the number of parameters when \code{adj="OLS"} or zero when \code{adj="ML"}. 

}
\references{
Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford: Oxford University Press (1990).
}

\author{Matthieu Stigler}
\examples{
#Lynx model as in Tong (1980, p. 387)
mod.setar <- setar(log10(lynx), mL=7,mH=2, thDelay=1, th=3.116)
summary(mod.setar)
#coefficients are same for lower regime but differ for higer

resVar(mod.setar, adj="ML")
#variance or the residuals is same for lower regime but differ for higer regime and hence for total

#Lynx model as in Tong (1980, p. 405)
mod.setar2 <- setar(log10(lynx), mL=1,mM=7,mH=2, thDelay=1, nthresh=2,th=c(2.373, 3.154))
round(coefficients(mod.setar2),3)

resVar(mod.setar2, adj="ML")


}

\keyword{ ts }
