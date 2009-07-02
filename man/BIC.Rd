\name{BIC}
\alias{BIC}
\alias{bic}
\title{Bayesian (Schwarz) information criterion}
\description{
Function calling AIC witk k= log(n) as default instead of k=2. 
}
\usage{
BIC(object, ...,k)
}
\arguments{
	\item{object}{ object of class \code{lm} or \code{nlar}}
	\item{...}{further arguments passed to \code{AIC}}
	\item{k}{penalty term}
}

\details{
This is just a direct function to use the BIC criterion from the AIC. The penalty term is: k=log(N). N is the number of observations in the initial series and hence does not correspond to the number of fitted observations with time series due to the lagged ariables. }
\value{
Numeric value with the corresponding BIC
}
\author{ Matthieu Stigler }
\examples{

mod.linear<-linear(log10(lynx), m=2)
mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)

AIC(mod.linear)
AIC(mod.setar)

BIC(mod.linear)
BIC(mod.setar)


}
\keyword{ ts }

