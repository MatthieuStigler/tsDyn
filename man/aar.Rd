\name{aar}
\alias{AAR}
\alias{aar}
\alias{summary.aar}
\alias{print.aar}
\alias{plot.aar}
\title{Additive nonlinear autoregressive model}
\description{
Additive nonlinear autoregressive model.
}
\usage{
aar(x, m, d=1, steps=d, series)
}
\arguments{
	\item{x}{ time series }
	\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
	\item{series}{ time series name (optional) }
}

\details{
  Nonparametric additive autoregressive model of the form:
  \deqn{
        x_{t+s} = \mu + \sum_{j=1}^{m} s_j(x_{t-(j-1)d})
  }{
        x[t+steps] = mu + sum_j s_j(x[t-(j-1)d])
  }

  where \eqn{s_j} are nonparametric univariate functions of lagged time
  series values. They are represented by cubic regression splines.
  \eqn{s_j} are estimated together with their level of
  smoothing using routines in the \pkg{mgcv} package (see references).
}
\value{
  An object of class \code{nlar}, subclass \code{aar}, i.e. a list
  with mostly internal structures for the fitted \code{\link[mgcv]{gam}} object.
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit an AAR model:
mod <- aar(log(lynx), m=3)
#Summary informations:
summary(mod)
#Diagnostic plots:
plot(mod)
}
\keyword{ ts }
\references{
Wood, mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
(2001)

Wood and Augustin, GAMs with integrated model selection using penalized regression splines and applications to environmental modelling. Ecological Modelling 157:157-177 (2002)
}
