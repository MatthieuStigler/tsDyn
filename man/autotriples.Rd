\name{autotriples}
\alias{autotriples}
\title{Trivariate time series plots}
\description{
Trivariate time series plots: kernel autoregression using functions in the 'sm' package
}
\usage{
autotriples(x, lags=1:2, h,
type=c("levels","persp","image", "lines", "points"),
GUI=interactive())
}
\arguments{
\item{x}{ time series}
\item{lags}{ vector of regressors lags }
\item{h}{ kernel window }
\item{type}{ type of plot: countour levels, perspective plots, image }
\item{GUI}{ should a GUI be displayed? }
}
\details{
  This function displays trivariate time series plots, i.e. kernel
  regression of \eqn{x[t-lags[1]], x[t-lags[2]]}{x_{t-l_1}, x_{t-l_2}}
  against \eqn{x[t]}{x_t} using functions in the package 'sm'.
  In particular, \code{\link[sm]{sm.regression}} is used, with smoothing
  parameter  defaulting to \code{\link[sm]{hnorm}(x)}.
  If requested, a simple GUI is displayed, to change interactively
  functions parameters and watching corresponding outputs.
}
\value{
None. Plots are produced on the default graphical device.
}
\seealso{
For finer control on kernel regression, consider using directly \code{\link[sm]{sm.regression}} and, especially, \code{\link[sm]{sm.autoregression}} in package \code{\link[sm]{sm}}.
}
\examples{
autotriples(log(lynx))
autotriples(log(lynx), type="persp")
autotriples(log(lynx), type="image")
}
\author{ Wrappers to 'sm' and GUI by Antonio, Fabio Di Narzo  }
\keyword{ ts }
