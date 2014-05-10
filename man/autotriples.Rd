\name{autotriples}
\alias{autotriples}
\title{Trivariate time series plots}
\usage{
autotriples(x, lags = 1:2, h, type = c("levels", "persp", "image", "lines",
  "points"))
}
\arguments{
  \item{x}{time series}

  \item{lags}{vector of regressors lags}

  \item{h}{kernel window}

  \item{type}{type of plot: contour levels, perspective
  plots, image}
}
\value{
None. Plots are produced on the default graphical device.
}
\description{
Trivariate time series plots: kernel autoregression using
functions in the \pkg{sm} package
}
\details{
This function displays trivariate time series plots, i.e.
kernel regression of \eqn{x[t-lags[1]],
x[t-lags[2]]}{x_{t-l_1}, x_{t-l_2}} against \eqn{x[t]}{x_t}
using functions in the package \pkg{sm}.  In particular,
\code{\link[sm]{sm.regression}} is used, with smoothing
parameter defaulting to \code{\link[sm]{hnorm}(x)}.
}
\examples{
autotriples(log(lynx))
autotriples(log(lynx), type="persp")
autotriples(log(lynx), type="image")
}
\author{
Wrappers to \pkg{sm} by Antonio, Fabio Di Narzo
}
\seealso{
For finer control on kernel regression, consider using
directly \code{\link[sm]{sm.regression}} and, especially,
\code{\link[sm]{sm.autoregression}} in package
\code{\link[sm]{sm}}.
}
\keyword{ts}

