\name{autopairs}
\alias{autopairs}
\title{Bivariate time series plots}
\description{
Bivariate time series plots: scatterplots, directed lines and kernel density estimations using functions in the \pkg{sm} package.
}
\usage{
autopairs(x, lag=1, h,
type=c("levels","persp","image","lines","points","regression"),
GUI=interactive())
}
\arguments{
\item{x}{ time series}
\item{lag}{ time lag }
\item{h}{ kernel window (useful only for kernel estimations) }
\item{type}{ type of plot: contour levels, perspective plots, image, directed lines, points or points with superposed kernel regression }
\item{GUI}{ should a GUI be displayed? }
}
\details{
  Bivariate time series plots: scatterplots, directed lines and kernel
  density and regression functions estimations using functions in the
  package \pkg{sm}.
  In particular, for kernel density estimation
  \code{\link[sm]{sm.density}} is used, with smoothing parameter \code{h}
  defaulting to \code{\link[sm]{hnorm}}. For kernel regression, \code{\link[sm]{sm.regression}} is used.

  If \code{GUI==TRUE}, a simple graphical user interface is displayed
  to control graphical parameters.
}
\value{
None. Plots are produced on the default graphical device.
}
\seealso{
For finer control on density estimation, consider using directly 
\code{\link[sm]{sm.density}} and, especially, 
\code{\link[sm]{sm.ts.pdf}} from package \pkg{sm}.
}
\examples{
x <- log10(lynx)
autopairs(x, lag=2, type="lines")
}
\author{ Wrappers to \pkg{sm} and GUI by Antonio, Fabio Di Narzo }
\keyword{ ts }
