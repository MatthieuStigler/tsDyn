\name{llar}
\alias{llar}
\alias{llar.predict}
\alias{llar.fitted}
\alias{print.llar}
\alias{as.data.frame.llar}
\alias{plot.llar}
\title{Locally linear model}
\description{
Casdagli test of nonlinearity via locally linear forecasts
}
\usage{
llar(x, m, d = 1, steps = d, series, eps.min = sd(x)/2,
	eps.max = diff(range(x)), neps = 30, trace = 0)

llar.predict(x, m, d=1, steps=d, series, n.ahead=1, 
eps=stop("you must specify a window value"),
onvoid=c("fail","enlarge"), r = 20, trace=1)

llar.fitted(x, m, d=1, steps=d, series, eps, trace=0)
}

\arguments{
  \item{x}{ time series }
	\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
	\item{series}{time series name (optional)}
  \item{n.ahead}{ n. of steps ahead to forecast }
  \item{eps.min, eps.max}{ min and max neighborhood size }
  \item{neps}{number of neghborhood levels along which iterate}
  \item{eps}{neighborhood size }
  \item{onvoid}{what to do in case of an isolated point: stop or enlarge neghborhood size by an r\%}
  \item{r}{if an isolated point is found, enlarge neighborhood window by r\%}
  \item{trace}{tracing level: 0, 1 or more than 1 for \code{llar}, 0 or
    1 for \code{llar.forecast}}
}
\details{
  \code{llar} does the Casdagli test of non-linearity. Given the embedding state-space (of dimension \code{m} and time delay \code{d}) obtained from time series \code{series}, for a sequence of distance values \code{eps}, the relative error made by forecasting time series values with a linear autoregressive model estimated on points closer than \code{eps} is computed. 
If minimum error is reached at relatively small length scales, a global linear model may be inappropriate (using current embedding parameters).
This was suggested by Casdagli(1991) as a test for non-linearity.

\code{llar.predict} tries to extend the given time series by
\code{n.ahead} points by iteratively 
fitting locally (in the embedding space of dimension m and time delay d)
a linear model. If the spatial neighborhood window is too small, your
time series last point would be probably isolated. You can ask to
automatically enlarge the window \code{eps} by a factor of r\%
sequentially, until enough neighbours are found for fitting the linear
model.

\code{llar.fitted} gives out-of-sample fitted values from locally linear
models.
}
\note{
The C implementation was re-adapted from that in the TISEAN package ("ll-ar" routine, see references). However, here the euclidean norm is used, in place of the max-norm.
}
\section{Warning}{
For long time series, this can be slow, especially for relatively big neighborhood sizes.
}
\value{
\code{llar} gives an object of class 'llar'. I.e., a list of components: 
\tabular{ll}{
 RMSE \tab vector of relative errors \cr
 eps \tab vector of neighborhood sizes (in the same order of RMSE) \cr
 frac \tab vector of fractions of the time series used for RMSE computation \cr
 avfound \tab vector of average number of neighbours for each point in the time series
}

which can be plotted using the \code{plot} method, and transformed to a regular 
\code{data.frame} with the \code{\link{as.data.frame}} function.

\code{llar.forecast} gives the vector of n steps ahead locally linear iterated
forecasts.

\code{llar.fitted} gives out-of-sample fitted values from locally linear
models.
}

\author{ Antonio, Fabio Di Narzo }
\examples{
res <- llar(log(lynx), m=3, neps=7)
plot(res)

x.new <-  llar.predict(log(lynx),n.ahead=100, m=3, eps=1,
onvoid="enlarge", r=5)
lag.plot(x.new, labels=FALSE)

x.fitted <-  llar.fitted(log(lynx), m=3, eps=1)
lag.plot(x.fitted, labels=FALSE)
}
\references{
M. Casdagli, Chaos and deterministic versus stochastic nonlinear modelling, J. Roy. Stat. Soc. 54, 303 (1991)

Hegger, R., Kantz, H., Schreiber, T., Practical implementation of nonlinear time series methods: The TISEAN package; CHAOS 9, 413-435 (1999)
}
\keyword{ ts }
