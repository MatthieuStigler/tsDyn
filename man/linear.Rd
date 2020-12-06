\name{LINEAR}
\alias{LINEAR}
\alias{linear}
\alias{print.linear}
\alias{summary.linear}
\alias{print.summary.linear}
\title{Linear AutoRegressive models}
\description{
AR(m) model
}
\usage{
	linear(x, m, d=1, steps=d, series, include = c( "const", "trend","none", "both"), 
        type=c("level", "diff", "ADF"), warn_root=TRUE)
}
\value{
	A \code{\link{nlar}} object, \code{linear} subclass.
}
\arguments{
	\item{x}{ time series }
	\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
	\item{series}{ time series name (optional) }
	\item{include}{Type of deterministic regressors to include}
	\item{type}{Whether the variable is taken is level, difference or a mix (diff y= y-1, diff lags) as in the ADF test}
	\item{warn_root}{Whether to check (and warn) for roots outside the unit circle?}
}
\details{
AR(m) model:
\deqn{  x_{t+s} = \phi_0 + \phi_1 x_t + \phi_2 x_{t-d} + \dots + \phi_m x_{t -
    (m-1)d} + \epsilon_{t+s}}{
  x[t+steps] = phi[0] + phi[1] x[t] + phi[2] x[t-d] + \dots + phi[m]
  x[t - (m-1)d] + eps[t+steps]}
}
\seealso{
\code{\link{nlar}} for fitting this and other models to time series data
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit an AR(2) model
mod.linear <- linear(log(lynx), m=2)
mod.linear
summary(mod.linear)
}
\keyword{ ts }
