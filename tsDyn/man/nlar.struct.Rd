\name{nlar.struct}
\alias{nlar.struct}
\alias{nlar.struct-methods}
\alias{getXXYY}
\alias{getXX}
\alias{getYY}
\alias{getNUsed}
\title{NLAR common structure}
\description{
NLAR common structure
}
\usage{
nlar.struct(x, m, d=1, steps=d, series)
getXXYY(obj, ...)
getXX(obj, ...)
getYY(obj, ...)
getNUsed(obj, ...)
}
\value{
	\code{nlar.struct} builds an object of class \code{nlar.struct}, from which common NLAR models infos can be extracted using:
	\item{getXXYY}{design matrix}
	\item{getXX}{regressor matrix}
	\item{getYY}{response matrix}
	\item{getNUsed}{original time series length}
}
\arguments{
\item{x}{time series}
\item{m,d,steps}{embedding dimension, time delay and forecasting steps}
\item{series}{(optional) time series name}
\item{obj}{\code{nlar.struct} object}
\item{...}{arguments to be passed to and from other methods}
}
\author{ Antonio, Fabio Di Narzo }
\keyword{ ts }
\keyword{ internal }
