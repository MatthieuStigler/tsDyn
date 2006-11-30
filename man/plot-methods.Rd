\name{plot methods}
\alias{plot-methods}
\alias{plot.setar}
\alias{plot.lstar}
\title{Plotting methods for setar and lstar subclasses}
\description{
Plotting methods \sQuote{setar} and \sQuote{lstar} subclasses
}
\usage{
\method{plot}{setar}(x, ask=interactive(), legend=FALSE, regSwStart, regSwStop, \dots)
\method{plot}{lstar}(x, ask=interactive(), legend=FALSE, regSwStart, regSwStop, \dots)
}

\arguments{
  \item{x}{ fitted \sQuote{setar} or \sQuote{lstar} object }
  \item{ask}{ graphical option. See \code{\link{par}}}
  \item{legend}{ Should a legend be plotted? (logical)}
  \item{regSwStart, regSwStop}{ optional starting and stopping time indices for regime
  switching plot }
  \item{\dots}{ further arguments to be passed to and from other methods }
}
\details{
These plot methods produce a plot which gives to you an idea of the behaviour of the fitted model.

Firstly, if embedding dimension is, say, \var{m}, \var{m} scatterplots are
produced. On the x axis you have the lagged time series values. On the y
axis the \sQuote{response} time series values. Observed points are
represented with different colors-symbols depending on the level of the
threshold variable. Specifically, for the setar model, black means
\sQuote{low regime}, red means \sQuote{high regime}. For the lstar
model, where the self-threshold variable is continuous, threshold values
are grouped in 5 different zones with the same number of points in
each. Note that if more than 300 points are to be plotted, they all
share the same symbol, and regimes can be distinguished only by
color. If you want, by specifying \code{legend=TRUE} a legend is added
at the upper-left corner of each scatterplot. To each scatterplot, a
dashed line is superposed,
which links subsequent fitted values.

Finally, a new time series plot is produced, with lines segments
coloured depending on the regime (colors meanings are the same of those
in the preceedings scatterplots). Optionally, you can specify a starting
and ending time indices, for zooming on a particular segment of the
time series.
}
\author{ Antonio, Fabio Di Narzo }
\examples{
##
##See 'setar' examples
##
}
\seealso{
\code{\link{setar}}, \code{\link{lstar}}

\link{nlar-methods} for other generic available methods for this kind of objects.
}
\keyword{ ts }
