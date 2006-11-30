\name{delta}
\alias{delta}
\alias{delta.test}
\title{delta test of conditional indipendence}
\description{
delta statistic of conditional indipendence and associated bootstrap test
}
\usage{
delta(x, m, d=1, eps)
delta.test(x, m=2:3, d=1, eps=seq(0.5*sd(x),2*sd(x),length=4), B=49)
}
\arguments{
  \item{x}{time series}
  \item{m}{vector of embedding dimensions}
  \item{d}{time delay}
  \item{eps}{vector of length scales}
  \item{B}{number of bootstrap replications}
}
\details{
delta statistic of conditional indipendence and associated bootstrap test. For details, see Manzan(2003).
}
\section{Warning}{
Results are sensible to the choice of the window \code{eps}. So, try the test for a grid of \code{m} and \code{eps} values. Also, be aware of the course of dimensionality: m can't be too high for relatively small time series. See references for further details.
}
\value{
\code{delta} returns the computed delta statistic. \code{delta.test} returns the bootstrap based 1-sided p-value.
}
\examples{
delta(log10(lynx), m=3, eps=sd(log10(lynx)))
}
\seealso{
BDS marginal indipendence test: \code{\link[tseries]{bds.test}} in package \code{tseries}

Teraesvirta's neural network test for nonlinearity: \code{\link[tseries]{terasvirta.test}} in package \code{tseries}

delta test for nonlinearity: \code{\link{delta.lin.test}}
}
\references{
Sebastiano Manzan, Essays in Nonlinear Economic Dynamics, Thela Thesis (2003)
}
\author{Antonio, Fabio Di Narzo}
\keyword{ts}
