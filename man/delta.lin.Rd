\name{delta.lin}
\alias{delta.lin}
\alias{delta.lin.test}
\title{delta test of linearity}
\usage{
  delta.lin(x, m, d = 1)

  delta.lin.test(x, m = 2:3, d = 1,
    eps = seq(0.5 * sd(x), 2 * sd(x), length = 4), B = 49)
}
\arguments{
  \item{x}{time series}

  \item{m}{vector of embedding dimensions}

  \item{d}{time delay}

  \item{eps}{vector of length scales}

  \item{B}{number of bootstrap replications}
}
\value{
  \code{delta.lin} returns the parametrically estimated
  delta statistic for the given time series (assuming
  linearity). \code{delta.lin.test} returns the bootstrap
  based 1-sided p-value. The test statistic is the
  difference between the parametric and nonparametric delta
  estimators.
}
\description{
  delta test of linearity based on conditional mutual
  information
}
\details{
  delta test of linearity based on conditional mutual
  information
}
\examples{
delta.lin(log10(lynx), m=3)
}
\author{
  Antonio, Fabio Di Narzo
}
\references{
  Sebastiano Manzan, Essays in Nonlinear Economic Dynamics,
  Thela Thesis (2003)
}
\keyword{ts}

