\name{NNET}
\alias{NNET}
\alias{nnetTs}
\title{Neural Network nonlinear autoregressive model}
\description{
Neural Network nonlinear autoregressive model.
}
\usage{
nnetTs(x, m, d = 1, steps = d, series, size, 
	control = list(trace = FALSE))
}
\arguments{
	\item{x}{ time series }
	\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
	\item{series}{ time series name (optional) }
	\item{size}{number of hidden units in the neural network}
  \item{control}{control list to be passed to \code{nnet::nnet} optimizer}
}

\details{
  Neural network model with 1 hidden layer and linear output:

  \deqn{x_{t+s} = \beta_0 + \sum_{j=1}^D \beta_j g( \gamma_{0j} +
    \sum_{i=1}^{m} \gamma_{ij} x_{t-(i-1) d} )}{
    x[t+steps] = beta[0] + sum_j beta[j] g( gamma[0,j] +
    sum_i gamma[i,j] x[t-(i-1) d] )
  }

Model is estimated using the \code{\link[nnet]{nnet}} function in \pkg{nnet}
package. Optimization is done via the BFGS method of
\code{\link{optim}}. Note that for this model, no additional
model-specific summary and plot methods are made available from this package.
}
\value{
  An object of class \code{nlar}, subclass \code{nnetTs}, i.e. a list
  with mostly \code{nnet::nnet} internal structures.
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit a Neural Network model
mod.nnet <- nnetTs(log(lynx), m=2, size=3)
mod.nnet
}
\keyword{ ts }
\references{
Non-linear time series models in empirical finance, Philip Hans Franses and Dick van Dijk, Cambridge: Cambridge University Press (2000).

Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford: Oxford University Press (1990).

Chaos: A Statistical Perspective, Chan, K., Tong, H., New York: Springer-Verlag (2001).
}
