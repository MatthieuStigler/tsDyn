\name{TARCH}
\alias{TARCH}
\alias{tarch}
\alias{summary.tarch}
\alias{print.summary.tarch}
\title{Treshold-ARCH model}
\description{
  Treshold AutoRegressive Conditionally Heteroschedastic model
}
\usage{
  tarch(x, m, d=1, steps=d, series, coef, thDelay=0, control=list(), ...)
}
\value{
  An object of class \code{tarch}.
}
\arguments{
	\item{x}{ time series }
	\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
	\item{series}{ time series name (optional) }
  \item{coef}{vector of starting coefficients values. If missing, they
    are randomly generated from the log-normal distribution}
  \item{thDelay}{ time delay value for thresholding }
  \item{control, ...}{additional parameters to be passed to \code{\link{optim}}}
}
\details{
  Treshold-ARCH model:
  \deqn{x_t = \sigma_t \epsilon_t}{x[t] = sigma[t] eps[t]}
  with \eqn{\epsilon_t}{eps[t]} standard white noise, and \eqn{\sigma_t}{sigma[t]} conditional
  standard deviation which takes the form:
  \deqn{
    \sigma^2_{t+s} = [b_{0,0} + \sum_{j=1}^m b_{0,j} \sigma^2_{t-(j-1)d}]
    I(Z_t \leq 0) + [b_{1,0} + \sum_{j=1}^m b_{1,j} \sigma^2_{t-(j-1)d}] I(Z_t > 0)
  }{
    sigma2[t+steps] = ( b[0,0] + sum_j b[0,j] sigma2[t-(j-1)d] ) * (Z[t] <= 0)
    + ( b[1,0] + sum_j b[1,j] sigma2[t-(j-1)d] ) * (Z[t] > 0)
  }

  and \eqn{Z_t}{Z[t]} threshold variable defined as
  \eqn{Z_t = x_{t-thD\cdot d}}{Z[t] = x[t-thDelay*d]}.
  The model is estimated by Conditional Maximum Likelihood, with
  positivity of parameters restriction (strict for
  \eqn{b_{0,0}}{b[0,0]} and \eqn{b_{1,0}}{b[1,0]}), using the L-BFGS-B
  provided by the \code{\link{optim}} function.

  Standard errors provided in the summary are asymptoticals.

  No model specific plots are produced by the \code{plot} method.
}
\seealso{
  \code{\link{setar}}, \code{\link{lstar}}
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#
#Taken from tseries::garch man page
#
n <- 1100
a <- c(0.1, 0.5, 0.2)  # ARCH(2) coefficients
e <- rnorm(n)
x <- double(n)
x[1:2] <- rnorm(2, sd = sqrt(a[1]/(1.0-a[2]-a[3])))
for(i in 3:n)  # Generate ARCH(2) process
{
   x[i] <- e[i]*sqrt(a[1]+a[2]*x[i-1]^2+a[3]*x[i-2]^2)
}
x <- ts(x[101:1100])

x.tarch <- tarch(x, m=2)
summary(x.tarch)
}
\keyword{ ts }
\references{
  Threshold Arch Models and asymmetries in volatility, R. Rabemanajara
  and J. M. Zakoian, Journal of Applied Econometrics, vol. 8 (1993)

  Threshold heteroschedastic models, J. M. Zakoian, D. P. INSEE (1991)
}
