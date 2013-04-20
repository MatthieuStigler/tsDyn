\name{LSTAR}
\alias{LSTAR}
\alias{lstar}
\title{Logistic Smooth Transition AutoRegressive model}
\description{
Logistic Smooth Transition AutoRegressive model.
}
\usage{
lstar(x, m, d=1, steps=d, series, mL, mH, mTh, thDelay, 
          thVar, th, gamma, trace=TRUE, include = c("const", "trend","none", "both"), control=list(), starting.control=list())
}
\value{
  An object of class \code{nlar}, subclass \code{lstar}, i.e. a list
  with fitted model informations.
%  Among others:
%  \describe{
%    \item{mL,mH}{low and high regimes autoregressive orders}
%    \item{externThVar}{is the threshold variable external?}
%    \item{thVar}{threshold variable values}
%  }
%  If the thresholding variable isn't external, there is an additional
%  component \code{mTh}, containing coefficients for lagged time series
%  used for obtaining the threshold variable.
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{mL}{autoregressive order for 'low' regime (default: m). Must be <=m}
\item{mH}{autoregressive order for 'high' regime (default: m). Must be
  <=m}
\item{thDelay}{'time delay' for the threshold variable (as multiple of
  embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the
  threshold variable}
\item{thVar}{external threshold variable}
\item{th, gamma}{starting values for coefficients in the LSTAR model. If missing, a grid search is performed}
\item{trace}{should additional infos be printed? (logical)}
\item{include}{Type of deterministic regressors to include}
\item{control}{further arguments to be passed as \code{control} list to
  \code{\link{optim}}}
\item{starting.control}{further arguments for the grid search (dimension, bounds). See details below.  }
}
\details{
  \deqn{ x_{t+s} = ( \phi_{1,0} + \phi_{1,1} x_t + \phi_{1,2} x_{t-d} + \dots +
  \phi_{1,mL} x_{t - (mL-1)d} ) G( z_t, th, \gamma ) +
  ( \phi_{2,0} + \phi_{2,1} x_t + \phi_{2,2} x_{t-d} + \dots + \phi_{2,mH}
  x_{t - (mH-1)d} ) (1 - G( z_t, th, \gamma ) ) + \epsilon_{t+steps}}{
  x[t+steps] = ( phi1[0] + phi1[1] x[t] + phi1[2] x[t-d] + \dots +
  phi1[mL] x[t - (mL-1)d] ) G( z[t], th, gamma )
  + ( phi2[0] + phi2[1] x[t] + phi2[2] x[t-d] + \dots + phi2[mH] x[t -
  (mH-1)d] ) (1 - G( z[t], th, gamma ) ) + eps[t+steps]
  }
with \var{z} the threshold variable, and \eqn{G} the logistic function,
computed as \code{plogis(q, location = th, scale = 1/gamma)}, so see
\code{\link{plogis}} documentation for details on the logistic function
formulation and parameters meanings.  
The threshold variable can alternatively be specified by:
\describe{
\item{mTh}{ \eqn{z[t] = x[t] mTh[1] + x[t-d] mTh[2] + \dots +
    x[t-(m-1)d] mTh[m]} } 
\item{thDelay}{ \eqn{z[t] = x[t - thDelay*d ]} }
\item{thVar}{ \eqn{z[t] = thVar[t]} }
}

Note that if starting values for phi1 and phi2 are provided, isn't
necessary to specify \code{mL} and \code{mH}. Further, the user has to specify only
one parameter between \code{mTh}, \code{thDelay} and \code{thVar} for indicating the
threshold variable. 

Estimation of the transition parameters \var{th} and \var{gamma}, as well as the regression parameters 
\var{phi1} and \var{phi2}, is done using concentrated least squares, as suggested in \cite{Leybourne et al. (1996)}. 

Given \var{th} and \var{gamma}, the model is linear, so regression coefficients
can be obtained as usual by OLS. So the nonlinear numerical search needs only to be done for \var{th} and \var{gamma}; 
the regression parameters are then recovered by OLS again from the optimal \var{th} and \var{gamma}. 

For the nonlinear estimation of the
parameters \var{th} and \var{gamma}, the program uses the
\code{\link{optim}} function, with optimization method BFGS using the analytical gradient. 
For the estimation of standard values, \code{\link{optim}} is re-run
using the complete Least Squares objective function, and the standard errors are obtained by inverting the hessian. 
You can pass further arguments to \link{optim} directly with the \code{control} list argument. For instance, the option \code{maxit} maybe useful when
there are convergence issues (see examples).

Starting parameters are obtained doing a simple two-dimensional grid-search over \var{th} and \var{gamma}. 
Parameters of the grid (interval for the values, dimension of the grid) can be passed to \code{starting.control}. 

  \describe{
    \item{\code{nTh}}{The number of threshold values (\var{th}) in the grid. Defaults to 200}
    \item{\code{nGamma}}{The number of smoothing values (\var{gamma}) in the grid. Defaults to 40}
    \item{\code{trim}}{The minimal percentage of observations in each regime. Defaults to 10\% (possible threshold values are between the 0.1 and 0.9 quantile)}
    \item{\code{gammaInt}}{The lower and higher smoothing values of the grid. Defaults to c(1,40)}
    \item{\code{thInt}}{The lower and higher threshold values of the grid. When not specified (default, i.e NA), the interval are the \code{trim} quantiles above. }
  }
}

\seealso{
\code{\link{plot.lstar}} for details on plots produced for this model
from the \code{plot} generic. 
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit a LSTAR model. Note 'maxit': slow convergence
mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar

#fit a LSTAR model without a constant in both regimes. 
mod.lstar2 <- lstar(log10(lynx), m=1,  include="none")
mod.lstar2

#Note in example below that the initial grid search seems to be to narrow. Extend it, and evaluate more values (slow!):
mod.lstar3 <- lstar(log10(lynx), m=1,  include="none", starting.control=list(gammaInt=c(1,2000), nGamma=50))
mod.lstar3

# a few methods for lstar:
summary(mod.lstar)
residuals(mod.lstar)
AIC(mod.lstar)
BIC(mod.lstar)
plot(mod.lstar)
predict(mod.lstar, n.ahead=5)
}
\keyword{ ts }
\references{
Non-linear time series models in empirical finance, Philip Hans Franses
and Dick van Dijk, Cambridge: Cambridge University Press (2000). 

Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford:
Oxford University Press (1990). 

Leybourne, S., Newbold, P., Vougas, D. (1998) Unit roots and smooth transitions, Journal of Time Series Analysis, 19: 83-97

}
