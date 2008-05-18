\name{SETAR}
\alias{SETAR}
\alias{setar}
\alias{summary.setar}
\title{Self Threshold Autoregressive model}
\description{
Self Exciting Threshold AutoRegressive model.
}
\synopsis{
setar(x, m, d=1, steps=d, series, mL, mM,mH, thDelay=0, mTh, thVar, th, trace=FALSE, nested=FALSE, include = c( "const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL), MM=seq_len(mM), MH=seq_len(mH),nthresh=1,trim=0.15)
}
\usage{
setar(x, m, d=1, steps=d, series, mL=m, mH=m, thDelay=0, th,
trace=FALSE,include = c( "const", "trend","none", "both"),
common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL), MH=seq_len(mH),nthresh=1,trim=0.15)

setar(x, m, d=1, steps=d, series, mL=m, mH=m, mTh, th,
trace=FALSE,include = c( "const", "trend","none", "both"), 
common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL), MH=seq_len(mH),nthresh=1,trim=0.15)

setar(x, m, d=1, steps=d, series, mL=m, mH=m, thVar, th,
trace=FALSE,include = c( "const", "trend","none", "both"), 
common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL), MH=seq_len(mH),nthresh=1,trim=0.15)
}
\value{
  An object of class \code{nlar}, subclass \code{setar}
%, i.e. a list  with mostly internal structures.
%  Among others:
%  \describe{
%    \item{fixedTh}{boolean. Was the threshold value fixed?}
%    \item{mL,mH}{low and high regimes autoregressive orders}
%    \item{externThVar}{is the threshold variable external?}
%    \item{thVar}{threshold variable values}
%    \item{lowRegProp}{proportion of points in low regime}
%  }
%  If the thresholding variable isn't external, there is an additional
%  component \code{mTh}, containing coefficients for lagged time series
%  used for obtaining the treshold variable. If \code{th} coefficient was estimated,
%  there is also a \code{th} component, which is a vector containing all
%  the tried threshold values.
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{mL, mH,mM}{autoregressive order for \sQuote{low} (mL) \sQuote{middle} (mM, only useful if nthresh=2) and \sQuote{high} (mH)regime (default values: m). Must be <=m. Alternatively, you can specify \code{ML}}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the threshold variable}
\item{thVar}{external threshold variable}
\item{th}{threshold value (if missing, a search over a resonable grid is tried)}
\item{trace}{should additional infos be printed? (logical)}
\item{include} {Type of deterministic regressors to include}
\item{common}{Whether the deterministic regressors are regime specific (common=FALSE) or not.}
\item{ML,MM,MH}{vector of lags for  order for \sQuote{low} (ML) \sQuote{middle} (MM, only useful if nthresh=2) and \sQuote{high} (MH)regime. Max must be <=m}
\item{model}{Currently not implemented}
\item{nthresh}{Number of threshold of the model}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime. Default to 0.15}
\item{...}{further arguments to be passed to \code{nlar}}
}
\details{
Self Exciting Threshold AutoRegressive model.

\deqn{X_{t+s} = 
  x_{t+s} = ( \phi_{1,0} + \phi_{1,1} x_t + \phi_{1,2} x_{t-d} + \dots +
  \phi_{1,mL} x_{t - (mL-1)d} ) I( z_t \leq th) +
  ( \phi_{2,0} + \phi_{2,1} x_t + \phi_{2,2} x_{t-d} + \dots + \phi_{2,mH}
  x_{t - (mH-1)d} ) I(z_t > th) + \epsilon_{t+steps}}{
  x[t+steps] = ( phi1[0] + phi1[1] x[t] + phi1[2] x[t-d] + \dots + phi1[mL] x[t - (mL-1)d] ) I( z[t] <= th)
  + ( phi2[0] + phi2[1] x[t] + phi2[2] x[t-d] + \dots + phi2[mH] x[t -
  (mH-1)d] ) I( z[t] > th) + eps[t+steps]  
}

with \var{z} the treshold variable. The threshold variable can alternatively be specified by (in that order):
\describe{
\item{thDelay}{ \code{z[t] = x[t - thDelay*d ]} }
\item{mTh}{ \code{z[t] = x[t] mTh[1] + x[t-d] mTh[2] + \dots + x[t-(m-1)d] mTh[m]} }
\item{thVar}{ \code{z[t] = thVar[t]} }
}
For fixed \code{th} and threshold variable, the model is linear, so
\code{phi1} and \code{phi2} estimation can be done directly by CLS
(Conditional Least Squares).
Standard errors for phi1 and phi2 coefficients provided by the
\code{summary} method for this model are taken from the linear
regression theory, and are to be considered asymptoticals.
}
\seealso{
\code{\link{plot.setar}} for details on plots produced for this model from the \code{plot} generic.
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit a SETAR model, with threshold as suggested in Tong(1990)
mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)
mod.setar
summary(mod.setar)

if(require(FinTS)) {
	data(m.unrate)
	setar(diff(m.unrate), ML=c(2,3,4,12), MH=c(2,4,12), th=0.1, include="none")
}
}
\keyword{ ts }
\references{
Non-linear time series models in empirical finance, Philip Hans Franses and Dick van Dijk, Cambridge: Cambridge University Press (2000).

Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford: Oxford University Press (1990).
}
