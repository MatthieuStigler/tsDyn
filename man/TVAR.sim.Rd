\name{TVAR.sim}
\alias{TVAR.sim}
\title{Simulation and bootstrap of multivariate Treshold Autoregressive model}
\description{
Estimate or bootstraps a multivariate Threshold VAR
}
\usage{
TVAR.sim(data,B,TVARobject, Thresh, nthresh=1, type=c("simul","boot", "check"), n=200, lag=1, include = c("const", "trend","none", "both"),  thDelay=1,  thVar=NULL, mTh=1, starting=NULL, sigma, round=FALSE)
}


\value{
a list with the simulated/bootstraped data and the parameter matrix used. 
}
\arguments{
\item{data}{matrix of parameter to simulate}
\item{B}{Matrix of coefficients to simulate}
\item{TVARobject}{Object computed by function TVAR()}
\item{Thresh}{The threshold value(s). Vector of length nthresh}
\item{nthresh}{number of threshold (see details)}
\item{type}{Whether a bootstrap or simulation is to employ. See details}
\item{n}{Number of observations to create when type="simul"}
\item{lag}{Number of lags to include in each regime}
\item{include}{Type of deterministic regressors to include. NOT WORKING PROPERLY CURRENTLY if not const}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d) PLEASE NOTE that the notation is currently different to univariate models in tsDyn. The left side variable is taken at time t, and not t+1 as in univariate cases.  }
\item{thVar}{external transition variable}
\item{mTh}{combination of variables with same lag order for the transition variable. Either a single value (indicating which variable to take) or a combination}
\item{starting}{Starting values when a simulation with given parameter matrix is made}
\item{sigma}{variance of the errors}
\item{round}{Rounds the series created to have the same digits (hopefully) as original series. Useful when t}
}
\details{
This function offers the possibility to generate series following a TVAR from two approaches: bootstrap or simulation. When the data is given, one can use a simple residual bootstrap or simulate a serie with the same parameter matrix and with normal distributed residuals (with variance pre-specified). The function calls internaly TVAR() to compute the parameter matrix. an alternative is to give a 
When the parameter matrix is given, there is only the possibility to simulate series. The starting values (of length equal to the number of lags) can be given. The user should take care for the choice of the starting values and parameters values, since it is not sure that the simulated values will cross the threshold even once.
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{
##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
sim<-TVAR.sim(B=B,nthresh=1,n=500, type="simul",mTh=1, Thresh=5, starting=c(5.2, 5.5), sigma=c(0.3,0.4))$serie

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)

##Bootstrap a TVAR with two threshold (three regimes)
data(zeroyld)
serie<-zeroyld
TVAR.sim(data=serie,nthresh=2, type="boot",mTh=1, Thresh=c(7,9))

##Check the bootstrap
cbind(TVAR.sim(data=serie,nthresh=2, type="check",mTh=1, Thresh=c(7,9))$serie,serie)
}
\keyword{ ts }
\references{
TODO
}

