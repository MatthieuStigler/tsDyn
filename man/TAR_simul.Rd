\name{TAR_simul}
\alias{TAR_simul}
\title{Simulation and bootstrap of Treshold Autoregressive model}
\description{
Estimate or bootstraps a Threshold VAR
}
\usage{
TAR_simul<-function(data,B, sigma,n=200, lag=1, trend=TRUE, nthresh=0, thDelay=1, Thresh, thVar=NULL, mTh=NULL, type=c("boot", "simul", "check"), starting=NULL)
}
\value{
a list with the simulated/bootstraped data and the parameter matrix used. 
}
\arguments{
\item{y}{ univariate time series }
\item{B}{vector of parameter to simulate}
\item{Thresh}{The threshold value(s). Vector of length nthresh}
\item{nthresh}{number of threshold (see details)}
\item{lag}{Number of lags to include in each regime}
\item{type}{Whether a bootstrap or simulation is to employ. See details}
\item{trend}{If a trend should be included in the model}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d). }
\item{thVar}{external transition variable}
\item{mTh}{Currently not implemented}
\item{starting}{Starting values when a simulation with given parameter matrix is made}
}
\details{
This function offers the possibility to generate series following a TAR from two approaches: bootstrap or simulation. When the data is given, one can use a simple residual bootstrap or simulate a serie from the parameter matrix and with normal distributed residuals (with variance pre-specified). The choice "check" is just there to check the function: one should obtain the same values as the given values. Please report if different. 
When the parameter matrix is given, there is only the possibility to simulate series. The starting values (of length equal to the number of lags) can be given. The user should take care for the choice of the starting values and parameters values, since it is not sure that the simulated values will cross the threshold even once.
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{
##Simulation of a TAR with 1 threshold
sim<-TAR_simul(B=c(2.9,-0.4,-0.1,-1.5, 0.2,0.3),lag=2, type="simul", nthresh=1, sigma=1.2,Thresh=2, starting=c(2.8,2.2))$serie
mean(ifelse(sim>2,1,0))	#approximation of values over the threshold

#check the result
selectSETAR(sim, m=2)

##Bootstrap a TAR with two threshold (three regimes)
sun<-(sqrt(sunspot.year+1)-1)*2
TAR_simul(data=sun,nthresh=2,n=500, type="boot", Thresh=c(7,9))$serie

##Check the bootstrap
cbind(TAR_simul(data=sun,nthresh=2,n=500, type="check", Thresh=c(7,9))$serie,sun)
}
\keyword{ ts }
\references{
TODO
}

