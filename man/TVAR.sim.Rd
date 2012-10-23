\name{TVAR.sim}
\alias{TVAR.sim}
\alias{VAR.sim}

\title{Simulation of a multivariate Threshold Autoregressive model (TVAR)}

\description{
Simulate a multivariate Threshold VAR (TVAR)
}
\usage{
TVAR.sim(B, Thresh, nthresh=1, n=200, lag=1, include = c("const", "trend","none", "both"),  thDelay=1,  thVar=NULL, mTh=1, starting=NULL, innov=rmnorm(n, mean=0, varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE, round=FALSE, seed, ...)
VAR.sim(B,  n=200, lag=1, include = c("const", "trend","none", "both"),  starting=NULL, innov=rmnorm(n, mean=0, varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE, round=FALSE, seed, ...)

}


\value{
A matrix with the simulated/bootstraped series.
}
\arguments{
\item{B}{Matrix of coefficients to simulate}
\item{Thresh}{The threshold value(s). Vector of length nthresh}
\item{nthresh}{number of threshold (see details)}
\item{n}{Number of observations to create when type="simul"}
\item{lag}{Number of lags to include in each regime}
\item{include}{Type of deterministic regressors to include. NOT WORKING PROPERLY CURRENTLY if not const}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d) PLEASE NOTE that the notation is currently different to univariate models in tsDyn. The left side variable is taken at time t, and not t+1 as in univariate cases.}
\item{thVar}{external transition variable}
\item{mTh}{combination of variables with same lag order for the transition variable. Either a single value (indicating which variable to take) or a combination}
\item{starting}{Starting values (matrix of dimension lag x k). If not given, set to zero. }
\item{innov}{Innovations used for simulation. Should be matrix of dim n x k. By default multivariate normal.}
\item{varcov}{Variance-covariance matrix for the innovations. By default identity matrix.}
\item{show.parMat}{Logical. Should the parameter matrix be shown? Useful to understand how to give right input}
\item{round}{Rounds the series created to have the same digits (hopefully) as original series. }
\item{seed}{Optional. Seed for the random number generation.}
\item{\dots}{Further arguments passed to the underlying (un-exported) \code{TVAR.gen} function}
}
\details{
This function offers the possibility to generate series following a TVAR.

By giving a matrix of coefficients, on can only simulate a VAR (nthresh=0) or TVAR (nthresh=1 or 2). One can have a specification with constant (default), trend, both or none (see arg include). Order in parameters is include/lags (VECM) and include/lags/include/lags for TVECM, hence, a matrix for a TVECM with 3 regimes, a const and a 2 lags would have 2 lines and 2*(1+4) columns. The innovations can be given by the user (a matrix of dim nxk, here n does not include the starting values!), by default it uses a multivariate normal distribution, with covariance matrix specified by varcov. The starting values (of dim lags x k) can be given. The user should take care for their choice, since it is not sure that the simulated values will cross the threshold even once.

The matrix \sQuote{B} has to be in the form: constant, trend, lags, then repeated if many regimes. In case of uncertainty, using \code{who.parMat=TRUE} will print the matrix as interpreted by the function, helping the user to feed the right input. 

}

\seealso{
\code{\link{TVAR}} to estimate a TVAR, \code{\link{VAR.sim}} to simulate/bootstrap a VAR, \code{\link{TVECM.sim}} to simulate/bootstrap a TVECM. 
}
\author{Matthieu Stigler}

\examples{
##simulate VAR as in Enders 2004, p 268
B1<-matrix(c(0.7, 0.2, 0.2, 0.7), 2)
var1<-TVAR.sim(B=B1,nthresh=0,n=100, include="none")
ts.plot(var1, type="l", col=c(1,2))


B2<-rbind(c(0.5, 0.5, 0.5), c(0, 0.5, 0.5))
varcov<-matrix(c(1,0.2, 0.3, 1),2)
var2<-TVAR.sim(B=B2,nthresh=0,n=100, include="const", varcov=varcov)
ts.plot(var2, type="l", col=c(1,2))


##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
colnames(B) <- paste(rep(c("Const", "Lag_1_var1", "Lag_1_var2"), 2), c("Low", "High"), sep="_")
sim<-TVAR.sim(B=B,nthresh=1,n=500, mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)


}
\keyword{ ts }
\keyword{ bootstrap }

