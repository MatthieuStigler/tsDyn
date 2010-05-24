\name{VAR.sim}
\alias{VAR.sim}
\title{Simulation of VAR}
\description{
Estimate or bootstraps a multivariate Threshold VAR
}
\usage{
VAR.sim(B, n=200, lag=1, include = c("const", "trend","none", "both"), starting=NULL,  innov=rmnorm(n, mean=0, varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE)
}


\value{
A matrix with the simulated/bootstraped series.
}
\arguments{
\item{B}{Matrix of coefficients to simulate}
\item{n}{Number of observations to create when type="simul"}
\item{lag}{Number of lags to include in each regime}
\item{include}{Type of deterministic regressors to include. NOT WORKING PROPERLY CURRENTLY if not const}
\item{starting}{Starting values when a simulation with given parameter matrix is made}
\item{innov}{Innovations used for simulation. Should be matrix of dim nxk. By default multivariate normal.}
\item{varcov}{Variance-covariance matrix for the inov when default multivariate normal is used.}
\item{show.parMat}{Logical. Should the parameter matrix be shown? Usefull to understand how to give right input}
}
\details{
Bootstrap a VAR model. K, the number of variables, is determined by the number of rows of B matrix.

}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{
##simulate VAR as in Enders 2004, p 268
B1<-matrix(c(0.7, 0.2, 0.2, 0.7), 2)
var1<-VAR.sim(B=B1,n=100,include="none")
ts.plot(var1, type="l", col=c(1,2))
}
\keyword{ ts }
\references{
TODO
}

