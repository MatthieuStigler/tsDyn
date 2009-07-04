\name{linear2}
\alias{linear2}
\title{Multivariate liner models: VAR and VECM}
\description{
Estimate either a VAR or a VEC
}
\usage{
linear2(data, lag, include = c( "const", "trend","none", "both"), model=c("VAR", "VECM"), I=c("level", "diff"),beta=NULL, estim=c("2OLS", "ML"),LRinclude=c("none", "const", "trend","both"))
}
\value{
Fitted model data
}
\arguments{
\item{data}{multivariate time series }
\item{lag}{Number of lags to include in each regime}
\item{include}{Type of deterministic regressors to include}
\item{model}{Model to estimate. Either a VAR or a VECM}
\item{I}{For VAR only: whether in the VAR the variables are to be taken in levels (original series) or in diference}
\item{beta}{for VECM only: cointegrating value. If null, will be estimated}
\item{LRinclude}{Possibility to include in the long-run relationship and the ECT trend, constant... Can also be a matrix with exogeneous regressors}
\item{estim}{Type of estimator for the VECM: '2OLS' for the two-step approach or 'ML' for Johansen MLE}

}
\details{
This function provides basic functionaities for VAR and VECM models. More comprehensive functions are in package \pkg{vars}. A few differences appear in the VECM estimation:
\itemize{
\item{1 cointegrating value}{\code{linear2} is able ton hanfdle with only 1 cointegrating relatioship.}
\item{Engle-Granger estimator}{The Engle-Granger estimator is available}
\item{Presentation}{Results are printed in a different ways, using a matrix form}
\item{lateX export}{The matrix of coefficients can be exported to latex, with or without standard-values and significance stars}
}
 Here, only one cointegrating relationship can be estimated. Two estimators are available: the Engle-Granger two step approach (\code{2OLS}) or the Johansen (\code{ML}). For the 2OLs, deterministics regressors (or external variables if LRinclude is of class numeric) can be added for the estimation of the cointegrating value and for the ECT. This is only working when the beta value is not pre-specified.

The arg beta is the cointegrating value,Cointegrating vector will be taken as: (1, -beta).

Note that 
}
\seealso{
\code{\link{TVAR}} and \code{\link{TVECM}} for the correspoding threshold models. \code{linear} for the univariate AR model.
}
\author{Matthieu Stigler}
\examples{
data(zeroyld)
data<-zeroyld

#Fit a VAR
VAR<-linear2(data, lag=1)
VAR
summary(VAR)

#compare results with package vars:
if(require(vars)) {
	a<-VAR(data, p=1)
	vaco1<-coef(a)$short.run[c(3,1,2),1]
	vaco2<-coef(a)$long.run[c(3,1,2),1]
	round(coef(VAR),8)==round(rbind(vaco1, vaco2),8)
}

###VECM
VECM.EG<-linear2(data, lag=2, model="VECM")
VECM.EG
summary(VECM.EG)

VECM.ML<-linear2(data, lag=2, model="VECM", estim="ML")
VECM.ML
summary(VECM.ML)


###Check Johansen MLE
myVECM<-linear2(data, lag=1, include="const", model="VECM", estim="ML")
summary(myVECM, digits=7) 
#comparing with vars package
if(require(vars)){
	a<-ca.jo(data, spec="trans")
	summary(a)
	#same answer also!
}

##export to Latex
toLatex(VECM.EG)
toLatex(summary(VECM.EG))
options("show.signif.stars"=FALSE)
toLatex(summary(VECM.EG), parenthese="Pvalue")
options("show.signif.stars"=TRUE)


}
\keyword{ ts }


