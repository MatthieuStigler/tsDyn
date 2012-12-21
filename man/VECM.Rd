\name{VECM}
\alias{VECM}
\title{Estimation of Vector error correction model (VECM)}
\description{
Estimate either a  VECM by Engle-Granger or Johansen (MLE) method. 
}
\usage{
VECM(data, lag,r=1, include = c( "const", "trend","none", "both"),beta=NULL, estim=c("2OLS", "ML"),LRinclude=c("none", "const", "trend","both"))
}
\value{
Fitted model data
}
\arguments{
\item{data}{multivariate time series (first row  being first=oldest value)}
\item{lag}{Number of lags to include in each regime}
\item{r}{Number of cointegrating relationships}
\item{include}{Type of deterministic regressors to include}
\item{beta}{for VECM only: cointegrating value. If null, will be estimated}
\item{LRinclude}{Possibility to include in the long-run relationship and the ECT trend, constant... Can also be a matrix with exogeneous regressors}
\item{estim}{Type of estimator for the VECM: '2OLS' for the two-step approach or 'ML' for Johansen MLE}

}
\details{
This function is just a wrapper for the \code{\link{lineVar}}, with model="VECM".

More comprehensive functions for VECM are in package \pkg{vars}. A few differences appear in the VECM estimation:
\itemize{
\item{Engle-Granger estimator}{The Engle-Granger estimator is available}
\item{Presentation}{Results are printed in a different ways, using a matrix form}
\item{lateX export}{The matrix of coefficients can be exported to latex, with or without standard-values and significance stars}
}
 Here, only one cointegrating relationship can be estimated. Two estimators are available: the Engle-Granger two step approach (\code{2OLS}) or the Johansen (\code{ML}). For the 2OLS, deterministics regressors (or external variables if LRinclude is of class numeric) can be added for the estimation of the cointegrating value and for the ECT. This is only working when the beta value is not pre-specified.

The arg beta is the cointegrating value, the cointegrating vector will be taken as: (1, -beta).

Note that 
}
\seealso{
\code{\link{lineVar}}
\code{\link{TVAR}} and \code{\link{TVECM}} for the correspoding threshold models. \code{\link{linear}} for the univariate AR model.
}
\author{Matthieu Stigler}
\examples{
data(zeroyld)
data<-zeroyld

#Fit a VECM with Engle-Granger 2OLS estimator:
vecm.eg<-VECM(zeroyld, lag=2)

#Fit a VECM with Johansen MLE estimator:
vecm.jo<-VECM(zeroyld, lag=2, estim="ML")

#compare results with package vars:
if(require(vars)) {
data(finland)
#check long coint values
all.equal(VECM(finland, lag=2, estim="ML", r=2)$model.specific$coint, cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)$beta, check.attr=FALSE)
# check OLS parameters
all.equal(t(coefficients(VECM(finland, lag=2, estim="ML", r=2))), coefficients(cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)$rlm), check.attr=FALSE)

}


##export to Latex
toLatex(vecm.eg)
toLatex(summary(vecm.eg))
options("show.signif.stars"=FALSE)
toLatex(summary(vecm.eg), parenthese="Pvalue")
options("show.signif.stars"=TRUE)


}
\keyword{ ts }


