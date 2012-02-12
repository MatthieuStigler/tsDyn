\name{TVECM}
\alias{TVECM}
\title{Threshold Vector Error Correction model (VECM)}
\description{
Estimate a Threshold Vector Error Correction model (VECM)
}
\usage{
TVECM(data,lag=1,nthresh=1, trim=0.05, ngridBeta=50, ngridTh=50, plot=TRUE,  th1=list(exact=NULL, int=c("from","to"), around="val"), th2=list(exact=NULL, int=c("from","to"), around="val"), beta=list(exact=NULL, int=c("from","to"), around=c("val","by")), restr=c("none", "equal", "signOp"), common=c("All", "only_ECT"), include = c( "const", "trend","none", "both"),dummyToBothRegimes=TRUE,beta0=0,methodMapply=FALSE, trace=TRUE )
}
\value{
Fitted model data
}
\arguments{
\item{data}{ time series }
\item{lag}{Number of lags to include in each regime}
\item{nthresh}{number of threshold (see details)}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{ngridBeta}{number of elements to search for the cointegrating value}
\item{ngridTh}{number of elements to search for the threshold value}
\item{plot}{Whether a grid with the SSR of each threshold should be printed}
\item{th1}{different possibilities to pre-specify an exact value, an interval or a central point for the search of gamma1}
\item{th2}{different possibilities to pre-specify an exact value or a central point for the search of gamma2}
\item{beta}{different possibilities to pre-specify an exact value, an interval or a central point for the search of the cointegrating value}
\item{restr}{Currently not avalaible}
\item{common}{Whether the regime-specific dynamics are only for the ECT or for the ECT and the lags}
\item{include}{Type of deterministic regressors to include}
\item{dummyToBothRegimes}{Whether the dummy in the one threshold model is applied to each regime or not. }
\item{beta0}{Additional regressors to include in the cointegrating relation}
\item{methodMapply}{only for programming. Is to make the choice between a for loop or \code{mapply} implementation}
\item{trace}{should additional infos be printed? (logical)}
}
\details{
For fixed threshold  and cointegrating vector, the model is linear, so estimation of the regression parameters can be done directly by CLS (Conditional Least Squares). The search of the threshold and cointegrating parameters values which minimize the residual sum of squares (SSR) is made on a grid of potential values. Potential values of the threshold parameter are taken from the existing values of the error correction term, while the parameters of the grid for cointegrating value are taken from the confidence interval obtained in the linear VECM. 

The function can estimate 1 as well as two threhsolds:

nthresh=1: estimation of one threshold model (two regimes) upon a grid of \var{ngridTh} values (default to ALL) possible thresholds and delays values. 

nthresh=2: estimation of two thresholds model (three regimes) Conditional on the threshold found in model where nthresh=1, the second threshold is searched. When both are found, a second grid search is made with 30 values around each threshold.

The model can be either with a threshold effect on all variables ("All") or only on the error correction term (ECT) (argument "only ECT"). In the second case, the value for the middle threshold is taken a null, as in Balke and Fomby (1997). 
}
\seealso{
\code{\link{VECM}} for the linear VECM, \code{\link{TVAR}} for the threshold VAR, \code{\link{TVECM.SeoTest}} to test for TVECM, \code{\link{TVECM.sim}} to simulate/bootstrap a TVECM.
}
\author{Matthieu Stigler}
\examples{

data(zeroyld)
data<-zeroyld

##Estimate a TVECM (we use here minimal grid, it should be usually much bigger!)

tv<-TVECM(data, nthresh=2,lag=1, ngridBeta=20, ngridTh=30, plot=TRUE,trim=0.05, common="All")

print(tv)
summary(tv)

#Obtain diverse infos:
AIC(tv)
BIC(tv)

res.tv<-residuals(tv)

#export the equations as Latex:
toLatex(tv)
}
\keyword{ ts }
\references{
Hansen, B. and Seo, B. (2002), Testing for two-regime threshold cointegration in vector error-correction models, Journal of Econometrics, 110, pages 293 - 318

Seo, M. H. (2009) Estimation of non linear error-correction models, Working paper

}

