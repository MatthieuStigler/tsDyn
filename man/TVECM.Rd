\name{TVECM}
\alias{TVECM}
\title{Treshold Vector Error Correction model (VECM)}
\description{
Estimate a Treshold Vector Error Correction model (VECM)
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
\item{methodMapply}{only for programming. Is to make the choice between a for loop or mapply implementation}
\item{trace}{should additional infos be printed? (logical)}
}
\details{
For fixed threshold  and cointegrating vector, the model is linear, so estimation can be done directly by CLS (Conditional Least Squares). The search of the parameters values which minimize the residual sum of squares (SSR) is made upon a grid of potential values. The search on the threshold values is made upon existing ones, whereas the search for the cointegrating value is made by default based on a confidence interval around the estimate form linear model. 

nthresh=1: estimation of one threshold model (two regimes) upon a grid of \var{ngridTh} values (default to ALL) possible thresholds and delays values. 

nthresh=2: estimation of two thresholds model (three regimes) Conditional on the threshold found in model where nthresh=1, the second threshold is searched. When both are found, a second grid search is made with 30 values around each threshold.

The model can be either with a threshold effect on all variables ("All") or only on the error correction term (ECT) (argument "only ECT"). In the second case, the value for the middle threshold is taken a null, as in Balke and Fomby (1997). 
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{

data(zeroyld)
data<-zeroyld


TVECM(data, nthresh=2,lag=1, ngridBeta=20, ngridTh=30, plot=TRUE,trim=0.05, common	="All")
}
\keyword{ ts }
\references{
TODO
}

