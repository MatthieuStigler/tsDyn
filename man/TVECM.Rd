\name{TVECM}
\alias{TVECM}
\title{Treshold Vector Error Correction model (VECM)}
\description{
Estimate a Treshold Vector Error Correction model (VECM)
}
\usage{
TVECM<-function(data,lag=1,trend=TRUE, bn=50, ngridG=50, trim=0.05, nthresh=1,plot=TRUE, dummyToBothRegimes=TRUE, methodMapply=FALSE, gamma1=list(exact=NULL, int=c("from","to"), around="val"),gamma2=list(exact=NULL, int=c("from","to"), around="val"), beta=list(exact=NULL, int=c("from","to"), around=c("val","by")), rest=c("none", "equal", "signOp"), model=c("All", "only_ECT") )
}
\value{
Fitted model data
}
\arguments{
\item{data}{ time series }
\item{lag}{Number of lags to include in each regime}
\item{trend}{If a trend should be included in the model}
\item{nthresh}{number of threshold (see details)}
\item{ngridG}{number of elements to search for the threshold value}
\item{bn}{number of elements to search for the cointegrating value}
\item{beta}{different possibilities to pre-specify an exact value, an interval or a central point for the search of the cointegrating value}
\item{gamma1}{different possibilities to pre-specify an exact value, an interval or a central point for the search of gamma1}
\item{gamma2}{different possibilities to pre-specify an exact value or a central point for the search of gamma2}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{dummyToBothRegimes} {Whether the dummy in the one threshold model is applied to each regime or not. }
\item{rest} {Whether restrictions are put on the search.  }
\item{methodMapply}{only for programming. Is to make the choice between a for loop or mapply implementation}
\item{rest}{Currently not avalaible}
\item{model}{The choice of the model. See details}
\item{plot}{Whether a grid with the SSR of each threshold should be printed}
}
\details{
For fixed threshold  and cointegrating vector, the model is linear, so estimation can be done directly by CLS (Conditional Least Squares). The search of the parameters values which minimize the residual sum of squares (SSR) is made upon a grid of potential values. The search on the threshold values is made upon existing ones, whereas the search for the cointegrating value is made ba default from a 

nthresh=1: estimation of one threshold model (two regimes) upon a grid of \var{ngrid} values (default to ALL) possible thresholds and delays values. 

nthresh=2: estimation of two thresholds model (three regimes) Conditional on the threshold found in model where nthresh=1, the second threshold is searched. When both are found, a second grid search is made with 30 values around each threshold.

The model can be either with a threshold effect on all variables ("All") or only on the error correction term (ECT) ("only_ECT"). In the second case, the value for the middle threshold is taken a null, as in Balke and Fomby (1997). 
}
\seealso{
TODO
}
\author{Matthieu Stigler}
\examples{

data(zeroyld)
data<-zeroyld[,c(36,19)]
colnames(data)<-c("short", "long")

TVECM(data, nthresh=2,lag=1, bn=20, ngridG=30, plot=TRUE,trim=0.05, model="All", trend=FALSE)
}
\keyword{ ts }
\references{
TODO
}

