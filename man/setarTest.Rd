\name{setarTest}
\alias{setarTest}
\alias{setartest}
\title{Test of linearity}
\description{Test of linearity against threshold of Hansen (1999) with bootstrap distribution
}
\usage{
setarTest(x, m, d = 1, steps = d, series, thDelay = 0:1, mL, mH,
    mTh, thVar, nboot=10, plot=FALSE, trim=0.1, test=c("1vs", "2vs3"))
}

\value{A list containing:

-The residual Sum of squares of model AR, 1 threshold TAR and 2 thresholds TAR

-The Ftest statistic for each test

-The bootstrap Pvalues and critical values for the test selected
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{mL}{autoregressive order for 'low' regime (dafult: m). Must be <=m}
\item{mH}{autoregressive order for 'high' regime (default: m). Must be <=m}
\item{thDelay}{'time delay' for the threshold variable (as multiple of embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the threshold variable}
\item{thVar}{external threshold variable}
\item{nboot}{numer of bootstrap replications}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{plot}{Whether a plot with the bootstrap distribution and the Ftest should be shown}
\item{test}{Type of test. See details}

}
\details{
Estimation of the first threshold parameter is made with CLS, a conditional search with one iteration is made for the second threshold. The Ftest comparing the residual sum of squares (SSR) of each model is computed.


 \deqn{
F_{ij}=T( (S_{i}-S_{j})/S_{j} )}

where \eqn{S_{i}} is the SSR of the model with i regimes (and so i-1 thresholds). 

Three test are avalaible. The both first can be seen as linearity test, whereas the third can be seen as a specification test: once the 1vs2 or/and 1vs3 rejected the linearity and henceforth accepted the presence of a threshold, is a model with one or two thresholds preferable?

Test \bold{1vs}2: Linear AR versus 1 threshold TAR

Test \bold{1vs}3: Linear AR versus 2 threshold2 TAR

Test \bold{2vs3}: 1 threshold TAR versus 2 threshold2 TAR  

The both first are computed together and avalaible with test="1vs". The third test is avalaible with test="2vs3". 

The homoskedastik bootstrap distribution is based on resampling the residuals from H0 model, estimating the threshold parameter and then computing the Ftest, so it involves many computations and is pretty slow. 
}
\seealso{
\code{\link{TVAR_LRtest}} for the multivariate version. \code{\link{SETAR}} for estimation of the model.
}

\author{ Matthieu Stigler}
\examples{
#Data used by Hansen
sun<-(sqrt(sunspot.year+1)-1)*2

#Test 1vs2 and 1vs3
setarTest(sun, m=11, thDelay=0:1, nboot=5, plot=TRUE, trim=0.1, test="1vs")
}
\keyword{ ts }
\references{
Hansen (1999) Testing for linearity, Journal of Economic Surveys, Volume 13, Number 5, December 1999 , pp. 551-576(26)
avalaible at: \url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}

}
