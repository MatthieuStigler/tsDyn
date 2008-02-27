\name{TVAR_LRtest}
\alias{TVAR_LRtest}
\alias{Tvar_LRtest}
\title{Test of linearity}
\description{Multivariate extension of the linearity against threshold test from Hansen (1999) with bootstrap distribution 
}
\usage{
TVAR_LRtest <- function (x, m=1, d = 1, steps = d, trend=TRUE, 
series, thDelay = 1:2,    mTh=1, thVar, nboot=10, plot=FALSE,
trim=0.1, test=c("1vs", "2vs3"))
}

\value{A list containing:

-The values of each LR test

-The bootstrap Pvalues and critical values for the test selected
}
\arguments{
\item{x}{ multivariate time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{thDelay}{'time delay' for the threshold variable IN THE MULTIVARIATE IMPLEMENTATION (so 0:1 in univariate correspond to 1:2)}
\item{mTh}{combination of variables with same lag order for the transition variable. Either a single value (indicating which variable to take) or a combination}
\item{thVar}{external threshold variable}
\item{nboot}{numer of bootstrap replications}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{plot}{Whether a plot with the bootstrap distribution and the Ftest should be shown}
\item{test}{Type of test. See details}

}
\details{
This test is just the multivariate extension proposed by Lo and Zivot of the linearity test of Hansen (1999). As in univariate case, estimation of the first threshold parameter is made with CLS, for the second threshold a conditional search with one iteration is made. Instead of a Ftest comparing the SSR for the univariate case, a  Likelihood Ratio (LR) test comparing the covariance matrix of each model is computed. 

\deqn{
LR_{ij}=T( ln(\det \hat \Sigma_{i}) -ln(\det \hat \Sigma_{j}))}
where \eqn{ \hat \Sigma_{i}} is the estimated covariance matrix of the model with i regimes (and so i-1 thresholds). 

Three test are avalaible. The both first can be seen as linearity test, whereas the third can be seen as a specification test: once the 1vs2 or/and 1vs3 rejected the linearity and henceforth accepted the presence of a threshold, is a model with one or two thresholds preferable?

Test \bold{1vs}2: Linear VAR versus 1 threshold TVAR

Test \bold{1vs}3: Linear VAR versus 2 threshold2 TVAR

Test \bold{2vs3}: 1 threshold TAR versus 2 threshold2 TAR  

The both first are computed together and avalaible with test="1vs". The third test is avalaible with test="2vs3". 

The homoskedastik bootstrap distribution is based on resampling the residuals from H0 model, estimating the threshold parameter and then computing the Ftest, so it involves many computations and is pretty slow. 
}
\seealso{
\code{\link{setarTest}} for the univariate version. \code{\link{OlsTVAR}} for estimation of the model.
}

\author{ Matthieu Stigler}
\examples{

data(zeroyld)
data<-zeroyld

TVAR_LRtest(data, m=2, mTh=1,thDelay=1:2, nboot=3, plot=FALSE, trim=0.1, test="1vs")}
\keyword{ ts }
\references{
Hansen (1999) Testing for linearity, Journal of Economic Surveys, Volume 13, Number 5, December 1999 , pp. 551-576(26)
avalaible at: \url{http://www.ssc.wisc.edu/~bhansen/papers/cv.htm}

Lo and Zivot(2001) "Threshold Cointegration and Nonlinear Adjustment to the Law of One Price," Macroeconomic Dynamics, Cambridge University Press, vol. 5(4), pages 533-76, September.


}
