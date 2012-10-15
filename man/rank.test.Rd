\name{rank.test}
\alias{rank.test}
\alias{print.rank.test}
\alias{summary.rank.test}

\title{Test of the cointegrating rank}

\description{
Maximum-likelihood test of the cointegrating rank.
}

\usage{
rank.test(vecm, type=c("trace", "eigen"), r_null, cval=0.05)
\method{print}{rank.test}(x,...)
\method{summary}{rank.test}(object,...)
}


\arguments{
\item{vecm}{\sQuote{VECM} object computed with the function \code{\link{VECM}}. }
\item{r_null}{Rank to test specifically. }
\item{cval}{Critical value level for the automatic test. }
\item{x}{The output from \code{rank.test} for the print method.}
\item{object}{The output from \code{rank.test} for the summary method.}
\item{...}{Unused.}
}

\details{
This function computes the two maximum-likelihood tests for the cointegration rank from Johansen (1996). Tests are:
  \describe{
    \item{trace}{Test the hypothesis of rank \sQuote{h} against rank \sQuote{K}, i.e. against the alternative that the system is stationary.}
    \item{eigenvalue}{Test the hypothesis of rank \sQuote{h} against rank \sQuote{h+1}.}
  }

Two testing procedure canbe used:
  \describe{
    \item{Specific test}{By specifying a value for \sQuote{r_null}. The \sQuote{pval} value returned gives the speciifc p-value.}
    \item{Automatic test}{If not value is specified for \sQuote{r_null}, the function makes a simple automatic test: 
	    returns the rank (slot \sQuote{r}) of the first test not rejected (level specified by arg \sQuote{cval}).}
}

The two statistics are given in the summary method, while the print method just gives information on the specific test chosen, indicating the first hypothesis to be not-rejected. 

P-values are obtained from the gamam distribution from Doornik (1998). Small sample values adjused for the sample sitze are also available in the summary method. 
Note that the \sQuote{effective sample size} for the these values is different from output in gretl for example. 
}

\value{
An object of class \sQuote{rank.test}, with \sQuote{print} and \sQuote{summary methods}. 
}


\seealso{
\code{\link{VECM}} for estimating a VECM. \code{\link{rank.select}} to estimate the rank based on information criteria. 

\code{\link[urca]{ca.jo}} in package \pkg{urca} for another implementation of Johansen cointegration test. 
}

\author{Matthieu Stigler}

\references{
- Doornik, J. A. (1998) Approximations to the Asymptotic Distributions of Cointegration Tests, Journal of Economic Surveys, 12, 573-93

- Doornik, J. A. (1999) Erratum [Approximations to the Asymptotic Distribution of Cointegration Tests], Journal of Economic Surveys, 13, i

- Johansen, S. (1996) Likelihood-based inference in cointegrated Vector Autoregresive Models, Oxford University Press

}

\examples{

data(barry)

## estimate the VECM with Johansen! 
ve <- VECM(barry, lag=1, estim="ML")

## specific test:
ve_test_spec <- rank.test(ve, r_null=1)
ve_test_spec_eig <- rank.test(ve, r_null=1, type="eigen")

ve_test_spec
ve_test_spec_eig

## No specific test: automatic method
ve_test_unspec <- rank.test(ve)
ve_test_unspec_eig <- rank.test(ve, type="eigen")

ve_test_unspec
ve_test_unspec_eig

## summary method: output will be same for all types/test procedure:
summary(ve_test_unspec_eig)

## The function works for many specification of the VECM(), try:
rank.test(VECM(barry, lag=3, estim="ML"))
rank.test(VECM(barry, lag=3, include="both",estim="ML"))
rank.test(VECM(barry, lag=3, LRinclude="const",estim="ML"))

## Note that the tests are simple likelihood ratio, and hence can be obtained also manually:
-2*(logLik(ve, r=1)-logLik(ve, r=2)) # eigen test, 1 against 2
-2*(logLik(ve, r=1)-logLik(ve, r=3)) # eigen test, 1 against 3
}

\keyword{ ts }
\keyword{ VECM }
\keyword{ cointegration }


