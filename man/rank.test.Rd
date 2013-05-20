\name{rank.test}
\alias{print.rank.test}
\alias{rank.test}
\alias{summary.rank.test}
\title{Test of the cointegrating rank}
\usage{
  rank.test(vecm, type = c("eigen", "trace"), r_null,
    cval = 0.05)

  \method{print}{rank.test} (x, ...)

  \method{summary}{rank.test} (object, ...)
}
\arguments{
  \item{vecm}{\sQuote{VECM} object computed with the
  function \code{\link{VECM}}.}

  \item{type}{Type of test, either 'trace' or 'eigenvalue'.
  See details below.}

  \item{r_null}{Rank to test specifically.}

  \item{cval}{Critical value level for the automatic test.}

  \item{x}{The output from \code{rank.test} for the print
  method.}

  \item{object}{The output from \code{rank.test} for the
  summary method.}

  \item{...}{Unused.}
}
\value{
  An object of class \sQuote{rank.test}, with
  \sQuote{print} and \sQuote{summary methods}.
}
\description{
  Maximum-likelihood test of the cointegrating rank.
}
\details{
  This function computes the two maximum-likelihood tests
  for the cointegration rank from Johansen (1996). Tests
  are: \describe{ \item{trace}{Test the hypothesis of rank
  \sQuote{h} against rank \sQuote{K}, i.e. against the
  alternative that the system is stationary.}
  \item{eigenvalue}{Test the hypothesis of rank \sQuote{h}
  against rank \sQuote{h+1}.} }

  The test works for five specifications of the
  deterministic terms as in Doornik et al (1998), to be
  specified in the previous call to \code{\link{VECM}}:
  \describe{ \item{H_ql}{Unrestricted constant and trend:
  use \code{include="both"} } \item{H_l}{Unrestricted
  constant and restricted trend: use
  \code{include="const"}} and \code{LRinclude="trend"}
  \item{H_lc}{Unrestricted constant and no trend: use
  \code{include="const"}} \item{H_c}{Restricted constant
  and no trend: use \code{LRinclude="const"}} \item{H_z}{No
  constant nor trend: use \code{include="none"}} }

  Two testing procedures can be used: \describe{
  \item{Specific test}{By specifying a value for
  \sQuote{r_null}. The \sQuote{pval} value returned gives
  the speciifc p-value.} \item{Automatic test}{If not value
  is specified for \sQuote{r_null}, the function makes a
  simple automatic test: returns the rank (slot \sQuote{r})
  of the first test not rejected (level specified by arg
  \sQuote{cval}) as recommend i.a. in Doornik et al (1998,
  p. 544).} }

  A full table with both test statistics ad their
  respective p-values is given in the summary method.

  P-values are obtained from the gamma aproximation from
  Doornik (1998, 1999). Small sample values adjusted for
  the sample site are also available in the summary method.
  Note that the \sQuote{effective sample size} for the
  these values is different from output in gretl for
  example.
}
\section{Comparison with urca}{
  While \code{\link[urca]{ca.jo}} in package \pkg{urca} and
  \code{rank.test} both implement Johansen tests, there are
  a few differences:

  \itemize{ \item \code{rank.test} gives p-values, while
  \code{ca.jo} gives only critical values.  \item
  \code{rank.test} allows for five different specifications
  of deterministic terms (see above), \code{ca.jo} for only
  three.  \item \code{ca.jo} allows for seasonal and
  exogenous regressors, which is not available in
  \code{rank.test}.  \item The lag is specified
  differently: \code{K} from \code{ca.jo} corresponds to
  \code{lag}+1 in \code{rank.test}.  }
}
\examples{
data(barry)

## estimate the VECM with Johansen!
ve <- VECM(barry, lag=1, estim="ML")

## specific test:
ve_test_spec <- rank.test(ve, r_null=1)
ve_test_spec_tr <- rank.test(ve, r_null=1, type="trace")

ve_test_spec
ve_test_spec_tr

## No specific test: automatic method
ve_test_unspec <- rank.test(ve)
ve_test_unspec_tr <- rank.test(ve, type="trace")

ve_test_unspec
ve_test_unspec_tr

## summary method: output will be same for all types/test procedure:
summary(ve_test_unspec_tr)

## The function works for many specification of the VECM(), try:
rank.test(VECM(barry, lag=3, estim="ML"))
rank.test(VECM(barry, lag=3, include="both",estim="ML"))
rank.test(VECM(barry, lag=3, LRinclude="const",estim="ML"))

## Note that the tests are simple likelihood ratio, and hence can be obtained also manually:
-2*(logLik(ve, r=1)-logLik(ve, r=2)) # eigen test, 1 against 2
-2*(logLik(ve, r=1)-logLik(ve, r=3)) # eigen test, 1 against 3
}
\author{
  Matthieu Stigler
}
\references{
  - Doornik, J. A. (1998) Approximations to the Asymptotic
  Distributions of Cointegration Tests, Journal of Economic
  Surveys, 12, 573-93

  - Doornik, J. A. (1999) Erratum [Approximations to the
  Asymptotic Distribution of Cointegration Tests], Journal
  of Economic Surveys, 13, i

  - Doornik, Hendry and Nielsen (1998) Inference in
  Cointegrating Models: UK M1 Revisited, Journal of
  Economic Surveys, 12, 533-72

  - Johansen, S. (1996) Likelihood-based inference in
  cointegrated Vector Autoregresive Models, Oxford
  University Press
}
\seealso{
  \code{\link{VECM}} for estimating a VECM.
  \code{\link{rank.select}} to estimate the rank based on
  information criteria.

  \code{\link[urca]{ca.jo}} in package \pkg{urca} for
  another implementation of Johansen cointegration test
  (see section \sQuote{Comparison with urca} for more
  infos).
}
\keyword{ts}

