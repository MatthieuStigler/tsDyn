\name{rank.select}
\alias{rank.select}
\title{Selection of the cointegrating rank with Information criterion.}

\description{
Selection of the cointegrating rank and the lags with Information criterion (AIC, BIC).
}

\usage{
rank.select(data, lag.max=10, r.max=ncol(data)-1, include="intercept", fitMeasure=c("SSR", "LL"))
\method{print}{rank.select}(x,...)
\method{summary}{rank.select}(x,...)
}


\arguments{
\item{data}{multivariate time series }
\item{lag.max}{Maximum number of lags to investigate}
\item{r.max}{Maximum rank to investigate}
\item{include}{Type of deterministic regressors to include. For now only intercept. }
\item{fitMeasure}{Whether the AIC/BIC should be based on the full likelihood, or just the SSR. See explanations in \code{\link{logLik.VECM}}}
}

\details{
This function estimates the AIC, BIC and Hannan-Quinn for each rank (up to \code{lags.max}) and lags (up to \code{lags.max}). This method has been shown to be useful to select simultaneously the rank and the lags, see references. 
}

\value{
An object of class \sQuote{rank.select}, with \sQuote{print} and \sQuote{summary methods}, containing among other the matrices of AIC/BIC/HQ, the Likelihood, and best ranks according to each criterion. 
}


\seealso{
\code{\link{VECM}} for estimating a VECM. 
\code{\link[urca]{ca.jo}} in package \pkg{urca} for the classical Johansen cointegration test. 
}

\author{Matthieu Stigler}

\references{
- Aznar A and Salvador M (2002). Selecting The Rank Of The Cointegration Space And The Form Of The Intercept Using An Information
Criterion. Econometric Theory, *18*(04), pp. 926-947. <URL: http://ideas.repec.org/a/cup/etheor/v18y2002i04p926-947_18.html>.

-Cheng X and Phillips PCB (2009). Semiparametric cointegrating rank selection. Econometrics Journal , *12*(s1), pp. S83-S104. <URL:
http://ideas.repec.org/a/ect/emjrnl/v12y2009is1ps83-s104.html>.

- Gonzalo J and Pitarakis J (1998). Specification via model selection in vector error correction models. Economics Letters, *60*(3),
pp. 321 - 328. ISSN 0165-1765, <URL: http://dx.doi.org/DOI: 10.1016/S0165-1765(98)00129-3>.

- Kapetanios G (2004). The Asymptotic Distribution Of The Cointegration Rank Estimator Under The Akaike Information Criterion.
Econometric Theory, *20*(04), pp. 735-742. <URL: http://ideas.repec.org/a/cup/etheor/v20y2004i04p735-742_20.html>.

- Wang Z and Bessler DA (2005). A Monte Carlo Study On The Selection Of Cointegrating Rank Using Information Criteria. Econometric
Theory, *21*(03), pp. 593-620. <URL: http://ideas.repec.org/a/cup/etheor/v21y2005i03p593-620_05.html>.
}

\examples{

data(barry)

# 
rk_sel <- rank.select(barry)
rk_sel
summary(rk_sel)

}

\keyword{ ts }
\keyword{ VECM }
\keyword{ cointegration }


