\name{VECM}
\alias{VECM}
\title{Estimation of Vector error correction model (VECM)}
\usage{
  VECM(data, lag, r = 1,
    include = c("const", "trend", "none", "both"),
    beta = NULL, estim = c("2OLS", "ML"),
    LRinclude = c("none", "const", "trend", "both"),
    exogen = NULL)
}
\arguments{
  \item{data}{multivariate time series (first row being
  first=oldest value)}

  \item{lag}{Number of lags (in the VECM representation,
  see Details)}

  \item{r}{Number of cointegrating relationships}

  \item{include}{Type of deterministic regressors to
  include}

  \item{beta}{Possibility to impose a cointegrating value.
  By default is null, so values will be estimated}

  \item{LRinclude}{Type of deterministic regressors to
  include in the long-term relationship. Can also be a
  matrix with exogeneous regressors (2OLS only).}

  \item{estim}{Type of estimator: \code{2OLS} for the
  two-step approach or \code{ML} for Johansen MLE}

  \item{exogen}{Inclusion of exogenous variables (first row
  being first=oldest value). Is either of same size than
  data (then automatically cut) or than end-sample.}
}
\value{
  An object of class \code{VECM} (and higher classes
  \code{VAR} and \code{nlVar}) with methods: \describe{
  \item{Usual methods}{Print, summary, plot, residuals,
  fitted, vcov} \item{Fit criteria}{AIC, BIC,
  \code{\link{MAPE}}, \code{\link{mse}},
  \code{\link[=logLik.VECM]{logLik}} (latter only for
  models estimated with MLE)} \item{Prediction}{Predict and
  \code{\link{predict_rolling}}} \item{VAR/VECM
  methods}{Impulse response function
  (\code{\link[=irf.nlVar]{irf}}) and forecast error
  variance decomposition (\code{\link[=fevd.nlVar]{fevd}})}
  \item{LaTeX}{toLatex} }
}
\description{
  Estimate either a VECM by Engle-Granger or Johansen (MLE)
  method.
}
\details{
  This function is just a wrapper for the
  \code{\link{lineVar}}, with model="VECM".

  More comprehensive functions for VECM are in package
  \pkg{vars}. A few differences appear in the VECM
  estimation: \describe{ \item{Engle-Granger estimator}{The
  Engle-Granger estimator is available}
  \item{Presentation}{Results are printed in a different
  ways, using a matrix form} \item{lateX export}{The matrix
  of coefficients can be exported to latex, with or without
  standard-values and significance stars}
  \item{Prediction}{The \code{predict} method contains a
  \code{newdata} argument allowing to compute rolling
  forecasts.} }

  Two estimators are available: the Engle-Granger two step
  approach (\code{2OLS}) or the Johansen (\code{ML}). For
  the 2OLS, deterministics regressors (or external
  variables if LRinclude is of class numeric) can be added
  for the estimation of the cointegrating value and for the
  ECT. This is only working when the beta value is not
  pre-specified.

  The arg beta is the cointegrating value, the
  cointegrating vector will be taken as: (1, -beta).

  Note that the lag specification corresponds to the lags
  in the VECM representation, not in the VAR (as is done in
  package vars or software GRETL). Basically, a VAR with 2
  lags corresponds here to a VECM with 1 lag. Lag 0 in the
  VECM is not allowed.
}
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
   all.equal(VECM(finland, lag=2, estim="ML", r=2)$model.specific$coint,
             cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)  $beta, check.attr=FALSE)
 # check OLS parameters
   all.equal(t(coefficients(VECM(finland, lag=2, estim="ML", r=2))),
             coefficients(cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)$rlm), check.attr=FALSE)

}


##export to Latex
toLatex(vecm.eg)
toLatex(summary(vecm.eg))
options("show.signif.stars"=FALSE)
toLatex(summary(vecm.eg), parenthese="Pvalue")
options("show.signif.stars"=TRUE)
}
\author{
  Matthieu Stigler
}
\seealso{
  \code{\link{lineVar}} \code{\link{TVAR}} and
  \code{\link{TVECM}} for the correspoding threshold
  models. \code{\link{linear}} for the univariate AR model.
}
\keyword{ts}

