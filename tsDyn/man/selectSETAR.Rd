\name{selectSETAR}
\alias{selectSETAR}
\alias{selectSetar}
\alias{selectsetar}
\title{Automatic selection of SETAR hyper-parameters}
\description{
Automatic selection of SETAR hyper-parameters
}
\usage{
selectSETAR(x, m, d=1, steps=d, series, mL, mH,mM, thDelay=0, mTh, thVar, 
	    th=MakeThSpec(), trace=TRUE, include = c("const", "trend","none", "both"), 
	    common=c("none", "include","lags", "both"), model=c("TAR", "MTAR"), 
	    ML=seq_len(mL), MH=seq_len(mH), MM=seq_len(mM), nthresh=1, trim=0.15,
	    criterion = c("pooled-AIC", "AIC","BIC", "SSR"),  plot=TRUE,
	    max.iter=2, type=c("level", "diff", "ADF"), same.lags=FALSE, 
	    restriction=c("none","OuterSymAll","OuterSymTh"),  hpc=c("none", "foreach") )
}

\arguments{
  \item{x}{ time series }
  \item{m, d, steps}{ embedding parameters. For their meanings, see help about \code{\link{nlar}} }
\item{series}{ time series name (optional) }
\item{mL, mH,mM}{autoregressive order for \sQuote{low} (\code{mL}) \sQuote{middle} (\code{mM}, only useful if \code{nthresh}=2) and \sQuote{high} (\code{mH}) regime (default values: m). Must be <=m. Alternatively, you can specify \code{ML}}
\item{thDelay}{ Vector of possible \sQuote{threshold delay} values to check for}
\item{mTh}{coefficients for the lagged time series, to obtain the threshold variable}
\item{thVar}{external threshold variable}
\item{th}{Different specifications of the grid search, to pre-specify a value or set the number of points to search. See\code{\link{MakeThSpec}}}
\item{trace}{should additional infos be printed? (logical)}
\item{include}{Type of deterministic regressors to include}
\item{common}{Indicates which elements are common to all regimes: no, only the \code{include} variables, the lags or both}
\item{model}{Currently not implemented}
\item{ML,MM,MH}{vector of lags for  order for \sQuote{low} (ML) \sQuote{middle} (MM, only useful if nthresh=2) and \sQuote{high} (MH)regime. Max must be <=m}
\item{nthresh}{Number of threshold of the model}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime. Default to 0.15}
\item{criterion}{ Model selection criterion}
\item{plot}{Should a plot showing the criterion values be printed? (logical)}
\item{max.iter}{Number of iterations for the algorithm}
\item{type}{Whether the variable is taken is level, difference or a mix (diff y= y-1, diff lags) as in the ADF test}
\item{same.lags}{Logical. When AIC or pooled-AIC is used and arg m is given, should it search for same number of lags in each regime (TRUE) or allow for different (FALSE) lags in each regime. Different lags involves more computation}
\item{restriction}{Restriction on the threshold. \code{OuterSymAll} will take a symmetric threshold and symmetric coefficients for outer regimes. OuterSymTh currently unavailable}
\item{hpc}{Possibility to run the bootstrap on parallel core. See details}
}
\details{
  Routine for automatic selection of SETAR models hyper parameters.
  
  An exhaustive search over all possible combinations of values of
  specified hyper-parameters is performed. Thus the threshold delay, the number of lags in each regime and the threshold value are computed. 
  
  Embedding parameters \code{d,steps} are kept fixed.

  Possible criteria are the usual SSR, AIC and a pooled AIC formula:
  \eqn{AIC(low regime model) + AIC(high regime model)}. The default
  criterion is the pooled AIC formula. SSR criterion can't be used to compare models with different lags. 

When two thresholds(\code{nthresh}=2) have to be computed, the search for the second is made conditional on results for first threshold as suggested in Gonzalo and Pittarakis (2002). Refinements can be obtained by using \code{max.iter} (first threshold being re-estimated based on the second one). If SSR is used, the number of lags in the inner regime is either the same if only arg m was given, otherwise it has to be pre-specified. Criterion AIC can be used to determine the number of lags in the nner regime, whereas pooled-AIC is currently not implemented for nthresh=2. 

By default, all threshold values excluding the upper and lower \code{trim} of the threshold values are taken as potential threshold. restriction can be made with arg \code{th}. See function \code{\link{MakeThSpec}}. 

With the argument \code{hpc}, the heavy grid search can be run on parallel cores, thus alleviating the time of computation. Preliminary results 
indicate however that the length of the series must be very considerable in order that the parallel code becomes advantageous. To use it, the user needs simply to choose a package (among doMC, doMPI, doSNOW or doRedis) and register the backend. See the vignette for more details.
}
\value{
An object of class \code{selectSETAR} (print and plot methods) with: 
  \item{res}{A data-frame, with columns giving hyper-parameter values and the
  computed AIC for each row (only the best 10/5s are returned)}
\item{res2}{Same as res, returned if nthresh=2 otherwise set to NULL}
\item{bests}{estimated hyper-parameters}
\item{th, firstBests, bests2th, ML, MM, MH}{estimated parameters, from first and conditional search}
\item{criterion, nthresh,same.lags}{returns args given by user}
\item{allTh}{all threshold values and corresponding criterion from first search}
}
\seealso{
  \code{\link{selectLSTAR}}, \code{\link{selectNNET}}, \link{MakeThSpec}
}
\references{Gonzalo, J. & Pitarakis, J. (2002) Estimation and model selection based inference in single and multiple threshold models, Journal of Econometrics, 110, 319 - 352}
\author{ Antonio, Fabio Di Narzo and Stigler, Matthieu}
\examples{
llynx <- log10(lynx)
selectSETAR(llynx, m=2)
#Suggested model is the following:
setar(llynx, m=2, thDelay=1, th=3.4)
}
\keyword{ ts }
