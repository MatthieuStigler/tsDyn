#'Extract Log-Likelihood
#'
#'Log-Likelihood method for VAR and VECM models.
#'
#'For a VAR, the Log-Likelihood is computed as in Luetkepohl (2006) equ. 3.4.5 (p. 89) and
#'Juselius (2006) p. 56:
#'
#'\deqn{ LL = -(TK/2) \log(2\pi) - (T/2) \log|\Sigma| - (1/2) \sum^{T} \left [
#'(y_t - A^{'}x_t)^{'} \Sigma^{-1} (y_t - A^{'}x_t) \right ] } Where
#'\eqn{\Sigma} is the Variance matrix of residuals, and \eqn{x_t} is the matrix
#'stacking the regressors (lags and deterministic).
#'
#'However, we use a computationally simpler version:
#'
#'\deqn{ LL = -(TK/2) \log(2\pi) - (T/2) \log|\Sigma| - (TK/2) }
#'
#'See Juselius (2006), p. 57.
#'
#'(Note that Hamilton (1994) 11.1.10, p. 293 gives \eqn{+ (T/2)
#'\log|\Sigma^{-1}|}, which is the same as \eqn{-(T/2) \log|\Sigma|)}.
#'
#'For VECM, the Log-Likelihood is computed in two different ways, depending on whether the
#'\code{VECM} was estimated with ML (Johansen) or 2OLS (Engle and Granger).
#'
#'When the model is estimated with ML, the LL is computed as in Hamilton (1994)
#'20.2.10 (p. 637):
#'
#'\deqn{ LL = -(TK/2) \log(2\pi) - (TK/2) -(T/2) \log|\hat\Sigma_{UU}| - (T/2)
#'\sum_{i=1}^{r} \log (1-\hat\lambda_{i}) } Where \eqn{\Sigma_{UU}} is the
#'variance matrix of residuals from the first auxiliary regression, i.e.
#'regressing \eqn{\Delta y_t} on a constant and lags, \eqn{\Delta y_{t-1},
#'\ldots, \Delta y_{t-p}}. \eqn{\lambda_{i}} are the eigenvalues from the
#'\eqn{\Sigma_{VV}^{-1}\Sigma_{VU}\Sigma_{UU}^{-1}\Sigma_{UV}}, see 20.2.9 in
#'Hamilton (1994).
#'
#'When the model is estimated with 2OLS, the LL is computed as: \deqn{ LL =
#'\log|\Sigma| }
#'
#'Where \eqn{\Sigma} is the variance matrix of residuals from the the VECM
#'model. There is hence no correspondence between the LL from the VECM computed
#'with 2OLS or ML.
#'
#'
#'@param object object of class \code{VAR} computed by \code{\link{lineVar}}, or 
#' class \code{VECM} computed by \code{\link{VECM}}.
#'@param \dots additional arguments to \code{logLik}.
#'@return Log-Likelihood value.
#'@aliases  logLik.VAR
#'@author Matthieu Stigler
#'@references Hamilton (1994) \emph{Time Series Analysis}, Princeton University
#'Press
#'
#'Juselius (2006) \emph{The Cointegrated VAR model: methodology and
#'Applications}, Oxford Univesity Press
#'
#'Luetkepohl (2006) \emph{New Introduction to Multiple Time Series Analysis},
#'Springer
#'@keywords ts
#'@examples
#'
#'data(zeroyld)
#'
#'#Fit a VAR
#'VAR <- lineVar(zeroyld, lag=1)
#'logLik(VAR)
#'
#'#'#Fit a VECM
#'vecm <- VECM(zeroyld, lag=1, r=1, estim="ML")
#'logLik(vecm)
#'
#' @export
logLik.nlVar <- function(object,...){
  resids<-object$residuals
  k<-object$k
  t<-object$t
  Sigma<-matrix(1/t*crossprod(resids),ncol=k)
  # 	res <- -(t*k/2)*log(2*pi) - (t/2)* log(det(Sigma)) -1/2 *sum(diag(resids %*% solve(Sigma) %*% t(resids)))
  res <- -(t*k/2)*log(2*pi) -t*k/2 - (t/2)* log(det(Sigma)) 
  return(res)
}



### logLik.VAR see: Luetkepohl, 3.4.5 (p. 89), Juselius (2006) p. 56. Hamilton 11.1.10, p. 293 gives -t/2 log(solve(S))
### logLik.VECM see Hamilton 20.2.10, p. 637

#'@rdname logLik.nlVar
#'@param r The cointegrating rank. By default the rank specified in the call to
#'\code{\link{VECM}}, but can be set differently by user.
#' @export
logLik.VECM <- function(object,r,...){
  t<-object$t
  k<-object$k
  
  if(object$model.specific$estim=="ML"){
    S00 <- object$model.specific$S00 / t
    lambda <- object$model.specific$lambda
    Rank <- if(missing(r)) object$model.specific$r else r
    seq<-if(Rank==0) 0 else if(Rank%in%1:k) 1:Rank else warning("r can't be greater than k (number of variables)")
    res <- -(t*k/2)*log(2*pi) - t*k/2 - (t/2)*log(det(S00)) - (t/2)*sum(log(1-lambda[seq]))
  } else {
    Sigmabest <- matrix(1 / t * crossprod(object$residuals), ncol = k)
    res <- log(det(Sigmabest))
    if(!missing(r)) warning("Note this is computing the LL from a model estimated by 2 OLS\n")
  }
  return(res)
}


getP <- function(object) UseMethod("getP")
getP.ca.jo <- function(object) object@P
getP.cajo.test <- function(object) ncol(object@Z0)

#' @export
logLik.ca.jo <- function(object,r,...){
  t<-nrow(object@Z0)
  k<-getP(object)
  lambda<-object@lambda
  
  ## compute S00:
  M00 <- crossprod(object@Z0)/t
  M11 <- crossprod(object@Z1)/t
  M01 <- crossprod(object@Z0, object@Z1)/t
  M10 <- crossprod(object@Z1, object@Z0)/t
  M11inv <- solve(M11)
  S00 <- M00 - M01 %*% M11inv %*% M10
  
  
  seq<-if(r==0) 0 else if(r%in%1:k) 1:r else warning("r cann't be greater than k (number of variables)")
  res <- -(t*k/2)*log(2*pi) - t*k/2 - (t/2)*log(det(S00)) - (t/2)*sum(log(1-lambda[seq]))
  return(res)
}

#' @export
logLik.cajo.test <- function(object,r,...) logLik.ca.jo(object=object, r=r,...)

#### Small function: get number of estimated parameters
npar  <- function (object, ...)  
  UseMethod("npar")

npar.default<-function(object, ...) 
  length(coef(object))

npar.nlar<-function(object, ...) 
  object$x

npar.nlVar<-function(object, ...) 
  object$npar+object$model.specific$nthresh

npar.VECM<-function(object, ..., r) {
  nVar<-object$k
  Rank<-if(missing(r)) object$model.specific$r else r
  slopePars <- prod(dim(coef(object)[,-grep("^ECT[0-9]*$", colnames(coef(object)))])) ## get numb of al params but the alpha (ECT)
  nPar <- slopePars+2*nVar*Rank- Rank^2 ## formula: 2mr-r^2 (Cheng Phillips 2009, equ 1.2)
  return(nPar)
}

#' @export
#### AIC criterions
AIC.nlVar<-function(object,..., k=2, fitMeasure=c("SSR", "LL")){
  fitMeasure <- match.arg(fitMeasure)
  t<-object$t
  fit <- if(fitMeasure=="LL") -2*logLik.nlVar(object) else t*log(det(crossprod(residuals(object))/t))
  fit+k*npar(object)
}

#' @export
AIC.VECM<-function(object,..., k=2,r, fitMeasure=c("SSR", "LL")){
  fitMeasure <- match.arg(fitMeasure)
  Rank<-if(missing(r)) object$model.specific$r else r
  t<-object$t
  fit <- if(fitMeasure=="LL") -2*logLik.VECM(object,r=Rank) else t*log(det(crossprod(residuals(object))/t))
  fit+k*npar(object, r=Rank)
}

#' @export
#### BIC criterions
BIC.nlVar<-function(object,..., k=log(object$t), fitMeasure=c("SSR", "LL")){
  fitMeasure <- match.arg(fitMeasure)
  t<-object$t
  fit <- if(fitMeasure=="LL") -2*logLik.nlVar(object) else t*log(det(crossprod(residuals(object))/t))
  fit+k*npar(object)
}

#' @export
BIC.VECM<-function(object,..., k=log(object$t),r, fitMeasure=c("SSR", "LL")){
  fitMeasure <- match.arg(fitMeasure)
  nVar<-object$k
  Rank<-if(missing(r)) object$model.specific$r else r
  t<-object$t
  fit <- if(fitMeasure=="LL") -2*logLik.VECM(object,r=Rank) else t*log(det(crossprod(residuals(object))/t))
  fit+k*npar(object, r=Rank)
}

#' @export
deviance.nlVar<-function(object, ...){
  as.numeric(crossprod(c(object$residuals)))
}

#' @export
residuals.nlVar<-function(object, initVal=FALSE, ...){
  resids <- object$residuals
  
  if(initVal) {
    n_init <- object$T - object$t
    resids <- rbind(matrix(NA,  nrow=n_init, ncol = object$k),
                    resids)
  } 
  resids
}

#'fitted method for objects of class nlVar, i.e. VAR and VECM models.
#'
#'Returns the fitted values of the model, either as computed in the model, or
#'back to the original series level.
#'
#'
#'In case of a VAR in differences, in ADF specification, or a VECM, the fitted
#'values are actually in differences. With the option \code{level="original"},
#'the function returns the series in the original level.
#'
#'For VAR in levels, the two arguments are evidently the same and hence it is
#'not taken into account, returning a warning.
#'
#'@aliases fitted fitted.nlVar
#'@param object An object of class \sQuote{nlVar}; generated by
#'\code{\link{VECM}} or \code{\link{lineVar}}.
#'@param level How to return the fitted values. See below.
#'@param \dots Currently not used.
#'@return A matrix.
#'@author Matthieu Stigler
#'@keywords regression
#'@examples
#'
#'
#'## estimate models
#'data(barry)
#'
#'ve <- VECM(barry, lag=2)
#'va <- lineVar(barry, lag=1)
#'va_diff <- lineVar(barry, lag=1, I="diff")
#'va_ADF <- lineVar(barry, lag=1, I="ADF")
#'
#'
#'## get fitted values:
#'tail(fitted(ve))
#'tail(fitted(ve, level="original"))
#'
#'tail(fitted(va))
#'tail(fitted(object=va, level="original"))
#'
#'tail(fitted(va_diff))
#'tail(fitted(object=va_diff, level="original"))
#'
#'tail(fitted(va_ADF))
#'tail(fitted(object=va_ADF, level="original"))
#'
#'

#' @export
fitted.nlVar <- function(object, level=c("model", "original"),...){
  
  level <- match.arg(level)
  mod <- ifelse(inherits(object, "VECM"), "VECM", "VAR")
  
  if(mod=="VAR"&&level=="original" &&attr(object, "varsLevel")=="level"){
    warning("level='original' has no effect for VAR models in levels")
    level <- "model"
  }
  
  if(level=="model"){
    res <- object$fitted
  } else {
    original.data <- object$model[-c(1:(object$T-object$t-1),object$T),1:object$k]
    series <- cbind(original.data, object$fitted)
    res<- original.data+ object$fitted
  }
  
  return(res)
}

## #' @param regime Optional, allows to extract regime-specific coefficients. 
#' @export
coef.nlVar <- function(object,  regime = c("all", "L", "M", "H"), ...){
  regime <-  match.arg(regime)
  co <- object$coefficients
  commonInter <- object$model.specific$nthresh >0 & !is.list(co)
  
  ## select only one regime if requested
  if(regime != "all") {
    if(regime=="M" & object$model.specific$nthresh!=2) stop("No M regime if nthresh %in% 1, 2")
    if(commonInter) {
      # inc_terms <- switch(object$include, 
      #                     "const" =1,
      #                     "both" = 2,
      #                     "trend" = 1,
      #                     "none" = 0)
      # get <- switch(regime, 
      #               "L" = 1:object$lag * object$k,
      #               "M" = object$lag * object$k)
      # co <- co[inc_terms + object$lag * object$k]
      # 
      stop("Not implemented yet!")
    } else {
      co <- switch(regime, 
                   "L" = co$Bdown, 
                   "H" = co$Bup, 
                   "M" = co$Bmiddle)
    }
  }
  co
}

### Method coefMat
coefMat <- function (object, ...)  
  UseMethod("coefMat")

coefMat.default<-function(object, ...)
  coefficients(object)

coefMat.nlVar<-function(object,...){
  if(inherits(object, "VAR"))
    return(object$coefficients)
  else
    return(object$coeffmat)
}


### Method coefVec
coefVec <- function (object, ...)  
  UseMethod("coefVec")

coefVec.default<-function(object, ...)
  coefficients(object)

coefVec.nlVar<-function(object,...){
  coef_M <- coefMat(object)
  coef_c <- matrix(t(coef_M), ncol=1, byrow = FALSE)
  rownames(coef_c) <- paste(colnames(coef_M), rep(rownames(coef_M), each=ncol(coef_M)))
  coef_c
}



###Method toMlm
toMlm<- function(x, ...) {
  UseMethod("toMlm")
}

toMlm.default <- function(x){
  lm(x$model)
}

toMlm.nlVar<-function(x){
  mod<-as.data.frame(x$model[-c(1:(x$T-x$t)),] )
  ix <- 1:x$k
  Yt<-as.matrix(mod[,ix])
  Ytminusi<-mod[,-ix]
  mlm<-lm(Yt ~.-1, Ytminusi)
  return(mlm)
}


### Method eqNames
eqNames <- function (object, ...)  
  UseMethod("eqNames")

eqNames.nlVar <- function(object)
  gsub("Equation ", "", rownames(coef(object)))


###Tolatex preliminary###
#########################
###Latex vector
TeXVec<-function(vec){
  d<-vec[1]
  for(i in 1:(length(vec)	-1))
    d<-paste(d,"slashslash",vec[i+1] )
  d
}

###LateX elements of R matrix
TeXMat<-function(mat, oneLine=FALSE){
  mat<-matrix(mat, ncol=ifelse(inherits(mat, "matrix"), ncol(mat), length(mat)))
  nr<-nrow(mat)
  nc<-ncol(mat)	
  d<-mat[,1]
  for(i in 1:(nc-1))
    d<-paste(d,"&",mat[,i+1])
  d[seq_len(nr-1)]<-paste(d[seq_len(nr-1)],"slashslash")
  d[nr]<-paste(d[nr], "")
  matrix(d, nrow=ifelse(oneLine,1,nr), ncol=1)
}
if(FALSE){
  a<-matrix(c(1,2,3,4,5,6), ncol=2)
  TeXMat(a)
}
###Function include
include<-function(x, res, coef, skip=0, mat="smatrix"){
  n<-length(res)
  res[(n+1):(n+5)]<-"blank"
  if(x$include=="const"){
    res[n+1]<-paste("\\begin{",mat, "}     %const", sep="")
    res[n+2]<-TeXVec(coef[,1+skip])
    res[n+3]<-paste("\\end{",mat,"}", sep="")}
  if(x$include=="trend"){
    res[n+1]<-paste("\\begin{",mat,"}     %trend", sep="")
    res[n+2]<-TeXVec(coef[,1+skip])
    res[n+3]<-paste("\\end{",mat,"}     %trend", sep="")}
  if(x$include=="both"){
    res[n+1]<-paste("\\begin{",mat, "}     %const", sep="")
    res[n+2]<-TeXVec(coef[,1+skip])
    res[n+3]<-paste("\\end{",mat,"}+\\begin{",mat,"}     %trend", sep="")
    res[n+4]<-TeXVec(coef[,2+skip])
    res[n+5]<-paste("\\end{",mat, "}t", sep="")
  }
  return(res)
}

###Function lag
LagTeX<-function(res, x, coef, skip,mat="smatrix"){
  if(attr(x, "varsLevel")=="diff")
    delta<-"slashDelta "
  else
    delta<-NULL
  for(j in 1:x$lag){
    nres<-length(res)
    res[nres+1]<-paste("+\\begin{",mat,"}      %Lag", j,sep="")
    for(i in 1:x$k){
      res[nres+i+1]<-TeXMat(coef[,seq_len(x$k)+(j-1)*x$k+skip])[i]}
    nres<-length(res)
    res[nres+1]<-paste("\\end{",mat,"}",sep="")
    res[nres+2]<-paste("\\begin{",mat,"}", sep="")
    res[nres+3]<-TeXVec(paste(delta,"X_{t-",j,"}^{",seq(1, x$k),"}", sep=""))
    res[nres+4]<-paste("\\end{",mat,"}", sep="")
  }
  res
}

