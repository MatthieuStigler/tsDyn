print.nlVar<-function(object,...){
	if(object$model.specific$nthresh==0) 
		cat("Linear VAR model\n")
	else
		cat("\n\nNon Linear Model\n")
}

### logLik.VAR see: Luetkepohl, 3.4.5 (p. 89), Juselius (2006) p. 56. Hamilton 11.1.10, p. 293 gives -t/2 log(solve(S))
logLik.nlVar <- function(object,...){
	resids<-object$residuals
	k<-object$k
	t<-object$t
	Sigma<-matrix(1/t*crossprod(resids),ncol=k)
# 	res <- -(t*k/2)*log(2*pi) - (t/2)* log(det(Sigma)) -1/2 *sum(diag(resids %*% solve(Sigma) %*% t(resids)))
	res <- -(t*k/2)*log(2*pi) -t*k/2 - (t/2)* log(det(Sigma)) 
	return(res)
}

### logLik.VECM see Hamilton 20.2.10, p. 637
logLik.VECM <- function(object,r,...){
  t<-object$t
  k<-object$k
  
  if(object$model.specific$estim=="ML"){
    S00<-object$model.specific$S00/t
    lambda<-object$model.specific$lambda
    Rank <- if(missing(r)) object$model.specific$r else r
    seq<-if(Rank==0) 0 else if(Rank%in%1:k) 1:Rank else warning("r cann't be greater than k (number of variables)")
    res <- -(t*k/2)*log(2*pi) - t*k/2 - (t/2)*log(det(S00)) - (t/2)*sum(log(1-lambda[seq]))
  } else {
    Sigmabest<-matrix(1/t*crossprod(object$residuals),ncol=k)
    res <- log(det(Sigmabest))
    if(!missing(r)) warning("Note this is computing the LL from a model estimated by 2 OLS\n")
  }
  return(res)
}

#### Small function: get number of estimated parameters
npar  <- function (object, ...)  
  UseMethod("npar")
 
npar.default<-function(object, ...) NULL

npar.nlVar<-function(object, ...) 
  object$npar+object$model.specific$nthresh

npar.VECM<-function(object, ..., r) {
  nVar<-object$k
  Rank<-if(missing(r)) object$model.specific$r else r
  slopePars <- prod(dim(coef(object)[,-grep("^ECT[0-9]*$", colnames(coef(object)))])) ## get numb of al params but the alpha (ECT)
  nPar <- slopePars+2*nVar*Rank- Rank^2
  return(nPar)
}


#### AIC criterions
AIC.nlVar<-function(object,..., k=2, fitMeasure=c("SSR", "LL")){
	fitMeasure <- match.arg(fitMeasure)
	t<-object$t
	fit <- if(fitMeasure=="LL") -2*logLik.nlVar(object) else log(det(crossprod(residuals(object))/t))
	t*fit+k*npar(object)
}

AIC.VECM<-function(object,..., k=2,r, fitMeasure=c("SSR", "LL")){
	fitMeasure <- match.arg(fitMeasure)
	Rank<-if(missing(r)) object$model.specific$r else r
	t<-object$t
	fit <- if(fitMeasure=="LL") -2*logLik.VECM(object,r=Rank) else log(det(crossprod(residuals(object))/t))
	t*fit+k*npar(object, r=Rank)
}

#### BIC criterions
BIC.nlVar<-function(object,..., k=log(object$t), fitMeasure=c("SSR", "LL")){
	fitMeasure <- match.arg(fitMeasure)
	t<-object$t
	fit <- if(fitMeasure=="LL") -2*logLik.nlVar(object) else log(det(crossprod(residuals(object))/t))
	t*fit+k*npar(object)
}

BIC.VECM<-function(object,..., k=log(object$t),r, fitMeasure=c("SSR", "LL")){
	fitMeasure <- match.arg(fitMeasure)
	nVar<-object$k
	Rank<-if(missing(r)) object$model.specific$r else r
	t<-object$t
	fit <- if(fitMeasure=="LL") -2*logLik.VECM(object,r=Rank) else log(det(crossprod(residuals(object))/t))
	t*fit+k*npar(object, r=Rank)
}

deviance.nlVar<-function(object,...){
	as.numeric(crossprod(c(object$residuals)))
}

residuals.nlVar<-function(object,...){
	object$residuals
}

fitted.nlVar<-function(object,...){
	object$fitted
}

coef.nlVar<-function(object,...){
	return(object$coefficients)
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

