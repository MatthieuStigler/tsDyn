print.nlVar<-function(object,...){
	if(object$model.specific$nthresh==0) 
		cat("Linear VAR model\n")
	else
		cat("\n\nNon Linear Model\n")
}

### logLik.VAR see: Luetkepohl, 3.4.5 (p. 89), Juselius (2006) p. 56. Hamilton 11.1.10, p. 293 gives -t/2 log(solve(S))
logLik.nlVar <- function(object,...){
	res<-object$residuals
	k<-object$k
	t<-object$t
	Sigma<-matrix(1/t*crossprod(res),ncol=k)
	res <- -(t*k/2)*log(2*pi) - (t/2)* log(det(Sigma)) -1/2 *sum(diag(res %*% solve(Sigma) %*% t(res)))
	return(res)
}

### logLik.VECM see Hamilton 20.2.10, p. 637
logLik.VECM <- function(object,r=1,...){
  t<-object$t
  k<-object$k
  
  if(object$model.specific$estim=="ML"){
    S00<-object$model.specific$S00
    lambda<-object$model.specific$lambda
    seq<-if(r==0) 0 else if(r%in%1:k) 1:r else warning("r cann't be greater than k (number of variables)")
    
    res <- -(t*k/2)*log(2*pi) - t*k/2 -(t/2)*log(det(S00))-(t/2)*sum(log(1-lambda[seq]))
  } else {
    Sigmabest<-matrix(1/t*crossprod(object$residuals),ncol=k)
    res <- log(det(Sigmabest))
  }
  return(res)
}

AIC.nlVar<-function(object,..., k=2){
	t<-object$t
	t*logLik.nlVar(object)+k*(object$npar+object$model.specific$nthresh)
}

AIC.VECM<-function(object,..., k=2,r){
	kVar<-object$k
	Rank<-if(missing(r)) object$model.specific$r else r
	t<-object$t
#formula from Gonzalo pitarakis is p^2(k-1) +2pr -r^2+constan, here npar already contains pr
	nParFree<- object$npar+kVar*Rank- Rank^2
	t*logLik.VECM(object,r=Rank)+k*nParFree 
}

BIC.nlVar<-function(object,..., k=log(object$t)){
	t<-object$t
	t*logLik.nlVar(object)+k*(object$nparB+object$model.specific$nthresh)
}

BIC.VECM<-function(object,..., k=log(object$t),r){
	kVar<-object$k
	Rank<-if(missing(r)) object$model.specific$r else r
	t<-object$t
	#formula from Gonzalo pitarakis is p^2(k-1) +2pr -r^2+constan, here npar already contains pr
	nParFree<- object$npar+kVar*Rank- Rank^2
	t*logLik.VECM(object,r=Rank)+k*nParFree 
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

