linear2<-function(data, lag, include = c( "const", "trend","none", "both"), model=c("VAR", "VECM"), I=c("level", "diff"),beta=NULL, coINT=FALSE, LRinclude=c("none", "const", "trend","both"), beta0=0)
{
y <- as.matrix(data)
Torigin <- nrow(y) 	#Size of original sample
T <- nrow(y) 		#Size of start sample

if(length(lag)==1){
  p <- lag
  notAllLags<-FALSE
  Lags<-1:p
}
else{
  notAllLags<-TRUE
  p<-max(lag)
  Lags<-lag
}

t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables
t<-T-p			#Size of end sample
ndig<-getndp(y)
if(is.null(colnames(data)))
	colnames(data)<-paste("Var", c(1:k), sep="")

include<-match.arg(include)
model<-match.arg(model)
I<-match.arg(I)
Y <- y[(p+1):T,] #
X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

if(notAllLags)
  X<-X[,sort(rep((Lags*k-k+1), k))+0:(k-1)]

DeltaY<-diff(y)[(p+1):(T-1),]
Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
DeltaX<-embed(diff(y),p+1)[,-(1:k)]

if(model=="VAR"){
	Z<-X
	Y<-Y}
if(model=="VECM"|I=="diff"){
	Z<-DeltaX
	Y<-DeltaY
	t<-t-1}

##Long-run relationship OLS estimation
if(model=="VECM"){
	#beta has to be estimated
	if(is.null(beta) ){
	  if(class(LRinclude)=="character"){
	    LRinclude<-match.arg(LRinclude)
	    if(LRinclude=="none")
	      LRplus<-rep(0,T)
	    else if(LRinclude=="const")
	      LRplus<- rep(1,T)
	    else if(LRinclude=="trend")
	      LRplus<- seq_len(T)
	    else if(LRinclude=="both")
	      LRplus<-cbind(rep(1,T),seq_len(T))
	  }
	   else if(class(LRinclude)%in%c("matrix", "numeric"))
	    LRplus<-LRinclude
	else
	   stop("Argument LRinclude badly given")
	if(class(LRinclude)=="character") {
	  if(LRinclude=="none"){
	    LRplusplus<-matrix(0, nrow=T, ncol=1)
	    coint<-lm(y[,1] ~  y[,-1]-1)}
	  else{
	    coint<-lm(y[,1] ~  y[,-1]-1+ LRplus)
	    LRplusplus<-as.matrix(LRplus)%*%coint$coef[-1]}}
	else{
	  coint<-lm(y[,1] ~  y[,-1]-1+ LRplus)
	  LRplusplus<-as.matrix(LRplus)%*%coint$coef[-1]}
	
	betaLT<-coint$coef[[1]]
	betaLT_std <- sqrt(diag(summary(coint)$sigma*summary(coint)$cov))
	}
	#beta is given by user
	else{
	  coint<-c(1, -beta)
	  betaLT<-beta
	}

# ECT<-y%*%c(1,-betaLT)
# ECT<-round(ECT,ndig)
#ECTminus1<-ECT[-c(1:p,T)]

	#ECTminus1<-Xminus1%*%c(1,-betaLT)
	ECTminus1<-residuals(coint)[-c(seq_len(p), T)]
	Z<-cbind(ECTminus1,Z)

#check if other way works
#print(cbind(residuals(coint)[-c(seq_len(p),T)], (y%*%c(1,-betaLT)-(LRplusplus))[-c(seq_len(p),T),]))
#print(all.equal(residuals(coint)[-c(seq_len(p),T)], (y%*%c(1,-betaLT)-(LRplusplus))[-c(seq_len(p),T),]))
}

###Regressors matrix


if(include=="const")
	Z<-cbind(1, Z)
else if(include=="trend")
	Z<-cbind(seq_len(t), Z)
else if(include=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)
##################
###Linear model
#################
 B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar


npar<-ncol(B)*nrow(B)

rownames(B)<-paste("Equation",colnames(data))
LagNames<-c(paste(rep(colnames(data),length(Lags)), -rep(Lags, each=k)))


if(model=="VECM")
	ECT<-"ECT"
else
	ECT<-NULL


if(include=="const")
	Bnames<-c("Intercept",ECT, LagNames)
else if(include=="trend")
	Bnames<-c("Trend",ECT, LagNames)
else if(include=="both")
	Bnames<-c("Intercept","Trend",ECT, LagNames)
else
	Bnames<-c(ECT,LagNames)

colnames(B)<-Bnames
fitted<-Z%*%t(B)
res<-Y-fitted

model.specific<-list()
model.specific$nthresh<-0
if(model=="VECM"){
	model.specific$betaLT<-betaLT
	model.specific$betaLT_std<-betaLT_std}


z<-list(residuals=res,  coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", fitted.values=fitted, model.x=Z, include=include,lag=lag, model=model, model.specific=model.specific)
if(model=="VAR")
  class(z)<-c("VAR","nlVar")
if(model=="VECM")
  class(z)<-c("VAR","VECM", "nlVar")
return(z)
}


if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
environment(linear2)<-environment(star)
data(zeroyld)
dat<-zeroyld

aVAR<-linear2(dat[1:100,], lag=c(1,2), include="both", model="VAR")
#lag2, 2 thresh, trim00.05: 561.46
class(aVAR)
aVAR
print(aVAR)
logLik(aVAR)
AIC(aVAR)
BIC(aVAR)
deviance(aVAR)
coef(aVAR)
summary(aVAR)
toLatex(aVAR)
toLatex(summary(aVAR))
}




print.VAR<-function(x,...){
	print(coef(x))
}

summary.VAR<-function(object, digits=4,...){
	x<-object
	r<-4
	t<-x$t
	k<-x$k

	Sigma<-matrix(1/t*crossprod(x$residuals),ncol=k)
	VarCovB<-solve(crossprod(x$model.x))%x%Sigma
	StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)

	Tvalue<-x$coefficients/StDevB

	Pval<-pt(abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=TRUE)
	#Pval<-round(Pval,4)
	symp <- symnum(Pval, corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","** ","*  ",".  ","    "))
	stars<-matrix(symp, nrow=nrow(Pval))
	ab<-matrix(paste(myformat(x$coefficients,digits),"(", myformat(StDevB,digits),")",stars,sep=""), nrow=nrow(Pval))
	dimnames(ab)<-dimnames(x$coefficients)		

	x$bigcoefficients<-ab
	x$Sigma<-Sigma
	x$StDev<-StDevB
	x$Pvalues<-Pval
	x$stars<-stars
	x$starslegend<-symp
	x$aic<-AIC.nlVar(x)
	x$bic<-BIC.nlVar(x)
	x$SSR<-deviance.nlVar(x)
	class(x)<-c("summary.VAR", "VAR")
	return(x)
}



print.summary.VAR<-function(x,...){
	cat("#############\n###Model", x$model,"\n#############")
	cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
	cat("\nNumber of variables:", x$k,"\tNumber of estimated slope parameters", x$npar)
	cat("\nAIC",x$aic , "\tBIC", x$bic, "\tSSR", x$SSR)
	if(x$model=="VECM")
		cat("\nCointegrating vector:", x$model.specific$betaLT)
	cat("\n\n")
	print(noquote(x$bigcoefficients))

}


toLatex.VAR<-function(object,..., digits=4, parenthese=c("StDev","Pvalue")){
	x<-object
	parenthese<-match.arg(parenthese)
	if(inherits(x,"summary.VAR")){
		a<-myformat(x$coefficients,digits, toLatex=TRUE)
		if(parenthese=="StDev")
			b<-myformat(x$StDev,digits,toLatex=TRUE)
		else if(parenthese=="Pvalue")
			b<-myformat(x$Pvalues,digits,toLatex=TRUE)
		if(getOption("show.signif.stars"))

			stars<-paste("^{",x$stars,"}", sep="")
		else
			stars<-NULL
		coeftoprint<-matrix(paste(a,"(",b,")",stars, sep=""),ncol=ncol(a), nrow=nrow(a))
	}#end if x is of class summary

	else{
		coeftoprint <-myformat(x$coefficients, digits, toLatex=TRUE)}
	varNames<-rownames(x$coefficients)
	res<-character()
	res[1]<-"%insert in the preamble and uncomment the line you want for usual /medium /small matrix"
	res[2]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\begin{pmatrix}}{\\end{pmatrix}} %USUAL"
	res[3]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\left(\\begin{smallmatrix}}{\\end{smallmatrix}\\right)} %SMALL"
	res[4]<-"%\\usepackage{nccmath} \\newenvironment{smatrix}{\\left(\\begin{mmatrix}}{\\end{mmatrix}\\right)} %MEDIUM"
	res[5]<-"\\begin{equation}"
	res[6]<- "\\begin{smatrix} %explained vector"
	if(inherits(x, "VECM"))
	  res[7]<-TeXVec(paste("slashDelta X_{t}^{",seq(1, x$k),"}", sep=""))
	else
	  res[7]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
	res[8]<- "\\end{smatrix}="
 	res<-include(x, res, coeftoprint)
	ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
	if(x$model=="VECM"){
		len<-length(res)
		res[len+1]<-"+\\begin{smatrix}  %ECT"
		res[len+2]<-TeXVec(coeftoprint[,ninc+1])
		res[len+3]<-"\\end{smatrix}ECT_{-1}"
		ninc<-ninc+1}
	res<-LagTeX(res, x, coeftoprint, ninc)
	res[length(res)+1]<-"\\end{equation}"
	res<-gsub("slash", "\\", res, fixed=TRUE)
	res<-res[res!="blank"]
	
	return(structure(res, class="Latex"))
}



