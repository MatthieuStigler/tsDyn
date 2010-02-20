TiVaryVECM<-function(data, lag,m=2, r=1, include=c("const","none") , ecdet=c("none", "const"))
{

include<-match.arg(include)
ecdet<-match.arg(ecdet)
if(ecdet=="const")
  include<-"none"

  y <- as.matrix(data)
  Torigin <- nrow(y) 	#Size of original sample
  T <- nrow(y) 		#Size of start sample
if(m<1 | m>T) stop("M Should be greater than 0 and less than T\n")  
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
  t<-T-p -1		#Size of end sample
  
  if(is.null(colnames(data)))
    colnames(data)<-paste("Var", c(1:k), sep="")

  Y <- y[(p+1):T,] #
  X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix
  
#Set up of dependant and independant variables matrices
  DeltaY<-diff(y)[(p+1):(T-1),]
  DeltaX<-embed(diff(y),p+1)[,-(1:k)]

  cheb<-function(i,T){
    t<-1:T
    if(i==0) 1
    else sqrt(2)*cos(i*pi*(t-0.5)/T)
  }
  
  
  Z<-DeltaX
  Y<-DeltaY
  Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
if(include=="const")
	Z<-cbind(1, Z)
  
##VECM: ML (Johansen ) estimation of cointegrating vector
#Auxiliary regression 1
  reg_res1<-lm.fit(Z,Y)
  u<-residuals(reg_res1)
#Auxiliary regression 2
  reg_res2<-lm.fit(Z,Xminus1)
  v<-residuals(reg_res2)
#Auxiliary regression 3
if(ecdet=="const"){
  reg_res3<-lm.fit(Z,matrix(1, nrow=nrow(Z)))
  v<-cbind(v,residuals(reg_res3))
}
#Moment matrices
  S00<-crossprod(u)
  S11<-crossprod(v)
  S01<-crossprod(u,v)
  SSSS<-solve(S11)%*%t(S01)%*%solve(S00)%*%S01
  eig<-eigen(SSSS)
  ve<-eig$vectors
  va<-eig$values
#normalize eigenvectors
  ve_no<-apply(ve,2, function(x) x/sqrt(t(x)%*%S11%*%x))
  ve_2<-t(t(ve_no)/diag(ve_no)) 
  ve_3<-ve_2[,1:r, drop=FALSE]
  C2 <- matrix(0, nrow = nrow(ve_2) - r, ncol = r)
  C <- rbind(diag(r), C2)
  ve_4 <- ve_3 %*% solve(t(C) %*% ve_3)

#compute A (speed adjustment)
  z0<-t(u)%*%v%*%ve_no[,1:r]%*%t(ve_no[,1:r])

  ###Slope parameters
if(ecdet=="const"){
  ECTminus1<-cbind(Xminus1,1)%*%ve_4
}else{
  ECTminus1<-Xminus1%*%ve_4}
Z_final<-cbind(ECTminus1,Z)
B<-t(Y)%*%Z_final%*%solve(t(Z_final)%*%Z_final)		#B: OLS parameters, dim 2 x npar

# print(B)
# print(-S01%*%ve_no)

###naming of variables and parameters
rownames(ve_4)<-c(colnames(data), if(ecdet=="const") "const" else NULL)
  npar<-ncol(B)*nrow(B)
  rownames(B)<-paste("Equation",colnames(data))
  LagNames<-c(paste(rep(colnames(data),length(Lags)), -rep(Lags, each=k)))
  ECT<-paste("ECT", 1:r, sep="")

if(include=="const")
	Bnames<-c(ECT,"Intercept", LagNames)
else
	Bnames<-c(ECT,LagNames)
  colnames(B)<-Bnames

### Smoothed VECM
if(FALSE){
  Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
  Xminus1t<-Xminus1
  for(i in 1:m)
    Xminus1t<-cbind(Xminus1t, cheb(i,t)*Xminus1)
  
  Sig_xDy<-crossprod(Z,Y)
  Sig_xx_inv<-solve(crossprod(Z))
  Sig_xym<-crossprod(Z,Xminus1t)
  Sig_ymym<-crossprod(Xminus1t)
  Sig_DyDy<-crossprod(Y)
  Sig_Dyym<-crossprod(Y, Xminus1t)
  
  S00t<- Sig_DyDy-t(Sig_xDy)%*%Sig_xx_inv%*%Sig_xDy
  S11t<- Sig_ymym-t(Sig_xym)%*%Sig_xx_inv%*%Sig_xym
  S01t<- Sig_Dyym-t(Sig_xDy)%*%Sig_xx_inv%*%Sig_xym
  S10t<-t(S01t)
  
  m_t<-solve(S11t)%*%S10t%*%solve(S00t)%*%S01t	#m2x2 #m=inv(v'*v)*(v'*u)*inv(uu)*(u'*v);
  ve_t<-eigen(m_t)$vectors 			#eigenvectors 2x2		#[ve,va]=eig(m);
  va_t<-eigen(m_t)$values				#Eigenvalues 2
  maa<-which.max(va_t)				#Selection of the biggest eigenvalue	#   [temp,maa]=max(va);
  h<-ve_t[,maa]					#Eigenvector of the Biggest eigenvalue		#h=ve(:,maa);
  betaLT<- -h[2]/h[1]				#Normalization		#b0= -h(2)/h(1);
  
#   ECTminus1<-Xminus1%*%c(1,-betaLT)
#   Z<-cbind(ECTminus1,Z)


}
###Y and regressors matrix
  fitted<-Z_final%*%t(B)
  res<-Y-fitted
  naX<-rbind(matrix(NA, ncol=ncol(Z_final), nrow=T-t), Z_final)
  rownames(naX)<-rownames(data)
  colnames(naX)<-Bnames
  YnaX<-cbind(data, naX)
  

###Return outputs
  model.specific<-list()
  model.specific$nthresh<-0
  model.specific$r<-r
  model.specific$S00<-S00
  model.specific$lambda<-eig$values
  

  z<-list(residuals=res,  coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", fitted.values=fitted, model.x=Z_final,lag=lag, model=YnaX, df.residual=t-npar/k, model.specific=model.specific, coint=ve_4)
  class(z)<-c("VaryVECM","VECM", "nlVar")
  
  attr(z, "varsLevel")<-"diff"
  return(z)
}
####FIN
print.VaryVECM<-function(x) {
  cat("Cointegration parameters\n")
  print(x$coint)
  cat("\n VECM parameters\n")
  print(x$coefficients)}

logLik.VaryVECM<-function(object,...) {
  T<-object$T
  k<-object$k
  mod<-object$model.specific
-(T*k/2)*log(2*pi)-(T*k/2)-(T/2)*log(det(mod$S00))-(T/2)*sum(log(1-mod$lambda))
}

library(urca)
library(vars)
library(tsDyn)
data(zeroyld)
data(finland)

t<-TiVaryVECM(finland[,1:2], lag=2)

logLik(t)
###test


TiVaryVECM(finland[,1:2], lag=2)
TiVaryVECM(finland[,1:2], lag=2, include="none")
TiVaryVECM(finland[,1:2], lag=2, ecdet="const")
coefficients(cajorls(ca.jo(finland[,1:2], K=3, spec="transitory", ecdet="const"))$rlm)

#3 var, r=1
all.equal(TiVaryVECM(finland[,1:3], lag=2)$coint, cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"))$beta, check.attr=F)
logLik(TiVaryVECM(finland[,1:3], lag=2))
logLik(vec2var(ca.jo(finland[,1:3], K=3, spec="transitory")))





all.equal(t(coefficients(TiVaryVECM(finland[,1:3], lag=2))), coefficients(cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"))$rlm), check.attr=F)

#3 var, r=2
all.equal(TiVaryVECM(finland[,1:3], lag=2, r=2)$coint, cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"),r=2)$beta, check.attr=F)
all.equal(t(coefficients(TiVaryVECM(finland[,1:3], lag=2,r=2))), coefficients(cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"),r=2)$rlm), check.attr=F)

#4 var, r=3
all.equal(TiVaryVECM(finland, lag=2, r=3)$coint, cajorls(ca.jo(finland, K=3, spec="transitory"),r=3)$beta, check.attr=F)
all.equal(t(coefficients(TiVaryVECM(finland, lag=2,r=3))), coefficients(cajorls(ca.jo(finland, K=3, spec="transitory"),r=3)$rlm), check.attr=F)

#4 var, no int
all.equal(TiVaryVECM(finland, lag=2, r=3,ecdet="const")$coint, cajorls(ca.jo(finland, K=3, spec="transitory", ecdet="const"),r=3)$beta, check.attr=F)
all.equal(t(coefficients(TiVaryVECM(finland, lag=2,r=3,ecdet="const"))), coefficients(cajorls(ca.jo(finland, K=3, spec="transitory",ecdet="const"),r=3)$rlm), check.attr=F)

TiVaryVECM(finland[,1:3], lag=2)$coint
cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"))$beta
coefficients(TiVaryVECM(finland[,1:3], lag=2))
coefficients(cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"))$rlm)


TiVaryVECM(finland[,1:3], lag=2, r=2)$coint
cajorls(ca.jo(finland[,1:3], K=3, spec="transitory"), r=2)$beta


dec<-function(fun, T, m){
  cheb<-function(i,T){
    t<-1:T
    if(i==0) 1
    else sqrt(2)*cos(i*pi*(t-0.5)/T)
  }
  xi<-function(fun,m) sum(fun(T)*cheb(m,T))/T
  a<-0
  m<-if(missing(m)) T else m
  for(k in 0:(m-1)) a<-xi(gt,k )*cheb(k,T) +a
  a
}


mycheb<-function(y,p,chebdim){
  T<-nrow(y)
  if(chebdim==0){
    ret<-y
  }else {
    mat<-y[(p+1):(T-1),]
    n<-nrow(y)-p-1
    s<-seq(p+2, length.out=n)
    for(i in 1:chebdim)    mat<-cbind(mat, sqrt(2)*cos(i*pi*(s-0.5)/n)*y[(p+1):(T-1),])
    ret<-mat
  }
  return(ret)
}













if(FALSE){
  gt<-function(t) 1:t
  dec(gt,100)
  dec(gt,100, 90)
}



if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
environment(lineVar)<-environment(star)
environment(summary.VAR)<-environment(star)
data(zeroyld)
dat<-zeroyld

#tests
aVAR<-lineVar(dat[1:100,], lag=c(1,2), include="both", model="VAR")
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
	Sigma<-matrix(1/(object$df.residual)*crossprod(x$residuals),ncol=k)
	cov.unscaled<-solve(crossprod(x$model.x))
	VarCovB<-cov.unscaled%x%Sigma
	StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)
	Tvalue<-x$coefficients/StDevB

	Pval<-pt(abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=TRUE)
	#Pval<-round(Pval,4)
	symp <- symnum(Pval, corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","** ","*  ",".  ","    "))
	stars<-matrix(symp, nrow=nrow(Pval))
	ab<-matrix(paste(myformat(x$coefficients,digits),"(", myformat(StDevB,digits),")",stars,sep=""), nrow=nrow(Pval))
	dimnames(ab)<-dimnames(x$coefficients)		

	x$bigcoefficients<-ab
	x$cov.unscaled<-cov.unscaled
	x$sigma<-Sigma
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
	cat("#############\n###Model", attr(x,"model"),"\n#############")
	cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
	cat("\nNumber of variables:", x$k,"\tNumber of estimated slope parameters", x$npar)
	cat("\nAIC",x$aic , "\tBIC", x$bic, "\tSSR", x$SSR)
	if(attr(x,"model")=="VECM")
		cat("\nCointegrating vector: (1, -", x$model.specific$beta, ")")
	cat("\n\n")
	print(noquote(x$bigcoefficients))

}

vcov.VAR<-function(object, ...){
    sum<-summary.VAR(object)
    so<-sum$cov.unscaled%x%sum$sigma
    co.names<-gsub(" ", "", colnames(coef(object)))
    eq.names<-gsub("Equation ", "",rownames(coef(object)))
    together.names<-paste(rep(eq.names,each= length(co.names)), co.names, sep=":")
    dimnames(so)<-list(together.names, together.names)
    so
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
	###explained vector
	if(attr(x, "varsLevel")=="diff")
	  res[7]<-TeXVec(paste("slashDelta X_{t}^{",seq(1, x$k),"}", sep=""))
	else
	  res[7]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
	res[8]<- "\\end{smatrix}="
	###ECT
	ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
	if(attr(x,"model")=="VECM"){
		len<-length(res)
		res[len+1]<-"+\\begin{smatrix}  %ECT"
		res[len+2]<-TeXVec(coeftoprint[,1]) #see nlVar-methods.R
		res[len+3]<-"\\end{smatrix}ECT_{-1}"
	  }
	###Const
 	res<-include(x, res, coeftoprint)	#nlVar-methods.R
	###Lags
	a<-if(attr(x,"model")=="VECM") 1 else 0
	res<-LagTeX(res, x, coeftoprint, ninc+a)	#nlVar-methods.R
	res[length(res)+1]<-"\\end{equation}"
	res<-gsub("slash", "\\", res, fixed=TRUE)
	res<-res[res!="blank"]
	
	return(structure(res, class="Latex"))
}


if(FALSE){
###TODO
#check if const/trend/both in LR rel and VECM makes sense!
#check for standaard deviation of coint vector whith ML estim!
#consistency between ML and OLS coint estimator?
}


if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
data(zeroyld)
dat<-zeroyld
environment(lineVar)<-environment(star)
environment(summary.VAR)<-environment(star)

aVAR<-lineVar(dat, lag=1, include="both", model="VAR")
aVAR<-lineVar(dat, lag=1, include="const", model="VECM", estim="ML", beta=0.98)
#lag2, 2 thresh, trim00.05: 561.46
aVAR
summary(aVAR)
sqrt(diag(summary(aVAR, cov=0)$sigma))
vcov.VAR(aVAR)
vcovHC.VAR(aVAR)
logLik(aVAR)
AIC(aVAR)
BIC(aVAR)
deviance(aVAR)
coef(aVAR)
environment(toLatex.VAR)<-environment(star)
toLatex(aVAR)
toLatex(summary(aVAR))

###Check VAR: comparing with vars
myVAR<-lineVar(dat, lag=1)

library(vars)
var<-VAR(dat, lag=1)

vaco1<-coef(var)$short.run[c(3,1,2),1]
vaco2<-coef(var)$long.run[c(3,1,2),1]
round(coef(myVAR),8)==round(rbind(vaco1, vaco2),8)

###Check Johansen MLE
myVECM<-lineVar(dat, lag=1, include="const", model="VECM", estim="ML")
summary(myVECM, digits=7) 
#comparing with Hansen paper:reported in Gauss procedure is:
#coint vector: 1.02206: ok!
#coeff: 
#comparing with vars package
a<-ca.jo(dat, spec="trans")
summary(a)
#same answer also!
}

if(FALSE){
#alternative as in Hamilton
#compute other parameters
if(FALSE){
if(ecdet=="const"){
  last<-ncol(z0)-1
print(z0)
  z1<-z0[,last]
  z0<-z0[,-last]
# ve_4<-head(ve_4,-1)
}
  beta_reg1<-t(coefficients(reg_res1))
  beta_reg2<-t(coefficients(reg_res2))
  temp<-matrix(NA, nrow=nrow(beta_reg1), ncol=ncol(beta_reg1))
a<-if(include=="const") 1 else NULL
  temp[,-1]<-beta_reg1[,-1]- z0%*%beta_reg2[,-1]
  temp[,1]<-beta_reg1[,1]- z0%*%beta_reg2[,1]
#All parameters matrix
  B<-cbind(z0[,1:r], temp)
}
}