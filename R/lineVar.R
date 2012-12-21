lineVar<-function(data, lag, r=1,include = c( "const", "trend","none", "both"), model=c("VAR", "VECM"), 
		  I=c("level", "diff", "ADF"),beta=NULL, estim=c("2OLS", "ML"),
		  LRinclude=c("none", "const", "trend","both"), exogen=NULL){

  y <- as.matrix(data)
  Torigin <- nrow(y) 	#Size of original sample
  T <- nrow(y) 		#Size of start sample

  if(length(lag)==1){
    p <- lag
    notAllLags<-FALSE
    Lags<-1:p
  } else {
    notAllLags<-TRUE
    p<-max(lag)
    Lags<-lag
  }

  t <- T-p 		#Size of end sample
  k <- ncol(y) 		#Number of variables
  t<-T-p			#Size of end sample

  if(is.null(colnames(data)))
	  colnames(data)<-paste("Var", c(1:k), sep="")

###Check args
  include<-match.arg(include)
  LRinclude<-match.arg(LRinclude)
  if(lag<1) stop("Lag is minimum 1")
  if(LRinclude%in%c("const", "both"))  include<-"none"
  ninclude<-switch(include, "const"=1, "trend"=1,"none"=0, "both"=2)
  model<-match.arg(model)
  estim<-match.arg(estim)
  I<-match.arg(I)
  if(!r%in%1:(k-1)) stop("Arg r, the number of cointegrating relationships, should be between 1 and K-1\n")
  if(model=="VECM"&estim=="2OLS"&r>1){
    warning("Estimation of more than 1 coint relationship is not possible with estim '2OLS'. Switched to Johansen 'ML'\n")
    estim<-"ML"
  }

  minPara<-p*k+ninclude
  if(!t>minPara) stop("Not enough observations. Try reducing lag number\n")

###Construct variables
  Y <- y[(p+1):T,] #
  X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

  #Set up of dependant and independant variables matrices
  if(notAllLags)
    X<-X[,sort(rep((Lags*k-k+1), k))+0:(k-1)]

  DeltaY<-diff(y)[(p+1):(T-1),]
  Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
  DeltaX<-embed(diff(y),p+1)[,-(1:k)]

  if(model=="VAR"){
    if(I=="level"){
      Z<-X
      Y<-Y
    } else if(I=="diff"){
      Z<-DeltaX
      Y<-DeltaY
      t<-t-1
    } else if(I=="ADF"){
      Z<-cbind(Xminus1, DeltaX)
      Y<-DeltaY
      t<-t-1
    }
  } else if(model=="VECM"){
    Z<-DeltaX
    Y<-DeltaY
    t<-t-1
  }

###Regressors matrix
  if(include=="const")
    Z<-cbind(1, Z)
  else if(include=="trend")
    Z<-cbind(seq_len(t), Z)
  else if(include=="both")
    Z<-cbind(rep(1,t),seq_len(t), Z)


##VECM: Long-run relationship OLS estimation
  if(model=="VECM"&estim=="2OLS"){
	  #beta has to be estimated
    beta.estimated<-if(is.null(beta)) TRUE else FALSE
    if(is.null(beta) ){

    ## build LRplus: deterministic/exogeneous regressor in coint
      if(class(LRinclude)=="character"){
	LRplus <-switch(LRinclude, "none"=NULL,"const"=rep(1,T),"trend"=seq_len(T),"both"=cbind(rep(1,T),seq_len(T)))
	LRinc_name <- switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend"), "none"=NULL)
	LRinc_dim <- switch(LRinclude, "const"=1, "trend"=1, "both"=2, "none"=0)
      } else if(class(LRinclude)%in%c("matrix", "numeric")) {
	LRplus<-LRinclude
      } else{
	stop("Argument LRinclude badly given")
      }
    ## run coint regression
      if(LRinclude=="none"){
	cointLM<-lm(y[,1] ~  y[,-1]-1)
      } else {
	cointLM<-lm(y[,1] ~  y[,-1]-1+ LRplus)
	Xminus1 <- cbind(Xminus1, tail(LRplus,nrow(Xminus1)))
      }
      
      betaLT<-coint<-c(1,-cointLM$coef)
      betaLT_std <- c(1,summary(cointLM)$coef[,2])
      names(betaLT_std)<-c(colnames(data), LRinc_name)

## case beta pre-estimated
    } else {
      if(length(beta)!=k-1) stop("Arg 'beta' should be of length k-1")
      if(LRinclude!="none")
	warning("Arg LRinclude not taken into account when beta is given by user")
	LRinc_name <- NULL
	LRinc_dim <- 0
      coint<-c(1, -beta)
      betaLT<-c(1,-beta)
    }

    coint_export<-matrix(coint, nrow=k+LRinc_dim , dimnames=list(c(colnames(data),LRinc_name), "r1"))
    betaLT<-matrix(betaLT, nrow=k+LRinc_dim , dimnames=list(c(colnames(data),LRinc_name),"r1"))
    ECTminus1<-Xminus1%*%betaLT
    Z<-cbind(ECTminus1,Z)
  }

##VECM: ML (Johansen ) estimation of cointegrating vector
  else if(model=="VECM"&estim=="ML"){
  beta.estimated<-if(is.null(beta)) TRUE else FALSE
    if (is.null(beta)){
      #Auxiliary regression 1
      reg_res1<-lm.fit(Z,Y)
      u<-residuals(reg_res1)
      #Auxiliary regression 2
      reg_res2<-lm.fit(Z,Xminus1)
      v<-residuals(reg_res2)
      #Auxiliary regression 3
      if(LRinclude!="none"){
	add <- switch(LRinclude, "const"=matrix(1, nrow=nrow(Z)), "trend"=matrix(1:nrow(Z), nrow=nrow(Z)), "both"=cbind(1,1:nrow(Z)))
	reg_res3<-lm.fit(Z,add)
	v<-cbind(v,residuals(reg_res3)) # equ 20.2.46 in Hamilton 
      }
      #Moment matrices
      S00<-crossprod(u)
      S11<-crossprod(v)
      S01<-crossprod(u,v)
      SSSS<-solve(S11)%*%t(S01)%*%solve(S00)%*%S01
      eig<-eigen(SSSS)
      ve<-Re(eig$vectors)
      va<-Re(eig$values)
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
      if(LRinclude!="none"){
	ECTminus1<-cbind(Xminus1,add)%*%ve_4
      }else{
	ECTminus1<-Xminus1%*%ve_4
      }
      Z<-cbind(ECTminus1,Z)
      coin_ve_names <- switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend"), "none"=NULL)
      dimnames(ve_4)<-list(c(colnames(data), coin_ve_names), paste("r", 1:r, sep=""))
      betaLT<-ve_4
    }else{  #end beta to be estimated
      betaLT<-beta
      ECTminus1<-Xminus1%*%c(1,-betaLT)
      Z<-cbind(ECTminus1,Z)
    }
  }#end model=="VECM"&estim=="ML"


###Slope parameters, residuals and fitted
  B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
  fitted<-Z%*%t(B)
  res<-Y-fitted

###naming of variables and parameters
  npar<-ncol(B)*nrow(B)
  rownames(B)<-paste("Equation",colnames(data))
  LagNames<-c(paste(rep(colnames(data),length(Lags)), -rep(Lags, each=k)))
  if(I=="ADF") LagNames <- paste("D", LagNames,sep="_")
  ECT<- if(model=="VECM") paste("ECT", if(r>1) 1:r else NULL, sep="") else NULL
  Xminus1Names<- if(I=="ADF") paste(colnames(data),"-1",sep="") else NULL
  BnamesInter<-switch(include,"const"="Intercept","none"=NULL,"trend"="Trend","both"=c("Intercept","Trend"))
  Bnames<-c(ECT,BnamesInter,Xminus1Names, LagNames)
  colnames(B)<-Bnames

###Y and regressors matrix to be returned
  naX<-rbind(matrix(NA, ncol=ncol(Z), nrow=T-t), Z)
  rownames(naX)<-rownames(data)
  YnaX<-cbind(data, naX)
  colnames(YnaX)<-c(colnames(data),Bnames)

###Return outputs
  model.specific<-list()
  model.specific$nthresh<-0

  if(model=="VECM"){
    model.specific$beta<- betaLT
    model.specific$r<-r
    model.specific$estim<-estim
    model.specific$LRinclude<-LRinclude
    model.specific$beta.estimated<-beta.estimated

    if(estim=="ML"){
      model.specific$S00<-S00
      model.specific$lambda<-va
      model.specific$coint<-ve_4
      model.specific$estim<-"ML"
    } else {
      model.specific$r<-1
      model.specific$coint<-coint_export
      model.specific$estim<-"OLS"
    }
  }


  z<-list(residuals=res,  
	  coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", 
	  fitted.values=fitted, 
	  model.x=Z, 
	  include=include,
	  lag=lag, 
	  model=YnaX, 
	  df.residual=t-npar/k, 
	  model.specific=model.specific)
  if(model=="VAR"){
    class(z)<-c("VAR","nlVar")
  } else {
    class(z)<-c("VECM","VAR", "nlVar")
    I<-"diff"
  }

  attr(z, "varsLevel")<-I
  attr(z, "model")<-model
  return(z)
}


#### VECM function: wrapper to lineVar
VECM<-function(data, lag,r=1, include = c( "const", "trend","none", "both"), beta=NULL, estim=c("2OLS", "ML"),LRinclude=c("none", "const", "trend","both"))
  lineVar(data, lag, r=r,include = include, model="VECM" ,beta=beta, estim=estim,LRinclude=LRinclude)



###Testing
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
  if(attr(x,"model")=="VECM"){
    cat("\nCointegrating vector ", if(x$model.specific$beta.estimated) paste("(estimated by ", x$model.specific$estim, "):\n", sep="") else "(fixed by user):\n", sep="")
    print(t(x$model.specific$beta))
  }
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
  if(attr(x,"model")=="VECM"&&x$model.specific$LRinclude!="none") stop("toLatex not implemented now for models with arg 'LRinclude' different from 'none'") 
  parenthese<-match.arg(parenthese)
  if(inherits(x,"summary.VAR")){
    a<-myformat(x$coefficients,digits, toLatex=TRUE)
    inp<-switch(parenthese, "StDev"= x$StDev, "Pvalue"= x$Pvalues )
    b<-myformat(inp ,digits,toLatex=TRUE)
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
    r<-x$model.specific$r
    len<-length(res)
    if(r==1){
      res[len+1]<-"+\\begin{smatrix}  %ECT"
      res[len+2]<-TeXVec(coeftoprint[,1]) #see nlVar-methods.R
      res[len+3]<-"\\end{smatrix}ECT_{-1}"
    }else{
      res[len+1]<-"+\\begin{smatrix}  %ECT"
      res[(len+2):(len+x$k+1)]<-TeXMat(coeftoprint[,1:r]) #see nlVar-methods.R
      res[len+x$k+2]<-"\\end{smatrix}ECT_{-1}"
    }
  }
###Const
  a<-if(attr(x,"model")=="VECM") r else 0
  res<-include(x, res, coeftoprint, skip=a)	#nlVar-methods.R
###Lags
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
