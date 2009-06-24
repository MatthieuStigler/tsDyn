TVECM.sim<-function(data,B,TVECMobject, nthresh=1, Thresh, beta, n=200, lag=1, type=c("simul","boot", "check"),  include = c("const", "trend","none", "both"), starting=NULL, sigma){


if(!missing(data)&!missing(B))
	stop("You have to provide either B or y, but not both")
p<-lag
type<-match.arg(type)
include<-match.arg(include)

###check correct arguments
if(!nthresh%in%c(0,1,2))
  stop("Arg nthresh should  be either 0, 1 or 2")
if(!missing(n)&any(!missing(data), !missing(TVECMobject)))
  stop("arg n should not be given with arg data or TVECMobject")
if(!missing(TVECMobject)&any(!missing(Thresh), !missing(nthresh), !missing(lag)))
  warning("When object TVECMobject is given, only args 'type' and 'round' are relevant, others are not considered")
##include term
  if(include=="none")
    ninc<-0
  else if(include=="const")
    ninc<-1
  else if(include=="trend")
    ninc<-1
  else if(include=="both")
    ninc<-2

### possibility 1: only parameters matrix is given
if(!missing(B)){
  if(missing(beta))
    stop("please provide arg beta (cointegrating value)")
  if(type!="simul"){
    type<-"simul"
    warning("Type check or boot are only avalaible with pre specified data. The type simul was used")
    }
  if(missing(sigma)){
	warning("sigma is missing, the values taken are 0.25 for each variable")
	sigma<-rep(0.25,nrow(B))
	}
  nB<-nrow(B)
  ndig<-4
  esp<-p*nB+1+ninc #number of lags +ecm
  if(esp*(nthresh+1)!=ncol(B))
    stop(paste("Matrix B bad specified: should have ", esp*(nthresh+1), "columns, but has", ncol(B), "\n"))
  y<-matrix(0,ncol=nB, nrow=n)
  if(!is.null(starting)){
    if(all(dim(as.matrix(starting))==c(nB,p)))
      y[seq_len(p),]<-starting
    else
	stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables")
  }
    Bmat<-B
    k <- ncol(y) 		#Number of variables
    T <- nrow(y) 		#Size of start sample
}

### possibility 2: only data is given: compute it with linear or selectSETAR
else if(!missing(data)){
  if(nthresh==0){
    TVECMobject<-linear2(data, lag=p, include=include, model="VECM")
  }
  else{ 
    if(!missing(Thresh))
      TVECMobject<-TVECM(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE, gamma1=Thresh)
    else
      TVECMobject<-TVECM(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE)
  }
}
### possibility 3: setarobject is given by user (or by poss 2)
if(!missing(TVECMobject)){
  k<-TVECMobject$k
  T<-TVECMobject$T
  p<-TVECMobject$lag
  include<-TVECMobject$include
  if(include %in% c("trend", "both"))
    warning(paste("Accuracy of function (tested with arg type=check) is not good when arg include=",include," is given\n"))
  modSpe<-TVECMobject$model.specific
  beta<-modSpe$beta
  res<-residuals(TVECMobject)
  Bmat<-coefMat(TVECMobject)
  Sigma<- matrix(1/T*crossprod(res),ncol=k)
  y<-as.matrix(TVECMobject$model)[,1:k]
  ndig<-getndp(y[,1])
  if(missing(sigma)&type=="simul")
      sigma<-diag(Sigma)
  if(nthresh>0){
    Thresh<-modSpe$Thresh
    nthresh<-modSpe$nthresh
  }
}

t <- T-p 		#Size of end sample
npar<-k*(p+ninc+1)

  ##### put coefficients vector in right form according to arg include (arg both need no modif)
  a<-NULL
  if(include=="none")
    for(i in 0:nthresh) a<-c(a, i*(p*k+3)+c(2,3))
  else if(include=="const")
    for(i in 0:nthresh) a<-c(a, i*(p*k+3)+c(3))
  else if(include=="trend")
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(2))
    #if (include=="both"): correction useless
  Bmat<-myInsertCol(Bmat, c=a ,0)
  nparBmat<-p*k+2+1
  
##############################
###Reconstitution boot/simul
##############################
#initial values

#initial values
Yb<-matrix(0, nrow=nrow(y), ncol=k)
Yb[1:(p+1),]<-y[1:(p+1),]		

trend<-c(rep(NA, T-t),1:t)
BETA<-matrix(c(1, -beta), nrow=1)


#resampling/ simulation of residual/innovations
innov<-switch(type, "boot"=res[sample(seq_len(t), replace=TRUE),], "simul"=t(sqrt(sigma)*t(matrix(rnorm(k*t), ncol=k))), "check"=res)
resb<-rbind(matrix(0,nrow=p+1, ncol=k),innov)	


if(nthresh==0){
  for(i in (p+2):T){
    ECT<-Bmat[,1]%*%BETA%*%matrix(Yb[i-1,], ncol=1)
    Yb[i,]<-rowSums(cbind(Yb[i-1,],Bmat[,2], Bmat[,3]*trend[i], ECT,Bmat[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
  }
}

else if(nthresh==1){
  BD<-Bmat[,seq_len(nparBmat)]
  BU<-Bmat[,-seq_len(nparBmat)]
  
  for(i in (p+2):(nrow(y))){
   ECT<-BETA%*%matrix(Yb[i-1,], ncol=1)
   if(round(ECT,ndig)<=Thresh){
      Yb[i,]<-rowSums(cbind(Yb[i-1,],BD[,1]%*%ECT, BD[,2], BD[,3]*trend[i],BD[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
      }
    else{
      Yb[i,]<-rowSums(cbind(Yb[i-1,],BU[,1]%*%ECT, BU[,2], BU[,3]*trend[i],BU[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
      }
  }
}


else if(nthresh==2){
  BD <- Bmat[,seq_len(nparBmat)]
  BM <- Bmat[,seq_len(nparBmat)+nparBmat]
  BU <- Bmat[,seq_len(nparBmat)+2*nparBmat]
  for(i in (p+2):(nrow(y))){
  ECT<-BETA%*%matrix(Yb[i-1,], ncol=1)
    if(round(ECT,ndig)<=Thresh[1]){ 
      Yb[i,]<-rowSums(cbind(Yb[i-1,],BD[,1]%*%ECT,BD[,2], BD[,3]*trend[i], BD[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
      }
    else if(round(ECT,ndig)>Thresh[2]) {
      Yb[i,]<-rowSums(cbind(Yb[i-1,],BU[,1]%*%ECT,BU[,2], BU[,3]*trend[i],BU[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
      }
    else{
      Yb[i,]<-rowSums(cbind(Yb[i-1,],BM[,1]%*%ECT,BM[,2], BM[,3]*trend[i],BM[,-c(1,2,3)]%*%matrix(t(Yb[i-c(1:p),]-Yb[i-c(2:(p+1)),]), ncol=1),resb[i,]))
      }
  }
}


list(B=Bmat,serie=round(Yb,ndig))
}

if(FALSE){
library(tsDyn)
environment(TVECM.sim)<-environment(star)

##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.2, 0.11928245, 1.00880447, -0.009974585, 0.3, -0.089316, 0.95425564, 0.02592617),c( -0.1, 0.25283578, 0.09182279,  0.914763741, 0.35,-0.0530613, 0.02248586, 0.94309347))
sim<-TVECM.sim(B=B,beta=1, nthresh=1,n=500, type="simul",Thresh=5, starting=c(5.2, 5.5), sigma=c(0.3,0.4))$serie

#estimate the new serie
TVECM(sim, lag=1)

##Bootstrap a TVAR with two threshold (three regimes)
data(zeroyld)
dat<-zeroyld
TVECM.sim(data=dat,nthresh=2, type="boot", Thresh=c(7,9))

##Check the bootstrap
linObject<-linear2(dat, lag=1, model="VECM")
all(TVECM.sim(TVECMobject=linObject,type="check")$serie==dat)
all(TVECM.sim(TVECMobject=linear2(dat, lag=1, model="VECM", include="none"),type="check")$serie==dat)

#not working: (probably trend coefficients too small so digits errors)
all(TVECM.sim(TVECMobject=linear2(dat, lag=1, model="VECM", include="trend"),type="check")$serie==dat)
all(TVECM.sim(TVECMobject=linear2(dat, lag=1, model="VECM", include="both"),type="check")$serie==dat)

#nthresh=1
TVECMobject<-TVECM(dat, nthresh=1, lag=1, bn=20, ngridG=20, plot=FALSE)
all(TVECM.sim(TVECMobject=TVECMobject,type="check")$serie==dat)

all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=1, lag=2, bn=20, ngridG=20, plot=FALSE),type="check")$serie==dat)
all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=1, lag=1, bn=20, ngridG=20, plot=FALSE, include="none"),type="check")$serie==dat)
all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=1, lag=2, bn=20, ngridG=20, plot=FALSE, include="none"),type="check")$serie==dat)

#nthresh=2
TVECMobject2<-TVECM(dat, nthresh=2, lag=1, bn=20, ngridG=20, plot=FALSE)
all(TVECM.sim(TVECMobject=TVECMobject2,type="check")$serie==dat)
all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=2, lag=2, bn=20, ngridG=20, plot=FALSE),type="check")$serie==dat)

all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=2, lag=1, bn=20, ngridG=20, plot=FALSE, include="none"),type="check")$serie==dat) 
#famous rounding problem...
all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=2, lag=2, bn=20, ngridG=20, plot=FALSE, include="none"),type="check")$serie==dat)

###TODO:
#improve trend/both case
#TVECM: possibility to give args!
TVECM(dat, nthresh=1, lag=2, bn=20, ngridG=20, plot=FALSE, gamma1=list(exact=-1.4),include="none")
TVECM(dat, nthresh=1, lag=2, bn=20, ngridG=20, plot=FALSE, gamma1=list(exact=-1.4),beta=list(exact=1),include="none")
TVECM(dat, nthresh=2, lag=2, bn=20, ngridG=20, plot=FALSE, gamma1=list(exact=-1.4),gamma2=list(exact=0.5),include="none")
}
