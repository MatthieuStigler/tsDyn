TVAR.sim<-function(data,B,TVARobject, Thresh, nthresh=1, type=c("simul","boot", "check"), n=200, lag=1, include = c("const", "trend","none", "both"),  thDelay=1,  thVar=NULL, mTh=1, starting=NULL, innov=rmnorm(n, mean=0, varcov=varcov), varcov=diag(1,k), show.parMat=FALSE, round=FALSE){
if(!missing(data)&!missing(B))
	stop("You have to provide either B or y, but not both")
p<-lag
type<-match.arg(type)
include<-match.arg(include)

###check correct arguments
if(!nthresh%in%c(0,1,2))
  stop("Arg nthresh should  be either 0, 1 or 2")
if(!missing(n)&any(!missing(data), !missing(TVARobject)))
  stop("arg n should not be given with arg data or TVARobject")
if(max(thDelay)>p)
	stop("Max of thDelay should be smaller or equal to the number of lags")
if(!missing(TVARobject)&any(!missing(Thresh), !missing(nthresh), !missing(lag), !missing(thDelay), !missing(mTh)))
  warning("When object TVARobject is given, only args 'type' and 'round' are relevant, others are not considered")
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
  if(type!="simul"){
    type<-"simul"
    warning("Type check or boot are only avalaible with pre specified data. The type simul was used")}
  nB<-nrow(B)
  ndig<-4
  esp<-p*nB+ninc 
  if(esp*(nthresh+1)!=ncol(B))
    stop("Matrix B bad specified")
  y<-matrix(0,ncol=nB, nrow=n)
  if(!is.null(starting)){
    if(all(dim(as.matrix(starting))==c(nB,p)))
      y[seq_len(p),]<-starting
    else
      stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables")
  }
  Bmat<-B
  temp<-TVAR_thresh(mTh=mTh,thDelay=thDelay,thVar=thVar,y=y, p=p) #Stored in misc.R
  trans<-temp$trans
  combin<-temp$combin
  
  if(nthresh==1){
    if(missing(Thresh))
      Thresh<-mean(trans)
    if(length(Thresh)!=1){
      warning("Please only one Thresh value if you choose nthresh=1. First one was chosen")
      Thresh<-Thresh[1]}
  }
  if(nthresh==2){
    if(missing(Thresh))
      Thresh<-quantile(trans, probs=c(0.25, 0.75))
    if(length(Thresh)!=2)
      stop("please give two Thresh values if you choose nthresh=2")
  }
  if(nthresh%in%c(1,2))
    Thresh<-round(Thresh,ndig)
  k <- ncol(y) 		#Number of variables
  T <- nrow(y) 		#Size of start sample
}

### possibility 2: only data is given: compute it with linear or selectSETAR
else if(!missing(data)){
  if(nthresh==0){
    TVARobject<-lineVar(data, lag=p, include=include, model="VAR")
  }
  else{ 
    if(!missing(Thresh))
      TVARobject<-TVAR(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE, gamma=Thresh)
    else
      TVARobject<-TVAR(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE)
  }
}
### possibility 3: setarobject is given by user (or by poss 2)
if(!missing(TVARobject)){
  k<-TVARobject$k
  T<-TVARobject$T
  p<-TVARobject$lag
  modSpe<-TVARobject$model.specific
  res<-residuals(TVARobject)
  Bmat<-coefMat(TVARobject)
  y<-as.matrix(TVARobject$model)[,1:k]
  ndig<-getndp(y[,1])
  if(nthresh>0){
    if(modSpe$oneMatrix)
      stop("arg commoninter in TVAR currently not implemented in TVAR.sim")
    if(attr(TVARobject, "levelTransVar")=="MTAR")
      stop("arg model=MTAR in TVAR currently not implemented in TVAR.sim")
    if(attr(TVARobject, "transVar")=="external")
      stop("arg thVaR in TVAR currently not implemented in TVAR.sim")
    Thresh<-modSpe$Thresh
    thDelay<-modSpe$thDelay
    combin<-modSpe$transCombin
    nthresh<-modSpe$nthresh
  }
}

t <- T-p 		#Size of end sample
npar<-k*(p+ninc)

  ##### put coefficients vector in right form according to arg include (arg both need no modif)
  a<-NULL
  if(include=="none")
    for(i in 0:nthresh) a<-c(a, i*(p*k)+c(1,2))
  else if(include=="const")
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(2))
  else if(include=="trend")
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(1))
    #if (include=="both"): correction useless
  Bmat<-myInsertCol(Bmat, c=a ,0)
  nparBmat<-p*k+2
  
  
##############################
###Reconstitution boot/simul
##############################

#initial values
Yb<-matrix(0, nrow=T, ncol=k)
Yb[1:p,]<-y[1:p,]		

if(nthresh>0){
  z2<-vector("numeric", length=nrow(y))
  z2[1:p]<-y[1:p,]%*%combin
  }

trend<-1:T

##resampling
if(type=="simul"&&dim(innov)!=c(n,k))
  stop(paste("input innov is not of right dim, should be matrix with", n,"rows and ", k, "cols\n"))
innov<-switch(type, "boot"=res[sample(seq_len(t), replace=TRUE),], "simul"=innov, "check"=res)
resb<-rbind(matrix(0,nrow=p, ncol=k),innov)	

if(nthresh==0){
	for(i in (p+1):T){
		Yb[i,]<-rowSums(cbind(Bmat[,1], Bmat[,2]*trend[i], Bmat[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
	}
}

else if(nthresh==1){
	BDown<-Bmat[,seq_len(nparBmat)]
	BUp<-Bmat[,-seq_len(nparBmat)]

	for(i in (p+1):(nrow(y))){
		if(round(z2[i-thDelay],ndig)<=Thresh) {
			Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,2]*trend[i], BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
			if(round)
			  Yb[i,]<-round(Yb[i,], ndig)
		}
		else{
			Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
			if(round)
			  Yb[i,]<-round(Yb[i,], ndig)
		 }
		z2[i]<-Yb[i,]%*%combin
	}
}


else if(nthresh==2){
	BDown <- Bmat[,seq_len(nparBmat)]
	BMiddle <- Bmat[,seq_len(nparBmat)+nparBmat]
	BUp <- Bmat[,seq_len(nparBmat)+2*nparBmat]
	for(i in (p+1):(nrow(y))){
		if(round(z2[i-thDelay],ndig)<=Thresh[1]) 
			Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,2]*trend[i], BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		else if(round(z2[i-thDelay],ndig)>Thresh[2]) 
			Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		else
			Yb[i,]<-rowSums(cbind(BMiddle[,1],BMiddle[,2]*trend[i], BMiddle[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		z2[i]<-Yb[i,]%*%combin
	}
}


if(FALSE){
	new<-TVAR_thresh(y=Yb,mTh=mTh,thDelay=thDelay,thVar=thVar, p=p)$trans
	if(mean(ifelse(new[,thDelay]<Thresh, 1,0))>0.05)
		list(B=Bmat,Yb=Yb)
	else
		Recall(B=Bmat, n=n, lag=lag, nthresh=nthresh, thDelay=thDelay, Thresh, thVar=thVar, mTh=mTh, type=type, starting=starting)
}
# print(cbind(y, round(Yb,3)))

if(show.parMat)
  print(Bmat)
res<-round(Yb, ndig)
return(res)
}

if(FALSE){
library(tsDyn)
environment(TVAR.sim)<-environment(star)

##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
sim<-TVAR.sim(B=B,nthresh=1,n=500, type="simul",mTh=1, Thresh=5, starting=c(5.2, 5.5))$serie

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)


##Bootstrap a TVAR 
data(zeroyld)
serie<-zeroyld

TVAR.sim(data=serie,nthresh=0, type="sim")
all(TVAR.sim(data=serie,nthresh=0, type="check", lag=1)$serie==serie)

##with two threshold (three regimes)
TVAR.sim(data=serie,nthresh=2,type="boot",mTh=1, Thresh=c(7,9))

##Check the bootstrap: ok!
environment(TVAR.sim)<-environment(star)
all(TVAR.sim(data=serie,nthresh=0, type="check",mTh=1)$serie==serie) #TRUE
all(TVAR.sim(data=serie,nthresh=1, type="check",mTh=1)$serie==serie)#TRUE
all(TVAR.sim(data=serie,nthresh=2, type="check",mTh=1)$serie==serie) #TRUE

all(TVAR.sim(data=serie,nthresh=2, type="check",mTh=2)$serie==serie) #TRUE
all(TVAR.sim(data=serie,nthresh=0,lag=3, type="check",mTh=2)$serie==serie) #TRUE
all(TVAR.sim(data=serie,nthresh=1,lag=2, type="check",mTh=2)$serie==serie) #TRUE


###with TVARobject
all(TVAR.sim(TVARobject=TVAR(serie, nthresh=2, lag=1),type="check")$serie==serie) #TRUE
all(TVAR.sim(TVARobject=TVAR(serie, nthresh=1, lag=1),type="check")$serie==serie) #TRUE
all(TVAR.sim(TVARobject=TVAR(serie, nthresh=1, lag=2),type="check")$serie==serie) #TRUE
all(TVAR.sim(TVARobject=TVAR(serie, nthresh=1, lag=2),type="check")$serie==serie) #TRUE

all(TVAR.sim(TVARobject=lineVar(serie, lag=1),type="check")$serie==serie) #TRUE

##Check the bootstrap: no! prob with trend... both.. none...
TVAR.sim(data=serie,nthresh=1, type="check",mTh=1, include="trend", round=TRUE)$serie==serie
TVAR.sim(data=serie,nthresh=1, type="check",mTh=1, include="trend")$serie==serie
all(TVAR.sim(data=serie,nthresh=2, type="check",mTh=1, include="both")$serie==serie)
TVAR.sim(data=serie,nthresh=2, type="check",mTh=1, include="none", round=TRUE)$serie==serie
}
