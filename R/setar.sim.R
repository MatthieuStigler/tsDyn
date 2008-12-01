setar.sim<-function(data,B, setarObject, sigma,n=200, lag=1, trend=TRUE, nthresh=0, thDelay=0, Thresh, thVar=NULL, mTh=NULL, type=c("boot", "simul", "check"), starting=NULL){

if(!missing(data)&!missing(B))
	stop("You have to provide either B or data, but not both")
p<-lag
m<-lag
type<-match.arg(type)

if(!missing(B)){
	if(type!="simul"){
		type<-"simul"
		warning("Type check or boot are only avalaible with pre-specified data. The type simul was used")}
	if(missing(sigma)){
		warning("sigma is missing, the value taken is 0.25 ")
		sigma<-0.25}
	esp<-p
	if(trend)
		esp<-p+1
	if(esp*(nthresh+1)!=length(B))
		stop("Matrix B bad specified")
	y<-vector("numeric", length=n)
	if(!is.null(starting))
		if(length(starting)==p)
			y[seq_len(p)]<-starting
		else
			stop("Bad specification of starting values. Should have so many values as the number of lags")
}
if(!missing(data))
  y<-data

  str <- nlar.struct(x = y, m = m, d = 1, steps = 1, series = deparse(substitute(x)))
  xx <- str$xx
  yy <- str$yy

if(!missing(setarObject)){
  y<-setarObject$str$x
  xx <- getXX(setarObject$str)
  yy <- getYY(setarObject$str)
  Thresh<-getTh(coef(setarObject))
  nthresh<-setarObject$model.specific$nthresh
  incNames<-setarObject$model.specific$incNames
  thDelay<-setarObject$model.specific$thDelay
  m<-max(c(setarObject$model.specific$ML, setarObject$model.specific$MH))
  p<-m
  if(incNames%in%c("none", "trend"))
    stop("Arg include = none or trend currently not implemented")
  if(incNames=="trend")
    trend<-TRUE
  B<-coef(setarObject)
  res<-na.omit(residuals(setarObject))
  BUp<-getSetarXRegimeCoefs(setarObject, regime="H", returnFirstSlopeCoef=FALSE)
  BDown<-getSetarXRegimeCoefs(setarObject, regime="L", returnFirstSlopeCoef=FALSE)
  if(nthresh==2)
    BMiddle<-getSetarXRegimeCoefs(setarObject, regime="M", returnFirstSlopeCoef=FALSE)
  sigma<-NULL #should be removed after
}

if(type=="simul")
	ndig<-4
else
	ndig<-getndp(y)
if(ndig>.Options$digits){
	ndig<-.Options$digits
	y<-round(y,ndig)
}
npar<-ncol(xx)
if(trend)
	npar<-ncol(xx)+1



###Threshold transition variable
print(missing(thDelay))
    if (!missing(thDelay)) {  
      if (max(thDelay) >= m)             stop(paste("thDelay too high: should be < m (=",  m, ")"))
        z <- xx[, seq_len(max(thDelay)+1)]
    }
    else if (!missing(mTh)) { stop("mTh argument currenly not implemented")
	    if (length(mTh) != m) stop("length of 'mTh' should be equal to 'm'")
        z <- xx %*% mTh
        dim(z) <- NULL
    }
    else if (!missing(thVar)) { stop("thVar argument currenly not implemented")
        if (length(thVar) > nrow(x)) {thVar <- thVar[seq_len(nrow(x))]}
	if (length(thVar) < nrow(x)) {stop("The threshold variable should not be smaller than the serie") }
	z<-nlar.struct(x = thVar, m = m, d = 1, steps = 1, series = deparse(substitute(z)))[, thDelay + 1]
    }
    else {z <- xx[, thDelay]
    }
trans<-z


if(nthresh==1&missing(Thresh))
	Thresh<-mean(trans)
if(nthresh==2&missing(Thresh))
	Thresh<-quantile(trans, probs=c(0.25, 0.75))
if(nthresh!=0)
	Thresh<-round(Thresh, ndig)
##################
###Model
#################
if(!missing(data)){
  if(trend==TRUE)
	  xx <- cbind(1,xx)
  
  if(nthresh==1){
	  isL <- ifelse(round(trans,ndig)<= Thresh,1,0)	### isL: dummy 
	  xx <- cbind(xx * isL,xx * (1 - isL))	### Lower matrix
  }
  
  if(nthresh==2){
	  dummydown <- ifelse(round(trans,ndig)<=Thresh[1], 1, 0)
	  regimedown <- dummydown*xx
	  dummyup <- ifelse(round(trans,ndig)>Thresh[2], 1, 0)
	  regimeup <- dummyup*xx
	  xx <- cbind(regimedown, (1-dummydown-dummyup)*xx, regimeup)
  }
  
  
  model<-lm(yy~xx-1)
  B<-model$coeff
  Sigma<-summary(model)$sigma
  res<-model$residuals
  
  if(missing(sigma))
	  sigma<-Sigma
}

###Verification of stability
# is<-is.InUnitCircle(B,ninc=1, m=m, nthresh=nthresh)
# if(is$warn==TRUE){
# 	warning("The AR coefficients of one regime lie inside the unit circle, thus the serie can be non-stationnary")
# 	cat("\nUnit roots\n")
# 	print(is$root)}
##############################
###Bootstrap
##############################
thDelay<-thDelay+1
#initial values
Yb<-vector("numeric", length=length(y))		#Delta Y term
Yb[1:p]<-y[1:p]


z2<-vector("numeric", length=length(y))
z2[1:p]<-y[1:p]



innov<-switch(type, "boot"=sample(res, replace=TRUE), "simul"=rnorm(length(y)-p,sd=sigma), "check"=res)
resb<-c(rep(0,p),innov)	

if(missing(setarObject)){
  if(nthresh==1){
    BDown<-B[seq_len(npar)]
    BUp<-B[-seq_len(npar)]
  }
  if(nthresh==2){
    BDown <- B[seq_len(npar)]
    BMiddle <- B[seq_len(npar)+npar]
    BUp <- B[seq_len(npar)+2*npar]
  }
}
print(BUp)
print(BDown)
cat("\nnthresh", nthresh, "\n")
print(c(m,p))
print(Thresh)
cat("\nThdelay", thDelay, "\n")
if(nthresh==0){
for(i in (p+1):length(y)){
	Yb[i]<-sum(B[1],B[-1]*Yb[i-c(1:m)],resb[i])
	}
}

else if(nthresh==1){
	for(i in (m+1):length(y)){
		if(round(z2[i-thDelay],ndig)<=Thresh) 
			Yb[i]<-sum(BDown[1],BDown[-1]*Yb[i-c(1:m)],resb[i])
		else 
			Yb[i]<-sum(BUp[1],BUp[-1]*Yb[i-c(1:m)],resb[i])
		z2[i]<-Yb[i]
		}
}

else if(nthresh==2){
	for(i in (m+1):length(y)){
		if(round(z2[i-thDelay],ndig)<=Thresh[1]) 
			Yb[i]<-sum(BDown[1],BDown[-1]*Yb[i-c(1:m)],resb[i])
		else if(round(z2[i-thDelay],ndig)>Thresh[2]) 
			Yb[i]<-sum(BUp[1],BUp[-1]*Yb[i-c(1:m)],resb[i])
		else
			Yb[i]<-sum(BMiddle[1],BMiddle[-1]*Yb[i-c(1:m)],resb[i])
		z2[i]<-Yb[i]
		}
}

if(FALSE){
	while(mean(ifelse(Yb<Thresh, 1,0))>0.05){
		cat("not enough")
		if(!missing(thVar)) 
			Recall(B=B, sigma=sigma,n=n, lag=lag, trend=trend, nthresh=nthresh, thDelay=thDelay,thVar=thVar, type="simul", starting=starting)
		else
			Recall(B=B, sigma=sigma,n=n, lag=lag, trend=trend, nthresh=nthresh, thDelay=thDelay,mTh=mTh, type="simul", starting=starting)
		y<-NULL
		
		}
}


list(B=B, serie=round(Yb,ndig),Sigma=sigma)
}

if(FALSE){
library(tsDyn)
environment(setar.sim)<-environment(star)

##Simulation of a TAR with 1 threshold
sim<-setar.sim(B=c(2.9,-0.4,-0.1,-1.5, 0.2,0.3),lag=2, type="simul", nthresh=1, sigma=1.2,Thresh=2, starting=c(2.8,2.2))$serie
mean(ifelse(sim>2,1,0))	#approximation of values over the threshold

#check the result
selectSETAR(sim, m=2, criterion="SSR")
selectSETAR(sim, m=2, th=list(around=2), ngrid=20)
##Bootstrap a TAR with two threshold (three regimes)
sun<-(sqrt(sunspot.year+1)-1)*2
setar.sim(data=sun,nthresh=2, type="boot", Thresh=c(7,9))$serie

##Check the bootstrap
checkBoot<-setar.sim(data=sun,nthresh=1, type="check", Thresh=6.14)$serie
cbind(checkBoot,sun)

###setar object
environment(setar.sim)<-environment(star)
setarSun<-setar(sun, m=1, nthresh=1)
checkBoot2<-setar.sim(setarObject=setarSun, type="check")$serie
cbind(round(checkBoot2,6),sun)

setarSun<-setar(sun, m=3, nthresh=2)
checkBoot3<-setar.sim(setarObject=setarSun, type="check")$serie
cbind(round(checkBoot3,6),sun)

}


