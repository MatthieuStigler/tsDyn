library(tsDyn)

################################################################
######### From man file:
################################################################
#see that:
a<-matrix(c(-0.2, 0.2), ncol=1)
b<-matrix(c(1,-1), nrow=1)
a%*%b

innov<-rmnorm(100, varcov=diag(2))
vecm1<-TVECM.sim(B=rbind(c(-0.2, 0,0), c(0.2, 0,0)), nthresh=0, beta=1,n=100, lag=1,include="none", innov=innov)
ECT<-vecm1[,1]-vecm1[,2]

#add an intercept as in panel B
b<-TVECM.sim(B=rbind(c(-0.2, 0.1,0,0), c(0.2,0.4, 0,0)), nthresh=0, n=100,beta=1, lag=1,include="const", innov=innov)



##Bootstrap a TVAR with 1 threshold (two regimes)
data(zeroyld)
dat<-zeroyld
TVECMobject<-TVECM(dat, nthresh=1, lag=1, ngridBeta=20, ngridTh=20, plot=FALSE)
b<-TVECM.sim(TVECMobject=TVECMobject,type="boot")

##Check the bootstrap
all(TVECM.sim(TVECMobject=TVECM(dat, nthresh=1, lag=2, ngridBeta=20, ngridTh=20, plot=FALSE, include="none"),type="check")==dat)

################################################################
######### Check error message when matrix badly specified:
################################################################

B<-matrix(rnorm(14), byrow=TRUE,ncol=7)

## 0 thresh
try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=1,show.parMat=TRUE, include="none"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=1,show.parMat=TRUE, include="const"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=1,show.parMat=TRUE, include="both"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=1,show.parMat=TRUE, include="trend"))

try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=2,show.parMat=TRUE, include="none"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=0, n=100, lag=2,show.parMat=TRUE, include="const"))

## 1 thresh  
try(a<-TVECM.sim(B=B, beta=1, nthresh=1, Thresh=0, n=100, lag=1,show.parMat=TRUE, include="none"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=1, Thresh=0, n=100, lag=1,show.parMat=TRUE, include="const"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=1, Thresh=0, n=100, lag=1,show.parMat=TRUE, include="both"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=1, Thresh=0, n=100, lag=1,show.parMat=TRUE, include="trend"))

try(a<-TVECM.sim(B=B, beta=1, nthresh=1, Thresh=0, n=100, lag=2,show.parMat=TRUE, include="const"))

## 2 thresh
try(a<-TVECM.sim(B=B, beta=1, nthresh=2, Thresh=0, n=100, lag=1,show.parMat=TRUE, include="none"))
try(a<-TVECM.sim(B=B, beta=1, nthresh=2, Thresh=0, n=100, lag=2,show.parMat=TRUE, include="const"))

