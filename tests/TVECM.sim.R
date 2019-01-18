library(tsDyn)
library(mnormt)


TVECM.boot.check <- tsDyn:::TVECM.boot.check
options(useFancyQuotes=FALSE) # useful for all.equal comparison

################################################################
######### From man file:
################################################################
#see that:
a<-matrix(c(-0.2, 0.2), ncol=1)
b<-matrix(c(1,-1), nrow=1)
a%*%b

set.seed(123)
innov<-rmnorm(100, varcov=diag(2))
vecm1<-TVECM.sim(B=rbind(c(-0.2, 0,0), c(0.2, 0,0)), nthresh=0, beta=1,n=100, lag=1,include="none", innov=innov)
ECT<-vecm1[,1]-vecm1[,2]
ECT[1:5]

#add an intercept as in panel B
B2 <- rbind(c(-0.2, 0.1,0,0), c(0.2,0.4, 0,0))
b<- TVECM.sim(B=B2, nthresh=0, n=100,beta=1, lag=1,include="const", innov=innov, show.parMat=TRUE)
b_2 <- VECM.sim(B=B2,  n=100,beta=1, lag=1,include="const", innov=innov, show.parMat=TRUE)
b[1:5,]
all.equal(b, b_2)

## other ways to input beta:
b_beta_vec <- TVECM.sim(B=B2, nthresh=0, n=100,beta=c(1,-1), lag=1,include="const", innov=innov, show.parMat=TRUE)
b_beta_vec_2 <- VECM.sim(B=B2, n=100,beta=c(1,-1), lag=1,include="const", innov=innov, show.parMat=TRUE)
b_beta_vec[1:5,]
all.equal(b,b_beta_vec)
all.equal(b_beta_vec, b_beta_vec_2)

beta_mat <- matrix(c(1,-1), ncol=1)
b_beta_mat <- TVECM.sim(B=B2, nthresh=0, n=100,beta=beta_mat, lag=1,include="const", innov=innov, show.parMat=TRUE)
b_beta_mat[1:5,]
all.equal(b,b_beta_mat)


########################
######## TVECM
########################

##Bootstrap a TVECM with 1 threshold (two regimes)
data(zeroyld)
dat<-zeroyld
TVECMobject<-TVECM(dat, nthresh=1, lag=1, ngridBeta=20, th1=list(exact=-1.195), plot=FALSE)
tv_1_boot <-TVECM.boot(TVECMobject, seed=123, show.parMat=TRUE)
head(tv_1_boot)

##Check the bootstrap
TVECM.boot.check(TVECMobject)

## check correspondance bootstrap/simul:
tv_1_sim <-TVECM.sim(B=tsDyn:::coefMat.nlVar(TVECMobject),beta=TVECMobject$model.specific$beta,
                        Thresh=getTh(TVECMobject), show.parMat=TRUE, innov=matrix(0,200,2))
head(tv_1_boot)

tv_1_sim <-TVECM.sim(B=tsDyn:::coefMat.nlVar(TVECMobject), 
                     beta=TVECMobject$model.specific$beta,
                     Thresh=getTh(TVECMobject), show.parMat=TRUE)
head(tv_1_boot)

##Bootstrap a TVAR with two threshold (three regimes)
tv_2_const <- TVECM(dat, lag=1, nthresh=2, plot=FALSE, trace=FALSE, th1=list(exact=-1.312),
                   th2=list(exact=0.774))
tv_2_none <- TVECM(dat, lag=1, nthresh=2, plot=FALSE, trace=FALSE, th1=list(exact=-1.312),
                th2=list(exact=0.774), include="none")
tv_2_trend <- TVECM(dat, lag=1, nthresh=2, plot=FALSE, trace=FALSE, th1=list(exact=-1.312),
                   th2=list(exact=0.774), include="trend")
tv_2_both <- TVECM(dat, lag=1, nthresh=2, plot=FALSE, trace=FALSE, th1=list(exact=-1.312),
                    th2=list(exact=0.774), ngridBeta=5, include="both")

tv_2_const_common <- TVECM(dat, lag=1, nthresh=2, plot=FALSE, trace=FALSE, th1=list(exact=-1.451),
                    th2=list(exact=0.754), include="const", common="only_ECT")
tv_2_const_l2 <- TVECM(dat, nthresh=2, lag=2, ngridBeta=5,  plot=FALSE, include="none",
                       th1=list(exact=-1.312),th2=list(exact=0.774), trace=FALSE)

TVECM.boot(tv_2_const, show.parMat=TRUE, seed=456)[1:5,]
TVECM.boot(tv_2_none, show.parMat=TRUE, seed=456)[1:5,]
TVECM.boot(tv_2_trend,show.parMat=TRUE, seed=456)[1:5,]
TVECM.boot(tv_2_both, show.parMat=TRUE, seed=456)[1:5,]
try(TVECM.boot(tv_2_const_common, show.parMat=TRUE, seed=456)[1:5,], silent=TRUE)

TVECM.boot.check(tv_2_const)
TVECM.boot.check(tv_2_none)
TVECM.boot.check(tv_2_const_l2)
TVECM.boot.check(tv_2_both)

## does not work:
TVECM.boot.check(tv_2_trend)

###############
#### p>2
###############

data(barry)

ve_r1_l1 <- VECM(barry, lag=1)
ve_r1_l3 <- VECM(barry, lag=3)
ve_r2_l3 <- VECM(barry, lag=3, estim="ML", r=2)

TVECM.boot.check(ve_r1_l1)
TVECM.boot.check(ve_r1_l3)
TVECM.boot.check(ve_r2_l3)

VECM.boot(ve_r1_l1, show.parMat=TRUE, seed=234)[1:5,]
VECM.boot(ve_r1_l3, show.parMat=TRUE, seed=234)[1:5,]
VECM.boot(ve_r2_l3, show.parMat=TRUE, seed=234)[1:5,]

################################################################
######### Check error message when matrix badly specified:
################################################################

B <- matrix(rnorm(14), byrow=TRUE,ncol=7)

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

