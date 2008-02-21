TVAR_simul<-function(data,B,Thresh, nthresh=1, type=c("boot", "simul", "check"), sigma,n=200, lag=1, trend=TRUE,  thDelay=1,  thVar=NULL, mTh=1, starting=NULL){
if(!missing(data)&!missing(B))
	stop("You have to provide either B or y, but not both")
p<-lag

if(!missing(B)){
	if(type!="simul"){
		type<-"simul"
		warning("Type check or boot are only avalaible with pre specified data. The type simul was used")}
	if(missing(sigma)){
		warning("sigma is missing, the values taken are 0.25 for each variable")
		sigma<-rep(0.25,k)}
	nB<-nrow(B)
	esp<-p*nB
	if(trend)
		esp<-p*nB+1
	if(esp*(nthresh+1)!=ncol(B))
		stop("Matrix B bad specified")
	y<-matrix(0,ncol=nB, nrow=n)
	if(!is.null(starting))
		if(all(dim(as.matrix(starting))==c(nB,p)))
			y[seq_len(p),]<-starting
		else
			stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables")
}
else
	y<-as.matrix(data)
if(!missing(B))
	ndig<-4
else
	ndig<-getndp(y[,1])
T <- nrow(y) 		#Size of start sample
t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables

Y <- y[(p+1):T,]
Z <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix
if(trend==TRUE)
	Z <- cbind(1,Z)
npar<-ncol(Z)
if(max(thDelay)>p)
	stop("Max of thDelay should be smaller or equal to the number of lags")

type<-match.arg(type)

########################
### Threshold variable
########################

temp<-TVAR_thresh(mTh=mTh,thDelay=thDelay,thVar=thVar,y=y, p=p) #Stored in misc.R
trans<-temp$trans
combin<-temp$combin

if(nthresh==1&missing(Thresh))
	Thresh<-mean(trans)
if(nthresh==2&missing(Thresh))
	Thresh<-quantile(trans, probs=c(0.25, 0.75))
Thresh<-round(Thresh,ndig)
##################
###Model
#################
if(!missing(data)){

if(nthresh==1){
	dummydown <- ifelse(trans[,thDelay]<=Thresh, 1, 0)
	ndown <- mean(dummydown)
	regimeDown <- dummydown*Z
	regimeUp<-(1-dummydown)*Z
	Z <- cbind(regimeDown,regimeUp)		# dim k(p+1) x t
}
if(nthresh==2){
	dummydown <- ifelse(trans[,thDelay]<=Thresh[1], 1,0)
	ndown <- mean(dummydown)
	regimedown <- dummydown*Z
	dummyup <- ifelse(trans[,thDelay]>=Thresh[2], 1,0)
	nup <- mean(dummyup)
	regimeup <- dummyup*Z
	Z <- cbind(regimedown,(1-dummydown-dummyup)*Z, regimeup)	# dim k(p+1) x t
}

B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)
res<-Y-Z%*%t(B)
Sigma<- matrix(1/T*crossprod(res),ncol=k,dimnames=list(colnames(data), colnames(data)))
print(Sigma)
if(missing(sigma)&type=="simul")
	sigma<-diag(Sigma)
}

##############################
###Bootstrap for the F test
##############################
#initial values

Yb<-matrix(0, nrow=nrow(y), ncol=k)
Yb[1:p,]<-y[1:p,]		

z2<-vector("numeric", length=nrow(y))
z2[1:p]<-y[1:p,]%*%combin			



innov<-switch(type, "boot"=res[sample(seq_len(t), replace=TRUE),], "simul"=t(sqrt(sigma)*t(matrix(rnorm(k*t), ncol=k))), "check"=res)
resb<-rbind(matrix(0,nrow=p, ncol=k),innov)	


if(nthresh==0){
	for(i in (p+1):T){
		Yb[i,]<-rowSums(cbind(B[,1], B[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
	}
}

else if(nthresh==1){
	BDown<-B[,seq_len(npar)]
	BUp<-B[,-seq_len(npar)]

	for(i in (p+1):(nrow(y))){
		if(round(z2[i-thDelay],ndig)<=Thresh) 
			Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		else
			Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		z2[i]<-Yb[i,]%*%combin
		}
}


else if(nthresh==2){
	BDown <- B[,seq_len(npar)]
	BMiddle <- B[,seq_len(npar)+npar]
	BUp <- B[,seq_len(npar)+2*npar]
	for(i in (p+1):(nrow(y))){
		if(round(z2[i-thDelay],ndig)<=Thresh[1]) 
			Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		else if(round(z2[i-thDelay],ndig)>=Thresh[2]) 
			Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		else
			Yb[i,]<-rowSums(cbind(BMiddle[,1], BMiddle[,-1]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
		z2[i]<-Yb[i,]%*%combin
	}
}

if(FALSE){
	new<-TVAR_thresh(y=Yb,mTh=mTh,thDelay=thDelay,thVar=thVar, p=p)$trans
	if(mean(ifelse(new[,thDelay]<Thresh, 1,0))>0.05)
		list(B=B,Yb=Yb)
	else
		Recall(B=B, sigma=sigma,n=n, lag=lag, trend=trend, nthresh=nthresh, thDelay=thDelay, Thresh, thVar=thVar, mTh=mTh, type=type, starting=starting)
}
# print(cbind(y, round(Yb,3)))

list(B=B,serie=round(Yb,ndig))
}

if(FALSE){
library(tsDyn)
environment(TVAR_simul)<-environment(star)

##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
sim<-TVAR_simul(B=B,nthresh=1,n=500, type="simul",mTh=1, Thresh=5, starting=c(5.2, 5.5), sigma=c(0.3,0.4))$serie

#estimate the new serie
OlsTVAR(sim, lag=1, dummyToBothRegimes=TRUE)

##Bootstrap a TVAR with two threshold (three regimes)
data(zeroyld)
serie<-zeroyld
TVAR_simul(data=serie,nthresh=2,n=500, type="boot",mTh=1, Thresh=c(7,9))

##Check the bootstrap
TVAR_simul(data=serie,nthresh=2,n=500, type="check",mTh=1, Thresh=c(7,9))$serie==serie

}