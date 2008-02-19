TVAR_LRtest <- function (data, m=1, d = 1, steps = d, trend=TRUE, series, thDelay = 1:2, mTh=1, thVar, nboot=10, plot=FALSE, trim=0.1, test=c("1vs", "2vs3")) {
    if (missing(series))  series <- deparse(substitute(x))
y <- as.matrix(data) 
Torigin <- nrow(y) 	#Size of original sample
#if(type=="difference") {y<-diff(y)}
T <- nrow(y) 		#Size of start sample
t <- T-m 		#Size of end sample
k <- ncol(y) 		#Number of variables
p<-m
if(is.null(colnames(data)))
	colnames(data)<-paste("Var", seq_len(k), sep="")
if(max(thDelay)>m)
	stop("Max of thDelay should be smaller or equal to the number of lags")
ndig<-getndp(y)
Y <- y[(m+1):T,] #
Z <- embed(y, m+1)[, -seq_len(k)]	#Lags matrix
a<-0
if(trend==TRUE){
	Z <- cbind(1,Z)
	a<-1}
else 
	warning("The test was currently implemented for a model with trend. Results could be altered without trend")
npar <- ncol(Z)		

cat("Warning: the thDelay values do not correspond to the univariate implementation in tsdyn\n")

##################
###Linear model
#################
B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
res<-Y-Z%*%t(B)
Sigma<- matrix(1/T*crossprod(res),ncol=k)


Y<-t(Y)
########################
### Threshold variable
########################

###External threshold variable
if (!missing(thVar)) {		
        if (length(thVar) > Torigin) {
		z <- thVar[seq_len(Torigin)]
		warning("The external threshold variable is not of same length as the original variable")
        }
        else
		z <- thVar
	z<-embed(z,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
} 

###Combination (or single value indicating position) of contemporaneous variables

else {
	if (length(mTh) > k)
		stop("length of 'mTh' should be equal to the number of variables, or just one")
	if(length(mTh)==1) {
		if(mTh>p)
			stop("mTh too big, should be smaller or equal to the number of variables")
		combin <- matrix(0,ncol=1, nrow=k)
		combin[mTh,]<-1
	}
	else 
		combin<-matrix(mTh,ncol=1, nrow=k)
	zcombin <- y %*% combin
	z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
}




###############################
###Grid for transition variable
###############################

if(length(thDelay)==1) 
	b<-1
else b<-thDelay
	
allgammas<-sort(unique(z[,b]))
ng<-length(allgammas)
ninter<-round(trim*ng)
gammas<-allgammas[ceiling(trim*ng+1):floor((1-trim)*ng-1)]


###################
###Search function
##################

SSR_1thresh<- function(grid,Z,Y, trans){
	d<-grid[1]
	gam1<-grid[2]
	##Threshold dummies
	d1<-ifelse(trans[,d]<gam1, 1,0)
	regimeDown <- d1 * Z
	regimeUp<-(1-d1)*Z
	##SSR
	Z1 <- t(cbind(regimeDown, regimeUp))		# dim k(p+1) x t
	crossprod(c(Y - tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))%*% Z1))
} #end of the function

Sigma_1thresh<- function(gam1, d,Z,Y, trans){
	##Threshold dummies
	d1<-ifelse(trans[,d]<gam1, 1,0)
	regimeDown <- d1 * Z
	regimeUp<-(1-d1)*Z
	##SSR
	Z1 <- t(cbind(regimeDown, regimeUp))		# dim k(p+1) x t
	B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
	matrix(1/T*tcrossprod(Y - B1 %*% Z1),ncol=k)
} #end of the function

SSR_2thresh <- function(gam1,gam2,d,Z,Y,trans){
	##Threshold dummies
	dummydown <- ifelse(trans[,d]<gam1, 1, 0)
	regimedown <- dummydown*Z
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,d]>=gam2, 1, 0)
	regimeup <- dummyup*Z
	nup <- mean(dummyup)
	##SSR from TVAR(3)
	#print(c(ndown,1-nup-ndown,nup))
	if(min(nup, ndown, 1-nup-ndown)>trim){
		Z2 <- t(cbind(regimedown, (1-dummydown-dummyup)*Z, regimeup))		# dim k(p+1) x t	
		resid <- crossprod(c( Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2))	#SSR
	}
	else
		resid <- NA
	return(resid)
}

Sigma_2thresh <- function(gam1,gam2,d,Z,Y,trans){
	##Threshold dummies
	dummydown <- ifelse(trans[,d]<gam1, 1, 0)
	regimedown <- dummydown*Z
	dummyup <- ifelse(trans[,d]>=gam2, 1, 0)
	regimeup <- dummyup*Z
	##SSR from TVAR(3)
	#print(c(ndown,1-nup-ndown,nup))
	Z2 <- t(cbind(regimedown, (1-dummydown-dummyup)*Z, regimeup))		# dim k(p+1) x t	
	matrix(1/T*tcrossprod(Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2),ncol=k)
}

##################
###One threshold
##################

IDS<-as.matrix(expand.grid(thDelay, gammas)) 
result <- apply(IDS, 1, SSR_1thresh,Z=Z, Y=Y,trans=z)
bestDelay<-IDS[which.min(result),1]
bestThresh<-IDS[which.min(result),2]
cat("Best unique threshold", bestThresh, "\t\t\t\t SSR", min(result), "\n")

Sigma_mod1thresh<-Sigma_1thresh(gam1=bestThresh, d=bestDelay,Z=Z,Y=Y, trans=z)
##################
###Two thresholds
##################

###Function for conditional search
condiStep<-function(allgammas, threshRef, delayRef,ninter, fun,MoreArgs=NULL, target=NULL){

wh.thresh <- which.min(abs(allgammas-threshRef))
Thr2<-which.min(abs(allgammas-target))
#search for a second threshold smaller than the first

if(wh.thresh>2*ninter){
	gammaMinus<-allgammas[seq(from=max(ninter, Thr2-20), to=min(wh.thresh-ninter, Thr2+20))]
	storeMinus <- mapply(fun, gam1=gammaMinus,gam2=threshRef,MoreArgs=MoreArgs)	
}
else
	storeMinus <- NA

#search for a second threshold higher than the first
if(length(wh.thresh<length(allgammas)-2*ninter)){
	gammaPlus<-allgammas[seq(from=max(wh.thresh+ninter,Thr2-20), to=min(length(allgammas)-ninter, Thr2+20))]
	storePlus <- mapply(fun,gam1=threshRef,gam2=gammaPlus,  MoreArgs=MoreArgs)
}
else
	storePlus <- NA

#results
store2 <- c(storeMinus, storePlus)

positionSecond <- which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)[1]
if(positionSecond<=length(storeMinus))
	newThresh<-gammaMinus[positionSecond]
else
	newThresh<-gammaPlus[positionSecond-length(storeMinus)]

# cat("Second best: ",newThresh, " (conditionnal on ",threshRef, " ) \t SSR", min(store2, na.rm=TRUE), "\n")
list(newThresh=newThresh, SSR=min(store2, na.rm=TRUE))
}	#end function condistep


###Applying the function for conditional search to original data
More<-list(d=bestDelay, Z=Z, Y=Y,trans=z)
Thresh2<-condiStep(allgammas, bestThresh,ninter=ninter, fun=SSR_2thresh, MoreArgs=More)$newThresh
Thresh3<-condiStep(allgammas, Thresh2,ninter=ninter, fun=SSR_2thresh, MoreArgs=More)
smallThresh<-min(Thresh2, Thresh3$newThresh)
bigThresh<-max(Thresh2, Thresh3$newThresh)

cat("Second best: ",Thresh2, " (conditionnal on ",bestThresh, ")\n")
cat("Iterative best: ",Thresh3$newThresh, " (conditionnal on ",Thresh2, ")\n")

Sigma_mod2thresh<-Sigma_2thresh(gam1=smallThresh,gam2=bigThresh,d=bestDelay, Z=Z, Y=Y,trans=z)

###F test for original data
LRtest12<-as.numeric(t*(log(det(Sigma))-log(det(Sigma_mod1thresh))))
LRtest13<-as.numeric(t*(log(det(Sigma))-log(det(Sigma_mod2thresh))))
LRtest23<-as.numeric(t*(log(det(Sigma_mod1thresh))-log(det(Sigma_mod2thresh))))
LRs<-c(LRtest12, LRtest13, LRtest23)

##############################
###Bootstrap for the F test
##############################

### Reconstruction of series from linear model for LRtest 1vs2 and 1vs3
#initial data	
res_lin<-res
Yb<-matrix(0, nrow=nrow(y), ncol=k)		#Delta Y term
Yb[1:m,]<-y[1:m,]			

bootlinear<-function(x){
resi<-rbind(matrix(0,nrow=m, ncol=k),res_lin[sample(seq_len(nrow(res_lin)), replace=TRUE),])
resi<-rbind(matrix(0,nrow=m, ncol=k),res_lin)		#Uncomment this line to check the bootstrap

for(i in (m+1):(nrow(y))){
	Yb[i,]<-rowSums(cbind(B[,1], B[,-1]%*%matrix(t(Yb[i-c(1:m),]), ncol=1),resi[i,]))
}
return(Yb)
}#end bootlinear
# print(cbind(y, Yb))
### Reconstruction of series from 1 thresh model for LRtest 2vs3

dummydown <- ifelse(z[,bestDelay]<=bestThresh, 1, 0)
regimedown <- dummydown*Z
Z2 <- t(cbind(regimedown, (1-dummydown)*Z))		# dim k(p+1) x t	
B1thresh<-tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))
res_thresh <- t(Y - B1thresh%*%Z2)	#SSR

B1tDown<-B1thresh[,seq_len(k*m+1)]
B1tUp<-B1thresh[,-seq_len(k*m+1)]

#initial data	
Yb2<-matrix(0, nrow=nrow(y), ncol=k)		#Delta Y term
Yb2[1:m,]<-y[1:m,]			
z2<-vector("numeric", length=nrow(y))
z2[1:m]<-y[1:m,]%*%combin			

boot1thresh<-function(x){
resiT<-rbind(matrix(0,nrow=m, ncol=k),res_thresh)	
 resiT<-rbind(matrix(0,nrow=m, ncol=k),res_thresh[sample(seq_len(nrow(res_thresh)), replace=TRUE),])

for(i in (m+1):(nrow(y))){
	if(round(z2[i-bestDelay],ndig)<=bestThresh) 
		Yb2[i,]<-rowSums(cbind(B1tDown[,1], B1tDown[,-1]%*%matrix(t(Yb2[i-c(1:m),]), ncol=1),resiT[i,]))
	else
		Yb2[i,]<-rowSums(cbind(B1tUp[,1], B1tUp[,-1]%*%matrix(t(Yb2[i-c(1:m),]), ncol=1),resiT[i,]))
	z2[i]<-Yb2[i,]%*%combin
	}
# print(cbind(x,z2))
return(Yb2)
}#end boot1thresh



#####Bootstrap loop
model<-match.arg(test)
bootModel<-switch(model, "1vs"=bootlinear, "2vs3"=boot1thresh)

bootstraploop<-function(x, thVar=NULL){

xboot<-round(bootModel(x=x),ndig)


# Sigma of linear boot model
string<-embed(xboot,m+1)
Yboot <- string[,seq_len(k)] 	#
Zb <- string[, -seq_len(k)]	#Lags matrix
if(trend==TRUE)
	Zb <- cbind(1,Zb)
Bboot<-t(Yboot)%*%Zb%*%solve(t(Zb)%*%Zb)		#B: OLS parameters, dim 2 x npar
resboot<-Yboot-Zb%*%t(Bboot)
Sigmab<- matrix(1/T*crossprod(resboot),ncol=k)

#grid for threshold boot model

if(!missing(thVar)) 
	{stop("thVar currently badly implemented"); zcombin<-thVar}
else 
	zcombin<-xboot%*%combin
zb <- embed(zcombin,m+1)[,seq_len(max(thDelay))+1]


#  print(cbind(z,zb))

allgammasb<-sort(unique(zb[,b]))
ng<-length(allgammasb)
gammasb<-allgammasb[(ceiling(trim*ng)+1):floor((1-trim)*ng-1)]


###One threshold Search on bootstrap data

IDSb<-as.matrix(expand.grid(thDelay, gammasb))
resultb <- apply(IDSb, 1, SSR_1thresh,Z=Zb, Y=t(Yboot),trans=zb)
bestDelayb<-IDSb[which.min(resultb),1]
bestThreshb<-IDSb[which.min(resultb),2]

Sigma_mod1threshb<-Sigma_1thresh(gam1=bestThreshb, d=bestDelayb,Z=Zb,Y=t(Yboot), trans=zb)

###Two threshold Search (conditional and 1 iteration) on bootstrap data
Moreb<-list(d=bestDelayb, Z=Zb, Y=t(Yboot),trans=zb)
Thresh2b<-condiStep(allgammasb, bestThreshb,ninter=ninter, fun=SSR_2thresh, MoreArgs=Moreb)$newThresh
Thresh3b<-condiStep(allgammasb, Thresh2b,ninter=ninter, fun=SSR_2thresh, MoreArgs=Moreb)
smallThreshb<-min(Thresh2b, Thresh3b$newThresh)
bigThreshb<-max(Thresh2b, Thresh3b$newThresh)


Sigma_mod2threshb<-Sigma_2thresh(gam1=smallThreshb,gam2=bigThreshb,d=bestDelayb, Z=Zb, Y=t(Yboot),trans=zb)

###Test statistic on bootstrap data
LRtest12b<-as.numeric(t*(log(det(Sigmab))-log(det(Sigma_mod1threshb))))
LRtest13b<-as.numeric(t*(log(det(Sigmab))-log(det(Sigma_mod2threshb))))
LRtest23b<-as.numeric(t*(log(det(Sigma_mod1threshb))-log(det(Sigma_mod2threshb))))

list(LRtest12b, LRtest13b, LRtest23b)
 }#end of bootstraploop



LRtestboot<-replicate(n=nboot,bootstraploop(x=data))
LRtestboot12<-unlist(LRtestboot[1,])
LRtestboot13<-unlist(LRtestboot[2,])
LRtestboot23<-unlist(LRtestboot[3,])

PvalBoot12<-mean(ifelse(LRtestboot12>LRtest12,1,0))
CriticalValBoot12<-quantile(LRtestboot12, probs=c(0.9, 0.95, 0.975,0.99))
PvalBoot13<-mean(ifelse(LRtestboot13>LRtest12,1,0))
CriticalValBoot13<-quantile(LRtestboot13, probs=c(0.9, 0.95, 0.975,0.99))
PvalBoot23<-mean(ifelse(LRtestboot23>LRtest12,1,0))
CriticalValBoot23<-quantile(LRtestboot23, probs=c(0.9, 0.95, 0.975,0.99))

if(test=="1vs"){
	CriticalValBoot<-list(Test1vs2=CriticalValBoot12, Test1vs3=CriticalValBoot13)
	PvalBoot<-list(Test1vs2=PvalBoot12,Test1vs3=PvalBoot13)
	}
else{
	CriticalValBoot<-CriticalValBoot23
	PvalBoot<-PvalBoot23
	}


###Grahical output
if(plot==TRUE&nboot>0){
	if(test=="1vs"){
	layout(c(1,2))
	plot(density(LRtestboot12), xlab="LRtest12", xlim=c(0,max(LRtest12+1,max(LRtestboot12))),ylim=c(0,max(density(LRtestboot12)$y,dchisq(0:LRtest12, df=1+m))), main="Test linear VAR vs 1 threshold TVAR")
	abline(v=LRtest12, lty=2, col=2)
	curve(dchisq(x, df=1+k*m, ncp=0), from=0, n=LRtest12+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))


	plot(density(LRtestboot13), xlab="LRtest13", xlim=c(0,max(LRtest13+1,max(LRtestboot13))),ylim=c(0,max(density(LRtestboot13)$y,dchisq(0:LRtest12, df=2*(1+m)))),main="Test linear VAR vs 2 thresholds TVAR")
	abline(v=LRtest13, lty=2, col=2)
	curve(dchisq(x, df=2*(1+k*m), ncp=0), from=0, n=LRtest13+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	}
	else {
	plot(density(LRtestboot23), xlab="LRtest23", xlim=c(0,max(LRtest23+1,LRtestboot23)), ylim=c(0,max(density(LRtestboot23)$y,dchisq(0:LRtest12, df=1+m))), main="Test 1 threshold TVAR vs 2 thresholds TVAR")
	abline(v=LRtest23, lty=2, col=2)
	curve(dchisq(x, df=1+k*m, ncp=0), from=0, n=LRtest23+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	}
}


#nlar=extend(nlar(str, coef = res$coef, fit = res$fitted.values, res = res$residuals, k = res$k,
#list( model.specific = res),"setar")
return(list(bestDelay=bestDelay, LRtest.val=LRs, Pvalueboot=PvalBoot, CriticalValBoot=CriticalValBoot))
}#End of thw whole function




if(FALSE){ #usage example
environment(TVAR_LRtest)<-environment(star)
data(zeroyld)
data<-zeroyld


TVAR_LRtest(data, m=2, mTh=1,thDelay=1:2, nboot=2, plot=TRUE, trim=0.1, test="1vs")
}

