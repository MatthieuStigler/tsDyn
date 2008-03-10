OlsTVAR <- function(data, lag, demean = c( "const", "trend","none", "both"), model=c("TAR", "MTAR"), commonInter=FALSE, nthresh=1,thDelay=1, mTh=1,thVar, trim=0.1,ngrid, gamma=NULL,  around, plot=TRUE, dummyToBothRegimes=TRUE){
y <- as.matrix(data)
Torigin <- nrow(y) 	#Size of original sample
T <- nrow(y) 		#Size of start sample
p <- lag
t <- T-p 		#Size of end sample
k <- ncol(y) 		#Number of variables
t<-T-p			#Size of end sample
if(is.null(colnames(data)))
	colnames(data)<-paste("Var", c(1:k), sep="")
if(max(thDelay)>p)
	stop("Max of thDelay should be smaller or equal to the number of lags")
if(dummyToBothRegimes==FALSE&nthresh!= 1) 
	warning("The 'dummyToBothRegimes' argument is only relevant for one threshold models")
model<-match.arg(model)
demean<-match.arg(demean)

Y <- y[(p+1):T,] #
Z <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

if(demean=="const")
	Z<-cbind(1, Z)
else if(demean=="trend")
	Z<-cbind(seq_len(t), Z)
else if(demean=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)
if(commonInter & demean!="const")
	stop("commonInter argument only avalaible with demean = const")
npar <- ncol(Z)			#Number of parameters

##################
###Linear model
#################
 B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar


allpar<-ncol(B)*nrow(B)

rownames(B)<-paste("Equation",colnames(data))
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))

if(demean=="const")
	Bnames<-c("Intercept",LagNames)
else if(demean=="trend")
	Bnames<-c("Trend",LagNames)
else if(demean=="both")
	Bnames<-c("Intercept","Trend",LagNames)
else 
	Bnames<-c(LagNames)
colnames(B)<-Bnames

res<-Y-Z%*%t(B)

Sigma<- matrix(1/t*crossprod(res),ncol=k,dimnames=list(colnames(data), colnames(data)))
nlike<-log(det(Sigma))		#	nlike=(t/2)*log(det(sige));
aic<-t*nlike+2*(allpar)
bic<-t*nlike+log(t)*(allpar)	#bic #=nlike+log10(t)*4*(1+k); ###BIC
info_Lin<-c(aic, bic)


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
} ###Combination (or single value indicating position) of contemporaneous variables
else {
	if (!length(mTh)%in%c(1,k))
		stop("length of 'mTh' should be equal to the number of variables, or just one")
	if(length(mTh)==1) {
		if(mTh>k)
			stop("Unable to select the ",mTh, "variable for the threshold. Please see again mTh ")
		combin <- matrix(0,ncol=1, nrow=k)
		combin[mTh,]<-1
	}
	else 
		combin<-matrix(mTh,ncol=1, nrow=k)
	zcombin <- y %*% combin
	if(model=="MTAR"){
		if(max(thDelay)<p)
			z<-embed(diff(zcombin),p)[,seq_len(max(thDelay))+1]
		else if(max(thDelay)==p){
			z<-embed(diff(zcombin),p+1)[,seq_len(max(thDelay))+1]
			z<-rbind(0,as.matrix(z))}
	}
	else
		z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
}

trans<-as.matrix(z)

###############################
###Grid for transition variable
###############################

allgammas <- sort(unique(trans[,1]))
nga <- length(allgammas)
ninter <- round(trim*nga)
gammas <- allgammas[(trim*nga):((1-trim)*nga)]
if(!missing(ngrid)){
	gammas <- allgammas[seq(from=ceiling(trim*nga), to=floor((1-trim)*nga), length.out=ngrid)]
}
if(!missing(gamma)){
	gammas<-gamma
	plot<-FALSE}
Y<-t(Y)					#dim k x t

aroundGrid <- function(around,allvalues,ngrid,trim){
	ng <- length(allvalues)
	wh.around <- which.min(abs(allvalues-around))
	if(length(which(allvalues==around))==0)
		stop("\nThe value ", around, " did not match to existing ones", allvalues[wh.around], "was taken instead")
	if(length(wh.around)>1){
		warning("\nThere were", length(wh.around)," values corresponding to the around argument. The first one was taken")
		wh.around<-wh.around[1]}
	ar <- c((wh.around-round(ngrid/2)): (wh.around+round(ngrid/2)))		#Values around the point
	ar2 <- ar[ar>=round(trim*ng)&ar<=round((1-trim)*ng)]			#Bounding with trim 
	gammas <- allvalues[ar2]
	return(gammas)
}#end if missing around

if(!missing(around)){
	if(missing(ngrid)) ngrid<-20
	if(length(around)==1)
		gammas <- aroundGrid(around, allgammas,ngrid,trim)
	if(length(around)==2) {
		gammas <- aroundGrid(around[1], allgammas,ngrid,trim)
		gammas2 <- aroundGrid(around[2], allgammas,ngrid,trim)
	}
}

######################
###One threshold model					
######################

#Model with dummy applied to only one regime
loop1_onedummy <- function(gam1, d){
	##Threshold dummies
	dummyDown <- ifelse(trans[,d]<gam1, 1,0) * Z
	ndown<-mean(dummyDown)
	regimeDown<-dummyDown*Z
	##SSR
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- t(cbind(regimeDown, Z))		# dim k(p+1) x t
		B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
		res<-crossprod(c( Y - B1 %*% Z1))}
	else
		res<-NA
	return(res)
} #end of the function


#Model with dummy applied to both regimes
loop1_twodummy <- function(gam1, d){
	##Threshold dummies
	d1<-ifelse(trans[,d]<gam1, 1,0)
	ndown<-mean(d1)
	##SSR
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- t(cbind(d1 * Z, (1-d1)*Z))		# dim k(p+1) x t
		B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
		res<-crossprod(c( Y - B1 %*% Z1))}
	else
		res<-NA
	return(res)
} #end of the function

#Model with dummy applied to both regimes and a common intercept
loop1_twodummy_oneIntercept <- function(gam1, d){
	##Threshold dummies
	d1<-ifelse(trans[,d]<gam1, 1,0)
	ndown<-mean(d1)
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- t(cbind(1,d1 * Z[,-1], (1-d1)*Z[,-1]))		# dim k(p+1) x t
		B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
		res<-crossprod(c( Y - B1 %*% Z1))}
	else
		res<-NA
	return(res)
} #end of the function

onesearch <- function(thDelay,gammas){
	grid1 <- expand.grid(thDelay,gammas)				#grid with delay and gammas
	if(dummyToBothRegimes){
		if(commonInter)
			store<-mapply(loop1_twodummy_oneIntercept,d=grid1[,1],gam1=grid1[,2])
		else
			store <- mapply(loop1_twodummy,d=grid1[,1],gam1=grid1[,2])	#search for values of grid
	}
	else	
		store <- mapply(loop1_onedummy,d=grid1[,1],gam1=grid1[,2])
	posBestThresh <- which(store==min(store, na.rm=TRUE), arr.ind=TRUE)[1]

	if(plot&is.null(gamma)){
		result1 <- cbind(grid1,store)
		col <- rep(thDelay,length.out=nrow(result1))
		plot(result1[,2], result1[,3], col=col,xlab="Treshold Value",ylab="SSR", main="Results of the grid search")
		legend("topleft", pch=1, legend=paste("Threshold Delay", thDelay), col=thDelay)
	}
	cat("Best unique threshold", grid1[posBestThresh,2],"\n")
	if(length(thDelay)>1)
		cat("Best Delay", grid1[posBestThresh,1],"\n")
	list(bestThresh=grid1[posBestThresh,2],bestDelay=grid1[posBestThresh,1])
}#end of function one search

#######################
###Two thresholds model
#######################

loop2 <- function(gam1, gam2,d){
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
		res <- crossprod(c( Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2))	#SSR
	}
	else
		res <- NA
	return(res)
}

loop2_oneIntercept <- function(gam1, gam2,d){
	##Threshold dummies
	dummydown <- ifelse(trans[,d]<gam1, 1, 0)
	regimedown <- dummydown*Z[,-1]
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,d]>=gam2, 1, 0)
	regimeup <- dummyup*Z[,-1]
	nup <- mean(dummyup)
	##SSR from TVAR(3)
	#print(c(ndown,1-nup-ndown,nup))
	if(min(nup, ndown, 1-nup-ndown)>trim){
		Z2 <- t(cbind(1,regimedown, (1-dummydown-dummyup)*Z, regimeup))		# dim k(p+1) x t	
		res <- crossprod(c( Y - tcrossprod(Y,Z2) %*% solve(tcrossprod(Z2))%*%Z2))	#SSR
	}
	else
		res <- NA
	return(res)
}
############################
###Search for one threshold
############################
if(nthresh==1){
if(!missing(around))
	gammas <- aroundGrid(around,allgammas,ngrid,trim)
bestone <- onesearch(thDelay,gammas)
bestThresh <- bestone$bestThresh
bestDelay <- bestone$bestDelay
}

############################
###Search for two threshold
############################

###Conditional search
if(nthresh==2){

###Conditionnal step
bestone <- onesearch(thDelay, gammas)
bestThresh <- bestone$bestThresh
bestDelay <- bestone$bestDelay


condiStep<-function(allgammas, threshRef, delayRef,ninter, fun){

wh.thresh <- which.min(abs(allgammas-threshRef))

#search for a second threshold smaller than the first
if(wh.thresh>2*ninter){
	gammaMinus<-allgammas[seq(from=ninter, to=wh.thresh-ninter)]
	storeMinus <- mapply(fun,gam1=gammaMinus,gam2=threshRef, d=delayRef)	
}
else
	storeMinus <- NA

#search for a second threshold higher than the first
if(length(wh.thresh<length(allgammas)-2*ninter)){
	gammaPlus<-allgammas[seq(from=wh.thresh+ninter, to=length(allgammas)-ninter)]
	storePlus <- mapply(fun,gam1=threshRef,gam2=gammaPlus, d=delayRef)
}
else
	storePlus <- NA

#results
store2 <- c(storeMinus, storePlus)

positionSecond <- which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)
if(positionSecond<=length(storeMinus))
	newThresh<-gammaMinus[positionSecond]
else
	newThresh<-gammaPlus[positionSecond-length(storeMinus)]

cat("Second best: ",newThresh, " (conditionnal on ",threshRef, " ) \t SSR", min(store2, na.rm=TRUE), "\n")
list(threshRef=threshRef, newThresh=newThresh)
}

secondBestThresh<-condiStep(allgammas, threshRef=bestThresh, delayRef=bestDelay,ninter=ninter, fun=loop2)$newThresh

###Iterative step
if(commonInter)
	func<-loop2_oneIntercept
else
	func<-loop2
condiStep(allgammas, threshRef=secondBestThresh, delayRef=bestDelay,ninter=ninter, fun=func)

###Alternative step: grid around the points from first step
smallThresh <- min(bestThresh,secondBestThresh)
gammasDown <- aroundGrid(around=smallThresh,allgammas,ngrid=30, trim=trim)

bigThresh <- max(bestThresh,secondBestThresh)
gammasUp <- aroundGrid(around=bigThresh,allgammas,ngrid=30, trim=trim)

storeIter <- matrix(NA,ncol=length(gammasUp), nrow=length(gammasDown))

#Grid search
for(i in seq_len(length(gammasDown))){
	gam1 <- gammasDown[i]
	for(j in 1: length(gammasUp)){
		gam2 <- gammasUp[j]
		storeIter[i,j] <- func(gam1=gam1, gam2=gam2, d=bestDelay)
	}
}

#Finding the best result
positionIter <- which(storeIter==min(storeIter, na.rm=TRUE), arr.ind=TRUE)
rIter <- positionIter[1]
cIter <- positionIter[2]

gamma1Iter <- gammasDown[rIter]
gamma2Iter <- gammasUp[cIter]

bestThresh <- c(gamma1Iter, gamma2Iter)

cat("\nSecond step best thresholds", bestThresh, "\t\t\t SSR", min(storeIter, na.rm=TRUE), "\n")
}#end if nthresh=2

###Search both thresholds with d given

if(nthresh==3){
bestDelay <- thDelay
if(missing(gamma)==FALSE){
	gammas <- gamma[1]
	gammas2 <- gamma[2]
	ninter<- 2
}
if(missing(around)==FALSE){
	if(length(around)!=2)
		stop("Please give two thresholds possible values to search around")
	gammas <- aroundGrid(min(around), allgammas, ngrid=ngrid, trim=trim)
	gammas2 <- aroundGrid(max(around), allgammas, ngrid=ngrid, trim=trim)
}
else {
	gammas2 <- gammas
	if(length (gammas) * length(gammas2)/2>10000)
		cat("The function will compute about", length(gammas)*length(gammas2)/2, "operations. Take a coffee and come back\n")
}
if(length(thDelay)>1) stop("length of thDelay should not be bigger than 1. The whole search is made only upon the thresholds with given delay")

store3 <- matrix(NA,ncol=length(gammas2), nrow=length(gammas))

###Loop
for(i in seq_len(length(gammas))){
	gam1 <- gammas[i]
	for (j in seq_len(length(gammas))){
		gam2 <- gammas2[j]
		store3[i,j] <- loop2(gam1, gam2, d=bestDelay)
	}
}

position <- which(store3==min(store3, na.rm=TRUE), arr.ind=TRUE)
r <- position[1]
c <- position[2]

gamma1 <- gammas[r]
gamma2 <- gammas2[c]
bestThresh <- c(gamma1, gamma2)

}#end n

#############
###Best Model
#############
if(commonInter)
	val<- 1
else
	val<- -(seq_len(ncol(Z)))

if(nthresh==1){
	dummydown <- ifelse(trans[,bestDelay]<bestThresh, 1, 0)
	ndown <- mean(dummydown)
	regimeDown <- dummydown*Z[,-val]
	if(dummyToBothRegimes) 
		regimeUp<-(1-dummydown)*Z[,-val]
	else regimeUp<-Z
	if(commonInter)
		Zbest<-t(cbind(1,regimeDown,regimeUp))
	else
		Zbest <- t(cbind(regimeDown,regimeUp))		# dim k(p+1) x t
}

if(nthresh==2|nthresh==3){
	dummydown <- ifelse(trans[,bestDelay]<bestThresh[1], 1,0)
	ndown <- mean(dummydown)
	regimedown <- dummydown*Z[,-val]
	dummyup <- ifelse(trans[,bestDelay]>=bestThresh[2], 1,0)
	nup <- mean(dummyup)
	regimeup <- dummyup*Z[,-val]
	if(commonInter)
		Zbest <- t(cbind(1,regimedown,(1-dummydown-dummyup)*Z[,-1], regimeup))	# dim k(p+1) x t
	else
		Zbest <- t(cbind(regimedown,(1-dummydown-dummyup)*Z, regimeup))	# dim k(p+1) x t
}

Bbest <- Y %*% t(Zbest) %*% solve(Zbest %*% t(Zbest))
resbest <- t(Y - Bbest %*% Zbest)
SSRbest <- as.numeric(crossprod(c(resbest)))
nparbest<-nrow(Bbest)*ncol(Bbest)

Sigmabest<-matrix(1/t*crossprod(resbest),ncol=k,dimnames=list(colnames(data), colnames(data)))
SigmabestOls<-Sigmabest*(t/(t-ncol(Bbest)))

VarCovB<-solve(tcrossprod(Zbest))%x%SigmabestOls
StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)

Tvalue<-Bbest/StDevB
Pval<-pt(abs(Tvalue), df=(ncol(Zbest)-nrow(Zbest)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(ncol(Zbest)-nrow(Zbest)), lower.tail=TRUE)
Pval<-round(Pval,4)


### Info criteria
nlikebest<-log(det(Sigmabest))			#nlike=(t/2)*log(det(sige));
aicbest<-t*nlikebest+2*(nparbest)
bicbest<-t*nlikebest+log(t)*(nparbest)	#bic #=nlike+log10(t)*4*(1+k); ###BIC
info_Thresh<-c(aicbest, bicbest)
aicbest2<-t*nlikebest+2*(nparbest+nthresh)
bicbest2<-t*nlikebest+log(t)*(nparbest+nthresh)	#bic #=nlike+log10(t)*4*(1+k); ###BIC
info_Thresh2<-c(aicbest2, bicbest2)

nobsdown<-sum(dummydown)
nobsup<-length(dummydown)-nobsdown

SigmabestUnder<-matrix(1/nobsdown*crossprod(resbest*dummydown),ncol=k)
SigmabestOver<-matrix(1/nobsup*crossprod(resbest*(1-dummydown)),ncol=k)


nlikebest<-nobsdown*log(det(SigmabestUnder))+nobsup*log(det(SigmabestOver))
pooled_aicbest<-nlikebest+2*(nparbest)
pooled_bicbest<-nlikebest+(log(nobsup)+log(nobsdown))*(nparbest)	#bic #=nlike+log10(t)*4*(1+k); ###BIC
pooled_info_Thresh<-c(pooled_aicbest, pooled_bicbest)




###naming and dividing B
rownames(Bbest) <- paste("Equation", colnames(data))
Bcolnames <- c("Trend", c(paste(rep(colnames(data),p),"t", -rep(1:p, each=k))))
sBnames<-Bnames[-which(Bnames=="Intercept")]

if(nthresh==1){
	if(commonInter){
		colnames(Bbest)<-c("Intercept",paste("Dn",sBnames), paste("Up",sBnames))
		colnames(Pval)<-c("Intercept",paste("Dn",sBnames), paste("Up",sBnames))
		Blist<-Bbest
		Plist<-Pval}
	else{
		colnames(Bbest) <- rep(Bnames,2)
		colnames(Pval) <- rep(Bnames,2)
		Bdown <- Bbest[,c(1:npar)]
		Bup <- Bbest[,-c(1:npar)]
		Pdown <- Pval[,c(1:npar)]
		Pup <- Pval[,-c(1:npar)]
		Blist <- list(Bdown=Bdown, Bup=Bup)
		Plist <- list(Pdown=Pdown, Pup=Pup)}
	nobs <- c(ndown=ndown, nup=1-ndown)
}
else{
	if(commonInter){
		colnames(Bbest)<-c("Intercept",paste("Dn",sBnames), paste("Mi",sBnames), paste("Up",sBnames))
		colnames(Pval)<-c("Intercept",paste("Dn",sBnames), paste("Mi",sBnames), paste("Up",sBnames))
		Blist<-Bbest
		Plist<-Pval}
	else{
		colnames(Bbest)<-rep(Bnames,3)
		colnames(Pval)<-rep(Bnames,3)
		Bdown <- Bbest[,c(1:npar)]
		Bmiddle <- Bbest[,c(1:npar)+npar]
		Bup <- Bbest[,c(1:npar)+2*npar]
		Pdown <- Pval[,c(1:npar)]
		Pmiddle <- Pval[,c(1:npar)+npar]
		Pup <- Pval[,c(1:npar)+2*npar]			
		colnames(Bmiddle) <- Bnames
		colnames(Pmiddle) <- Bnames
		Blist <- list(Bdown=Bdown, Bmiddle=Bmiddle,Bup=Bup)
		Plist <- list(Pdown=Pdown, Pmiddle=Pmiddle,Pup=Pup)}
	nobs <- c(ndown=ndown, nmiddle=1-nup-ndown,nup=nup)
}


list(resids=list(lin=res, Thresh=resbest), VAR=VarCovB, Thresh=bestThresh, Parameters=Blist, Pvalues=Plist,SSR=SSRbest, aic_bic=list(info_Lin=info_Lin,info_Thresh=info_Thresh,info_Thresh2=info_Thresh2), nobs_regimes=nobs)
}	#end of the whole function

if(FALSE) { #usage example
###Hansen Seo data
data(zeroyld)
data<-zeroyld

OlsTVAR(data[1:100,], lag=2, nthresh=2, thDelay=1:2,trim=0.1, plot=FALSE, commonInter=FALSE, demean="const")
#lag2, 2 thresh, trim00.05: 561.46
}

