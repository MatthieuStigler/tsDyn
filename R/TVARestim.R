TVAR <- function(data, lag, include = c( "const", "trend","none", "both"), model=c("TAR", "MTAR"), commonInter=FALSE, nthresh=1,thDelay=1, mTh=1,thVar, trim=0.1,ngrid, gamma=NULL,  around, plot=TRUE, dummyToBothRegimes=TRUE, trace=TRUE, trick="for", max.iter=2){
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
include<-match.arg(include)

Y <- y[(p+1):T,] #
Z <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

if(include=="const")
	Z<-cbind(1, Z)
else if(include=="trend")
	Z<-cbind(seq_len(t), Z)
else if(include=="both")
	Z<-cbind(rep(1,t),seq_len(t), Z)
if(commonInter & include!="const")
	stop("commonInter argument only avalaible with include = const")
npar <- ncol(Z)			#Number of parameters


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
###One threshold functions					
######################

#Model with dummy applied to only one regime
loop1_onedummy <- function(gam1, thDelay){
	##Threshold dummies
	dummyDown <- ifelse(trans[,thDelay]<gam1, 1,0) * Z
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
loop1_twodummy <- function(gam1, thDelay){
	##Threshold dummies
	d1<-ifelse(trans[,thDelay]<gam1, 1,0)
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
loop1_twodummy_oneIntercept <- function(gam1, thDelay){
	##Threshold dummies
	d1<-ifelse(trans[,thDelay]<gam1, 1,0)
	ndown<-mean(d1)
	if(min(ndown, 1-ndown)>=trim){
		Z1 <- t(cbind(1,d1 * Z[,-1], (1-d1)*Z[,-1]))		# dim k(p+1) x t
		B1 <- tcrossprod(Y,Z1) %*% solve(tcrossprod(Z1))
		res<-crossprod(c( Y - B1 %*% Z1))}
	else
		res<-NA
	return(res)
} #end of the function


#######################
###Two thresholds functions
#######################

loop2 <- function(gam1, gam2,thDelay){
	##Threshold dummies
	dummydown <- ifelse(trans[,thDelay]<gam1, 1, 0)
	regimedown <- dummydown*Z
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,thDelay]>=gam2, 1, 0)
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

loop2_oneIntercept <- function(gam1, gam2,thDelay){
	##Threshold dummies
	dummydown <- ifelse(trans[,thDelay]<gam1, 1, 0)
	regimedown <- dummydown*Z[,-1]
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[,thDelay]>=gam2, 1, 0)
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
if(!missing(around))
	gammas <- aroundGrid(around,allgammas,ngrid,trim)
if(dummyToBothRegimes){
	if(commonInter)
		func<-loop1_twodummy_oneIntercept
	else
		func <-loop1_twodummy}
else	
	func <- loop1_onedummy

if(nthresh==1){
bestone <- onesearch(thDelay,gammas, fun=func, trace=trace, gamma=gamma, plot=plot)
bestThresh <- bestone$bestThresh
bestDelay <- bestone$bestDelay
}

############################
###Search for two threshold
############################

###Conditional search
if(nthresh==2){

###Conditionnal step
bestone <- onesearch(thDelay, gammas, fun=func, trace=trace, gamma=gamma, plot=plot)
bestThresh <- bestone$bestThresh
bestDelay <- bestone$bestDelay

if(commonInter)
	func2<-loop2_oneIntercept
else
	func2<-loop2
# secondBestThresh<-condiStep(allgammas, threshRef=bestThresh, delayRef=bestDelay,ninter=ninter, fun=func2)$newThresh
# step2FirstBest<-condiStep(allgammas, threshRef=secondBestThresh, delayRef=bestDelay,ninter=ninter, fun=func2)

last<-condiStep(allgammas, threshRef=bestThresh, delayRef=bestDelay,ninter=ninter, fun=func2)
i<-1
while(i<max.iter){
	b<-condiStep(allgammas, threshRef=last$newThresh, delayRef=bestDelay,ninter=ninter, fun=func2)
	if(b$SSR<last$SSR){	#minimum still not reached
		i<-i+1
		last<-b}
	else{			#minimum reached
		i<-max.iter
		last<-b}
}

bests<-c(last$threshRef, last$newThresh)

###Alternative step: grid around the points from first step
smallThresh <- min(bests)		#bestThresh,secondBestThresh)
gammasDown <- aroundGrid(around=smallThresh,allgammas,ngrid=30, trim=trim)

bigThresh <- max(bests)			#bestThresh,secondBestThresh)
gammasUp <- aroundGrid(around=bigThresh,allgammas,ngrid=30, trim=trim)

bestThresh<-grid(gammasUp, gammasDown, bestDelay, fun=func2, method=trick)

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
		store3[i,j] <- loop2(gam1, gam2, thDelay=bestDelay)
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
fitted<-Bbest %*% Zbest
resbest <- t(Y - fitted)
SSRbest <- as.numeric(crossprod(c(resbest)))
nparbest<-nrow(Bbest)*ncol(Bbest)

Sigmabest<-matrix(1/t*crossprod(resbest),ncol=k,dimnames=list(colnames(data), colnames(data)))
SigmabestOls<-Sigmabest*(t/(t-ncol(Bbest)))

VarCovB<-solve(tcrossprod(Zbest))%x%SigmabestOls
StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)

Tvalue<-Bbest/StDevB
Pval<-pt(abs(Tvalue), df=(ncol(Zbest)-nrow(Zbest)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(ncol(Zbest)-nrow(Zbest)), lower.tail=TRUE)
Pval<-round(Pval,4)




###Pooled AIC
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
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))
if(include=="const")
	Bnames<-c("Intercept",LagNames)
else if(include=="trend")
	Bnames<-c("Trend",LagNames)
else if(include=="both")
	Bnames<-c("Intercept","Trend",LagNames)
else 
	Bnames<-c(LagNames)

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

specific<-list()
specific$Thresh<-bestThresh
specific$nthresh<-nthresh
specific$nreg<-nthresh+1
specific$rowB<-npar
specific$nobs<-nobs
specific$oneMatrix<-commonInter
specific$threshEstim<-ifelse(is.null(gamma), TRUE, FALSE)

z<-list(coefficients=Blist, residuals=resbest, VAR=VarCovB,   Pvalues=Plist, nobs_regimes=nobs, k=k, t=t, nparB=nparbest, fitted.values=fitted, lag=lag, include=include,model.specific=specific)
class(z)<-c("TVAR","nlVar")
return(z)
}	#end of the whole function

if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
data(zeroyld)
dat<-zeroyld

VAR<-TVAR(dat[1:100,], lag=2, nthresh=2,thDelay=1,trim=0.1, plot=FALSE, commonInter=TRUE, include="const")
#lag2, 2 thresh, trim00.05: 561.46
class(VAR)
VAR
print(VAR)
logLik(VAR)
AIC(VAR)
BIC(VAR)
coef(VAR)
deviance(VAR)
summary(VAR)
toLatex(VAR)
toLatex(summary(VAR))
VAR[["model.specific"]][["oneMatrix"]]
###TODO
#pre specified gamma
}

print.TVAR<-function(x,...){
# 	NextMethod(...)
	cat("Model TVAR with ", x$model.specific$nthresh, " thresholds\n\n")
	print(x$coefficients)
	cat("\nThreshold value")
	print(x$model.specific$Thresh)
}

summary.TVAR<-function(object, ...){
	x<-object
	if(x$model.specific$oneMatrix) {
		x$coefficients<-list(x$coefficients)
		x$Pvalues<-list(x$Pvalues)}
	ab<-list()
	symp<-list()
	stars<-list()
	for(i in 1:length(x$Pvalues)){
		symp[[i]] <- symnum(x$Pvalues[[i]], corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
		stars[[i]]<-matrix(symp[[i]], nrow=nrow(x$Pval[[i]]))
		ab[[i]]<-matrix(paste(round(x$coefficients[[i]],4),"(", x$Pval[[i]],")",stars[[i]], sep=""), nrow=nrow(x$Pvalues[[1]]))
		dimnames(ab[[i]])<-dimnames(x$coefficients[[1]])
	}
	attributes(ab)<-attributes(x$coefficients)
	x$bigcoefficients<-ab
	x$aic<-AIC.nlVar(x)
	x$bic<-BIC.nlVar(x)
	x$SSR<-deviance(x)
	class(x)<-c("summary.TVAR", "TVAR")
	return(x)
}

print.summary.TVAR<-function(x,...){
	cat("Model TVAR with ", x$model.specific$nthresh, " thresholds\n")
	cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
	cat("\nNumber of variables:", x$k,"\tNumber of estimated parameters", x$npar)
	cat("\nAIC",x$aic , "\tBIC", x$bic, "\t SSR", x$SSR,"\n\n")
	print(noquote(x$bigcoefficients))	
	cat("\nThreshold value",x$model.specific$Thresh)
	if(!x$model.specific$threshEstim)
		cat(" (user specified)")
	cat("\nPercentage of Observations in each regime", x$model.specific$nobs)
}



toLatex.TVAR<-function(object,..., digits=4){
	x<-object
	if(x$model.specific$oneMatrix){
		x$coefficients<-list(x$coefficients)}
	if(inherits(x,"summary.TVAR")){
		coef<-x$bigcoefficients
		for(i in 1:length(coef)){
			coef[[i]]<-gsub(")", ")^{",coef[[i]], extended=FALSE)
			coef[[i]]<-matrix(paste(coef[[i]], "}"), ncol=ncol(coef[[i]]), nrow=nrow(coef[[i]]))
		}
	}
	else{
		coef<-as.list(x$coefficients)
		coef<-rapply(coef,round, digits=digits, how="list")}
	varNames<-rownames(coef[[1]])
	res<-character()
	res[1]<-"%This needs package amsmath. Write \\usepackage{amsmath}"
	res[2]<-"\\begin{equation}"
	res[3]<- "\\begin{pmatrix} %explained vector"
	res[4]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
	res[5]<- "\\end{pmatrix}="
	res[6]<- "\\left\\{"
 	res[7]<-"\\begin{array}{rl}"
	Th<-x$model.specific$Thresh
	nthresh<-length(Th)
	###Condition for the threshold
	if(nthresh%in%c(1,2)){
		cond<-paste(c("& \\text{if Th}<","& \\text{if Th}>"), Th)}
	if(nthresh==2){
		cond[3]<-cond[2]
		cond[2]<-paste("& \\text{if }",Th[1], "< \\text{Th} <", Th[2])	
 		}
	###Adds the const/trend and lags
	for(i in 1:(nthresh+1)){
		if(x$model.specific$oneMatrix){
			regimei<-coef[[1]]
			j<-i}
		else{
			regimei<-coef[[i]]
			j<-1}
 		res<-include(x, res, regimei)
		ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
		res<-LagTeX(res, x, regimei, skip=ninc+x$lag*x$k*(j-1))
		res[length(res)+1]<- paste(cond[i], "\\\\")
	}
	res[length(res)+1]<-"\\end{array}"
	res[length(res)+1]<-"\\right."
	res[length(res)+1]<-"\\end{equation}"
	res<-gsub("kik", "\\\\", res, fixed=TRUE)
	res<-res[res!="blank"]
	
	return(structure(res, class="Latex"))

}

onesearch <- function(thDelay,gammas, fun, trace, gamma, plot){
	grid1 <- expand.grid(thDelay,gammas)				#grid with delay and gammas
	store<-mapply(fun,thDelay=grid1[,1],gam1=grid1[,2])
	posBestThresh <- which(store==min(store, na.rm=TRUE), arr.ind=TRUE)[1]

	if(plot&is.null(gamma)){
		result1 <- cbind(grid1,store)
		col <- rep(thDelay,length.out=nrow(result1))
		plot(result1[,2], result1[,3], col=col,xlab="Treshold Value",ylab="SSR", main="Results of the grid search")
		legend("topleft", pch=1, legend=paste("Threshold Delay", thDelay), col=thDelay)
	}
	if(trace)
		cat("Best unique threshold", grid1[posBestThresh,2],"\n")
	if(length(thDelay)>1&trace)
		cat("Best Delay", grid1[posBestThresh,1],"\n")
	list(bestThresh=grid1[posBestThresh,2],bestDelay=grid1[posBestThresh,1])
}#end of function one search

condiStep<-function(allgammas, threshRef, delayRef,ninter, fun, trace=TRUE, More=NULL){

	wh.thresh <- which.min(abs(allgammas-threshRef))
	
	#search for a second threshold smaller than the first
	if(wh.thresh>2*ninter){
		gammaMinus<-allgammas[seq(from=ninter, to=wh.thresh-ninter)]
# print(gammaMinus)
		storeMinus <- mapply(fun,gam1=gammaMinus,gam2=threshRef, thDelay=delayRef, MoreArgs=More)	
	}
	else
		storeMinus <- NA

	#search for a second threshold higher than the first
	if(wh.thresh<length(allgammas)-2*ninter){
		gammaPlus<-allgammas[seq(from=wh.thresh+ninter, to=length(allgammas)-ninter)]
		storePlus <- mapply(fun,gam1=threshRef,gam2=gammaPlus, thDelay=delayRef,MoreArgs=More)
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
	SSR<-min(store2, na.rm=TRUE)
	if(trace)
		cat("Second best: ",newThresh, " (conditionnal on ",threshRef, " ) \t SSR:", SSR, "\n")
	list(threshRef=threshRef, newThresh=newThresh, SSR=SSR)
}

grid<-function(gammasUp, gammasDown, bestDelay, fun, trace=TRUE, method=c("for", "apply", "mapply")){
	store <- matrix(NA,ncol=length(gammasUp), nrow=length(gammasDown))
	method<-match.arg(method)
	if(method=="for"){
		#Grid search
		for(i in seq_len(length(gammasDown))){
			gam1 <- gammasDown[i]
			for(j in 1: length(gammasUp)){
				gam2 <- gammasUp[j]
				store[i,j] <- fun(gam1=gam1, gam2=gam2, thDelay=bestDelay)
			}
		}

		#Finding the best result
		positionIter <- which(store==min(store, na.rm=TRUE), arr.ind=TRUE)
		rIter <- positionIter[1]
		cIter <- positionIter[2]
	
		gamma1Iter <- gammasDown[rIter]
		gamma2Iter <- gammasUp[cIter]

		bestThresh <- c(gamma1Iter, gamma2Iter)
	}
	else if(method=="apply"){
		grid<-expand.grid(gammasDown,gammasUp)
# 		fun(gam1=gam1, gam2=gam2, d=bestDelay)
		temp<-function(a) fun(gam1=c(a)[1],gam2=c(a)[2],thDelay=bestDelay)
		store<-apply(grid,1,temp)
		bests<-which(store==min(store, na.rm=TRUE))
		if(length(bests)>1) {
			warning("There were ",length(bests), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
			bests<-bests[1]}
# 		beta_grid<-grid[bests,1]
# 		bestGamma1<-grid[bests,2]
		bestThresh <- c(grid[bests,1], grid[bests,2])
	}
	else if(method=="mapply"){
		grid<-expand.grid(gammasDown,gammasUp)	
		store<-mapply(fun, gam1=grid[,1],gam2=grid[,2], MoreArgs=list(thDelay=bestDelay))
		bests<-which(store==min(store, na.rm=TRUE))
		if(length(bests)>1) {
			warning("There were ",length(bests), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
			bests<-bests[1]}
# 		beta_grid<-grid[bests,1]
# 		bestGamma1<-grid[bests,2]
		bestThresh <- c(grid[bests,1], grid[bests,2])
	}
	if(trace)
		cat("\nSecond step best thresholds", bestThresh, "\t\t SSR:", min(store, na.rm=TRUE), "\n")
	return(bestThresh)
}#end of grid function
#  VAR<-TVAR(dat[1:300,], lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply", max.iter=5)

#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="apply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="mapply"))
#  system.time(TVAR(dat, lag=2, nthresh=2,thDelay=1, plot=FALSE, commonInter=TRUE, include="const", trick="for"))