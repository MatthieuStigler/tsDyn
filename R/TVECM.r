TVECM<-function(data,lag=1,trend=TRUE, bn=50, ngridG=50, trim=0.05, nthresh=1,plot=TRUE, dummyToBothRegimes=TRUE, methodMapply=FALSE, gamma1=list(exact=NULL, int=c("from","to"), around="val"),gamma2=list(exact=NULL, int=c("from","to"), around="val"), beta=list(exact=NULL, int=c("from","to"), around=c("val","by")), rest=c("none", "equal", "signOp"), model=c("All", "only_ECT") ) {
y<-as.matrix(data)
T<-nrow(y)		#T: number of observations
p<-lag 		#p: Number of lags
k<-ncol(y) 		#k: Number of equations
if(k>2 & is.null(beta$exact)) {stop("Sorry, the search is only possible with two variables. If more, please provide pre-specified beta values")}
if(is.null(colnames(data))==TRUE){colnames(data)<-paste("Var", c(1:k), sep="")}

model<-match.arg(model)
model<-switch(model, "All"="All", "only_ECT"="only_ECT")

ysmall<-y[(p+1):T,]
DeltaY<-diff(y)[(p+1):(T-1),]
Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
DeltaX<-embed(diff(y),p+1)[,-(1:k)]
if(trend==TRUE){DeltaX<-cbind(rep(1,T-p-1), DeltaX)}

##Long-run relationship OLS estimation
coint<-lm(y[,1]~ -1 +y[,2]) #Relation cointégrante estimée par OLS

betaLT<-coint$coef
betaLT_std <- sqrt(diag(summary(coint)$sigma*summary(coint)$cov))

ECT<-y%*%c(1,-betaLT)
#ECTminus1<-ECT[-c(1:p,T)]
ECTminus1<-Xminus1%*%c(1,-betaLT)

ECM<-residuals(coint)
ECMminus1<-ECM[-c(1:p,T)]

##Linear VECM estimation (Engle-Granger two step approach)
Z<-cbind(ECTminus1,DeltaX)		#Z: All regressors,ECT, trend and lags dim: t x npar
Y<-DeltaY				#Y: Delta Y, dim t x 2

B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
npar<-ncol(B)

rownames(B)<-paste("Equation",colnames(data))
if(trend==TRUE){colnames(B)<-c("ECM","Intercept",c(paste(rep(colnames(data),p), -rep(1:p, each=k))))}
else {colnames(B)<-c("ECM",c(paste(rep(colnames(data),p), -rep(1:p, each=k)))) }

res<-Y-Z%*%t(B)

Sigma<- matrix(crossprod(res), nrow=k,dimnames=list(colnames(data), colnames(data)))
nlike<-log(det(Sigma))/(T-p-1)		#	nlike=(t/2)*log(det(sige));
aic<-nlike+2*(npar)/(T-p-1)	
bic<-nlike+log(T-p-1)*(npar)/(T-p-1)	#bic #=nlike+log10(t)*4*(1+k); ###BIC

###Output for linear VECM
cat("Size of full sample:", T)
cat("\nSize of end sample:", T-p-1)
cat("\nNumber of lags:", p)
cat("\nNumber of variables:", k)
cat("\nVariables:", c(colnames(data)))

cat("\n############################ \n###Linear VECM estimates (by OLS)\n############################\n \n")
cat("Cointegrating vector",c(1, -betaLT), "\n")
cat("Standard error of the cointegrating value", betaLT_std, "\n")
cat("Parameters \n")
print(B)
cat("\nNegative LL \t", nlike, "\n")
cat("AIC \t\t", aic, "\n")
cat("BIC \t\t", bic, "\n")


#########################
###Set up of the grid
#########################

#Function to select values around a given point
aroundGrid <- function(around,allvalues,ngrid,trim){
	ng <- length(allvalues)
	wh.around <- which(allvalues==around)
	if(length(wh.around)==0)
		stop("\nSorry, the value you gave for the around argument did not match")
	if(length(wh.around)>1){
		warning("\nThere were", length(wh.around)," values corresponding to the around argument. The first one was taken")
		wh.around<-wh.around[1]}
	ar <- seq(from=wh.around-round(ngrid/2), to=(wh.around+round(ngrid/2)))		#Values around the point
	ar2 <- ar[ar>=round(trim*ng)&ar<=round((1-trim)*ng)]			#Bounding with trim 
	values <- allvalues[ar2]
	return(values)
}


###grid for gamma1
allgammas<-sort(unique(ECM))
ng<-length(allgammas)

#Default method: grid from lower to higher point
gammas<-allgammas[round(seq(from=trim, to=1-trim, length.out=ngridG)*ng)]
#gamma pre-specified
if(is.null(gamma1$exact)==FALSE){
	if(any(allgammas==gamma1$exact)==FALSE)
		warning("The value you gave for gamma does not correspond to an existing value. This causes problems currently")
	gammas<-gamma1$exact
	ngridG<-1
	}
#interval to search between given by user
if(is.numeric(gamma1$int)){
	intDown<-which(round(allgammas,3)==round(gamma1$int[1],3))
	intUp<-which(round(allgammas,3)==round(gamma1$int[2],3))
	gammas<-allgammas[seq(from=intDown, to=intUp, length.out=ngridG)]
	}
#value to search around	given by user
if(is.numeric(gamma1$around))		
	gammas<-aroundGrid(gamma$around,allvalues=allgammas,ngridG,trim)


###Grid for beta
#Default method: interval to search based on confidnce interval from linear model
betas<- seq(from=betaLT -2*betaLT_std, to=betaLT +2*betaLT_std, length.out=bn)
#beta pre-specified
if(is.null(beta$exact)==FALSE)			
	{betas<-beta$exact; bn<-1}
#interval to search between given by user
if(is.numeric(beta$int))			
	betas<-seq(from=beta$int[1], to=beta$int[2], length.out=bn)
#value to search around	given by user
if(is.numeric(beta$around)){
	by<-beta$around[2]
	betas<-seq(from=beta$around[1]-bn*by/2, to=beta$around[1]+bn*by/2, by=by)
	}

################
####One threshold model
################

oneSearch<-function(betas, gammas){

###Search function

#New operator if the dummy apply to both regimes. Is made outside the function so is evaluated only once
#dummy<-match.arg(dummy, "both regimes"="both regimes", "only one regime"="only one regime")
if(dummyToBothRegimes==TRUE){"%a%"<-function(matrix,dummy) matrix*dummy}
else {"%a%"<-function(matrix,dummy) matrix}

oneThresh<-function(betai, gam, Y, Xminus1,DeltaX){ 
	ECTi<-Xminus1%*%c(1,-betai)	#ECT in a column
	zi<-cbind(ECTi,DeltaX)		#All variables: ECT and lag, of dim t x kp+1+1
	d1<-ifelse(ECTi<=gam, 1,0)	#Dummy vector 		#d1=(w<=gam);
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	if(is.na(n1)==TRUE) n1<-0
	if (min(n1,1-n1)>trim) {
		zigamma<-c(d1)*zi
		zi<-zi%a%c(1-d1)		#new operator for choice between set up of first matrix
		Z<-cbind(zigamma,zi)
		LS<-crossprod(c(Y-tcrossprod(Z,crossprod(Y,Z)%*%solve(crossprod(Z)))))
		}
	else LS<-NA
	return(LS)
}#end function oneThresh


###Grid search


####Method with for
if(methodMapply==FALSE){

store<-matrix(NA,nrow=length(gammas), ncol=length(betas))	

for (i in seq_len(length(gammas))){
   gam<-gammas[i]
    for (j in seq_len(length(betas))){
	betai<-betas[j]
	store[i,j]<-oneThresh(betai=betai, gam=gam, DeltaX=DeltaX, Xminus1=Xminus1,Y=Y)
	}
}

#m<-min(store, na.rm=TRUE)
na<-sum(ifelse(is.na(store),1,0))
if(na>0) cat("\n",na, " points of the grid lead to regimes with percentage of observations < trim and were not computed\n")


pos<-which(store==min(store, na.rm=TRUE), arr.ind=TRUE)		#Best gamma
if(nrow(pos)>1) {
	warning("There were ",nrow(pos), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
	pos<-pos[,1]}

bestGamma1<-gammas[pos[1]]
beta_grid<-betas[pos[2]]

} #end methodMapply false
###Method with mapply

if(methodMapply==TRUE){
grid<-expand.grid(betas,gammas)
oneThreshTemp<-function(betai,gam) oneThresh(betai=betai, gam=gam, DeltaX=DeltaX,Xminus1=Xminus1, Y=Y)
storemap<-mapply(oneThreshTemp, betai=grid[,1], gam=grid[,2])
bests<-which(storemap==min(storemap, na.rm=TRUE))
if(length(bests)>1) {
	warning("There were ",length(bests), " thresholds values which minimize the SSR in the first search, the first one was taken") 	
	bests<-bests[1]}
beta_grid<-grid[bests,1]
bestGamma1<-grid[bests,2]
}



gammaMLE<-0.02321329
betaMLE<-0.8916303

#bestGamma1<-gammaMLE
#beta_grid<-betaMLE

###Plot results of grid search
if(is.null(gamma1$exact)==FALSE&is.null(beta$exact)==FALSE){plot<-FALSE}

if(plot==TRUE){
if(is.null(beta$exact)==FALSE&is.null(gamma1$exact)==TRUE){
	plot(gammas,store, type="l", xlab="Threshold parameter gamma", ylab="Residual Sum of Squares", main="Grid Search")}
if(is.null(beta$exact)==TRUE&is.null(gamma1$exact)==FALSE){
	plot(betas,store, type="l", xlab="Cointegrating parameter beta", ylab="Residual Sum of Squares", main="Grid Search")}
if(is.null(beta$exact)==TRUE&is.null(gamma1$exact)==TRUE){
	betaRSS<-apply(store,2,FUN=min, na.rm=TRUE)
	gammaRSS<-apply(store,1,FUN=min, na.rm=TRUE)
	layout(c(1,2))
	plot(gammas,gammaRSS, type="l", xlab="Threshold parameter gamma", ylab="Residual Sum of Squares", main="Grid Search")
	plot(betas,betaRSS, type="l", xlab="Cointegrating parameter beta", ylab="Residual Sum of Squares")
	abline(v=betaLT, lty=3)
	legend("topright", "OLS estimate from linear VECM", lty=3)}
}#end of the plot

#result of the whole function to search for one threshold
list(beta=beta_grid, gamma=bestGamma1)
}	#end of function oneSearch

if(nthresh==1){
results<-oneSearch(betas, gammas)
bestBeta<-results$beta
bestThresh<-results$gamma
}#end of if nthresh=1
############################
###Search for two thresholds
############################

###Two thresholds model: all variables change, three ECT

if(nthresh==2){


two_Thresh<-function(betai,gam1,gam2){
	ECTi<-Xminus1%*%c(1,-betai)	#ECT in a column
	zi<-cbind(ECTi,DeltaX)		#All variables: ECT and lag
	d1<-ifelse(ECTi<=gam1, 1,0)	#Dummy vector 
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	d2<-ifelse(ECTi>gam2,1,0)
	n2<-mean(d2)
	if(is.na(n1)==TRUE) {n1<-0; n2<-0}
	if (min(n1,n2,1-n1-n2)>trim) {
		ziUnder<-c(d1)*zi
		ziOver<-c(d2)*zi
		ziMiddle<-c(1-d1-d2)*zi
		Z<-cbind(ziUnder,ziMiddle,ziOver)
		LS<-crossprod(c(Y-tcrossprod(Z,crossprod(Y,Z)%*%solve(crossprod(Z)))))
		}
	else LS<-NA

	return(LS)
}

###Two thresholds model: only ECT change (and is zero in the middle regime)

two_partial_Thresh<-function(betai,gam1,gam2){
	ECTi<-Xminus1%*%c(1,-betai)	#ECT in a column
	d1<-ifelse(ECTi<=gam1, 1,0)	#Dummy vector 
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	d2<-ifelse(ECTi>gam2,1,0)
	n2<-mean(d2)
	if(is.na(n1)==TRUE) {n1<-0; n2<-0}
	if (min(n1,n2,1-n1-n2)>trim) {
		ectUnder<-c(d1)*ECTi
		ectOver<-c(d2)*ECTi
		Z<-cbind(ectUnder,ectOver,DeltaX)
		result<-crossprod(c(Y-tcrossprod(Z,crossprod(Y,Z)%*%solve(crossprod(Z)))))
		}
	else result<-NA

	return(result)
}#end two partial thresh

###Conditional search


###Conditionnal step
bestone <- oneSearch(betas, gammas)
bestThresh <- bestone$gamma
bestBeta <- bestone$beta
func<-switch(model, "All"=two_Thresh, "only_ECT"=two_partial_Thresh)

cat("\nBest threshold from first search", bestThresh)

if(!is.null(gamma2$exact))
	secondBestThresh<-gamma2$exact

if(is.null(gamma2$exact) & !is.numeric(gamma2$around)){


wh.thresh <- which(allgammas==bestThresh)
ninter<-round(trim*nrow(Xminus1))

#search for a second threshold smaller than the first
if(wh.thresh>2*ninter){
	gammaMinus<-allgammas[seq(from=ninter, to=wh.thresh-ninter)]
	storeMinus <- mapply(func,betai=bestBeta, gam1=gammaMinus,gam2=bestThresh)	
}
else	storeMinus <- NA

#search for a second threshold higher than the first
if(length(wh.thresh<length(allgammas)-2*ninter)){
	gammaPlus<-allgammas[seq(from=wh.thresh+ninter, to=length(allgammas)-ninter)]
	storePlus <- mapply(func,betai=bestBeta, gam1=bestThresh,gam2=gammaPlus)
}
else	storePlus <- NA

#results
store2 <- c(storeMinus, storePlus)
positionSecond <- which(store2==min(store2, na.rm=TRUE))
if(length(positionSecond)>1) warning("There were ",length(positionSecond), " thresholds values which minimize the SSR in the conditional step, the first one was taken") 
positionSecond<-positionSecond[1]
if(positionSecond<=length(storeMinus)){secondBestThresh<-gammaMinus[positionSecond]}
else {secondBestThresh<-gammaPlus[positionSecond-length(storeMinus)]}

cat("\nSecond best (conditionnal on the first one)", c(bestThresh,secondBestThresh), "\t SSR", min(store2, na.rm=TRUE))
}#end if ...many conditions


###Iterative step

bestThresh<-bestThresh
secondBestThresh<-secondBestThresh

if(is.numeric(gamma2$around))
	secondBestThresh<-gamma2$around

smallThresh <- min(bestThresh,secondBestThresh)
gammasDown <- aroundGrid(around=smallThresh,allgammas,ngrid=30, trim=trim)

bigThresh <- max(bestThresh,secondBestThresh)
gammasUp <- aroundGrid(around=bigThresh,allgammas,ngrid=30, trim=trim)

storeIter <- matrix(NA,ncol=length(gammasUp), nrow=length(gammasDown))

if(!is.null(gamma2$exact)){
	if(gamma2$exact<bestThresh)
		gammasDown<-gamma2$exact
	if(gamma2$exact>bestThresh)
		gammasUp<-gamma2$exact
}

if(!is.null(gamma1$exact)){
	if(gamma1$exact<secondBestThresh)
		gammasDown<-gamma2$exact
	if(gamma1$exact>secondBestThresh)
		gammasUp<-gamma2$exact
}

#Grid search
for(i in seq_len(length(gammasDown))){
	gam1 <- gammasDown[i]
	for(j in 1: length(gammasUp)){
		gam2 <- gammasUp[j]
		storeIter[i,j] <- func(gam1=gam1, gam2=gam2, beta=bestBeta)
	}
}

#Finding the best result
positionIter <- which(storeIter==min(storeIter, na.rm=TRUE), arr.ind=TRUE)
if(nrow(positionIter)>1) 
	{warning("There were ",length(positionIter), " thresholds values which minimize the SSR in the iterative step, the first one was taken") 
	positionIter<-positionIter[1,]}
rIter <- positionIter[1]
cIter <- positionIter[2]

bestThresh1Iter <- gammasDown[rIter]
bestThresh2Iter <- gammasUp[cIter]

bestThresh <- c(bestThresh1Iter, bestThresh2Iter)

cat("\nSecond step best thresholds", bestThresh, "\t\t\t SSR", min(storeIter, na.rm=TRUE), "\n")
}#end if nthresh=2



############################
###Details of the Best model
############################

if(nthresh==1){
	Xminus1<-Xminus1
	ECT_best<-Xminus1%*%c(1,-bestBeta)	#ECT
	Z_temp<-cbind(ECT_best,DeltaX)		#All variables: ECT and lag
	d1<-ifelse(ECT_best<=bestThresh, 1,0)	#Dummy vector 		#d1=(w<=gam);
	ndown<-mean(d1)				#Number of elements of the ECT under the threshold
	nup<-1-ndown
	Zunder<-c(d1)*Z_temp
	if(dummyToBothRegimes==TRUE){Zover<-c(1-d1)*Z_temp}
	else Zover<-Z_temp
	Zbest<-cbind(Zunder, Zover)
}

if(nthresh==2){
	Xminus1<-Xminus1
	ECT_best<-Xminus1%*%c(1,-bestBeta)	#ECT
	Z_temp<-cbind(ECT_best,DeltaX)		#All variables: ECT and lag
	d1<-ifelse(ECT_best<=bestThresh1Iter, 1,0)	#Dummy vector 		#d1=(w<=gam);
	ndown<-mean(d1)				#Number of elements of the ECT under the threshold
	d2<-ifelse(ECT_best>bestThresh2Iter,1,0)
	nup<-mean(d2)
	if(model=="All"){
		Zunder<-c(d1)*Z_temp
		Zover<-c(d2)*Z_temp
		Zmiddle<-(1-c(d1)-c(d2))*Z_temp
		Zbest<-cbind(Zunder,Zmiddle, Zover)
		}
	if(model=="only_ECT"){
		Zunder<-Zunder<-c(d1)*ECT_best
		Zover<-c(d2)*ECT_best
		Zbest<-cbind(Zunder,Zover, DeltaX)
		}
}


###Parameters, SSR AIC, BIC...
Bbest<-t(Y)%*%Zbest%*%solve(t(Zbest)%*%Zbest)	
npar<-ncol(Z)
nparTot<-npar+1+nthresh			#addition of threshold and cointegrating vector


resbest <- Y - Zbest%*% t(Bbest)
SSRbest <- as.numeric(crossprod(c(resbest)))
Sigma<- matrix(1/T*crossprod(resbest), nrow=k)
nlike_thresh<-log(det(Sigma))/(T-p-1)		#	nlike=(t/2)*log(det(sige));
aic_thresh<-nlike_thresh+2*(nparTot)/(T-p-1)
bic_thresh<-nlike_thresh+log(T-p-1)*(nparTot)/(T-p-1)		#bic #=nlike+log10(t)*4*(1+k); ###BIC




###naming the parameter matrix
rownames(Bbest) <- paste("Equation", colnames(data))


Bcolnames<-c("ECT", paste(rep(colnames(data),p), "t",-rep(1:p, each=k)))
if(trend) Bcolnames <- c("Trend", Bcolnames)

if(nthresh==1)
	colnames(Bbest) <- rep(Bcolnames,2)
else{
	if(model=="All")
		colnames(Bbest)<-rep(Bcolnames,3)
	if(model=="only_ECT"){
		print(dim(Bbest))
		print(length(Bcolnames))
		colnames(Bbest)<-c("ECT-","ECT+", Bcolnames[-2])
	}
}

###partitionning the matrix following the regimes

if(nthresh==1){
	Bdown <- Bbest[,c(1:npar)]
	Bup <- Bbest[,-c(1:npar)]
	Blist <- list(Bdown=Bdown, Bup=Bup)
	nobs <- c(ndown=ndown, nup=1-ndown)
} else {
	if(model=="All"){
		Bdown <- Bbest[,c(1:npar)]
		Bmiddle <- Bbest[,c(1:npar)+npar]
		Bup <- Bbest[,c(1:npar)+2*npar]
		colnames(Bmiddle) <- Bcolnames
		Blist <- list(Bdown=Bdown, Bmiddle=Bmiddle,Bup=Bup)
		nobs <- c(ndown=ndown, nmiddle=1-nup-ndown,nup=nup)
	}
	if(model=="only_ECT"){
		Blist<-Bbest
		Bdown<-Bbest
		Bup<-NA
		nobs <- c(ndown=ndown, nmiddle=1-nup-ndown,nup=nup)
	}
}

##Print
cat("\n############################ \n###Threshold VECM estimates \n############################\n \n")
if(is.null(gamma1$exact)==TRUE){cat("Threshold estimate \t\t", bestThresh)}
else {cat("User specified threshold \t\t", bestThresh)}
if(is.null(beta$exact)==TRUE){cat("\nCointegrating vector Estimate: \t", c(1,-bestBeta))}
else{cat("\nUser specified Cointegrating vector: \t", c(1,-bestBeta))}
cat("\nNegative Log-Like: \t\t", nlike_thresh)
cat("\nAIC \t\t\t\t", aic_thresh,"\nBIC \t\t\t\t", bic_thresh,"\nSSR \t\t\t\t", SSRbest)

cat("\n###Under regime \n \n Percentage of observations \t",  ndown)
cat("\n Parameters \n")
print(Bdown)

if(nthresh==2&model=="All"){
	cat("\n##Middle regime \n \n Percentage of observations \t",  1-ndown-nup)
	cat("\n Parameters \n")
	print(Bmiddle)
	}

if(!model=="only_ECT"){
	cat("\n###Over regime \n \n Percentage of observations \t",nup)
	cat("\n Parameters \n")
	print(Bup)
	}
}




if(FALSE) {
data(zeroyld)
data<-zeroyld[,c(36,19)]
colnames(data)<-c("short", "long")

TVECM(data, nthresh=2,lag=1, bn=20, ngridG=30, plot=TRUE,trim=0.05, model="All", trend=FALSE)
}
