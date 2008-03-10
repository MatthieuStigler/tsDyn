TVECM<-function(data,lag=1,trend=TRUE, bn=50, ngridG=50, trim=0.05, nthresh=1,plot=TRUE, dummyToBothRegimes=TRUE, methodMapply=FALSE, gamma1=list(exact=NULL, int=c("from","to"), around="val"),gamma2=list(exact=NULL, int=c("from","to"), around="val"), beta=list(exact=NULL, int=c("from","to"), around=c("val","by")), restr=c("none", "equal", "signOp"), model=c("All", "only_ECT"), demean = c( "const", "trend","none", "both"),beta0=0 ) {
y<-as.matrix(data)
T<-nrow(y)		#T: number of observations
p<-lag 		#p: Number of lags
t <- T-p-1 		#Size of end sample
k<-ncol(y) 		#k: Number of equations
if(k>2 & is.null(beta$exact)) 
	stop("Sorry, the search is only possible with 2 variables. If more, please provide pre-specified beta values")
if(is.null(colnames(data))==TRUE)
	colnames(data)<-paste("Var", c(1:k), sep="")
ndig<-getndp(y)
restr<-match.arg(restr)
demean<-match.arg(demean)
model<-match.arg(model)
model<-switch(model, "All"="All", "only_ECT"="only_ECT")

ysmall<-y[(p+1):T,]
DeltaY<-diff(y)[(p+1):(T-1),]
Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
DeltaX<-embed(diff(y),p+1)[,-(1:k)]

if(demean=="const")
	DeltaX<-cbind(rep(1,t), DeltaX)
else if(demean=="trend")
	DeltaX<-cbind(seq_len(t), DeltaX)
else if(demean=="both")
	DeltaX<-cbind(rep(1,t),seq_len(t), DeltaX)


##Long-run relationship OLS estimation
beta0<-as.matrix(beta0)

if(beta0[1]!=0){
	if(nrow(beta0)!=nrow(y))
		stop("Length of beta0 should be ", nrow(y), "\n")
	coint<-lm(y[,1]~ y[,2]+beta0-1) #Relation cointÃ©grante estimÃ©e par OLS
	beta0<-(beta0%*%coint$coef[-1])[-c(1:p,T),]}
else{
	coint<-lm(y[,1]~ y[,2]-1) #Relation cointÃ©grante estimÃ©e par OLS
	beta0<-rep(0, t)}

betaLT<-coint$coef[1]
betaLT_std <- sqrt(diag(summary(coint)$sigma*summary(coint)$cov))[1]


ECT<-y%*%c(1,-betaLT)
ECT<-round(ECT,ndig)
#ECTminus1<-ECT[-c(1:p,T)]

ECTminus1<-round(Xminus1%*%c(1,-betaLT),ndig)



# ECM<-residuals(coint)
# ECMminus1<-ECM[-c(1:p,T)]

##Linear VECM estimation (Engle-Granger second step approach)
Z<-cbind(ECTminus1-beta0,DeltaX)		#Z: All regressors,ECT, trend and lags dim: t x npar
Y<-DeltaY				#Y: Delta Y, dim t x 2

B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
npar<-ncol(B)
allpar<-ncol(B)*nrow(B)


rownames(B)<-paste("Equation",colnames(data))
LagNames<-c(paste(rep(colnames(data),p), -rep(seq_len(p), each=k)))
if(demean=="const")
	colnames(B)<-c("ECT","Intercept",LagNames)
else if(demean=="trend")
	colnames(B)<-c("ECT","Trend",LagNames)
else if(demean=="both")
	colnames(B)<-c("ECT","Intercept","Trend",LagNames)
else 
	colnames(B)<-c("ECT",LagNames)
res<-Y-Z%*%t(B)

Sigma<- matrix(1/t*crossprod(res),ncol=k,dimnames=list(colnames(data), colnames(data)))
VarCov<-solve(crossprod(Z))%x%Sigma
StDev<-matrix(diag(VarCov)^0.5, nrow=k)

Tvalue<-B/StDev
Pval<-pt(abs(Tvalue), df=(t-ncol(Z)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(t-ncol(Z)), lower.tail=TRUE)
Pval<-round(Pval,4)
colnames(Pval)<-colnames(B)


nlike<-log(det(Sigma))		#	nlike=(t/2)*log(det(sige));
aic<-t*nlike+2*(allpar)	
bic<-t*nlike+log(t)*(allpar)	#bic #=nlike+log10(t)*4*(1+k); ###BIC

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
cat("P-Values\n")
print(Pval)
cat("\nNegative LL \t", nlike, "\n")
cat("AIC \t\t", aic, "\n")
cat("BIC \t\t", bic, "\n")


#########################
###Set up of the grid
#########################

#Function to select values around a given point
aroundGrid <- function(around,allvalues,ngrid,trim){
	ng <- length(allvalues)
	wh.around <- which.min(abs(allvalues-around))
	if(length(which(allvalues==around))==0)
		warning("\nThe value ", around, " did not match to existing ones", allvalues[wh.around], "was taken instead")
	if(length(wh.around)>1){
		warning("\nThere were", length(wh.around)," values corresponding to the around argument. The first one was taken")
		wh.around<-wh.around[1]}
	ar <- seq(from=wh.around-round(ngrid/2), to=(wh.around+round(ngrid/2)))		#Values around the point
	ar2 <- ar[ar>=round(trim*ng)&ar<=round((1-trim)*ng)]			#Bounding with trim 
	values <- allvalues[ar2]
	return(values)
}


###grid for gamma1
allgammas<-sort(unique(ECTminus1-beta0))
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
	intDown<-which.min(abs(allgammas-gamma1$int[1]))
	intUp<-which.min(abs(allgammas-gamma1$int[2]))
	gammas<-allgammas[seq(from=intDown, to=intUp, length.out=min(ngridG,intUp-intDown))]
	}
#value to search around	given by user
if(is.numeric(gamma1$around))		
	gammas<-aroundGrid(gamma$around,allvalues=allgammas,ngridG,trim)
gammas<-round(gammas, ndig)

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
	ECTi<-Xminus1%*%c(1,-betai)-beta0	#ECT in a column
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

###threshold effect only in ECT
one_partial_Thresh<-function(betai, gam, Y, Xminus1,DeltaX){ 
	ECTi<-Xminus1%*%c(1,-betai)-beta0	#ECT in a column
	d1<-ifelse(ECTi<=gam, 1,0)	#Dummy vector 		#d1=(w<=gam);
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	if (min(n1,1-n1)>trim) {
		Z<-cbind(ECTi*d1, ECTi*(1-d1),DeltaX)
		LS<-crossprod(c(Y-tcrossprod(Z,crossprod(Y,Z)%*%solve(crossprod(Z)))))
		}
	else LS<-NA
	return(LS)
}#end function oneThresh

###Grid search


####Method with for
func_onethresh<-switch(model, "All"=oneThresh, "only_ECT"=one_partial_Thresh)
if(methodMapply==FALSE){

store<-matrix(NA,nrow=length(gammas), ncol=length(betas), dimnames=list(round(gammas,3), betas))	

for (i in seq_len(length(gammas))){
   gam<-gammas[i]
    for (j in seq_len(length(betas))){
	betai<-betas[j]
	store[i,j]<-func_onethresh(betai=betai, gam=gam, DeltaX=DeltaX, Xminus1=Xminus1,Y=Y)
	}
}

#m<-min(store, na.rm=TRUE)
na<-sum(ifelse(is.na(store),1,0))
if(na>0) cat("\n",na,"(", na/(nrow(store)*ncol(store)), ") points of the grid lead to regimes with percentage of observations < trim and were not computed\n")


pos<-which(store==min(store, na.rm=TRUE), arr.ind=TRUE)		#Best gamma
if(nrow(pos)>1) {
	warning("There were ",nrow(pos), " thresholds/cointegrating combinations (",paste(gammas[pos[,1]],"/",betas[pos[,2]],", "), ") \nwhich minimize the SSR in the first search, the first one  ", round(gammas[pos[1,1]],ndig), " ",round(betas[pos[1,2]],ndig)," was taken") 	
	pos<-pos[1,]}


bestGamma1<-gammas[pos[1]]
beta_grid<-betas[pos[2]]

} #end methodMapply false
###Method with mapply

if(methodMapply==TRUE){
grid<-expand.grid(betas,gammas)
oneThreshTemp<-function(betai,gam) func_onethresh(betai=betai, gam=gam, DeltaX=DeltaX,Xminus1=Xminus1, Y=Y)
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
	#mat[!is.na(apply(mat,1,sum)),]
	options(warn=-1)
	betaRSS<-apply(store,2,FUN=min, na.rm=TRUE)
	gammaRSS<-apply(store,1,FUN=min, na.rm=TRUE)
	options(warn=0)
	gammaRSS[is.infinite(gammaRSS)]<-NA
	betaRSS[is.infinite(betaRSS)]<-NA
	layout(c(1,2))
	plot(gammas,gammaRSS, type="l", xlab="Threshold parameter gamma", ylab="Residual Sum of Squares", main="Grid Search")
	points(x=bestGamma1, y=min(store, na.rm=TRUE), col=2, cex=2)
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
	ECTi<-Xminus1%*%c(1,-betai)-beta0	#ECT in a column
	zi<-cbind(ECTi,DeltaX)		#All variables: ECT and lag
	d1<-ifelse(ECTi<=gam1, 1,0)	#Dummy vector 
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	d2<-ifelse(ECTi>gam2,1,0)
	n2<-mean(d2)
	if(is.na(n1)==TRUE) {n1<-0; n2<-0}
# print(c(n1, 1-n1-n2,n2))
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
	ECTi<-Xminus1%*%c(1,-betai)-beta0	#ECT in a column
	d1<-ifelse(ECTi<=gam1, 1,0)	#Dummy vector 
	n1<-mean(d1)			#Number of elements of the ECT under the threshold
	d2<-ifelse(ECTi>gam2,1,0)
	n2<-mean(d2)
# print(c(n1, 1-n1-n2,n2))
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


wh.thresh <- which.min(abs(allgammas-bestThresh))
ninter<-round(trim*nrow(Xminus1))
print(restr)
if(restr=="none"){
	#search for a second threshold smaller than the first
	if(wh.thresh>2*ninter){
		gammaMinus<-allgammas[seq(from=ninter, to=wh.thresh-ninter)]
		storeMinus <- mapply(func,betai=bestBeta, gam1=gammaMinus,gam2=bestThresh)	
	}
	else	storeMinus <- NA

	#search for a second threshold higher than the first
	if(length(wh.thresh<length(allgammas)-2*ninter)){
		gammaPlus<-allgammas[seq(from=wh.thresh+ninter, to=length(allgammas)-ninter)]
		storePlus <- mapply(func,gam2=gammaPlus, MoreArgs=list(betai=bestBeta,gam1=bestThresh) )
	}
	else	storePlus <- NA
}

else if(restr=="signOp"){
	zero<-which.min(abs(allgammas))
	if(sign(bestThresh)>0){
		gammaMinus<-allgammas[seq(from=ninter, to=min(wh.thresh-ninter,zero))]
		storeMinus <- mapply(func,betai=bestBeta, gam1=gammaMinus,gam2=bestThresh)
		storePlus <- NA}
	else{	
		gammaPlus<-allgammas[seq(from=max(wh.thresh+ninter,zero), to=length(allgammas)-ninter)]
		storePlus <- mapply(func,betai=bestBeta, gam1=bestThresh,gam2=gammaPlus)
		storeMinus<-NA}
}
#results
store2 <- c(storeMinus, storePlus)

positionSecond <- which(store2==min(store2, na.rm=TRUE))
if(length(positionSecond)>1) warning("There were ",length(positionSecond), " thresholds values which minimize the SSR in the conditional step, the first one was taken") 

positionSecond<-positionSecond[1]
if(positionSecond<=length(storeMinus)){
	secondBestThresh<-gammaMinus[positionSecond]}
else {
	secondBestThresh<-gammaPlus[positionSecond-length(storeMinus)]}

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
	ECT_best<-Xminus1%*%c(1,-bestBeta)-beta0	#ECT
	Z_temp<-cbind(ECT_best,DeltaX)		#All variables: ECT and lag
	d1<-ifelse(ECT_best<=bestThresh, 1,0)	#Dummy vector
	ndown<-mean(d1)				#Number of elements of the ECT under the threshold
	nup<-1-ndown
	if(model=="All"){
		Zunder<-c(d1)*Z_temp
		if(dummyToBothRegimes==TRUE){Zover<-c(1-d1)*Z_temp}
		else Zover<-Z_temp
		Zbest<-cbind(Zunder, Zover)}
	else{
		Zbest<-cbind(d1*ECT_best, (1-d1)*ECT_best, DeltaX)}
}

if(nthresh==2){
	ECT_best<-Xminus1%*%c(1,-bestBeta)-beta0	#ECT
	d1<-ifelse(ECT_best<=bestThresh1Iter, 1,0)	#Dummy vector 		#d1=(w<=gam);
	ndown<-mean(d1)				#Number of elements of the ECT under the threshold
	d2<-ifelse(ECT_best>bestThresh2Iter,1,0)
	nup<-mean(d2)
	if(model=="All"){
		Z_temp<-cbind(ECT_best,DeltaX)		#All variables: ECT and lag
		Zunder<-c(d1)*Z_temp
		Zover<-c(d2)*Z_temp
		Zmiddle<-(1-c(d1)-c(d2))*Z_temp
		Zbest<-cbind(Zunder,Zmiddle, Zover)
		}
	if(model=="only_ECT"){
		Zunder<-c(d1)*ECT_best
		Zover<-c(d2)*ECT_best
		Zbest<-cbind(Zunder,Zover, DeltaX)
		}
}


###Parameters, SSR AIC, BIC...
Bbest<-t(Y)%*%Zbest%*%solve(t(Zbest)%*%Zbest)	
# npar<-ncol(Zbest)
allpar<-ncol(Bbest)*nrow(Bbest)

# nparTot<-npar+1+nthresh			#addition of threshold and cointegrating vector


resbest <- Y - Zbest%*% t(Bbest)
SSRbest <- as.numeric(crossprod(c(resbest)))
Sigmathresh<- matrix(1/t*crossprod(resbest), ncol=k)
nlike_thresh<-log(det(Sigmathresh))		#	nlike=(t/2)*log(det(sige));
aic_thresh<-t*nlike_thresh+2*(allpar+1+nthresh)
bic_thresh<-t*nlike_thresh+log(T-p-1)*(allpar+1+nthresh)		#bic #=nlike+log10(t)*4*(1+k); ###BIC

VarCovB<-solve(crossprod(Zbest))%x%Sigmathresh
StDevB<-matrix(diag(VarCovB)^0.5, nrow=k)

Tvalue<-Bbest/StDevB
Pval<-pt(abs(Tvalue), df=(t-ncol(Zbest)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(t-ncol(Zbest)), lower.tail=TRUE)
Pval<-round(Pval,4)


###naming the parameter matrix
rownames(Bbest) <- paste("Equation", colnames(data))
rownames(Pval) <- paste("Equation", colnames(data))


DeltaXnames<-c(paste(rep(colnames(data),p), "t",-rep(1:p, each=k)))


if(demean=="const") 
	Bcolnames <- c("ECT","Const", DeltaXnames)
else if(demean=="trend") 
	Bcolnames <- c("ECT","Trend", DeltaXnames)
else if(demean=="both") 
	Bcolnames <- c("ECT","Const","Trend", DeltaXnames)
else if(demean=="none") 
	Bcolnames <- c("ECT",DeltaXnames)

if(nthresh==1){
	if(model=="All")
		colnames(Bbest) <- rep(Bcolnames,2)
	if(model=="only_ECT")
		colnames(Bbest)<-c("ECT-","ECT+", Bcolnames[-which(Bcolnames=="ECT")])
}
else{
	if(model=="All")
		colnames(Bbest)<-rep(Bcolnames,3)
	if(model=="only_ECT")
		colnames(Bbest)<-c("ECT-","ECT+", Bcolnames[-which(Bcolnames=="ECT")])
}
colnames(Pval)<-colnames(Bbest)
###partitionning the matrix following the regimes

if(nthresh==1){
	if(model=="All"){
		Bdown <- Bbest[,c(1:npar)]
		Bup <- Bbest[,-c(1:npar)]
		Pdown <- Pval[,c(1:npar)]
		Pup <- Pval[,-c(1:npar)]
		Blist <- list(Bdown=Bdown, Bup=Bup)
		Plist <- list(Pdown=Pdown, Pup=Pup)
		nobs <- c(ndown=ndown,nup=nup)
	}
	if(model=="only_ECT"){
		Bdown<-Bbest
		Blist<-Bbest
		Plist<-Pval
		Pdown<-Pval
 		Bup<-NA
		nobs <- c(ndown=ndown, nup=nup)
	}
}
else if(nthresh==2){
	if(model=="All"){
		Bdown <- Bbest[,c(1:npar)]
		Bmiddle <- Bbest[,c(1:npar)+npar]
		Bup <- Bbest[,c(1:npar)+2*npar]
		colnames(Bmiddle) <- Bcolnames
		Blist <- list(Bdown=Bdown, Bmiddle=Bmiddle,Bup=Bup)
		Pdown <- Pval[,c(1:npar)]
		Pmiddle <- Pval[,c(1:npar)+npar]
		Pup <- Pval[,c(1:npar)+2*npar]		
		Plist <- list(Pdown=Pdown, Pmiddle=Pmiddle,Pup=Pup)	
		colnames(Pmiddle) <- Bcolnames
		nobs <- c(ndown=ndown, nmiddle=1-nup-ndown,nup=nup)
	}
	if(model=="only_ECT"){
		Blist<-Bbest
		Plist<-Pval
		Bdown<-Bbest
		Pdown<-Pval
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
cat("\nAIC \t\t\t\t", aic_thresh, "\nBIC \t\t\t\t", bic_thresh, "\nSSR \t\t\t\t", SSRbest)

cat("\n###Under regime \n \n Percentage of observations \t",  ndown)
cat("\n Parameters \n")
print(Bdown)
cat("\n P-values\n")
print(Pdown)

if(nthresh==2&model=="All"){
	cat("\n##Middle regime \n \n Percentage of observations \t",  1-ndown-nup)
	cat("\n Parameters \n")
	print(Bmiddle)
	cat("\n P-values\n")
	print(Pmiddle)
	}

if(!model=="only_ECT"){
	cat("\n###Over regime \n \n Percentage of observations \t",nup)
	cat("\n Parameters \n")
	print(Bup)
	cat("\n P-values\n")
	print(Pup)
	}
list(resids=list(res=res, resbest=resbest), VAR=VarCovB, Parameters=Blist, nobs_regime=nobs)
}



if(FALSE) {
data(zeroyld)
data<-zeroyld

environment(TVECM)<-environment(star)

summary(lm(zeroyld[,1]~zeroyld[,2]-1))
summary(lm(zeroyld[,1]~zeroyld[,2]))

TVECM(data, nthresh=1,lag=1, bn=80, ngridG=300, plot=TRUE,trim=0.05, model="All", beta=list(int=c(0.7,1.2)))
beta0<-rep(1.12,480)
TVECM(data, nthresh=1,lag=1, bn=20, ngridG=20, plot=FALSE,trim=0.05, model="only_ECT", beta0=beta0)
}

