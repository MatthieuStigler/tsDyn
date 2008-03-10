TVECM_SeoTest<-function(data,lag, beta, trim=0.1,nboot, plot=FALSE, check=FALSE) {
y<-as.matrix(data)
T<-nrow(y)
k<-ncol(y)
p<-lag
if(is.null(colnames(data))==TRUE){colnames(data)<-paste("Var", c(1:k), sep="")}

DeltaY<-t(diff(y))[,(p+1):(T-1)] 		#Matrix Delta Y(t) of dim 2 x t
DeltaX<-rbind(1,t(embed(diff(y),p+1)[,-(1:k)]))	#trend and lags matrixs DeltaX(t-1,..) of dim: pk x t

Y<-DeltaY
M<-diag(1,T-p-1)-t(DeltaX)%*%solve(DeltaX%*%t(DeltaX))%*%DeltaX		#Projection matrix: all variables but ECT

##ECT
if(missing(beta)){
	warning("The article is constructed with given beta. Estimating beta can distort the test")
	coint<-lm(y[,1]~ -1 +y[,-1])
	beta<-c(1,-coint$coef)
	print(coint)}	 		#OLS estimation of beta
else 	{betas<-c(1,-beta)}			#Pre-specified beta

ECTfull<-y%*%betas
ECT<-ECTfull[-c(1:p,T)]				#ECT
ECT2<-embed(ECT,p)
Z<-rbind(t(ECT),DeltaX)

#########
###Grid
#########
allgammas<-sort(unique(ECT))
ng<-length(allgammas)
ninter<-round(trim*ng)		#minimal number of obs in each regime
inf<-ceiling(trim*ng+1)
sup<-floor((1-trim)*ng-1)
gammas<-allgammas[inf:sup]

store<-matrix(0, nrow=ng,ncol=ng)
store2<-matrix(NA, nrow=ng,ncol=ng)

###Function which gives Wald test and detSigma from model with two different thresholds	
loop<-function(gam1,gam2, ECT, DeltaX,Y,M){
	##Threshold dummies
	ECTminus <-ifelse(ECT<gam1,1,0)*ECT			
	ECTplus <-ifelse(ECT>gam2, 1,0)*ECT
	ThreshECT<-cbind(ECTplus,ECTminus) 

	##Total and partial parameter matrix from TVECM
	Z2<-rbind(t(ThreshECT),DeltaX)
	alpha2<-t(solve(crossprod(ThreshECT,M)%*%ThreshECT)%*%crossprod(ThreshECT,M)%*%t(DeltaY)) #Parameter of TECM
	res<-Y-tcrossprod(Y,Z2)%*%chol2inv(chol(tcrossprod(Z2)))%*%Z2 #All parameters
	Sigma<-(1/T)*tcrossprod(res,res)
	detSigma<-det(Sigma)

	###Wald Test
	Wald<-t(matrix(alpha2, ncol=1))%*%solve(solve(crossprod(ThreshECT,M)%*%ThreshECT)%x%Sigma)%*%matrix(alpha2, ncol=1)

	list(Wald=Wald, detSigma=detSigma)
}#end  loop

loop2<-function(gam1,gam2, ECT, DeltaX,Y,M){
	##Threshold dummies
	ECTminus <-ifelse(ECT<gam1,1,0)*ECT			
	ECTplus <-ifelse(ECT>gam2, 1,0)*ECT
	ThreshECT<-cbind(ECTplus,ECTminus) 

	##Total and partial parameter matrix from TVECM
	Z2<-rbind(t(ThreshECT),DeltaX)
	alpha2<-t(solve(crossprod(ThreshECT,M)%*%ThreshECT)%*%crossprod(ThreshECT,M)%*%t(DeltaY)) #Parameter of TECM
	res<-Y-tcrossprod(Y,Z2)%*%chol2inv(chol(tcrossprod(Z2)))%*%Z2 #All parameters
	###Wald Test
	t(matrix(alpha2, ncol=1))%*% chol2inv(chol(solve(crossprod(ThreshECT,M)%*%ThreshECT)%x%(1/T)*tcrossprod(res,res)))%*%matrix(alpha2, ncol=1)

}#end  loop

###Loop for values of the grid
for(i in 1:length(gammas)){
	gam1<-gammas[i]
	for (j in 1:length(gammas)){
		if(j>i+ninter){		
			gam2<-gammas[j]
			res<-loop(gam1=gam1,gam2=gam2,ECT=ECT, DeltaX=DeltaX,Y=Y,M=M)
			store[i,j]<-res$Wald
			store2[i,j]<-res$detSigma
		} #End if
	}	#End for j

}		#end for i


#################
###Sup Wald model
#################
supWald<-max(store,na.rm=TRUE)

position<-which(store==max(store, na.rm=TRUE), arr.ind=TRUE)
row<-position[1]
col<-position[2]

gamma1<-gammas[row]
gamma2<-gammas[col]

resultw<-loop(gamma1,gamma2,ECT=ECT, DeltaX=DeltaX,Y=Y,M=M)
cat("SupWald Model \n Wald", resultw$Wald, "Sigma",resultw$detSigma,"\n")

##################
###Bootstrap test
##################

###Initial parameters taken from model which minimise Sigma
minSigma<-min(store2,na.rm=TRUE)
r2<-which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)[1]
c2<-which(store2==min(store2, na.rm=TRUE), arr.ind=TRUE)[2]

gamma1Sigma<-gammas[r2]
gamma2Sigma<-gammas[c2]


ECTminussig <-ifelse(ECT<gamma1Sigma, 1,0)*ECT	
ECTplussig <-ifelse(ECT>gamma2Sigma, 1,0)*ECT
Zsig<-rbind(ECTminussig,ECTplussig,DeltaX)
Y<-DeltaY

Bsig<-tcrossprod(Y,Zsig)%*%chol2inv(chol(tcrossprod(Zsig)))
resSig<-t(Y-Bsig%*%Zsig)


resultSig<-loop(gamma1Sigma, gamma2Sigma,ECT=ECT, DeltaX=DeltaX,Y=Y,M=M)

cat("Min Sigma Model\n Wald", resultSig$Wald, "Sigma",resultSig$detSigma,"\n")

colnames(Bsig)<-c("ECTunder","ECTover","Trend",c(paste(rep(colnames(data),p), -rep(1:p, each=k))))
rownames(Bsig)<-colnames(data)

###Initial values
Yb<-matrix(0, nrow=nrow(y)-1, ncol=k)		#Delta Y term
Yb[1:(lag+1),]<-diff(y)[1:(lag+1),]			
Xminus1<-matrix(0,nrow=nrow(y), ncol=k)		#Xminus1 term
Xminus1[1:(lag+1),]<-y[(1):(lag+1),]
ECTtminus1<-matrix(0,nrow=nrow(y), ncol=1)		#ECT term


###Boostrap the residuals
bootstraploop<-function(vec_beta){


resb<-rbind(matrix(0,nrow=lag, ncol=k),apply(resSig,2,sample, replace=TRUE))
if(check)
	resb<-rbind(matrix(0,nrow=lag, ncol=k),resSig)		#uncomment this line to check the adequacy
for(i in (lag+2):(nrow(y)-1)){
	Xminus1[i,]<-Xminus1[i-1,]+Yb[i-1,]
	ECTtminus1[i]<-Xminus1[i,]%*%vec_beta
	Yb[i,]<-apply(cbind(Bsig[,3], Bsig[,-c(1:3)]%*%matrix(t(Yb[i-c(1:lag),]),ncol=1),resb[i,]),1,sum)
	if(ECTtminus1[i]<gamma1Sigma) {
		Yb[i,]<-apply(cbind(Bsig[,1]*ECTtminus1[i,],Yb[i,]),1,sum)}
	if(ECTtminus1[i]>gamma2Sigma) {
		Yb[i,]<-apply(cbind(Bsig[,2]*ECTtminus1[i,],Yb[i,]),1,sum)}
}

yboot<-apply(rbind(y[1,],Yb),2,cumsum)			#same as diffinv but it did not work

### Regression on the new series

ECTboot<-(y%*%c(1,-beta))[-c(1:p,T)]
DeltaYboot<-t(diff(yboot))[,(p+1):(T-1)] 		
DeltaXboot<-rbind(1,t(embed(diff(yboot),p+1)[,-(1:k)]))
Mboot<-diag(1,T-p-1)-t(DeltaXboot)%*%chol2inv(chol(tcrossprod(DeltaXboot)))%*%DeltaXboot

gammasb<-sort(unique(ECTboot))[inf:sup]
storeb<-matrix(0, nrow=ng,ncol=ng)

###Loop for values of the grid

for(i in 1:length(gammasb)){
	gam1<-gammasb[i]
	for (j in 1:length(gammasb)){
		if(j>i+ninter){		
			gam2<-gammasb[j]
			res<-loop(gam1,gam2,ECT=ECTboot, DeltaX=DeltaXboot,Y=DeltaYboot,M=Mboot)
			storeb[i,j]<-res$Wald
		} #End if
	}	#End for j
}		#end for i

supWaldboot<-max(storeb)
return(supWaldboot)
}#end of the bootstrap loop

Waldboots<-replicate(nboot,bootstraploop(c(1,-beta)))
PvalBoot<-mean(ifelse(Waldboots>supWald,1,0))
CriticalValBoot<-quantile(Waldboots, probs=c(0.9, 0.95, 0.975,0.99))

###Graphical output
if(plot==TRUE){
	plot(density(Waldboots))
	abline(v=c(supWald, CriticalValBoot[c(1,2)]), lty=c(1,2,3), col=c(2,3,4))
	legend("topright", legend=c("SupWald", "Boot 10%", "Boot 5%"), col=c(2,3,4), lty=c(1,2,3))
}


###Output
cat("\n Size f full sample", T)
cat("\n Size of end sample", ncol(DeltaY), "\n")
list(supWald=supWald,gamma1=gamma1, gamma2=gamma2, B=Bsig,PvalBoot=PvalBoot,CriticalValBoot=CriticalValBoot)
}

if(FALSE) {#usage example

data(zeroyld)
data<-zeroyld


TVECM_SeoTest(data[1:100,],lag=1, beta=1.1, trim=0.15, nboot=1, plot=FALSE, check=TRUE)
}
