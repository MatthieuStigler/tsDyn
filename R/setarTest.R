setarTest <- function (x, m, d = 1, steps = d, series, thDelay = 0:1, mL, mH,
    mTh, thVar, nboot=10, plot=FALSE, trim=0.1, test=c("1vs", "2vs3")) {
    if (missing(m)) m <- max(mL, mH, max(thDelay) + 1)
    if (missing(series))  series <- deparse(substitute(x))
    str <- nlar.struct(x = x, m = m, d = d, steps = steps, series = series)
    xx <- str$xx
    yy <- str$yy
	n<-nrow(xx)
	if(!missing(mL) |!missing(mH)) warning("The function may have some errors if you specify mL or mH")
    if (missing(mL)) { mL <- m    }
    if (missing(mH)) { mH <- m    }
	

###Linear model
xxlin<-cbind(1,xx)
linear <- lm.fit(xxlin, yy)
SSR<-as.numeric(crossprod(linear$residuals))


###Threshold transition variable

    if (!missing(thDelay)) {  
      if (max(thDelay) >= m)             stop(paste("thDelay too high: should be < m (=",  m, ")"))
        z <- xx[, seq_len(max(thDelay)+1)]
    }
    else if (!missing(mTh)) { stop("Currenly not implemented")
	    if (length(mTh) != m) stop("length of 'mTh' should be equal to 'm'")
        z <- xx %*% mTh
        dim(z) <- NULL
    }
    else if (!missing(thVar)) {
        if (length(thVar) > nrow(x)) {thVar <- thVar[seq_len(nrow(x))]}
	if (length(thVar) < nrow(x)) {stop("The threshold variable should not be smaller than the serie") }
	z<-nlar.struct(x = thVar, m = m, d = d, steps = steps, series = deparse(substitute(z)))[, thDelay + 1]
    }
    else {z <- xx[, thDelay]
    }


z<-as.matrix(z)
if(length(thDelay)==1) b<-thDelay
else b<-1
allgammas<-sort(unique(z[,b]))
ng<-length(allgammas)
ninter<-round(trim*ng)
gammas<-allgammas[ceiling(trim*ng+1):floor((1-trim)*ng)]


### Threshold model

xxl <- cbind(1, xx[, seq_len(mL)])
xxh <- cbind(1, xx[, seq_len(mH)])


###################
###Search function
##################
SSRestim<-function(parameters,yy, xxl,xxh,z) {#
	thDelay<-parameters[1]
	gammai<-parameters[2]
        isL <- ifelse(z[, thDelay + 1]<= gammai,1,0)	### isL: dummy 
        xxthresh <- cbind(xxl * isL,xxh * (1 - isL))	### Lower matrix
	crossprod(yy - xxthresh %*% chol2inv(chol(crossprod(xxthresh)))%*%crossprod(xxthresh,yy))
        }

SSRestim2 <- function(gam1,gam2,thDelay, yy, xx,z){
	##Threshold dummies
	dummydown <- ifelse(z[, thDelay + 1]<=gam1, 1, 0)
	regimedown <- dummydown*xx
	ndown <- mean(dummydown)
	dummyup <- ifelse(z[, thDelay + 1]>gam2, 1, 0)
	regimeup <- dummyup*xx
	nup <- mean(dummyup)
	##SSR from TAR(3)
	if(min(nup, ndown, 1-nup-ndown)>trim){
		XX <- cbind(regimedown, (1-dummydown-dummyup)*xx, regimeup)		# dim k(p+1) x t	
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSR
	}
	else
		res <- NA
# 	print(c(ndown,1-nup-ndown,nup, gam1, gam2,res))
	return(res)
}

##################
###One thresholds
##################

IDS<-as.matrix(expand.grid(thDelay, gammas))
result <- apply(IDS, 1, SSRestim,yy=yy, xxl=xxl,xxh=xxh,z=z)#
bestDelay<-IDS[which.min(result),1]
bestThresh<-IDS[which.min(result),2]
cat("Best unique threshold", bestThresh, "\t\t\t\t SSR", min(result), "\n")

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
	storeMinus <- mapply(fun,gam1=gammaMinus,gam2=threshRef,MoreArgs=MoreArgs)	
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
More<-list(thDelay=bestDelay, yy=yy, xx=xxlin,z=z)
Thresh2<-condiStep(allgammas, bestThresh,ninter=ninter, fun=SSRestim2, MoreArgs=More)$newThresh
Thresh3<-condiStep(allgammas, Thresh2,ninter=ninter, fun=SSRestim2, MoreArgs=More)

SSR1thresh<-min(result)
SSR2thresh<-Thresh3$SSR
SSRs<-matrix(c(SSR, SSR1thresh, SSR2thresh), ncol=3, dimnames=list("SSR", c("AR", "TAR(1)", "TAR(2)")))

###F test for original data
Ftest12<-as.numeric(n*(SSR-SSR1thresh)/SSR1thresh)
Ftest13<-as.numeric(n*(SSR-SSR2thresh)/SSR2thresh)
Ftest23<-as.numeric(n*(SSR1thresh-SSR2thresh)/SSR2thresh)

##############################
###Bootstrap for the F test
##############################

### Reconstruction of series from linear model for Ftest 1vs2 and 1vs3
#initial data	
xboot<-vector("numeric", length=length(x))		#Delta Y term
xboot[1:m]<-x[1:m]
reslin<-linear$residuals
B<-linear$coefficients
if(!missing(thVar)) zext<-embed(thVar,m+1)[round(trim*ng):round((1-trim)*ng),]

#loop
bootlinear<-function(x){
resib1<-c(rep(0,m),sample(reslin, replace=TRUE))	#residual sampling, 
#resib1<-c(rep(0,m),linear$residuals)				#uncomment this line to verify the bootstrap
for(i in (m+1):length(x)){
	xboot[i]<-sum(B[1],B[-1]*xboot[i-c(1:m)],resib1[i])
	}
return(xboot)
}
### Reconstruction of series from 1 thresh model for Ftest 2vs3
#initial data
isL <- ifelse(z[, bestDelay + 1]<= bestThresh,1,0)	### isL: dummy 
xxthresh <- cbind(xxl * isL,xxh * (1 - isL))	### Lower matrix

B1thresh<-chol2inv(chol(crossprod(xxthresh)))%*%crossprod(xxthresh,yy)
res1thresh<-yy - xxthresh %*%B1thresh
B1tDown<-B1thresh[seq_len(mL+1)]
B1tUp<-B1thresh[-seq_len(mL+1)]



xboot2<-vector("numeric", length=length(x))		#Delta Y term
xboot2[1:m]<-x[1:m]

z2<-vector("numeric", length=length(x))
z2[1:m]<-x[1:m]
#Loop
boot1thresh<-function(x){
resib2<-c(rep(0,m),sample(res1thresh, replace=TRUE))			#residual sampling, 
#resib2<-c(rep(0,m),res1thresh)					#uncomment this line to verify the bootstrap


for(i in (m+1):length(x)){
	if(round(z2[i-bestDelay-1],6)<=round(bestThresh,6)) 
		xboot2[i]<-sum(B1tDown[1],B1tDown[-1]*xboot2[i-c(1:m)],resib2[i])
	else 
		xboot2[i]<-sum(B1tUp[1],B1tUp[-1]*xboot2[i-c(1:m)],resib2[i])
	z2[i]<-xboot2[i]
	}
return(xboot2)
}


#####Bootstrap loop
model<-match.arg(test)
bootModel<-switch(model, "1vs"=bootlinear, "2vs3"=boot1thresh)

bootstraploop<-function(x, thVar=NULL){

xboot<-round(bootModel(x=x),6)

# SSR of linear boot model
string<-embed(xboot,m+1)
stringxx<-string[,-1]
xxb<-cbind(1,stringxx)
yyb <- string[,1]
SSRb<-crossprod(yyb-xxb%*%chol2inv(chol(crossprod(xxb)))%*%crossprod(xxb,yyb)) #SSR 


#SSR for threshold boot model
xxlb <- cbind(1, stringxx[, seq_len(mL)])
xxhb <- cbind(1, stringxx[, seq_len(mH)])

if(!is.null(thVar)) zb<-as.matrix(zext)
else zb<-as.matrix(stringxx[, seq_len(max(thDelay)+1)])

allgammasb<-sort(unique(zb[,b]))
ng<-length(allgammasb)
gammasb<-allgammasb[(ceiling(trim*ng)+1):floor((1-trim)*ng-1)]


###Search
IDSb<-as.matrix(expand.grid(thDelay, gammasb))			#Combinations of thresold and delays
storeb<-apply(IDSb, 1, SSRestim,yy=yyb,xxl=xxlb,xxh=xxhb,z=zb)				#Minimal SSR
SSR1threshb<-min(storeb, na.rm=TRUE)


#two thresh model
bestDelayb<-IDSb[which.min(storeb),1]
bestThreshb<-IDSb[which.min(storeb),2]

More<-list(thDelay=bestDelayb, yy=yyb, xx=xxb,z=zb)
Thresh2<-condiStep(allgammasb, bestThreshb, ninter=ninter, fun=SSRestim2, MoreArgs=More)$newThresh
Thresh3<-condiStep(allgammasb, Thresh2, ninter=ninter, fun=SSRestim2, MoreArgs=More, target=bestThreshb)
SSR2threshb<-Thresh3$SSR

#Test statistic

Ftest12b<-as.numeric(nrow(xx)*(SSRb-SSR1threshb)/SSR1threshb)
Ftest13b<-as.numeric(nrow(xx)*(SSRb-SSR2threshb)/SSR2threshb)
Ftest23b<-as.numeric(nrow(xx)*(SSR1threshb-SSR2threshb)/SSR2threshb)
list(Ftest12b, Ftest13b, Ftest23b)
}#end of bootstraploop



Ftestboot<-replicate(n=nboot,bootstraploop(x=x))
Ftestboot12<-unlist(Ftestboot[1,])
Ftestboot13<-unlist(Ftestboot[2,])
Ftestboot23<-unlist(Ftestboot[3,])

PvalBoot12<-mean(ifelse(Ftestboot12>Ftest12,1,0))
CriticalValBoot12<-quantile(Ftestboot12, probs=c(0.9, 0.95, 0.975,0.99))
PvalBoot13<-mean(ifelse(Ftestboot13>Ftest12,1,0))
CriticalValBoot13<-quantile(Ftestboot13, probs=c(0.9, 0.95, 0.975,0.99))
PvalBoot23<-mean(ifelse(Ftestboot23>Ftest12,1,0))
CriticalValBoot23<-quantile(Ftestboot23, probs=c(0.9, 0.95, 0.975,0.99))

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
	#Boundary kernel following Statistical computing with R, Rizzo
	#xx<-c(Ftestboot12,-Ftestboot12)
	#g<-density(xx, bw=bw.nrd0(x))
	#a<-seq(0,2*max(Ftestboot12),.01)
	#ghat<-approx(g$x, g$y, xout=a)
	#fhat<-2*ghat$y #density estimate along a
	#plot(a, fhat, type="l", xlim=c(0,max(Ftest12+1,max(Ftestboot12))),xlab="SSR value")
	
	#Usual kernel
	plot(density(Ftestboot12), xlab="Ftest12", xlim=c(0,max(Ftest12+1,max(Ftestboot12))),ylim=c(0,max(density(Ftestboot12)$y,dchisq(0:Ftest12, df=1+m))), main="Test linear AR vs 1 threshold TAR")
	abline(v=Ftest12, lty=2, col=2)
	curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest12+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))


	plot(density(Ftestboot13), xlab="Ftest13", xlim=c(0,max(Ftest13+1,max(Ftestboot13))),ylim=c(0,max(density(Ftestboot13)$y,dchisq(0:Ftest12, df=2*(1+m)))),main="Test linear AR vs 2 thresholds TAR")
	abline(v=Ftest13, lty=2, col=2)
	curve(dchisq(x, df=2*(1+m), ncp=0), from=0, n=Ftest13+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	}
	else {
	plot(density(Ftestboot23), xlab="Ftest23", xlim=c(0,max(Ftest23+1,Ftestboot23)), ylim=c(0,max(density(Ftestboot23)$y,dchisq(0:Ftest12, df=1+m))), main="Test 1 threshold TAR vs 2 thresholds TAR")
	abline(v=Ftest23, lty=2, col=2)
	curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest23+5, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(3,1,2), lty=c(1,1,2))
	}
}


#nlar=extend(nlar(str, coef = res$coef, fit = res$fitted.values, res = res$residuals, k = res$k,
#list( model.specific = res),"setar")
return(list(bestDelay=bestDelay,SSR=SSRs, test.val=c(Ftest12, Ftest13, Ftest23), Pvalueboot=PvalBoot, CriticalValBoot=CriticalValBoot))
}#End of thw whole function




if(FALSE){ #usage example
environment(setarTest)<-environment(setar)
#Transformation like in Hansen 1999
sun<-(sqrt(sunspot.year+1)-1)*2

setarTest(sun, m=11, thDelay=0:1, nboot=5, plot=TRUE, trim=0.1, test="1vs")
}

