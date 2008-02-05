setarTest <- function (x, m, d = 1, steps = d, series, thDelay = 0:1, mL, mH,
    mTh, thVar, th, nboot=10, plot=FALSE) {
    if (missing(m)) m <- max(mL, mH, max(thDelay) + 1)
    if (missing(series))  series <- deparse(substitute(x))
    str <- nlar.struct(x = x, m = m, d = d, steps = steps, series = series)
    xx <- str$xx
    yy <- str$yy
    externThVar <- FALSE

###Linear model
xxlin<-cbind(1,xx)
linear <- lm.fit(xxlin, yy)
SSR<-crossprod(linear$residuals)


###Threshold transition variable
    if (missing(mL)) { mL <- m    }
    if (missing(mH)) { mH <- m    }
    if (!missing(thDelay)) {  
      if (max(thDelay) >= m)             stop(paste("thDelay too high: should be < m (=",  m, ")"))
        z <- xx[, thDelay + 1]
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

### Threshold model
ng<-nrow(xx)

SSRestim<-function(parameters) {
	thDelay<-parameters[1]
	gammai<-parameters[2]
        isL <- ifelse(z[, thDelay + 1]< gammai,1,0)	### isL: dummy 
        xxL <- cbind(1, xx[, seq_len(mL)]) * isL	### Lower matrix
        xxH <- cbind(1, xx[, seq_len(mH)]) * (1 - isL)	### Upper matrix
	crossprod(lm.fit(cbind(xxL, xxH), yy)$residuals)### SSR of thresh model
        }

gammas<-apply(xx[-c(1:m),],2, sort)[round(trim*ng):round((1-trim)*ng),]
IDS<-as.matrix(expand.grid(thDelay, gammas))
result <- apply(IDS, 1, SSRestim )
SSRthresh<-min(result)

###F test
Ftest<-as.numeric(nrow(xx)*(SSR-SSRthresh)/SSRthresh)

###Bootstrap for the F test
#initial values
xboot<-vector("numeric", length=length(x))		#Delta Y term
xboot[1:m]<-x[1:m]			
resi<-c(rep(0,m),linear$residuals)
ng<-nrow(xx)

#loop
bootstraploop<-function(B){

#reconstruction of the serie 
resi<-sample(resi, replace=TRUE)			#residual sampling, comment this line to verify the function
for(i in (m+1):length(x)){
	xboot[i]<-sum(B[1],B[-1]*xboot[i-c(1:m)],resi[i])
	}
str <- nlar.struct(x = xboot, m = m, d = d, steps = steps)

#SSR for linear model
xxb <- cbind(1,str$xx)
yyb <- str$yy
SSRb<-crossprod(yyb-xxb%*%chol2inv(chol(crossprod(xxb)))%*%crossprod(xxb,yyb)) #SSR 

#SSR for threshold model
SSRestim<-function(parameters) {
	thDelay<-parameters[1]
	gammai<-parameters[2]
        isL <- ifelse(z[, thDelay + 1]< gammai,1,0)	### isL: dummy 
        xxbT <- cbind(xxb*isL,xxb * (1 - isL))		### Upper matrix
	#crossprod(yyb-xxbT%*%chol2inv(chol(crossprod(xxbT)))%*%crossprod(xxbT,yyb)) ### SSR of thresh model
        crossprod(qr.resid(qr(xxbT),yyb))
	#crossprod(lm.fit(xxbT,yyb)$residuals)
}

gammas<-sort(unique(xboot))[round(trim*ng):round((1-trim)*ng)]	#Treshold values
IDSb<-as.matrix(expand.grid(thDelay, gammas))			#Combinations of thresold and delays
SSRthreshb<-min(apply(IDSb, 1, SSRestim ))			#Minimal SSR

#Test statistic
ng*(SSRb-SSRthreshb)/SSRthreshb
}#end of bootstraploop

Ftestboot<<-replicate(n=nboot,bootstraploop(linear$coefficients))
Ftest<<-Ftest
PvalBoot<-sum(ifelse(Ftestboot>Ftest,1,0))/nboot
CriticalValBoot<-quantile(Ftestboot, probs=c(0.9, 0.95, 0.975,0.99))


###Grahical output
if(plot==TRUE){
	plot(density(Ftestboot), xlab="SSR value", ylim=c(0,0.10), xlim=c(-10,80))
	abline(v=Ftest, lty=2, col=2)
	curve(dchisq(x, df=1+m, ncp=0), from=0, n=100, add=TRUE, col=3)
	legend("topright", legend=c("Asymptotic Chi 2", "Bootstrap", "Test value"), col=c(2,1,3), lty=c(1,1,1))
}


#nlar=extend(nlar(str, coef = res$coef, fit = res$fitted.values, res = res$residuals, k = res$k,
#list( model.specific = res),"setar")
return(list(linear.SSR=SSR, thresh.SSR=SSRthresh, test.val=Ftest, Pvalueboot=PvalBoot, CriticalValBoot=CriticalValBoot))
}#End of thw whole function

if(FALSE){ #usage example
environment(setarTest)<-environment(setar)
#Transformation like in Hansen 1999
sun<-(sqrt(sunspot.year+1)-1)*2

setarTest(sun, m=11, thDelay=0:1, nboot=2, plot=TRUE)
}

