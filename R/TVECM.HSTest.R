TVECM.HSTest <- function(data, lag=1, ngridTh=300, trim=0.05, nboot=1000, fixed.beta=NULL,  intercept=TRUE, boot.type=c("FixedReg", "ResBoot")) {


## Check args:
boot.type<-match.arg(boot.type)
dir=FALSE #internal value, was used to try different implementation of lmtest


### Organize Data
data<-as.matrix(data)
if(ncol(data)>2) {warning("Please no more than two equations")}
if(is.null(colnames(data))){colnames(data)<-paste("Var", c(1:2), sep="")}

T<-nrow(data)
p<-lag
y<-diff(data)[(p+1):(T-1),]
DeltaX<-embed(diff(data),p+1)[,-(1:2)]
if(intercept)  DeltaX<-cbind(1, DeltaX)
x<-DeltaX
t<-nrow(y)

### Compute beta with VECM() and extract ect
if(is.null(fixed.beta)){
  ve<-VECM(data, lag=lag, include="const", estim="ML")
}else{
  ve<-VECM(data, lag=lag, include="const",  beta=fixed.beta, estim="2OLS")
}

ect<-ve$model[,"ECT"]
w0<-matrix(ect[!is.na(ect)], ncol=1)

### Deprecated: Standard error
if(FALSE){
b_like<-function(b){
  z<-cbind(xlag%*%c(1,-b),x)		# z=[xlag*[1;-b],x];
  sigma<-crossprod(residuals(lm.fit(z,y)))/t
  like<-(t/2)*log(det(sigma))
  like
}

nlike1<-b_like(b0+0.001) 		#nlike1=b_like(b0+.001);
nlike2<-b_like(b0-0.001)			#nlike2=b_like(b0-.001);
hp<-(nlike1+nlike2-2*nlike)/(0.001^2)		#hp=(nlike1+nlike2-2*nlike)/(.001^2);
seb<-1/sqrt(hp)				#seb=1/sqrt(hp);
}


###Set up of the grid
q<-sort(w0)
if(ngridTh>(1-2*trim)*T) {
  ngridTh<-round((1-2*trim)*T-1)
  warning("ngridTh bigger than number of potential threshold values, set to ", ngridTh, "\n")
}

gamma2<-q[round(seq(from=trim*T, to=(1-trim)*T,length.out=ngridTh))] 
gamma2<-unique(gamma2)
ngridTh<-length(gamma2)



###########
###Lm Test
###########
lmtest02<-function(y,x,w0,gammas,dir=dir){
#y: y var, x: intercept and lags matrix, w0: ECT term, gammas: potential thresholds

  X<-cbind(w0,x) 		#X: ECT and intercept and lags
  if(dir){
    q<-qr(X)
    res_restr<-qr.resid(q,y)
  } else{
    z0zz<-X%*%solve(t(X)%*%X)
    res_restr<-lm.fit(X,y)$residuals	#residuals from the linear VECM given b0
  }
  

  res_restr1<-res_restr[,1]
  res_restr2<-res_restr[,2]
  store<-rep(NA, ngridTh)
  
  Ttrim<-trim*t
  
  ngridTh<-min(t*(1-2*trim), length(gammas))
  
  for(j in 1:ngridTh){
    d1<-ifelse(w0<=gammas[j],1,0)	#d1: dummy variable		#=w0<=gammas(j);
    n1<-sum(d1)

    if (min(c(n1,(t-n1)))>Ttrim){
      z1<-c(d1)*X
      res_unrestr <-if(dir) qr.resid(q, z1) else z1-z0zz%*%(t(X)%*%z1) #z11: residuals from unrestricted model (with threhsold)
      zea<-res_restr1*res_unrestr
      zeb<-res_restr2*res_unrestr
      ze<-cbind(zea,zeb) #	[(z11.*(e1*ones(1,length(z11(1,:))))),(z11.*(e2*ones(1,length(z11(1,:)))))];
      v<-crossprod(ze)
      z11y<-crossprod(res_unrestr,y)
      s<-matrix(c(z11y), ncol=1)				#vectorization of the parameter matrix z11y
      store[j]<-t(s)%*%solve(t(v)%*%v)%*%t(v)%*%s 	
    } #end of the if	
  } #end of the whole loop
  return(store)
} #end of the function lmtest01

#### LM test for fixed regressor bootstrap: X'X^-1 evaluated only once
lmtest02_boot<-function(y,x,w0,gammas,dir=dir){
  X<-cbind(w0,x) 		#X: ECT and intercept and lags
  if(dir){
    res_restr<-qr.resid(q,y)
  } else{
    res_restr<-lm.fit(X,y)$residuals	#residuals from the linear VECM given b0
  }
  res_restr1<-res_restr[,1]
  res_restr2<-res_restr[,2]
  store<-rep(0, ngridTh)
  ngridTh<-min(t*(1-2*trim), length(gammas))
  
  for(j in 1:ngridTh){
  d1<-ifelse(w0<=gammas[j],1,0)	#d1: dummy variable		#=w0<=gammas(j);
  n1<-sum(d1)
  
  if (min(c(n1,(t-n1)))>Ttrim){
    z1<-c(d1)*X
    res_unrestr <-if(dir) qr.resid(q, z1) else z1-z0zz%*%(t(X)%*%z1) #z11: residuals from unrestricted model (with threhsold)
    zea<-res_restr1*res_unrestr
    zeb<-res_restr2*res_unrestr
    ze<-cbind(zea,zeb) #	[(z11.*(e1*ones(1,length(z11(1,:))))),(z11.*(e2*ones(1,length(z11(1,:)))))];
    v<-crossprod(ze)
    z11y<-crossprod(res_unrestr,y)
    s<-matrix(c(z11y), ncol=1)				#vectorization of the parameter matrix z11y
    store[j]<-t(s)%*%solve(t(v)%*%v)%*%t(v)%*%s 	
  } #end of the if	
} #end of the whole loop
  lm01<-max(store, na.rm=TRUE)
  lm01
} #end of the function lmtest01
###
lm01<-lmtest02(y,x,w0,gamma2, dir=dir)


teststat<-max(lm01, na.rm=TRUE)
##################################
### Bootstraps
##################################

if(nboot==0){
  CriticalValBoot<-NULL
  PvalBoot<-NULL
  boots.reps<-NULL
}else if (nboot>0){			
##################################
### Fixed Regressor Bootstrap %
##################################
  if(boot.type=="FixedReg"){
    X<-cbind(w0,x) 		#X: ECT and intercept and lags
    Ttrim<-trim*t
    if(dir){
      q<-qr(X)
    } else{
      z0zz<-X%*%solve(t(X)%*%X)
    }
    
    lmtest_withBoot<-function(e){
      yr<-rnorm(n=t,0,1)*e
      return(lmtest02_boot(yr,x,w0,gamma2,dir=dir))
    }
    boots.reps<-replicate(nboot, lmtest_withBoot(e=residuals(ve)))
    
##################################
###Parametric Bootstrap %
##################################
  } else{
      lmtest_with_resBoot<-function(ve){
      #bootstrap it
	data.boot<-TVECM.sim(TVECMobject=ve, type="boot")
      # estimate VECM
	if(is.null(fixed.beta)){
	  ve.boot<-VECM(data.boot, lag=lag, include="const", estim="ML")
	}else{
	  ve.boot<-VECM(data.boot, lag=lag, include="const",  beta=fixed.beta, estim="2OLS")
	}
      # extract w0, y and x
	ect.boot<-ve.boot$model[,"ECT"]
	w0.boot<-matrix(ect.boot[!is.na(ect.boot)], ncol=1)
	x.boot<-ve.boot$model[,-c(1:3)]
	y.boot<-ve.boot$model[,c(1:2)]
	y.boot<-diff(y.boot)[(p+1):(T-1),]
      # set-up grid
	w0.ord.boot<-sort(w0)
	gamma2.boot<-w0.ord.boot[round(seq(from=trim*T, to=(1-trim)*T,length.out=ngridTh))] 
	gamma2.boot<-unique(gamma2.boot)
	ngridTh.boot<-length(gamma2.boot)
	test.boot<-lmtest02(y.boot,x-boot,w0.boot,gamma2.boot,dir=dir)
	return(max(test.boot, na.rm=TRUE))
      }
    boots.reps<-replicate(nboot, lmtest_with_resBoot(ve))
  }#end if boot= ResBoot


## result: compute p values and critical values:
   PvalBoot<-mean(ifelse(boots.reps>teststat,1,0))
   CriticalValBoot<-quantile(boots.reps, probs= c(0.9, 0.95,0.99))
}#end if boot>0






####### Return args
args<-list()
args$nboot<-nboot
args$boot.type<-boot.type


ret<-list()
ret$args<-args
ret$stat<-teststat
ret$values<-lm01
ret$ths<-gamma2
ret$maxTh<-gamma2[which(lm01==ret$stat)]
ret$PvalBoot<-PvalBoot
ret$CriticalValBoot<-CriticalValBoot
ret$allBoots<-boots.reps


class(ret)<-"HanSeoTest"
return(ret)

}#End of the whole function
##################################################################################################################################################################
#### END OF FUNCTION
##################################################################################################################################################################


### Print method
print.HanSeoTest<-function(x,...){
  cat("## Test of linear versus threshold cointegration of Hansen and Seo (2002) ##\n\n", sep="")
  cat("Test Statistic:\t", x$stat)
  cat("\t(Maximized for threshold value:", x$maxTh, ")\n")

  if(x$args$nboot>0){
      boot.name<-switch(x$args$boot.type, "FixedReg"="Fixed regressor bootstrap", "ResBoot"="Residual Bootstrap")
      cat("P-Value:\t", x$PvalBoot, "\t\t(",boot.name, ")\n")
    }

}

### Summary method
summary.HanSeoTest<-function(object,...){
  print(object)

  if(object$args$nboot>0){
      cat("\nCritical values:\n")
      print(matrix(object$CriticalValBoot, ncol=3, dimnames=list("", c("0.90%", "0.95%", "0.99%"))))
      cat("Number of bootstrap replications:\t", object$args$nboot)
    }
}

### Plot method
plot.HanSeoTest<-function(x,which=c("LM values","Density"),...){
  
  if(x$args$nboot==0) which<- "LM values"
# set graphic parameters
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(length(which),1))

# plot of LM values
  if("LM values"%in% which){
    plot(x$ths,x$values, type="l", ylab="LM stats", xlab="ECT values")
    title("Hansen and Seo test of linear versus threshold cointegration")
    points(x$maxTh,x$stat, col=2)
    if(x$args$nboot>0){
      abline(h=x$CriticalValBoot, lty=2, col=3:5)
      boot.name<-switch(x$args$boot.type, "FixedReg"="Fixed regressor bootstrap", "ResBoot"="Residual Bootstrap")
      tit<-paste("Critical Values (", boot.name, ")")
      legend("bottomleft", lty=2, col=3:5, leg=c("0.90%", "0.95%", "0.99%"), horiz=TRUE, title=tit, bty="n")
    }
  }
# plot of density
  if("Density"%in%which){
    plot(density(na.omit(x$allBoots)), main="Density of bootstrap distribution")
    abline(v=x$stat, col=2)
    abline(v=x$CriticalValBoot, col=3:5, lty=2)
      legend("topleft", lty=c(1,2,2,2), col=2:5, leg=c("Test value", "0.90% cv", "0.95% cv", "0.99% cv"), bty="n")
  }
}


##################################################################################################################################################################
### ENd of function
##################################################################################################################################################################

if(FALSE) {#usage example
###Test
library(tsDyn)
data(zeroyld)
data<-zeroyld

## Test against paper:
all.equal(round(TVECM.HSTest(data, lag=1, intercept=TRUE, boot=0)$stat,4),20.5994)
all.equal(round(TVECM.HSTest(data, lag=2, intercept=TRUE, boot=0)$stat,4),28.2562 )
all.equal(round(TVECM.HSTest(data, lag=3, intercept=TRUE, boot=0)$stat,4), 29.9405 )


## prob:
all.equal(round(TVECM.HSTest(data, lag=2, intercept=TRUE, boot=0, fixed.beta=1)$stat,4),29.5295)
all.equal(round(TVECM.HSTest(data, lag=1, intercept=TRUE, boot=0, fixed.beta=1)$stat,4),21.5586 )
  
## Test: no boot
TVECM.HSTest(data, lag=1, intercept=TRUE, ngridTh=50, boot=0)
TVECM.HSTest(data, lag=1, intercept=FALSE, ngridTh=50, boot=0)
TVECM.HSTest(data, lag=1, intercept=TRUE, boot=0)
TVECM.HSTest(data, lag=1, intercept=FALSE, boot=0)


## Test: boot
t1<-TVECM.HSTest(data, lag=1, intercept=TRUE, ngridTh=50, boot=5)
t2<-TVECM.HSTest(data, lag=1, intercept=FALSE, ngridTh=50, boot=5)
t3<-TVECM.HSTest(data, lag=1, intercept=TRUE, ngridTh=50, boot=5, boot.type="ResBoot")
t4<-TVECM.HSTest(data, lag=1, intercept=FALSE, ngridTh=50, boot=5, boot.type="ResBoot")

## Test: methodst1
summary(t1)
plot(t1)
plot(t1, which="Density")
t2
summary(t2)
t3
summary(t3)
t4
summary(t4)
}

