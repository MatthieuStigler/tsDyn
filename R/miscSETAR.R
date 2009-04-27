###Build the xx matrix with 1 thresh and common=TRUE
buildXth1Common <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,ML]*isL,xx[,MH]*(1-isL))
}

buildXth1LagsIncCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,1]*isL,xx[,1]*(1-isL), xx[,-1])
}

buildXth1LagsCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const*isL,xx[,1]*isL,xx[,1]*(1-isL),const*isL, xx[,-1])
}

###Build the xx matrix with 1 thresh and common=FALSE
buildXth1NoCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
        isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
	xxL <- cbind(const,xx[,ML])*isL
	xxH <- cbind(const,xx[,MH])*(1-isL)
	xxLH<-cbind(xxL,xxH)
	nisL<-mean(isL)
	if(min(nisL, 1-nisL)<trim){
		cat("\n 1 T: Trim not respected: ", c(nisL, 1-nisL), "from th:", gam1)
	}
	return(xxLH)
}





###Build the xx matrix with 2 thresh and common=TRUE
buildXth2Common<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
# print(dummydown)
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
# print(dummyup)
	nup <- mean(dummyup)
	##Construction of the matrix
#  print(c(gam1, gam2,ndown, (1-ndown-nup),nup))
	xxLMH<-cbind(const,xx[,ML]*dummydown,xx[,MM]*(1-dummydown-dummyup),xx[,MH]*(dummyup))
##return result
	if(min(nup, ndown, 1-nup-ndown)>=trim){
		res <- xxLMH	#SSR
	}
	else
		res <- NA

	return(res)
}

buildXth2LagsIncCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
	##Construction of the matrix

	xxLMH<-cbind(const,xx[,1]*dummydown,xx[,1]*(1-dummydown-dummyup),xx[,1]*(dummyup), xx[,-1])
	return(xxLMH)
}


buildXth2LagsCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
	##Construction of the matrix

	xxLMH<-cbind(const*dummydown,xx[,1]*dummydown,const*(1-dummydown-dummyup), xx[,1]*(1-dummydown-dummyup),const*dummyup, xx[,1]*(dummyup), xx[,-1])
	return(xxLMH)
}


###Build the xx matrix with 2 thresh and common=FALSE
buildXth2NoCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
# print(dummydown)
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
# print(dummyup)
	nup <- mean(dummyup)
	##Construction of the matrix
#  print(c(gam1, gam2,ndown, (1-ndown-nup),nup))
	xxL <- cbind(const,xx[,ML])*dummydown
	xxM<-cbind(const, xx[,MM])*(1-dummydown-dummyup)
	xxH <- cbind(const,xx[,MH])*(dummyup)
	xxLMH<-cbind(xxL,xxM,xxH)

	##return result
	if(min(nup, ndown, 1-nup-ndown)>=trim-0.01){
		res <- xxLMH	#SSR
	}
	else{
		cat("\nTrim not respected: ", c( ndown,1-nup-ndown, nup), "from", c(gam1, gam2))
		res <- xxLMH
	}
	return(res)
}


SSR_1thresh<- function(gam1,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH,const=const,trim, fun=buildXth1Common){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH,const, trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#check if equal print(c(res,res2))
	}
	return(res)
}

SSR<-function(X,y){
	res<-crossprod(y- X %*%chol2inv(chol(crossprod(X)))%*%crossprod(X,y))
	#res2<-deviance(lm(y~X-1))
	return(res)
	}

AIC.matrices<-function(X,y, T, k=2){
	SSR<-crossprod(y- X %*%chol2inv(chol(crossprod(X)))%*%crossprod(X,y))
	res<-T*log(SSR/T)+k*(ncol(X)+2)
	return(res)
	}




SSR_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim){
	XX<-buildXth2Common(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}

# SSR_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2Common,trim=trim){
#   SSR_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2NoCommon,trim=trim)
# }

SSR_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim){
  	XX<-buildXth2NoCommon(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}

AIC_1thresh<-function(gam1,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH,const=const,trim=trim,fun=buildXth1Common , k=2, T, temp=FALSE){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH, const, trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA}
	else{
		SSR <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-T*log(SSR/T)+k*(ncol(XX)+1)
		#res2<-AIC(lm(yy~XX-1))}
	}
	return(res)
}

AIC_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2, T, temp=FALSE){
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA}
	else{
		SSR <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-T*log(SSR/T)+k*(ncol(XX)+2)
		#res2<-AIC(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}

AIC_2th <-function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2, T, temp=FALSE){
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA}
	else{
		SSR <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-T*log(SSR/T)+k*(ncol(XX)+2)
		print(SSR)
		print(T)
		print(k)
		print(ncol(XX))
		#res2<-AIC(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}

SSR_2th<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon){
  	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}



AIC_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2,T, temp=FALSE){
  AIC_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon, k=k,T, temp=temp)
}


buildConstants<-function(include=c("const", "trend","none", "both"), n){
  include<-match.arg(include)
  if(include=="const"){
	  const <- rep(1,n)
	  incNames<-"const"}
  else if(include=="trend"){
	  const<-seq_len(n)
	  incNames<-"trend"}
  else if(include=="both"){
	  const<-cbind(rep(1,n),seq_len(n))
	  incNames<-c("const","trend")} 
  else {
	  const<-NULL
	  incNames<-NULL
  }
  ninc<-length(incNames)
  res<-list(const=const, incNames=incNames, ninc=ninc)
  return(res)
}



MakeThSpec<-function(ngrid=c("All", "Half", "Third", "Quarter"), exact=NULL, int=c("from","to"), around="val",...){
  if(is.character(ngrid))
    ngrid<-match.arg(ngrid)
  exCheck<-ifelse(is.null(exact),0,1)
  inCheck<-ifelse(int==c("from","to"),0,1)
  aroundCheck<-ifelse(around=="val",0,1)
  if(sum(exCheck, inCheck, aroundCheck)>1)
    stop("Only one of the makeThSpec args should be specified")
  return(list(exact=exact, int=int, around=around, ngrid=ngrid))
} 

makeTh<-function(allTh, trim, th=list(exact=NULL, int=c("from","to"), around="val"), thSteps = 7,ngrid="ALL", trace=FALSE, nthresh=1){
  ng <- length(allTh)
  down<-ceiling(trim*ng)
  up<-floor(ng*(1-trim))
  allin<-up-down

#threshold pre-specified
if(!is.null(th$exact)){
  if(length(th$exact)==1){
    if(nthresh==2)
      stop("Please either provide two pre-specified threshold values or set nthresh to 1")
    th<-unique(allTh[which.min(abs(allTh-th$exact))])
    if(length(th)>1){
      th<-th[1]
      cat("Many values correspond to the threshold you gave. The first one",th, "was taken")
    }
    ngrid<-1
  }
  else if(length(th$exact)==2){
    th1<-unique(allTh[which.min(abs(allTh-th$exact[1]))])
    th2<-unique(allTh[which.min(abs(allTh-th$exact[2]))])
    if(length(c(th1, th2))>2){
      th1<-th1[1]
      th2<-th2[2]
      cat("Many values correspond to the threshold you gave. The first ones",c(th1, th2), "were taken")
    }
  }
  else
    warning("Too many threshold values given")
}	  
#interval to search inside given by user
else if(is.numeric(th$int)){
	if(missing(ngrid))
		ngrid<-20
	intDown<-which.min(abs(allTh-th$int[1]))
	intUp<-which.min(abs(allTh-th$int[2]))
	if(length(intDown)>1|length(intUp)>1)
		intDown<-intDown[1];intUp<-intUp[1];
	lengthInt<-min(ngrid,intUp-intDown)
	if(trace)
		cat("Searching within",lengthInt, "values between",allTh[intDown], "and", allTh[intUp],"\n")
	th<-allTh[seq(from=intDown, to=intUp, length.out=lengthInt)]
	}
#value to search around	given by user
else if(is.numeric(th$around)){
	if(missing(ngrid))
		ngrid<-20
	if(trace)
		cat("Searching within", ngrid, "values around", th$around,"\n")
	th<-aroundGrid(th$around,allvalues=allTh,ngrid=ngrid,trim=trim, trace=trace) #fun stored in TVECM.R
}

#Default method: grid from lower to higher point
else{
	if(ngrid=="ALL")
		ngrid<-allin
	else if(ngrid>allin)
		ngrid<-allin
	
	ths<-allTh[seq(from=down, to=up, length.out=ngrid)]
	th <-unique(ths)
	if(trace)
		cat("Searching on",length(th), "possible threshold values within regimes with sufficient (",percent(trim*100,2),") number of observations\n")
}
# th<-round(th, getndp(x)) bad idea, rather use format in print and summary
return(th)
}

if(FALSE){
environment(makeTh)<-environment(star)
a<-makeTh(unique(embed(lynx, 2)[,2, drop=FALSE]), trim=0.15)
a
length(a)
length(unique(a))
length(lynx)
length(unique(lynx))
  
  cat("ng", ng, "\n")
  cat("ng*(1-trim)", ng*(1-trim), "\n")
  cat("up", up, "\n")
  print(c(down, up, allin, ninter))
}


