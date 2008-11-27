###Build the xx matrix with 1 thresh and common=TRUE
buildXth1Common <- function(gam1, thDelay, xx,trans, ML, MH,const) {
  isL <- ifelse(trans[, thDelay + 1]< gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,ML]*isL,xx[,MH]*(1-isL))
}

###Build the xx matrix with 1 thresh and common=FALSE
buildXth1NoCommon <- function(gam1, thDelay, xx,trans, ML, MH,const) {
        isL <- ifelse(trans[, thDelay + 1]< gam1,1,0)	### isL: dummy variable
	xxL <- cbind(const,xx[,ML])*isL
	xxH <- cbind(const,xx[,MH])*(1-isL)
	xxLH<-cbind(xxL,xxH)
	return(xxLH)
}


###Build the xx matrix with 2 thresh and common=TRUE
buildXth2Common<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim){
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

###Build the xx matrix with 2 thresh and common=FALSE
buildXth2NoCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim){
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
	xxL <- cbind(const,xx[,ML])*dummydown
	xxM<-cbind(const, xx[,MM])*(1-dummydown-dummyup)
	xxH <- cbind(const,xx[,MH])*(dummyup)
	xxLMH<-cbind(xxL,xxM,xxH)

	##return result
	if(min(nup, ndown, 1-nup-ndown)>=trim){
		res <- xxLMH	#SSR
	}
	else
		res <- NA

	return(res)
}







SSR_1thresh<- function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH,const=const,fun=buildXth1Common){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH,const)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#check if equal print(c(res,res2))
	}
	return(res)
}



SSR_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2Common,trim=trim){
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
SSR_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2Common,trim=trim){
  SSR_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2NoCommon,trim=trim)
}
AIC_1thresh<-function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH,const=const,trim=trim,fun=buildXth1Common ){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH, const)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-length(yy)*log(res)+2*ncol(xx)
		#res2<-AIC(lm(yy~XX-1))}
	}
	return(res)
}

AIC_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common){
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-length(yy)*log(res)+2*ncol(xx)
		#res2<-AIC(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}

AIC_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common){
  AIC_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon)
}


buildConstants<-function(include=c("const", "trend","none", "both"), n){
  include<-match.arg(include)
  if(include=="const"){
	  const <- rep(1,n)
	  inc<-"const"}
  else if(include=="trend"){
	  const<-seq_len(n)
	  inc<-"trend"}
  else if(include=="both"){
	  const<-cbind(rep(1,n),seq_len(n))
	  inc<-c("const","trend")} 
  else {
	  const<-NULL
	  inc<-NULL
  }
  ninc<-length(inc)
  res<-list(const=const, inc=inc, ninc=ninc)
  return(res)
}