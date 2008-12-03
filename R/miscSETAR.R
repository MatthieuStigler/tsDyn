###Build the xx matrix with 1 thresh and common=TRUE
buildXth1Common <- function(gam1, thDelay, xx,trans, ML, MH,const) {
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,ML]*isL,xx[,MH]*(1-isL))
}

###Build the xx matrix with 1 thresh and common=FALSE
buildXth1NoCommon <- function(gam1, thDelay, xx,trans, ML, MH,const) {
        isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
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



makeTh<-function(allTh, trim, th=list(exact=NULL, int=c("from","to"), around="val"), thSteps = 7,ngrid="ALL", trace=FALSE){
  ng <- length(allTh)
  down<-ceiling(trim*ng)
  up<-floor(ng*(1-trim))
  allin<-up-down
  ninter<-max(down, ng-up)

#gamma pre-specified
if(!is.null(th$exact)){
	th<-allTh[which.min(abs(allTh-th$exact))]
	if(length(th)>1){
		th<-th[1]
		cat("Many values correspond to the one you gave. The first one",th, "was taken")
		}
	ngrid<-1
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
	
	th<-allTh[seq(from=down, to=up, length.out=ngrid)]
	if(trace)
		cat("Searching on",ngrid, "possible threshold values within regimes with sufficient (",percent(trim*100,2),") number of observations\n")
}
# th<-round(th, getndp(x)) bad idea, rather use format in print and summary
res<-list(th=th, ninter=ninter)
return(res)
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


