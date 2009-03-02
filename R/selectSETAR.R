## Copyright (C) 2005/2006/2009  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY;without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.


###Exhaustive search over a grid of model parameters
selectSETAR<- function (x, m, d=1, steps=d, series, mL, mH,mM, thDelay=0, mTh, thVar, th=MakeThSpec(), trace=TRUE, include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL),MH=seq_len(mH), MM=seq_len(mM),nthresh=1,trim=0.15,criterion = c("pooled-AIC", "AIC", "BIC","SSR"),thSteps = 7,ngrid="ALL",  plot=TRUE,max.iter=2, type=c("level", "diff", "ADF"), same.lags=FALSE) 
#This internal function called by setar() makes a grid search over the values given by setar or user
#1: Build the regressors matrix, cut Y to adequate (just copy paste of function setar() )
#2: establish the transition variable z (just copy paste of function setar() )
#3. set-up the grid following argument th or by default ngrid=="ALL" on all values
#4: Sets up functions to compute the SSR or AIC on model with 1 or 2 thresh, or Pooled AIC
#5: establish a grid with addition of thDelay (crit SSR) or thDelay, Ml and Mh  (crit AIC and pooledAIC)
#6: apply the function in 3 to grid in 4
#7: sorts the results to show the best values
#8: Computation for 2 thresh using condiStep (see TVARestim.r)
#9: plot of the results of the grid search
#10: return results
{
### SelectSETAR 0: preliminary checks
  include<-match.arg(include)
  model<-match.arg(model)
  type<-match.arg(type)
  criterion<-match.arg(criterion)
  if(missing(m))
    m <- max(ML, MH, ifelse(nthresh==2, max(MM),0),thDelay+1)
  if(missing(series))
      series <- deparse(substitute(x))

  if (criterion == "pooled-AIC"&nthresh==2){
      warning("\ncriterion pooled AIC currently not implemented for nthresh=2, criterion SSR is used\n")
      criterion<-"SSR"
    }
### SelectSETAR 1:  Build the regressors matrix, cut Y to adequate (just copy paste of function setar() )    
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	
	if(type=="level"){
	  xx <- getXX(str)
	  yy <- getYY(str)
	}
	else{ 
	  if(type=="diff"){
	    xx <- getdXX(str)
	    yy <- getdYY(str)
	  }
	  else if(type=="ADF"){
	    xx <- cbind(getdX1(str),getdXX(str))
	    yy <- getdYY(str)
	  }
	  str$xx<-xx
	  str$yy<-yy
	}
	
  externThVar <- FALSE
  
##SelectSETAR: Lags selection
# 
#   if(criterion=="SSR"&nthresh==2){
#     if(missing(m)&all(missing(MM),missing(mM)))
#       stop("With criterion SSR, for 2 threholds, number of lags in inner regime has to be given")
#   }
  
  ###ML
  if(missing(ML)) {		#ML: different lags
    if (missing(mL)) {	#mL: suit of lags
      if(missing(m))
	cat("arg m is missing")
      else
	mL <- m
      if (trace) 
	cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
    }
    ML <- seq_len(mL)
  }

  ###MH
  if(missing(MH)) {
    if (missing(mH)) {
    mH <- m
      if (trace) 
	cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
    }
  MH <- seq_len(mH)
  }
  
  ###MM
  if(missing(MM)) {
    if (missing(mM)) {
    mM <- m
      if (trace&nthresh==2) 
	cat("Using maximum autoregressive order for middle regime: mM =", m,"\n")
    }
  MM <- seq_len(mM)
  }

      
  if(!all(isTRUE(all.equal(ML,MM))&isTRUE(all.equal(ML,MH)))&same.lags==TRUE)
    stop("Arg same.lag set to true but different lags given. Choose either one of those")
    
  if(criterion=="SSR")
    same.lags<- if(isTRUE(all.equal(ML,MM))&isTRUE(all.equal(ML,MH))) TRUE else FALSE


###SelectSETAR: includes const, trend
  if(include=="none" && any(c(ML,MM,MH)==0))
    stop("you cannot have a regime without constant and lagged variable")
  constMatrix<-buildConstants(include=include, n=nrow(xx)) #stored in miscSETAR.R
  incNames<-constMatrix$incNames #vector of names
  const<-constMatrix$const #matrix of none, const, trend, both
  ninc<-constMatrix$ninc #number of terms (0,1, or 2)

### selectSETAR 2: Set-up of transition variable
#two models: TAR or MTAR (z is differenced)
#three possibilitiees for thVar:
#thDelay: scalar: prespecified lag, or vector: of lags to search for. Z is a matrix
#mTh: combination of lags. Z is one vector-matrix
#thVar: external variable Zis one vector-matrix.
# Default: thDelay=0
	if(!missing(thDelay)) {
		if(max(thDelay)>=m) 
			stop(paste("thDelay too high: should be < m (=",m,")"))
		if(model=="TAR"){
			z <- xx
			z2<-embedd(x, lags=c((0:(m-1))*(-d), steps) )[,1:m,drop=FALSE]
			z4<-embed(x,m+1)[,-1]
		}
		else{
			if(thDelay==m-1)
				stop("th Delay too high, should be <m-1 (=",m,")(because of differencing)")
 		z<-embed(diff(x),m)[,thDelay+2]
		#print(cbind(yy,xx,z))
		#z<-z[if(m-thDelay-2>0)-seq_len(m-thDelay-2>0) else seq_len(nrow(z)),thDelay+2]
		}
	}
 	else if(!missing(mTh)) {
		if(length(mTh) != m) 
			stop("length of 'mTh' should be equal to 'm'")
		z <- xx %*% mTh #threshold variable
		dim(z) <- NULL
	} else if(!missing(thVar)) {
		if(length(thVar)>nrow(xx)) {
			z <- thVar[seq_len(nrow(xx))]
			if(trace) 
				cat("Using only first", nrow(xx), "elements of thVar\n")
		}
		else 
			z <- thVar
		externThVar <- TRUE
	} else {
		z <- xx
	}


z<-as.matrix(z)

### selectSETAR 3: set-up of the grid
#Possibilities:
#gamma pre-specified
#interval to search inside given by user
#value to search around	given by user
#Default method: grid from lower to higher point
  allTh <- sort(z[,1])
  th<-makeTh(allTh=allTh, trim=trim, th=th,ngrid="ALL", trace=trace, nthresh=nthresh)
  
### selectSETAR 4: Sets up functions to compute the SSR/AIC/Pooled-AIC for th= 1 or 2
#have been moved to miscSETAR.R, call before functions to build the regressors matrix
### selectSETAR 4b:Function pooled AIC 
pooledAIC <- function(parms) {	
	nParms<-length(parms)
	thDelayVal <- parms[1] + 1
	mLVal <- parms[2]
	mHVal <- parms[nParms-1]
	m <- max(thDelayVal, mLVal, mHVal)
	lags <- c( (seq_len(m)-1) * (-d), steps)

	xxyy <- embedd(x, lags = lags)

	z <- xxyy[, thDelayVal]
	isLow <- (z <= parms[nParms])

	if ((sum(isLow) < mLVal) | (sum(!isLow) < mHVal)) 
	    return(NA)

	xx <- xxyy[isLow, seq_len(mLVal)]
	y <- xxyy[isLow, m + 1]
	AIC1 <- AIC(lm(y ~ xx))

	xx <- xxyy[!isLow, seq_len(mHVal)]
	y <- xxyy[!isLow, m + 1]
	AIC2 <- AIC(lm(y ~ xx))
	return(AIC1 + AIC2)
}


###selectSETAR 5: Grid of combinations of all parameters
  if(criterion=="SSR"){
    ncombin<-length(thDelay)*length(th)
    if(trace)
      cat("Searching on ", ncombin, " combinations of threshold and thDelay\n", sep="")
    IDS <- as.matrix(expand.grid(thDelay, th))
    colnames(IDS)<-c("thDelay", "th")
  }
  else{
    if(same.lags){
      ncombin<-length(thDelay)*length(th)*length(ML)
      if(trace)
	cat("Searching on ", ncombin, " combinations of threshold, thDelay and m \n")
      IDS <- as.matrix(expand.grid(thDelay,  ML, th))			
      colnames(IDS)<-c("thDelay", "m", "th")
    }
    else{
      ncombin<-length(thDelay)*length(th)*length(ML)*length(MH)
      if(trace)
	cat("Searching on ", ncombin, " combinations of threshold, thDelay, ML and MM\n")
      IDS <- as.matrix(expand.grid(thDelay,  ML, MH, th))			
      colnames(IDS)<-c("thDelay", "mL", "mH", "th")
   }
  }


###selectSETAR 6: Computation for 1 thresh
  if (criterion == "pooled-AIC") {
      computedCriterion <- apply(IDS, 1, pooledAIC)
  } 
  else if(criterion%in%c("AIC","BIC")){
    kaic<-switch(criterion, "AIC"=2, "BIC"=log(length(yy)))
    mHPos<-ifelse(same.lags, 2,3)
    if(common)
      computedCriterion <- mapply(AIC_1thresh, gam1=IDS[,"th"], thDelay=IDS[,1],ML=IDS[,2],MH=IDS[,mHPos], MoreArgs=list(xx=xx,yy=yy,trans=z,const=const,trim=trim,fun=buildXth1Common, k=kaic))
    else
      computedCriterion <- mapply(AIC_1thresh, gam1=IDS[,"th"], thDelay=IDS[,1],ML=IDS[,2],MH=IDS[,mHPos], MoreArgs=list(xx=xx,yy=yy,trans=z,const=const,trim=trim,fun=buildXth1NoCommon, k=kaic))
  } 
  else if(criterion=="SSR"){
    if(common)
      computedCriterion <- mapply(SSR_1thresh, gam1=IDS[,2],thDelay=IDS[,1],MoreArgs=list(xx=xx,yy=yy,trans=z, ML=ML, MH=MH, const=const,trim=trim, fun=buildXth1Common))
    else
      computedCriterion <- mapply(SSR_1thresh, gam1=IDS[,2],thDelay=IDS[,1],MoreArgs=list(xx=xx,yy=yy,trans=z, ML=ML, MH=MH, const=const,trim=trim, fun=buildXth1NoCommon))
  }

###selectSETAR 7: sorts the results to show the best values
  allres <- cbind(IDS, computedCriterion)
  colnames(allres) <- c(colnames(IDS), criterion)
  idSel <- sort(computedCriterion, index = TRUE)$ix
  idSel <- idSel[seq_len(min(ifelse(nthresh==1,10,5), length(idSel)))]
  res <- data.frame(allres[idSel, , drop=FALSE], row.names = NULL, check.names=FALSE)
  
  if(criterion=="SSR"){
    bests<-c(res[1,"thDelay"],res[1,"th"],res[1,criterion])
    names(bests)<-c("thDelay", "th",criterion)
  }
  else{
    if(same.lags){
      bests<-c(res[1,"thDelay"],res[1,"th"],res[1,criterion], res[1, "m"])
      names(bests)<-c("thDelay", "th",criterion, "m")
    }
    else{
      bests<-c(res[1,"thDelay"],res[1,"th"],res[1,criterion], res[1, "mL"], res[1,"mH"])
      names(bests)<-c("thDelay", "th",criterion, "mL", "mH")
    }
  }
  
  firstBests<-bests
###selectSETAR 8: Computation for 2 thresh 
#use function condistep (see TVARestim.r) to estimate the second th given the first
#iterate the algortihm: once the second is estimate, reestimate the first, and alternatively until convergence
  if(nthresh==2){
    ###trace infos
    if(trace){
      if(criterion=="SSR")
	cat("Result of the one threshold search:\n -Thresh: ",res[1,"th"],"\t-Delay: ",res[1,"thDelay"],"\t-",criterion, res[1,criterion], "\n" )
      else {
	if(same.lags)
	  cat("Result of the one threshold search:\n -Thresh: ",res[1,"th"],"\t-Delay: ",res[1,"thDelay"],"\t-m:",res[1,"m"], "\t-",criterion, res[1,criterion], "\n" )
	else
	  cat("Result of the one threshold search:\n -Thresh: ",res[1,"th"],"\t-Delay: ",res[1,"thDelay"],"\t-mL:",res[1,"mL"], "\t-mH:",res[1,"mH"],"\t-",criterion, res[1,criterion], "\n" )
      }
    }
    
    ###set lags according to first search
    potMM<-list()
    
    if(criterion%in%c("AIC", "BIC")){
      if(same.lags){
	ML<-firstBests["m"]
	MH<-ML
	potMM[[1]]<-ML
      }
      else{
	ML<-firstBests["mL"]
	MH<-firstBests["mH"]
	print(MM)
	for(i in 1:length(MM))
	  potMM[[i]]<-seq_len(MM[i])
      }
    }
    else
      potMM[[1]]<-ML
    
    if(common)
      func<-switch(criterion, "AIC" = AIC_2threshCommon,"BIC" = AIC_2threshCommon, "SSR" = SSR_2threshCommon)	
    else
      func<-switch(criterion, "AIC" = AIC_2threshNoCommon,"BIC" = AIC_2threshNoCommon, "SSR" = SSR_2threshNoCommon)
    
    Bests<-matrix(NA, nrow=length(potMM), ncol=4)
    colnames(Bests)<-c("th1", "th2", criterion, "pos")
    for(j in 1:length(potMM)){
  
      More<-list(yy=yy, xx=xx,trans=z, ML=ML, MH=MH,MM=potMM[[j]], const=const,trim=trim)
      #first conditional search
      last<-condiStep(allTh,threshRef=res[1,"th"], delayRef=res[1,"thDelay"], fun=func, trim=trim, trace=trace, More=More)
      
      #iterative loop for conditional search
      i<-1	#initialise the loop
      while(i<max.iter){
	      b<-condiStep(allTh,last$newThresh, delayRef=res[1,1], fun=func, trim=trim, trace=trace, More=More)
	      if(b$SSR<last$SSR){	#minimum still not reached
		      i<-i+1
		      last<-b}
	      else{			#minimum reached
		      i<-max.iter
		      last<-b}
      }
    Bests[j,]<-c(min(c(last$threshRef, last$newThresh)),max(c(last$threshRef, last$newThresh)), last$SSR, max(potMM[[j]]))  
    }  
    if(length(potMM)>1){
      index<-sort(Bests[,3], index=TRUE)$ix
      BestsOrdered<-Bests[index,,drop=FALSE]
      MM<-potMM[[BestsOrdered[1,"pos"]]]
      
      res2<-unlist(lapply(potMM, max))[index]
      res2<-cbind(res2,BestsOrdered[,-4])
      colnames(res2)<-c("mM", "th1", "th2", criterion)
      res2 <- data.frame(res2[1:min(nrow(res2),5),], row.names = NULL, check.names=FALSE)
    }
    else{
      res2 <- data.frame(Bests[,-4, drop=FALSE], row.names = NULL, check.names=FALSE)
    }
    
    if(length(potMM)>1){
      bests2th<-c(max(MM), res2[1,"th1"],res2[1,"th2"], res2[1,criterion])
      names(bests2th)<-c("mM", "th1","th2", criterion)
    }
    else{
      bests2th<-c(res2[1,"th1"],res2[1,"th2"], res2[1,criterion])
      names(bests2th)<-c("th1","th2", criterion)
    }
    
    bests<-c(res[1,"thDelay"], res2[1,"th1"],res2[1,"th2"], res2[1,criterion])
    names(bests)<-c("thDelay","th1","th2", criterion)
  }
  else {#ntresh!=2
    res2<-NULL
    bests2th<-NULL
  }
###selectSETAR 9: plot of the results of the grid search
if(plot==TRUE){
	allcol <- seq_len(max(thDelay+1)*max(ML)*max(MH))
	col <- switch(criterion, "AIC"=allcol, "BIC"=allcol,"pooled-AIC"=allcol,"SSR"=(thDelay+1) )
	big <- apply(expand.grid(thDelay,ML, MH),1,function(a) paste("Th:", a[1],"mL:", a[2], "mH:", a[3]))
	legend <- switch(criterion, "AIC"=big,"BIC"=big, "pooled-AIC"=big, "SSR"=paste("Threshold Delay", thDelay))

	plot(allres[,"th"], allres[,ncol(allres)], col=col, xlab="Treshold Value",ylab=criterion, main="Results of the grid search")
 	legend("topleft", pch=1, legend=legend, col=col, bg=0)
}

###selectSETAR 10: return results

ret<-list(res=res, res2=res2,bests=bests, th=getTh(bests), nthresh=nthresh, allTh=allres, criterion=criterion,firstBests=firstBests, bests2th=bests2th, ML=ML, MM=MM, MH=MH, same.lags=same.lags)
class(ret)<-"selectSETAR"
return(ret)
}

print.selectSETAR<-function(x,...){
  cat("Results of the grid search for 1 threshold\n")
  if(x$criterion=="SSR"){
    if(x$same.lags)
      cat("\t Conditional on m= ",max(x$ML), "\n")
    else
	cat("\t Conditional on ML= ",x$ML, "and MH= ", x$MH, "\n")
    }
  print(x$res)
  if(x$nthresh==2){
    if(x$criterion=="SSR"){
      if(x$same.lags)
	cat("\nResults of the grid search for 2 thresholds\n \t Conditional on thDelay = ",x$bests["thDelay"], " and m =",max(x$ML), "\n")
      else
	cat("\nResults of the grid search for 2 thresholds\n \t Conditional on thDelay = ",x$bests["thDelay"], ", and on ML=",x$ML, ", MM=", x$MM, ", MH=", x$MH, "\n")
    }
    else{
      if(x$same.lags)
	cat("\nResults of the grid search for 2 thresholds\n \t Conditional on thDelay = ",x$bests["thDelay"], " and m =",max(x$ML), "\n")
      else
	cat("\nResults of the grid search for 2 thresholds\n \t Conditional on thDelay = ",x$bests["thDelay"], " mL=",x$res[1,"mL"], ", mH =", x$res[1,"mL"], "\n")
    }
   print(x$res2)
  cat("\nOverall best results:\n")
  print(x$bests)
  cat("With lags:\n\t-ML:", x$ML, "\n")
  if(x$nthresh==2)
    cat("\t-MM:", x$MM, "\n")
  cat("\t-MH:", x$MH, "\n")
  }
}



plot.selectSETAR<-function(x,...){
  if(x$criterion!="SSR")
    stop("Currently only implemented for criterion SSR")
  plot3(th=x$th,nthresh=x$nthresh,allTh= x$allTh)
#  plot3a(th=x$th,nthresh=x$nthresh,allTh= x$allTh, bestDelay=x$bests[1])
}

###Try it
if(FALSE) { #usage example
library(tsDyn)
environment(selectSETARmat)<-environment(selectNNET)
#Transformation like in Hansen 1999
sun<-(sqrt(sunspot.year+1)-1)*2		

###Full grid search with OLS
selectSETARmat(sun, m=3, criterion="SSR", d=1, thDelay=0:2,model="TAR", trim=0.15, max.iter=10, plot=FALSE, nthresh=2)

selectSETARmat(sun, m=3, criterion="AIC", d=1, thDelay=0:2,model="TAR", trim=0.25, max.iter=10, plot=FALSE, nthresh=1)

###restricted search with AIC or AIC pooled around the max selected by OLS
selectSETARmat(sun, m=2, criterion="AIC", d=1, thDelay=0:1, around=7.444575)
}
