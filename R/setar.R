## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
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

#SETAR model contructor	(sequential conditional LS)
#	th: threshold. If not specified, a grid of reasonable values is tried
#	m: general autoregressive order (mL=mH)
#	mL: autoregressive order below the threshold ('Low')
#	mH: autoregressive order above the threshold ('High')
#	nested: is this a nested call? (useful for correcting final model df)
#	trace: should infos be printed?
setar <- function(x, m, d=1, steps=d, series, mL,mM,mH, thDelay=0, mTh, thVar, th, trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL),MM=seq_len(mM), MH=seq_len(mH), nthresh=1,trim=0.15){
# 1: preliminaries
# 2:  Build the regressors matrix and Y vector
# 3: Set-up of transition variable
# 4: Search of the treshold if th not specified by user
# 5: Build the threshold dummies and then the matrix of regressors
# 6: compute the model, extract and name the vec of coeff
# 7: return the infos



###SETAR 1: preliminaries
	include<-match.arg(include)

	model<-match.arg(model)
	if(missing(m))
		m <- max(ML, MH, thDelay+1)
	if(!missing(th)){
		if(length(th)==2)
			nthresh<-2
	}
  if(missing(series))
    series <- deparse(substitute(x))
### SETAR 2:  Build the regressors matrix and Y vector
        str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- getXX(str)
	yy <- getYY(str)
	externThVar <- FALSE
	##Lags selection
	if(missing(ML)) {		#ML: different lags
		if (missing(mL)) {	#mL: suit of lags
			mL <- m
			if (trace) 
				cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
		}
		ML <- seq_len(mL)
	}

	if(missing(MM)) {
		if (missing(mM)) {
			mM <- m
		    if (trace&nthresh==2) 
				cat("Using maximum autoregressive order for middle regime: mM =", m,"\n")
		}
		MM <- seq_len(mM)
	}
	if(missing(MH)) {
		if (missing(mH)) {
			mH <- m
			if (trace) 
				cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
		}
		MH <- seq_len(mH)
	}

	###includes const, trend
  if(include=="none" && any(c(ML,MM,MH)==0))
    stop("you cannot have a regime without constant and lagged variable")
	if(include=="const"){
		const <- rep(1,nrow(xx))
		inc<-"const"}
	else if(include=="trend"){
		const<-seq_len(nrow(xx))
		inc<-"trend"}
	else if(include=="both"){
		const<-cbind(rep(1,nrow(xx)),seq_len(nrow(xx)))
		inc<-c("const","trend")} 
	else {
		const<-NULL
		inc<-NULL
	}
	ninc<-length(inc)

### SETAR 3: Set-up of transition variable
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
			z <- xx		#xx <- getXX(str) 
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
		if(trace) 
			cat("Using default threshold variable: thDelay=0\n")
		z <- xx
		thDelay<-0
	}

z<-as.matrix(z)

### SETAR 4: Search of the treshold if th not specified by user
#if nthresh==1, try over a reasonable grid (30), if nthresh==2, whole values
#call the function selectSETAR
	if (missing(th)) { 
		ngrid<-ifelse(nthresh==1,30,"ALL") #if 1 thresh grid with 30 values, if 2 th all values
		search<-selectSETAR(x, m, d=d, steps=d, series, mL=mL, mH=mH,mM=mM, thDelay=thDelay, mTh, thVar, trace=trace, include = include, common=common, model=model, ML=ML,MH=MH, MM=MM,nthresh=nthresh,trim=trim,criterion = "SSR",thSteps = 7,ngrid=ngrid, plot=FALSE,max.iter=2)
		thDelay<-search$bests[1]
		th<-search$bests[2:(nthresh+1)]
		nested<-TRUE
		if(trace) {
			cat("Selected threshold: ", th,"\n")
			cat("Selected delay: ", thDelay,"\n")
		}
		#missing(th)<-FALSE
	}
### SETAR 5: Build the threshold dummies and then the matrix of regressors
	#if(!missing(th)) {
		#check number of observations)
		if(nthresh==1){
			isL <- ifelse(z[, thDelay + 1]<=th, 1, 0)
			isM<-NA
			isH <- 1-isL}	
		else{
			isL <- ifelse(z[, thDelay + 1]<=th[1], 1, 0)
			isH <- ifelse(z[, thDelay + 1]>th[2], 1, 0)
			isM <- 1-isL-isH
		}

		nobs<-na.omit(c(mean(isL),mean(isM),mean(isH)))	#N of obs in each regime
		if(min(nobs)<trim){
			if(trace)
				cat("\nWith the threshold you gave, there a not ",trim,"observations in each regime")
		}
		if(min(nobs)==0)
			stop("With the threshold you gave, there is a regime with no observations!")
		#build the X matrix
		if(nthresh==1){
		  if(common)
		    xxLH<-buildXth1Common(gam1=th, thDelay, xx=xx,trans=z, ML=ML, MH=MH,const)	
		  else
		    xxLH<-buildXth1NoCommon(gam1=th, thDelay, xx=xx,trans=z, ML=ML, MH=MH,const)	
			midCommon<-mid<-NA}
		else if(nthresh==2){
		  if(common)
		    xxLH<-buildXth2Common(gam1=th[1],gam2=th[2],thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const,trim=trim)
		  else
		    xxLH<-buildXth2NoCommon(gam1=th[1],gam2=th[2],thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const,trim=trim)
			#midCommon<-c(paste(inc,rep(3,ninc)),paste("phi3", MH, sep=".")) no more used: 
			#mid<-c(paste("phi3", MH, sep=".")) #no more used
		}

### SETAR 6: compute the model, extract and name the vec of coeff
		res <- lm.fit(xxLH, yy)
#Coefficients and names
		res$coefficients <- c(res$coefficients, th)

		if(FALSE){
			names(res$coefficients) <- na.omit(c(paste(inc, rep(1,ninc)), paste("phi1", ML, sep="."), paste(inc,rep(2,ninc)),paste("phi2", if(nthresh==1)MH else MM, sep="."),midCommon, rep("th",nthresh)))
		} 
if(FALSE){
			names(res$coefficients) <- na.omit(c(inc,paste("phi1", ML, sep="."),paste("phi2",if(nthresh==1)MH else MM, sep="."),mid, rep("th",nthresh)))
		}

if(!common){
  if(nthresh==1)
    co<-c(getIncNames(inc,ML), getArNames(ML), getIncNames(inc,MH), getArNames(MH),"th")
  else
    co<-c(getIncNames(inc,ML), getArNames(ML), getIncNames(inc,MM), getArNames(MM), getIncNames(inc,MH), getArNames(MH),"th1","th2")
}
else{
  if(nthresh==1)
    co<-c(inc, getArNames(ML), getArNames(MH),"th")
  else 
    co<-c(inc, getArNames(ML), getArNames(MM), getArNames(MH),"th1","th2")
}
  names(res$coefficients) <- na.omit(co)
### SETAR 7: return the infos
		res$k <- if(nested) (res$rank+nthresh) else res$rank	#If nested, 1 more fitted parameter: th
		res$thDelay<-thDelay
		res$fixedTh <- if(nested) FALSE else TRUE
		res$mL <- max(ML)
		res$mH <- max(MH)
		res$ML <- ML
		res$MH <- MH
		res$externThVar <- externThVar
		res$thVar <- z[, thDelay + 1]
		res$nconst<-inc
		res$common<-common  	#wheter arg common was given by user
		res$nthresh<-nthresh 	#n of threshold
		res$model<-model
		res$RegProp <- c(mean(isL),mean(isH))
		res$usedThVar<-z[,thDelay+1]
		res$trim<-trim
		if(nthresh==2)
			res$RegProp <- c(mean(isL),mean(isM),mean(isH))
		res$VAR<-as.numeric(crossprod(na.omit(res$residuals))/(nrow(xxLH)))*solve(crossprod(xxLH))
		if(!externThVar) {
			if(missing(mTh)) {
				mTh <- rep(0,m)
				mTh[thDelay+1] <- 1
			}
			res$mTh <- mTh
		}
		return(extend(nlar(str,	coef=res$coef,	fit=res$fitted.values,	res=res$residuals,
			k=res$k,model.specific=res), "setar"))
	#}
}

getSetarXRegimeCoefs <- function(x, regime=c("L","M","H")) {
	regime <- match.arg(regime)
	x <- x$coef
	x1 <- x[grep(paste("^phi", regime, "\\.", sep=""), names(x))]
	x2 <- x[grep(paste("^const ", regime, "$", sep=""), names(x))]
	x3 <- x[grep(paste("^trend ", regime, "$", sep=""), names(x))]
	return(c(x1, x2, x3))
}

getTh<-function(x){
	x[grep("th",names(x))]}

#gets a vector with names of the arg inc
getIncNames<-function(inc,ML){
  ninc<-length(inc)
    letter<-deparse(substitute(ML))
    letter<-sub("M","",letter)
    paste(inc, rep(letter,ninc))
}

#get a vector with names of the coefficients
getArNames<-function(ML){
  if(any(ML==0))
    return(NA)
  else{
    letter<-deparse(substitute(ML))
    letter<-sub("M","",letter)
    paste(paste("phi",letter, sep=""),".", ML, sep="")
  }
}
print.setar <- function(x, ...) {
	NextMethod(...)
	x.old <- x
	x <- x$model.specific
	order.L <- x$mL
	order2.L <- length(x$ML)
	order.H <- x$mH
	order2.H <- length(x$MH)
	common <- x$common
	nconst <- x$nconst
	nthresh<-x$nthresh
	externThVar <- x$externThVar
	cat("\nSETAR model (",nthresh+1,"regimes)\n")
	cat("Coefficients:\n")
	if(common==FALSE){
		lowCoef <- getSetarXRegimeCoefs(x.old, "L")
		highCoef<- getSetarXRegimeCoefs(x.old, "H")
		cat("Low regime:\n")
		print(lowCoef, ...)
		if(nthresh==2){
			midCoef<- getSetarXRegimeCoefs(x.old, "M")
			cat("\nMid regime:\n")
			print(midCoef, ...)}
		cat("\nHigh regime:\n")
		print(highCoef, ...)
	} else {
		print(x.old$coeff[-length(x.old$coeff)], ...)
	}
	thCoef<-getTh(coef(x.old))
	cat("\nThreshold:")
	if(x$model=="MTAR"){
		cat("\nMomentum Threshold (MTAR) Adjustment")
		D<-"Δ"}
	else
		D<-NULL
	cat("\n-Variable: ")
        if(externThVar)
          cat("external")
        else {
          cat('Z(t) = ')
          cat('+ (',format(x$mTh[1], digits=2), paste(")",D," X(t)", sep=""), sep="")
          if(length(x$mTh)>1)
            for(j in 1:(length(x$mTh) - nthresh)) {
              cat('+ (', format(x$mTh[j+1], digits=2), paste(")",D,"X(t-", j, ")", sep=""), sep="")
            }
          cat('\n')
        }
	cat("-Value:", format(thCoef, digits=4))
	if(x$fixedTh) cat(" (fixed)")
	cat("\n")
	cat("Proportion of points in ")
	if(nthresh==1)
		cat(paste(c("low regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)), "\n")
	else
		cat(paste(c("low regime:","\t Middle regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)), "\n")
	invisible(x)
}

summary.setar <- function(object, ...) {
	ans <- list()
	mod <- object$model.specific
	order.L <- mod$mL
	order2.L <- length(mod$ML)
	order.H <- mod$mH
	order2.H <- length(mod$MH)
	nthresh<-mod$nthresh		#numer of thresholds
	nconst<-mod$nconst		
	common<-mod$common
	ans$lowCoef <- getSetarXRegimeCoefs(object, "L")
	ans$highCoef<- getSetarXRegimeCoefs(object, "H")
	ans$thCoef <- getTh(coef(object))
	ans$fixedTh <- mod$fixedTh
	ans$externThVar <- mod$externThVar
	ans$lowRegProp <- mod$lowRegProp
	n <- getNUsed(object$str)
	coef <- object$coef[seq_len(length(object$coef)-nthresh)] #all coeffients except of the threshold
 	p <- length(coef)			#Number of slope coefficients
	resvar <- mse(object) * n / (n-p)
	Qr <- mod$qr
	p1 <- 1:p
	est <- coef[Qr$pivot[p1]]
	R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) #compute (X'X)^(-1) from the (R part) of the QR decomposition of X.
	se <- sqrt(diag(R) * resvar) #standard errors
	tval <- est/se			# t values
	coef <- cbind(est, se, tval, 2*pt(abs(tval), n-p, lower.tail = FALSE))
	dimnames(coef) <- list(names(est), c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
	ans$coef <- coef
	ans$mTh <- mod$mTh
	extend(summary.nlar(object), "summary.setar", listV=ans)
}

print.summary.setar <- function(x, digits=max(3, getOption("digits") - 2),
	signif.stars = getOption("show.signif.stars"), ...) {
	NextMethod(digits=digits, signif.stars=signif.stars, ...)
	cat("\nCoefficient(s):\n\n")
	printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)		
	cat("\nThreshold")
	cat("\nVariable: ")
        if(x$externThVar)
          cat("external")
        else {
          cat('Z(t) = ')
          cat('+ (',format(x$mTh[1], digits=2), ') X(t) ', sep="")
          if(length(x$mTh)>1)
            for(j in 1:(length(x$mTh) - 1)) {
              cat('+ (', format(x$mTh[j+1], digits=2), ') X(t-', j, ')', sep="")
            }
          cat('\n')
        }        
	cat("\nValue:", format(x$thCoef, digits=4))
	if(x$fixedTh) cat(" (fixed)")
	cat('\n')
	invisible(x)
}

plot.setar <- function(x, ask=interactive(), legend=TRUE, regSwStart, regSwStop, ...) {
	op <- par(no.readonly=TRUE)
	par(ask=ask)
	NextMethod(ask=ask, ...)
	str <- x$str
	xx <- getXX(str)
	yy <- getYY(str)
	nms <- colnames(xx)
	m <- str$m
	d <- str$d
	lags <- c((0:(m-1))*(-d), str$steps)
	xxyy <- getXXYY(str)
	x.old <- x
	x <- c(x, x$model.specific)
	series <- str$x
	z <- x$thVar
	th <- getTh(coef(x)) #x$coefficients["th"]
	nthresh<-x.old$model.specific$nthresh
	regime<-ifelse(z<=th[1],1,2)
	if(nthresh==2)
	  regime[regime==2]<-ifelse(z[regime==2]<th[2], 2,3)
	regime.id<-regime
	#regime <- factor(z <= th, levels=c(TRUE, FALSE), labels=c("low","high"))
	#regime.id <- as.numeric(regime) #1/2 vector indicating the regime
	if(length(regime)<=300) {
		pch <- c(20,23)[regime.id]
		  cex <- 1
	}
	else {
		pch <- '.'
		cex <- 4
	}
#Phase plot
	for(j in 1:m) {
		plot(xxyy[,j], xxyy[,m+1], xlab=paste("lag", -lags[j]), ylab=paste("lag", -lags[m+1]),
			col=regime.id, pch=pch, cex=cex, ...)
		lines.default(xxyy[,j], x.old$fitted, lty=2)
                if(nthresh==2)
                  title("Curently not implemented for nthresh=2!")
		if(legend)
			legend("topleft", legend=c("low","high"), pch=pch[c(1,1)], col=1:2, merge=FALSE, title="regime")
	}
#Regime switching plot
	sta <- 1
	sto <- length(regime.id)
	if(!missing(regSwStart))
		sta <- regSwStart
	if(!missing(regSwStop))
		sto <- regSwStop
	t <- sta:sto
	regime.id <- regime.id[t]
	series <- series[t+(m*d)]
	ylim <- range(series)
	l <- ylim[1] * 0.9
	h <- ylim[2] * 1.1
	ylim[1] <- ylim[1] * 0.8
	ylim[2] <- ylim[2] * 1.2
	x0 <- t
	x1 <- t+1
	y0 <- series[t]
	y1 <- series[t+1]
	par(mar=c(0,4,4,0))	#number of lines of margin to be specified on the 4 sides of the plot
	plot(t, series, type="p", ax=FALSE, ylab="time series values", main="Regime switching plot")
	if(legend)
	    legend("topright", legend=c("low",if(nthresh==2) "middle","high"), pch=pch[c(1,1)], col=1:(nthresh+1), merge=FALSE, title="regime")
	abline(h=th)
	axis(2)		#adds an axis on the left
	segments(x0,y0,x1,y1,col=regime.id)	#adds segments between the points with color depending on regime
	#plot for the transition variable
	par(mar=c(5,4,4,2))
	layout(matrix(1:2, ncol=1))
	plot1(th, nthresh,usedThVar=x$model.specific$usedThVar)
	plot2(th, nthresh,usedThVar=x$model.specific$usedThVar, trim=x$model.specific$trim)
	par(op)
	invisible(x)
}

###Exhaustive search over a grid of model parameters
selectSETAR<- function (x, m, d=1, steps=d, series, mL, mH,mM, thDelay=seq_len(m)-1, mTh, thVar, th=list(exact=NULL, int=c("from","to"), around="val"), trace=TRUE, include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL),MH=seq_len(mH), MM=seq_len(mM),nthresh=1,trim=0.15,criterion = c("pooled-AIC", "AIC", "SSR"),thSteps = 7,ngrid="ALL",  plot=TRUE,max.iter=2) 
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
### SelectSETAR 1:  Build the regressors matrix, cut Y to adequate (just copy paste of function setar() )
	include<-match.arg(include)
	model<-match.arg(model)
	if(missing(m))
		m <- max(ML, MH, thDelay+1)
  if(missing(series))
    series <- deparse(substitute(x))
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- getXX(str)
	yy <- getYY(str)
	externThVar <- FALSE
	##Lags selection
	if(missing(ML)) {		#ML: different lags
		if (missing(mL)) {	#mL: suit of lags
			mL <- m
			if (trace) 
				cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
		}
		ML <- seq_len(mL)
	}

	if(missing(MM)) {
		if (missing(mM)) {
			mM <- m
			if (trace&nthresh==2) 
				cat("Using maximum autoregressive order for middle regime: mM =", m,"\n")
		}
		MM <- seq_len(mM)
	}
	if(missing(MH)) {
		if (missing(mH)) {
			mH <- m
			if (trace) 
				cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
		}
		MH <- seq_len(mH)
	}

	###includes const, trend
	if(include=="const"){
		const <- rep(1,nrow(xx))
		inc<-"const"}
	else if(include=="trend"){
		const<-seq_len(nrow(xx))
		inc<-"trend"}
	else if(include=="both"){
		const<-cbind(rep(1,nrow(xx)),seq_len(nrow(xx)))
		inc<-c("const","trend")} 
	else {
		const<-NULL
		inc<-NULL
	}
	ninc<-length(inc)

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
allTh <- sort(unique(z[,1]))
ng <- length(allTh)
ninter<-round(trim*ng)
nmax<-ng-2*ninter

#gamma pre-specified
if(!is.null(th$exact)){
	th<-allTh[which.min(abs(allTh-th$exact))]
	if(length(th)>1){
		cat("Many values correspond to the one you gave. The first one was taken")
		th<-th[1]}
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
	if(trace)
		cat("Searching within",min(ngrid,intUp-intDown), "values between",allTh[intDown], "and", allTh[intUp],"\n")
	th<-allTh[seq(from=intDown, to=intUp, length.out=min(ngrid,intUp-intDown))]
	}
#value to search around	given by user
else if(is.numeric(th$around)){
	if(missing(ngrid))
		ngrid<-20
	if(trace)
		cat("Searching within", ngrid, "values around", th$around,"\n")
	th<-aroundGrid(th$around,allvalues=allTh,ngrid=ngrid,trim=trim)
}

#Default method: grid from lower to higher point
else{
	if(ngrid=="ALL")
		ngrid<-nmax
	else if(ngrid>nmax)
		ngrid<-nmax
	th<-allTh[round(seq(from=trim, to=1-trim, length.out=ngrid)*ng)]
	if(trace)
		cat("Searching on",ngrid, "possible threshold values within regimes with sufficient (",percent(trim*100,2),") number of observations\n")
}
# th<-round(th, getndp(x)) bad idea, rather use format in print and summary
gammas<-th

### selectSETAR 4: Sets up functions to compute the SSR/AIC/Pooled-AIC for th= 1 or 2
SSR_1thresh<- function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH,const=const,fun=buildXth1Common){
	XX<-fun(gam1,thDelay, xx,trans=z, ML=ML, MH=MH,const)
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
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
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
  SSR_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2NoCommon,trim=trim)
}
AIC_1thresh<-function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH,const=const,trim=trim,fun=buildXth1Common ){
	XX<-fun(gam1,thDelay, xx,trans=z, ML=ML, MH=MH, const)
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
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
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

AIC_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common){
  AIC_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon)
}
### selectSETAR 4b:Function pooled AIC 
pooledAIC <- function(parms) {	
	thDelayVal <- parms[1] + 1
	mLVal <- parms[2]
	mHVal <- parms[3]

	m <- max(thDelayVal, mLVal, mHVal)
	lags <- c( (seq_len(m)-1) * (-d), steps)

	xxyy <- embedd(x, lags = lags)

	z <- xxyy[, thDelayVal]
	isLow <- (z <= parms[4])

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
IDS <- as.matrix(expand.grid(thDelay,  ML, MH, th))			
colnames(IDS)<-c("thDelay", "mL", "mH", "th")
IDS2 <- as.matrix(expand.grid(thDelay, th))
colnames(IDS2)<-c("thDelay", "th")


###selectSETAR 6: Computation for 1 thresh
criterion <- match.arg(criterion)
IDS <- switch(criterion, "AIC" = IDS, "pooled-AIC" = IDS, "SSR" = IDS2)	###Selection of the grid

if (criterion == "pooled-AIC") {
    computedCriterion <- apply(IDS, 1, pooledAIC)} 
else if(criterion=="AIC"){
  if(common)
    computedCriterion <- mapply(AIC_1thresh, gam1=IDS[,4], thDelay=IDS[,1],ML=IDS[,2],MH=IDS[,3], MoreArgs=list(xx=xx,yy=yy,trans=z,const=const,trim=trim,fun=buildXth1Common))
  else
    computedCriterion <- mapply(AIC_1thresh, gam1=IDS[,4], thDelay=IDS[,1],ML=IDS[,2],MH=IDS[,3], MoreArgs=list(xx=xx,yy=yy,trans=z,const=const,trim=trim,fun=buildXth1NoCommon))
} 
else if(criterion=="SSR"){
  if(common)
    computedCriterion <- mapply(SSR_1thresh, gam1=IDS[,2],thDelay=IDS[,1],MoreArgs=list(xx=xx,yy=yy,trans=z, ML=ML, MH=MH, const=const,fun=buildXth1Common))
  else
    computedCriterion <- mapply(SSR_1thresh, gam1=IDS[,2],thDelay=IDS[,1],MoreArgs=list(xx=xx,yy=yy,trans=z, ML=ML, MH=MH, const=const,fun=buildXth1NoCommon))
}

###selectSETAR 7: sorts the results to show the best values
allres <- cbind(IDS, computedCriterion)
colnames(allres) <- c(colnames(IDS), criterion)
idSel <- sort(computedCriterion, index = TRUE)$ix
idSel <- idSel[seq_len(min(ifelse(nthresh==1,10,5), length(idSel)))]
res <- data.frame(allres[idSel, ], row.names = NULL)
bests<-c(res[1,"thDelay"],res[1,"th"])
names(bests)<-c("tDelay", "th")

###selectSETAR 8: Computation for 2 thresh 
#use function condistep (see TVARestim.r) to estimate the second th given the first
#iterate the algortihm: once the second is estimate, reestimate the first, and alternatively until convergence
if(nthresh==2){
  if(trace){
	cat("Result of the one threshold search: -Thresh: ",res[1,"th"],"\t-Delay: ",res[1,"thDelay"],"\n" )
  }
  More<-list(yy=yy, xx=xx,trans=z, ML=ML, MH=MH,MM=MM, const=const,trim=trim)
  if(criterion=="pooled-AIC")
    warning("\ncriterion pooled AIC currently not implemented for nthresh=2, use rather SSR\n")
  if(common)
    func<-switch(criterion, "AIC" = AIC_2threshCommon, "pooled-AIC" = IDS, "SSR" = SSR_2threshCommon)	
  else
    func<-switch(criterion, "AIC" = AIC_2threshNoCommon, "pooled-AIC" = IDS, "SSR" = SSR_2threshNoCommon)
#first conditional search
last<-condiStep(th,threshRef=res[1,"th"], delayRef=res[1,"thDelay"],ninter=ninter, fun=func, trace=trace, More=More)

#iterative loop for conditional search
i<-1	#initialise the loop
while(i<max.iter){
	b<-condiStep(th,last$newThresh, delayRef=res[1,1],ninter=ninter, fun=func, trace=trace, More=More)
	if(b$SSR<last$SSR){	#minimum still not reached
		i<-i+1
		last<-b}
	else{			#minimum reached
		i<-max.iter
		last<-b}
}

bests<-matrix(c(res[1,"thDelay"], min(c(last$threshRef, last$newThresh)),max(c(last$threshRef, last$newThresh))),nrow=1)
colnames(bests)<-c("tDelay", "th1","th2")


}
###selectSETAR 9: plot of the results of the grid search
if(plot==TRUE){
	allcol <- seq_len(max(thDelay+1)*max(mL)*max(mH))
	col <- switch(criterion, AIC=allcol, "pooled-AIC"=allcol,"SSR"=(thDelay+1) )
	big <- apply(expand.grid(thDelay,mL, mH),1,function(a) paste("Th:", a[1],"mL:", a[2], "mH:", a[3]))
	legend <- switch(criterion, "AIC"=big, "pooled-AIC"=big, "SSR"=paste("Threshold Delay", thDelay))

	plot(allres[,"th"], allres[,ncol(allres)], col=col, xlab="Treshold Value",ylab=criterion, main="Results of the grid search")
 	legend("topleft", pch=1, legend=legend, col=col, bg=0)
}

###selectSETAR 10: return results
if(trace)
	cat("\nNumber of combinations tested:", nrow(IDS),"\n")
return(list(res=res, bests=bests))

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

oneStep.setar <- function(object, newdata, itime, thVar, ...){
	mL <- object$model$mL
	mH <- object$model$mH
	phi1 <- object$coefficients[1:(mL+1)]
	phi2 <- object$coefficients[mL+1+ 1:(mH+1)]
	th <- object$coefficients[mL+mH+3]
	ext <- object$model$externThVar
	if(ext)	{
		z <- thVar[itime]
	}
	else {
		z <- newdata %*% object$model$mTh
		dim(z) <- NULL
	}
	z <- (z<=th)+0
	if(nrow(newdata)>1) {
	xL <- cbind(1,newdata[,1:mL])
	xH <- cbind(1,newdata[,1:mH])
	} else {
	xL <- c(1,newdata[,1:mL])
	xH <- c(1,newdata[,1:mH])
	}
	(xL %*% phi1) * z + (xH %*% phi2) * (1-z)
}

toLatex.setar <- function(object, digits=3, ...) {
	obj <- object
	res <- character()
  beta <- formatSignedNum(coefficients(obj),digits=digits,...)
  mL <- obj$model$mL
  mH <- obj$model$mH
  steps <- obj$str$steps
  d <- obj$str$d
  namesL <- paste("phi1",0:mL,sep=".")
  namesH <- paste("phi2",0:mH,sep=".")
  betaL <- beta[namesL]
  betaH <- beta[namesH]
	res[1] <- "\\["
  res[2] <- paste("X_{t+",steps,"} = \\left\\{\\begin{array}{lr}",sep="")
	res[3] <- betaL[1]
  for(j in seq_len(mL-1))
		res[3] <- paste(res[3],betaL[j+1]," X_{t-",  (j-1)*d,"}",sep="")
	res[3] <- paste(res[3],betaL[mL+1]," X_{t-",  (mL-1)*d,"}", "& Z_t \\leq ",beta["th"],"\\\\",sep="")
  res[4] <- betaH[1]
  for(j in seq_len(mH-1))
    res[4] <- paste(res[4], betaH[j+1]," X_{t-",  (j-1)*d,"} ",sep="")
	res[4] <- paste(res[4], betaH[mL+1]," X_{t-",  (mH-1)*d,"}", "& Z_t > ",beta["th"],"\\\\",sep="")
	res[5] <- "\\end{array}\\right."
	res[6] <- "\\]"
	res[7] <- ""

  if(!obj$model$externThVar) {
    mTh <- formatSignedNum(obj$model$mTh)
    m <- obj$str$m
		res[8] <- "\\["
		res[9] <- "Z_t = "
    for(j in seq_len(m)) {
      if(obj$model$mTh[j]==1)
        res[9] <- paste(res[9],"X_{t-",(j-1)*d,"} ",sep="")
      else if(obj$model$mTh[j]!=0)
        res[9] <- paste(res[9], obj$model$mTh[j]," X_{t-",(j-1)*d,"} ",sep="")
    }
		res[10] <- "\\]"
		res[11] <- ""
  }
	return(structure(res, class="Latex"))
}

showDialog.setar <- function(x, ...) {
	vD <- tclVar(1)
	vSteps <- tclVar(1)
  vML  <- tclVar(1)
  vMH  <- tclVar(1)
  vThDelay <- tclVar(0)
  vFixedTh <- tclVar(0)
  vTh <- tclVar(0)

	frTop <- Frame(opts=list(side="left"))
	frLeft <- Frame()
	add(frLeft,
		namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
		namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4)
	)
	frRight <- Frame()
	add(frRight,
		Widget(opts=list(type="checkbutton", text="fixed treshold", variable=vFixedTh)),
		namedEntry("treshold value", vTh, width=4),
		namedSpinbox("treshold delay", vThDelay, from=0, to=1000)
	)
	frExtra <- Frame()
	add(frExtra,
		namedSpinbox("time delay", vD, from=1, to=1000),
		namedSpinbox("forecasting steps", vSteps, from=1, to=1000)
	)
	add(frTop, frExtra, frLeft, frRight)
	frRoot <- Frame()	

  onFinish <- function() {
		d <- as.numeric(tclObj(vD))
		steps <- as.numeric(tclObj(vSteps))
    mL <- as.numeric(tclObj(vML))
    mH <- as.numeric(tclObj(vMH))
    fixedTh <- as.logical(tclObj(vFixedTh))
    th <- as.numeric(tclObj(vTh))
    thDelay <- as.numeric(tclObj(vThDelay))
    tkdestroy(frRoot$tkvar)
    if(fixedTh)
			res <- setar(x, d=d, steps=steps, mL=mL, mH=mH, th=th, thDelay=thDelay)
    else
			res <- setar(x, d=d, steps=steps, mL=mL, mH=mH, thDelay=thDelay)
    assign("nlarModel", res, .GlobalEnv)
  }
  onCancel <- function()
    tkdestroy(frRoot$tkvar)

	bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
	add(frRoot, frTop, bttnFrame)
	buildDialog(title="SETAR model", frRoot)
}
