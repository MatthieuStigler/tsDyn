## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#SETAR model contructor	(sequential conditional LS)
#	th: threshold. If not specified, a grid of reasonable values is tried
#	mL: autoregressive order below the threshold ('Low')
#	mH: autoregressive order above the threshold ('High')
#	nested: is this a nested call? (useful for correcting final model df)
#	trace: should infos be printed?
setar <- function(x, m, d=1, steps=d, series, mL,mM,mH, thDelay=0, mTh, thVar, th, trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL),MM=seq_len(mM), MH=seq_len(mH), nthresh=1,trim=0.15){
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

	###Set-up of transition variable
	if(!missing(thDelay)) {
		if(thDelay>=m) 
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
		if(trace) 
			cat("Using default threshold variable: thDelay=0\n")
		z <- xx
		thDelay<-0
	}

z<-as.matrix(z)
	###Threshold search
	if (missing(th)) { #if 'th' not specified, try over a reasonable grid
		ngrid<-ifelse(nthresh==1,30,"ALL")
		search<-selectSETAR(x, m, d=d, steps=d, series, mL=mL, mH=mH,mM=mM, thDelay=thDelay, mTh, thVar, trace=trace, nested=FALSE,include = include, common=common, model=model, ML=ML,MH=MH, MM=MM,nthresh=nthresh,trim=trim,criterion = "SSR",thSteps = 7,ngrid=ngrid, plot=FALSE,max.iter=2)
		if(nthresh==1) 
			th<-c(search[1,"th"])
		else
			th<-c(search)
		nested<-TRUE
		if(trace) {
			cat("Selected threshold: ", th,"\n")
		}
		#missing(th)<-FALSE
	}

	##specified threshold 
	if(!missing(th)) {
		#check number of observations
		if(nthresh==1){
			isL <- ifelse(z[, thDelay + 1]<=th, 1, 0)
			isM<-NA
			isH <- 1-isL}	
		else{
			isL <- ifelse(z[, thDelay + 1]<=th[1], 1, 0)
			isH <- ifelse(z[, thDelay + 1]>th[2], 1, 0)
			isM <- 1-isL-isH
		}
		nobs<-na.omit(c(mean(isL),mean(isM),mean(isH)))

		if(min(nobs)<trim){
			if(trace)
				cat("\nWith the threshold you gave, there a not ",trim,"observations in each regime")
		}
		if(min(nobs)==0)
			stop("With the threshold you gave, there is a regime with no observations!")
		#build the X matrix
		if(nthresh==1){
			xxLH<-thresh1(gam1=th, thDelay, xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const,common,trim=trim)	
			midCommon<-mid<-NA}
		else if(nthresh==2){
			xxLH<-thresh2(gam1=th[1],gam2=th[2],thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const,common,trim=trim)
			midCommon<-c(paste(inc,rep(3,ninc)),paste("phi3", MH, sep="."))
			mid<-c(paste("phi3", MH, sep="."))
		}	

		##compute model
		res <- lm.fit(xxLH, yy)

		##Coefficients and names
		res$coefficients <- c(res$coefficients, th)
		if(common==FALSE){
			names(res$coefficients) <- na.omit(c(paste(inc, rep(1,ninc)), paste("phi1", ML, sep="."), paste(inc,rep(2,ninc)),paste("phi2", if(nthresh==1)MH else MM, sep="."),midCommon, rep("th",nthresh)))
		} else{
			names(res$coefficients) <- na.omit(c(inc,paste("phi1", ML, sep="."),paste("phi2",if(nthresh==1)MH else MM, sep="."),mid, rep("th",nthresh)))
		}

		res$k <- if(nested) (res$rank+nthresh) else res$rank	#If nested, 1 more fitted parameter: th
		res$fixedTh <- if(nested) FALSE else TRUE
		res$mL <- max(ML)
		res$mH <- max(MH)
		res$ML <- ML
		res$MH <- MH
		res$externThVar <- externThVar
		res$thVar <- z
		res$nconst<-inc
		res$common<-common
		res$nthresh<-nthresh
		res$model<-model
		res$RegProp <- c(mean(isL),mean(isH))
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
	}
}

getSetarXRegimeCoefs <- function(x, regime=c("1","2","3")) {
	regime <- match.arg(regime)
	x <- x$coef
	x1 <- x[grep(paste("^phi", regime, "\\.", sep=""), names(x))]
	x2 <- x[grep(paste("^const ", regime, "$", sep=""), names(x))]
	x3 <- x[grep(paste("^trend ", regime, "$", sep=""), names(x))]
	return(c(x1, x2, x3))
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
		lowCoef <- getSetarXRegimeCoefs(x.old, "1")
		highCoef<- getSetarXRegimeCoefs(x.old, ifelse(nthresh==1,"2","3"))
		cat("Low regime:\n")
		print(lowCoef, ...)
		if(nthresh==2){
			midCoef<- getSetarXRegimeCoefs(x.old, "2")
			cat("\nMid regime:\n")
			print(midCoef, ...)}
		cat("\nHigh regime:\n")
		print(highCoef, ...)
	} else {
		print(x.old$coeff[-length(x.old$coeff)], ...)
	}
	thCoef<-coef(x.old)[which(names(coef(x.old))=="th")]
	cat("\nThreshold")
	if(x$model=="MTAR"){
		cat("\nMomentum Threshold (MTAR) Adjustment")
		D<-"Diff"}
	else
		D<-NULL
	cat("\nVariable: ")
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
	cat("Value:", format(thCoef, digits=4))
	if(x$fixedTh) cat(" (fixed)")
	cat("\n")
	cat("Proportion of points in ")
	if(nthresh==1)
		cat(paste(c("low regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)))
	else
		cat(paste(c("low regime:","\t Middle regime:","\t High regime:"), percent(x$RegProp, digits=4,by100=TRUE)))
	invisible(x)
}
summary.setar <- function(object, ...) {
	ans <- list()
	mod <- object$model.specific
	order.L <- mod$mL
	order2.L <- length(mod$ML)
	order.H <- mod$mH
	order2.H <- length(mod$MH)
	nconst<-mod$nconst
	common<-mod$common
	ans$lowCoef <- getSetarXRegimeCoefs(object, "1")
	ans$highCoef<- getSetarXRegimeCoefs(object, "2")
	ans$thCoef <- coef(object)["th"]
	ans$fixedTh <- mod$fixedTh
	ans$externThVar <- mod$externThVar
	ans$lowRegProp <- mod$lowRegProp
	n <- getNUsed(object$str)
	coef <- object$coef[-length(object$coef)]
	p <- length(coef)
	resvar <- mse(object) * n / (n-p)
	Qr <- mod$qr
	p1 <- 1:p
	est <- coef[Qr$pivot[p1]]
	R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
	se <- sqrt(diag(R) * resvar)
	tval <- est/se
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
	cat("\nValue:", format(x$thCoef[1], digits=4))
	if(x$fixedTh) cat(" (fixed)")
  cat("\nProportion of points in low regime: ", format(x$lowRegProp*100, digits=3), "%\n", sep="")
	invisible(x)
}

plot.setar <- function(x, ask=interactive(), legend=FALSE, regSwStart, regSwStop, ...) {
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
	th <- x$coefficients["th"]
	regime <- factor(z <= th, levels=c(TRUE, FALSE), labels=c("low","high"))
	regime.id <- as.numeric(regime)
	if(length(regime)<=300) {
		pch <- c(20,23)[regime.id]
		cex <- 1
	}
	else {
		pch <- '.'
		cex <- 4
	}
	for(j in 1:m) {
		plot(xxyy[,j], xxyy[,m+1], xlab=paste("lag", -lags[j]), ylab=paste("lag", -lags[m+1]),
			col=regime.id, pch=pch, cex=cex, ...)
		lines.default(xxyy[,j], x.old$fitted, lty=2)
		if(legend)
			legend("topleft", legend=c("low","high"), pch=pch[c(1,1)], col=1:2, merge=FALSE, title="regime")
	}
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
	par(mar=c(0,4,4,0))
	plot(t, series, type="n", ax=FALSE, ylab="time series values", main="Regime switching plot")
	axis(2)
	segments(x0,y0,x1,y1,col=regime.id)
	par(op)
	invisible(x)
}

#Exhaustive search over a grid of model parameters
selectSETAR<- function (x, m, d=1, steps=d, series, mL, mH,mM, thDelay=seq_len(m)-1, mTh, thVar, th=list(exact=NULL, int=c("from","to"), around="val"), trace=TRUE, include = c("const", "trend","none", "both"), common=FALSE, model=c("TAR", "MTAR"), ML=seq_len(mL),MH=seq_len(mH), MM=seq_len(mM),nthresh=1,trim=0.15,criterion = c("pooled-AIC", "AIC", "SSR"),thSteps = 7,ngrid="ALL",  plot=TRUE,max.iter=2) 

{
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
	ninc<-length(include)

	###Set-up of transition variable
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

##Grid
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
#interval to search between given by user
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
th<-round(th, getndp(x))
gammas<-th

###Computation
SSR_1thresh<- function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim){
	XX<-thresh1(gam1,thDelay, xx,trans=z, ML=ML, MH=MH, MM=MM,const,common,trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#check if equal print(c(res,res2))
	}
	return(res)
}


SSR_2thresh<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim){
	XX<-thresh2(gam1,gam2,thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	#SSRres <- NA
		#res2<-deviance(lm(yy~XX-1))
		#print(c(res,res2))
	}
	return(res)
}


AIC_1thresh<-function(gam1,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim){
	XX<-thresh1(gam1,thDelay, xx,trans=z, ML=ML, MH=MH, MM=MM,const,common,trim)
	if(any(is.na(XX))){
		res<-NA}
	else{
		res <- crossprod(yy- XX %*%chol2inv(chol(crossprod(XX)))%*%crossprod(XX,yy))	
		res<-length(yy)*log(res)+2*ncol(xx)
		#res2<-AIC(lm(yy~XX-1))}
	}
	return(res)
}

AIC_2thresh<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim){
	XX<-thresh2(gam1,gam2,thDelay,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,common=common,trim=trim)
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


###Function pooled AIC
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


###Grid of combinations of all parameters
    IDS <- as.matrix(expand.grid(thDelay,  ML, MH, th))			
	colnames(IDS)<-c("thDelay", "mL", "mH", "th")
    IDS2 <- as.matrix(expand.grid(thDelay, th))
	colnames(IDS2)<-c("thDelay", "th")


###Computation for 1 thresh
criterion <- match.arg(criterion)
IDS <- switch(criterion, "AIC" = IDS, "pooled-AIC" = IDS, "SSR" = IDS2)	###Selection of the grid

if (criterion == "pooled-AIC") {
    computedCriterion <- apply(IDS, 1, pooledAIC)} 
else if(criterion=="AIC"){
   computedCriterion <- mapply(AIC_1thresh, gam1=IDS[,4], thDelay=IDS[,1],ML=IDS[,2],MH=IDS[,3], MoreArgs=list(xx=xx,yy=yy,trans=z,const=const,common=common,trim=trim))} 
else if(criterion=="SSR"){
    computedCriterion <- mapply(SSR_1thresh, gam1=IDS[,2],thDelay=IDS[,1],MoreArgs=list(xx=xx,yy=yy,trans=z, ML=ML, MH=MH, const=const,common=common,trim=trim))}

###Results
allres <- cbind(IDS, computedCriterion)
colnames(allres) <- c(colnames(IDS), criterion)
idSel <- sort(computedCriterion, index = TRUE)$ix
idSel <- idSel[seq_len(min(ifelse(nthresh==1,10,5), length(idSel)))]
res <- data.frame(allres[idSel, ], row.names = NULL)

###Computation for 2 thresh
if(nthresh==2){
if(trace){
	print(res)
	cat("Number of combinations tested:", nrow(IDS),"\n")
	cat("Result of the one threshold search: -Thresh: ",res[1,"th"],"\t-Delay: ",res[1,"thDelay"],"\n" )
}
More<-list(yy=yy, xx=xx,trans=z, ML=ML, MH=MH,MM=MM, const=const,common=common, trim=trim)
func<-switch(criterion, "AIC" = AIC_2thresh, "pooled-AIC" = IDS, "SSR" = SSR_2thresh)	

last<-condiStep(th,threshRef=res[1,"th"], delayRef=res[1,"thDelay"],ninter=ninter, fun=func, trace=trace, More=More)

i<-1
while(i<max.iter){
	b<-condiStep(th,last$newThresh, delayRef=res[1,1],ninter=ninter, fun=func, trace=trace, More=More)
	if(b$SSR<last$SSR){	#minimum still not reached
		i<-i+1
		last<-b}
	else{			#minimum reached
		i<-max.iter
		last<-b}
}

bests<-matrix(c(min(c(last$threshRef, last$newThresh)),max(c(last$threshRef, last$newThresh))),nrow=1)
colnames(bests)<-c("th1","th2")


}
###Graphical output
if(plot==TRUE){
	allcol <- seq_len(max(thDelay+1)*max(mL)*max(mH))
	col <- switch(criterion, AIC=allcol, "pooled-AIC"=allcol,"SSR"=(thDelay+1) )
	big <- apply(expand.grid(thDelay,mL, mH),1,function(a) paste("Th:", a[1],"mL:", a[2], "mH:", a[3]))
	legend <- switch(criterion, "AIC"=big, "pooled-AIC"=big, "SSR"=paste("Threshold Delay", thDelay))

	plot(allres[,"th"], allres[,ncol(allres)], col=col, xlab="Treshold Value",ylab=criterion, main="Results of the grid search")
 	legend("topleft", pch=1, legend=legend, col=col, bg=0)
}

##Result
if(nthresh==1){
	if(trace)
		cat("\nNumber of combinations tested:", nrow(IDS),"\n")
	return(res)
	}
else{
	return(bests)}
}

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



thresh1 <- function(gam1, thDelay, xx,trans, ML, MH, MM,const,common,trim) {
        isL <- ifelse(trans[, thDelay + 1]< gam1,1,0)	### isL: dummy variable
	if(common==FALSE){
		xxL <- cbind(const,xx[,ML])*isL
		xxH <- cbind(const,xx[,MH])*(1-isL)
		xxLH<-cbind(xxL,xxH)}
	else
		xxLH<-cbind(const,xx[,ML]*isL,xx[,mH]*(1-isL))
	return(xxLH)
}




thresh2<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,common,trim){
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
	if(common==FALSE){
		xxL <- cbind(const,xx[,ML])*dummydown
		xxM<-cbind(const, xx[,MM])*(1-dummydown-dummyup)
		xxH <- cbind(const,xx[,MH])*(dummyup)
		xxLMH<-cbind(xxL,xxM,xxH)
	}
	else
		xxLMH<-cbind(const,xx[,ML]*dummydown,xx[,MM]*(1-dummydown-dummyup),xx[,MH]*(dummyup))
	##computation of SSR
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
