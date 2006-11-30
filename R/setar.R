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
setar <- function(x, m, d=1, steps=d, series, mL, mH, thDelay=0, mTh, thVar, th, trace=FALSE, nested=FALSE){
	if(missing(m))
		m <- max(mL, mH, thDelay+1)
  if(missing(series))
    series <- deparse(substitute(x))
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- str$xx
	yy <- str$yy
	externThVar <- FALSE
	if (missing(mL)) {
		mL <- m
		if (trace) 
			cat("Using maximum autoregressive order for low regime: mL =", m,"\n")
	}
	if (missing(mH)) {
		mH <- m
		if (trace) 
			cat("Using maximum autoregressive order for high regime: mH =", m,"\n")
	}
	if(!missing(thDelay)) {
		if(thDelay>=m) 
			stop(paste("thDelay too high: should be < m (=",m,")"))
		z <- xx[,thDelay+1]
	} else if(!missing(mTh)) {
		if(length(mTh) != m) 
			stop("length of 'mTh' should be equal to 'm'")
		z <- xx %*% mTh #threshold variable
		dim(z) <- NULL
	} else if(!missing(thVar)) {
		if(length(thVar)>nrow(xx)) {
			z <- thVar[1:nrow(xx)]
			if(trace) 
				cat("Using only first", nrow(xx), "elements of thVar\n")
		}
		else 
			z <- thVar
		externThVar <- TRUE
	} else {
		if(trace) 
			cat("Using default threshold variable: thDelay=0\n")
		z <- xx[,1]
	}
	if (missing(th)) { #if 'th' not specified, try over a reasonable grid
		th <- quantile(z, prob=c(0.15, 0.85))	#interval such that we have enough observations
		if(trace)
			cat("Searching inside threshold range: ", th[1], ", ",th[2],"\n")
		th <- seq(th[1], th[2], length=20)
		mses <- numeric(length(th))
		ress <- list()
		for(i in 1:length(th)) {
			if(externThVar)
				ress[[i]] <- Recall(x=x, m=m, d=d, steps=steps, 
					series=series, mL=mL, mH=mH, thVar=x, th=th[i], nested=TRUE)
			else if(missing(thDelay)){
				ress[[i]] <- Recall(x=x, m=m, d=d, steps=steps, 
					series=series, mL=mL, mH=mH, mTh=mTh, th=th[i], nested=TRUE)
			} else {
				ress[[i]] <- Recall(x=x, m=m, d=d, steps=steps, 
					series=series, mL=mL, mH=mH, thDelay=thDelay, th=th[i], nested=TRUE)
			}
			mses[i] <- var(ress[[i]]$residuals)
		}
		if(trace) 
			cat(" Selected threshold: ", th[which.min(mses)],"\n")
		res <- ress[[which.min(mses)]]
		res$model.specific$th <- th
		return(res)
	} else {	#else fit with the specified threshold
		isL <- 0+(z <= th)		#regime-switching indicator variable
		xxL <- cbind(1,xx[,1:mL])*isL
		xxH <- cbind(1,xx[,1:mH])*(1-isL)
		res <- lm.fit(cbind(xxL, xxH), yy)
		res$coefficients <- c(res$coefficients, th)
		names(res$coefficients) <- c(paste("phi1", 0:mL, sep="."), 
			paste("phi2", 0:mH, sep="."), "th")
		res$k <- if(nested) (res$rank+1) else res$rank	#If nested, 1 more fitted parameter: th
		res$fixedTh <- if(nested) FALSE else TRUE
		res$mL <- mL
		res$mH <- mH
		res$externThVar <- externThVar
		res$thVar <- z
    res$lowRegProp <- mean(isL)
		if(!externThVar) {
			if(missing(mTh)) {
				mTh <- rep(0,m)
				mTh[thDelay+1] <- 1
			}
			res$mTh <- mTh
		}
		return(extend(nlar(str,
			coef=res$coef,
			fit=res$fitted.values,
			res=res$residuals,
			k=res$k,
			model.specific=res), "setar"))
	}
}

print.setar <- function(x, ...) {
	NextMethod(...)
	cat("\nSETAR model (2 regimes)\n")
	x.old <- x
	x <- x$model.specific
	order.L <- x$mL
	order.H <- x$mH
	lowCoef <- x.old$coef[1:(order.L+1)]
	highCoef<- x.old$coef[(order.L+1)+1:(order.H+1)]
	thCoef <- x.old$coef[order.L+order.H+3]
	externThVar <- x$externThVar
	cat("Coefficients:\n")
	cat("Low regime:\n")
	print(lowCoef, ...)
	cat("\nHigh regime:\n")
	print(highCoef, ...)
	cat("\nThreshold")
	cat("\nVariable: ")
        if(externThVar)
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
	cat("Value:", format(thCoef, digits=4))
	if(x$fixedTh) cat(" (fixed)")
	cat("\n")
	cat("Proportion of points in low regime: ", format(x$lowRegProp*100, digits=3), "%\n", sep="")
	invisible(x)
}

summary.setar <- function(object, ...) {
	ans <- list()
	mod <- object$model.specific
	order.L <- mod$mL
	order.H <- mod$mH
	ans$lowCoef <- object$coef[1:(order.L+1)]
	ans$highCoef<- object$coef[(order.L+1)+1:(order.H+1)]
	ans$thCoef <- object$coef[order.L+order.H+3]
	ans$fixedTh <- mod$fixedTh
	ans$externThVar <- mod$externThVar
	ans$lowRegProp <- mod$lowRegProp
	n <- object$str$n.used
	coef <- object$coef[-length(object$coef)]
	p <- length(coef)
	resvar <- mse(object)*n/(n-p)
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
	xx <- str$xx
	yy <- str$yy
	nms <- colnames(xx)
	m <- str$m
	d <- str$d
	lags <- c((0:(m-1))*(-d), str$steps)
	xxyy <- cbind(xx,yy)	#construct design matrix
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
selectSETAR <- function (x, m, d=1, steps=d, thSteps = 7, mL = 1:m, mH = 1:m, 
    th = quantile(x, prob = seq(0.15, 0.85, length = thSteps)), 
    thDelay = 0:(m - 1), criterion = c("pooled-AIC","AIC")){
	str <- nlar.struct(x,m,d,steps)
	pooledAIC <- function(parms) {
		thDelayVal <- parms[1] + 1
		mLVal <- parms[3]
		mHVal <- parms[4]
		m <- max(thDelayVal, mLVal, mHVal)
		lags <- c((0:(m - 1)) * (-d), steps)
		xxyy <- embedd(x, lags = lags)
		z <- xxyy[, thDelayVal]
		isLow <- (z <= parms[2])
		if ((sum(isLow) < mLVal) | (sum(!isLow) < mHVal)) 
				return(NA)
		xx <- xxyy[isLow, 1:mLVal]
		y <- xxyy[isLow, m + 1]
		AIC1 <- AIC(lm(y ~ xx))
		xx <- xxyy[!isLow, 1:mHVal]
		y <- xxyy[!isLow, m + 1]
		AIC2 <- AIC(lm(y ~ xx))
		return(AIC1 + AIC2)
	}
	parsToModel <- function(parms) {
		thDelayVal <- parms[1]
		thVal <- parms[2]
		mLVal <- parms[3]
		mHVal <- parms[4]
		m <- max(thDelayVal+1, mLVal, mHVal)
		return(setar(x, m=m, d=d, steps=steps, mL=mLVal, mH=mHVal, th=thVal))
	}
	x <- str$x
	IDS <- as.matrix(expand.grid(thDelay, th, mL, mH))
	criterion <- match.arg(criterion)
	if(criterion=="pooled-AIC") {
		computedCriterion <- apply(IDS, 1, pooledAIC)
	} else {
		critFun <- switch(criterion, AIC=AIC)
		computedCriterion <- apply(IDS, 1, function(x) critFun(parsToModel(x)))
	}
	res <- cbind(IDS, computedCriterion)
	colnames(res) <- c("thDelay", "th", "mL", "mH", criterion)
	idSel <- sort(computedCriterion, index=TRUE)$ix
	idSel <- idSel[1:min(10, length(idSel))]
	res <- data.frame(res[idSel,], row.names=NULL)
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
  for(j in 1:(mL-1))
		res[3] <- paste(res[3],betaL[j+1]," X_{t-",  (j-1)*d,"}",sep="")
	res[3] <- paste(res[3],betaL[mL+1]," X_{t-",  (mL-1)*d,"}", "& Z_t \\leq ",beta["th"],"\\\\",sep="")
  res[4] <- betaH[1]
  for(j in 1:(mH-1))
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
    for(j in 1:m) {
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