## Copyright (C) 2005, 2006/2006  Antonio, Fabio Di Narzo
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

# STAR fitter
#	mTh, thDelay, thVar: as in setar
#       maxRegressors[i]: maximum number of autoregressors in regime i
#       noRegimes: number of regimes of the model
#       phi[i]: vector with the maxRegressors[i]+1 parameters of regime i
#       transFunct[i]: vector with the parameters of tr. function i.
#	trace: should infos be printed?
#	control: 'control' options to be passed to optim
star <- function(x, m, d=1, steps=d, series, noRegimes,
                  maxRegressors, noRegimes, phi, transFunct,
                  mTh, thDelay, thVar, trace=TRUE, control=list())
{

  if(missing(m))
    m <- max(maxRegressors, thDelay+1)

  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  xx <- str$xx
  yy <- str$yy
  externThVar <- FALSE
  T <- NROW(xx)

  if (missing(maxRegressors)) {
    maxRegressors <- rep(m, times = noRegimes)
    if (trace) 
      cat("Using maximum autoregressive order for all regimes: ", m,"\n")
  }

  if(!missing(thDelay)) {
    if(thDelay>=m) 
      stop(paste("thDelay too high: should be < m (=",m,")"))
    z <- xx[,thDelay+1]
  } else if(!missing(mTh)) {
    if(length(mTh) != m) 
      stop("length of 'mTh' should be equal to 'm'")
    z <- xx %*% mTh #threshold variable
    dim(x) <- NULL
  }
  else if(!missing(thVar)) {
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
  
#Automatic starting values####################
  if(missing(phi)) {
    phi <- array(rnorm(noRegimes * m+1), c(noRegimes, m+1))
    if(trace) {
      cat('Missing starting linear values. Using random values.\n')
    }
  }

  if(missing(transFunct)) {
    range <- range(z)
    transFunct[1:noRegimes,1] <-
      range[1] + ((range[2] - range[1]) / (noRegimes)) * (0:(noRegimes-1))
    transFunct[,2] <- 4
    if(trace) {
      cat('Missing starting transition values. Using uniform distribution.\n')
    }
  }
  
#############################################

  #Logistic transition function
  #y: variable
  #g: smoothing parameter
  #c: threshold value
  G <- function(y, g, c) {
    plogis(y, c, 1/g)
  }
  
  #Sum of squares function
  #p: vector of parameters
  SS <- function(p) {
    # Extract parameters from p into phi and transFunct again
    int_phi <- p[1:(noRegimes*(m+1))]
    dim(int_phi) <- c(noRegimes, m+1)
    
    int_transFunct <- p[(noRegimes*(m+1) + 1):length(p)]
    dim(int_transFunct) <- c(noRegimes,2)
    
    y.hat <- F(int_phi,int_transFunct)
    crossprod(yy - y.hat)
  }
 
  #Fitted values, given parameters
  # phi: vector of linear parameters
  # transFunct: vector of tr. functions' parameters
  F <- function(phi, transFunct) {
    local <- array(0, c(noRegimes, T))
    int_xx <- cbind(1, xx)
    for (i in 1:noRegimes) {
      local[i,] <- (int_xx %*% phi[i,]) * G(z, transFunct[i,2], transFunct[i,1])
    }

#    res <- 1:T
#    for (i in 1:T) 
#      res[i] <- sum(local[1:noRegimes,i])
    result <- apply(local, 2, sum)
    result

  }
  
#Numerical minimization##########
  p <- c(phi, transFunct)		#pack parameters in one vector
  res <- optim(p, SS, hessian = TRUE, control = control)

  if(trace)
    if(res$convergence!=0)
      cat("Convergence problem. Convergence code: ",res$convergence,"\n")
    else
      cat("Optimization algorithm converged\n")
################################
  
  #Results storing################
  res$coefficients <- res$par
#  names(res$coefficients) <- c(paste("phi1", 0:mL, sep="."),
#                               paste("phi2", 0:mH, sep="."),
#                               "gamma", "th")
  res$phi <- p[1:(noRegimes*(m+1))]
  dim(res$phi) <- c(noRegimes, m+1)
    
  res$transFunct <- p[(noRegimes*(m+1) + 1):length(p)]
  dim(res$transFunct) <- c(noRegimes,2)
    
  res$m <- m
  
  res$externThVar <- externThVar

  if(!externThVar) {
    if(missing(mTh)) {
      mTh <- rep(0,m)
      mTh[thDelay+1] <- 1
    }
    res$mTh <- mTh
  }

  res$thVar <- z
  res$fitted <- F(res$coef[1:(mL+1)], res$coef[(mL+2):(mL+mH+2)], 
                  res$coef[mL+mH+3], res$coef[mL+mH+4])
  res$residuals <- yy - res$fitted
  dim(res$residuals) <- NULL	#this should be a vector, not a matrix
  res$k <- length(res$coefficients)
  
################################

  return(extend(nlar(str, 
                     coef=res$coef,
                     fit =res$fitted,
                     res =res$residuals,
                     k   =res$k,
                     model.specific=res), "star"))
}

print.lstar <- function(x, ...) {
  NextMethod(...)
  cat("\nLSTAR model\n")
  x <- x$model.specific
  order.L <- x$mL
  order.H <- x$mH
  lowCoef <- x$coef[1:(order.L+1)]
  highCoef<- x$coef[(order.L+1)+1:(order.H+1)]
  gammaCoef <- x$coef[order.L + order.H + 3]
  thCoef <- x$coef[order.L+order.H+4]
  externThVar <- x$externThVar
  cat("Coefficients:\n")
  cat("Low regime:\n")
  print(lowCoef, ...)
  cat("\nHigh regime:\n")
  print(highCoef, ...)
  cat("\nSmoothing parameter: gamma =", format(gammaCoef, digits=4),"\n")
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
  cat("\nValue:", format(thCoef, digits=4), "\n")
  invisible(x)
}

summary.lstar <- function(object, ...) {
  ans <- list()
  # Taken from 'arma.R' in package tseries####
  coef <- object$coefficients
  n <- object$str$n.used
  rank <- qr(object$mod$hessian, 1e-07)$rank
  if(rank != length(coef)) {
    se <- rep(NA, length(coef))
    warning("singular Hessian\n")
  } else{
    di <- diag(2*object$mod$value/n*solve(object$mod$hessian))
    if(any(di < 0))
      warning("Hessian negative-semidefinite\n")
    se <- sqrt(di)
  }
  tval <- coef / se
  ans$coef <- cbind(coef, se, tval, 2 * ( 1 - pnorm( abs(tval) ) ) )
  dimnames(ans$coef) <- list(names(object$coef),
                          c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
  
############################################
  
  #Non-linearity test############
  xx <- object$str$xx
  sX <- object$mod$thVar
  dim(sX) <- NULL
  xx1<- xx*sX		#predictors set B (approximated non-linear component)
  yy <- object$str$yy
  modA <- lm(yy ~ xx)
  modB <- update(modA, . ~ . + xx1)
  a <- anova(modA, modB)
  ans$nlTest.value <- a[["F"]][2]
  ans$nlTest.pval  <- a[["Pr(>F)"]][2]
  
###############################

  order.L <- object$mod$mL
  order.H <- object$mod$mH
  ans$lowCoef <- object$coef[1:(order.L+1)]
  ans$highCoef<- object$coef[(order.L+1)+1:(order.H+1)]
  ans$thCoef <- object$coef[order.L+order.H+3]
  ans$externThVar <- object$mod$externThVar
  ans$mTh <- object$mod$mTh
  return(extend(summary.nlar(object), "summary.lstar", listV=ans))
}

print.summary.lstar <- function(x, digits=max(3, getOption("digits") - 2),
                       signif.stars = getOption("show.signif.stars"), ...)
{
  NextMethod(digits=digits, signif.stars=signif.stars, ...)
  cat("\nCoefficient(s):\n")
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
  cat("\nNon-linearity test of full-order LSTAR model against full-order AR model\n")
  cat(" F =", format(x$nlTest.value, digits=digits),"; p-value =", format(x$nlTest.pval, digits=digits),"\n")
  cat("\nThreshold ")
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
  invisible(x)
}

plot.lstar <- function(x, ask=interactive(), legend=FALSE,
                       regSwStart, regSwStop, ...) {
  
  op <- par(no.readonly=TRUE)
  par(ask=ask)
  NextMethod(ask=ask, ...)
  str <- x$str
  xx <- str$xx
  yy <- str$yy
  nms <- colnames(xx)
  z <- x$mod$thVar
  z <- plogis(z, x$coefficients["th"], 1/x$coefficients["gamma"])
  regime.id <- cut(z, breaks=quantile(z, 0:5/5), include.lowest=TRUE)
  regime.id <- as.numeric(regime.id)
  if(length(regime.id)<=300) {
    pch <- regime.id
    cex <- 1
  }
  else {
    pch <- '.'
    cex <- 4
  }
  palette <- rgb(0:4/4,0,0)
  for(j in 1:x$str$m) {
    plot(xx[,j], yy, xlab=nms[j], ylab=paste("lag",x$str$steps),
         col=palette[regime.id], pch=pch, cex=cex, ...)
    lines.default(xx[,j], x$mod$fitted, lty=2)
    if(legend) {
      labels <- c("[0;0.2]","(0.2,0.4]","(0.4;0.6]","(0.6;0.8]","(0.8;1]")
      legend("topleft", legend=labels, pch=sort(unique(regime.id)),
             col=palette[sort(unique(regime.id))], title="regime quantiles")
    }
  }
  sta <- 1
  sto <- length(regime.id)
  if(!missing(regSwStart))
    sta <- regSwStart
  if(!missing(regSwStop))
    sto <- regSwStop
  t <- sta:sto
  regime.id <- regime.id[t]
  m <- x$str$m
  d <- x$str$d
  series <- x$str$x[t+(m*d)]
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
  plot(t, series, type="n", ax=FALSE, ylab="time series values",
       main="Regime switching plot")
  axis(2)
  segments(x0,y0,x1,y1,col=palette[regime.id])
  par(op)
  invisible(x)
}

oneStep.lstar <- function(object, newdata, itime, thVar, ...){
  mL <- object$model$mL
  mH <- object$model$mH
  phi1 <- object$coefficients[1:(mL+1)]
  phi2 <- object$coefficients[mL+1+ 1:(mH+1)]
  gamma <- object$coefficients["gamma"]
  th <- object$coefficients["th"]
  ext <- object$model$externThVar
  if(ext) {
    z <- thVar[itime]
  }
  else {
    z <- newdata %*% object$model$mTh
    dim(z) <- NULL
  }
  z <- plogis(z, th, 1/gamma)
  if(nrow(newdata)>1) {
    xL <- cbind(1,newdata[,1:mL])
    xH <- cbind(1,newdata[,1:mH])
  } else {
    xL <- c(1,newdata[,1:mL])
    xH <- c(1,newdata[,1:mH])
  }
  (xL %*% phi1) * z + (xH %*% phi2) * (1-z)
}

#Exhaustive search over a grid of model parameters
#x: time series
#m: maximum autoregressive order
#d: time delay
#steps: steps ahead
selectLSTAR <- function(x, m, d=1, steps=d, mL = 1:m, mH = 1:m, thDelay=0:(m-1)) {
  op <- options(warn=-1)
  computeAIC <- function(parms) {
    mLVal <- parms[2]
    mHVal <- parms[3]
    thDelayVal <- parms[1]
    m <- max(mLVal,mHVal,thDelayVal+1)
    return(AIC( lstar(x, m=m, mL=mLVal, mH=mHVal, thDelay=thDelayVal,
                      control=list(maxit=1e2)) ) )
  }
  IDS <- as.matrix( expand.grid(thDelay, mL, mH) )
  colnames(IDS) <- c("thDelay","mL","mH")
  computedAIC <- apply(IDS, 1, computeAIC)
  options(op)
  res <- cbind(IDS, AIC = computedAIC)
  idSel <- sort(computedAIC, index=TRUE)$ix
  idSel <- idSel[1:min(10, length(idSel))]
  res <- data.frame(res[idSel,], row.names=NULL)
  return(res)
}

showDialog.lstar <- function(x, ...) {
  vML  <- tclVar(1)
  vMH  <- tclVar(1)
  vThDelay <- tclVar(0)
  vTh <- tclVar(0)
  vMaxit <- tclVar(3000)
  
  frTop <- Frame(opts=list(side="left"))
  frLeft <- Frame()
  add(frLeft,
      namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
      namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4)
      )
  frRight <- Frame()
  add(frRight,
      namedEntry("treshold value", vTh, width=4),
      namedSpinbox("treshold delay", vThDelay, from=0, to=1000)
      )
  frExtra <- Frame()
  add(frExtra,
      namedSpinbox("max iterations", vMaxit, from=100, to=1e5, increment=100, width=6)
      )
  add(frTop, frLeft, frRight, frExtra)
  frRoot <- Frame()	
  
  onFinish <- function() {
    mL <- as.numeric(tclObj(vML))
    mH <- as.numeric(tclObj(vMH))
    th <- as.numeric(tclObj(vTh))
    thDelay <- as.numeric(tclObj(vThDelay))
    maxit <- as.numeric(tclObj(vMaxit))
    tkdestroy(frRoot$tkvar)
    res <- lstar(x, mL=mL, mH=mH, th=th, thDelay=thDelay, control=list(maxit=maxit))
    assign("nlarModel", res, .GlobalEnv)
  }
  onCancel <- function()
    tkdestroy(frRoot$tkvar)
  
  bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
  add(frRoot, frTop, bttnFrame)
  buildDialog(title="LSTAR model", frRoot)

}
