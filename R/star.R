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
#       phi_1[i]: vector with the maxRegressors[i]+1 parameters of regime i
#       phi_2[i]: vector with the parameters of tr. function i.
#	trace: should infos be printed?
#	control: 'control' options to be passed to optim
star <- function(x, m, noRegimes, d=1, steps=d, series,
                  maxRegressors, phi_1, phi_2,
                  mTh, thDelay=1, thVar, trace=TRUE, control=list())
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
  if(missing(phi_1)) {
    phi_1 <- array(rnorm(noRegimes * m+1), c(noRegimes, m+1))
    if(trace) {
      cat('Missing starting linear values. Using random values.\n')
    }
  }

  if(missing(phi_2)) {
    phi_2 <- array(0, c(noRegimes, 2))
    range <- range(z)
    phi_2[1:noRegimes,1] <-
      range[1] + ((range[2] - range[1]) / (noRegimes)) * (0:(noRegimes-1))
    phi_2[,2] <- 4
    if(trace) {
      cat('Missing starting transition values. Using uniform distribution.\n')
    }
  }
  
#############################################

  #Logistic transition function
  #y: variable
  #g: smoothing parameter
  #c: threshold value
  G <- function(z, g, c) {
    if((length(c) > 1) && (length(g) > 1))
      t( apply( as.matrix(z) , 1, plogis, c, 1/g ) )
    else
      plogis(z, c, 1/g)
  }

  #Sum of squares function
  #p: vector of parameters
  SS <- function(phi_2) {
    dim(phi_2) <- c(noRegimes, 2)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 1:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, phi_2[i,2], phi_2[i,1])

    new_phi_1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(new_phi_1) <- c(noRegimes, m + 1)

    # Return the sum of squares
    y.hat <- F(new_phi_1, phi_2)
    crossprod(yy - y.hat)
  }
  
  #Fitted values, given parameters
  # phi_1: vector of linear parameters
  # phi_2: vector of tr. functions' parameters
  F <- function(phi_1, phi_2) {
    local <- array(0, c(noRegimes, T))
    int_xx <- cbind(1, xx)
    for (i in 1:noRegimes) 
      local[i,] <-
        (int_xx %*% phi_1[i,]) * G(z, phi_2[i,2], phi_2[i,1])

    result <- apply(local, 2, sum)
    result
  }
  
#Numerical optimization##########
  if(trace) 
    cat('Optimizing...')

  p <- as.vector(phi_2)
  res <- optim(p, SS, hessian = TRUE, control = control)
  phi_2 <- res$par
  dim(phi_2) <- c(noRegimes,2)

  # A last linear optimization
  tmp <- rep(cbind(1,xx), noRegimes)
  dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
  for (i in 1:noRegimes) 
    tmp[,,i] <- tmp[,,i] * G(z, phi_2[i,2], phi_2[i,1])
  
  dim(tmp) <- c(NROW(xx), (NCOL(xx) + 1) * noRegimes)
  phi_1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
  dim(phi_1) <- c(noRegimes, m+1)

  if(trace) 
    cat(' Done.\n')

  if(trace)
    if(res$convergence!=0)
      cat("Convergence problem. Convergence code: ",res$convergence,"\n")
    else
      cat("Optimization algorithm converged\n")
################################
  
  #Results storing################
  res$coefficients <- c(phi_1, phi_2)
  res$phi_1 <- phi_1
  res$phi_2 <- phi_2
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
  
  res$fitted <- F(res$phi_1, res$phi_2)
  res$residuals <- yy - res$fitted
  dim(res$residuals) <- NULL	#this should be a vector, not a matrix

  res$k <- length(res$coefficients)
  
  return(extend(nlar(str, 
                     coef=res$coef,
                     fit =res$fitted,
                     res =res$residuals,
                     k   =res$k,
                     model.specific=res), "star"))
}


oneStep.star <- function(object, newdata, itime, thVar, ...){
  m <- object$model.specific$m

  phi_1 <- object$model.specific$phi_1
  phi_2 <- object$model.specific$phi_2
  
  th <- object$coefficients["th"]
  ext <- object$model.specific$externThVar

  if(ext) {
    z <- thVar[itime]
  } else {
    z <- newdata %*% object$model$mTh
    dim(z) <- NULL
  }

  if(nrow(newdata) > 1) {
    accum <- array(0, c(noRegimes, nrow(newdata)))
    for (i in 1:noRegimes) 
      accum[i,] <-
        (cbind(1,newdata) %*% phi_1[i,]) * G(z, phi_2[i,2], phi_2[i,1])
    
    result <- apply(accum, 2, sum)
  }
  else {
    accum <- 0
    for (i in 1:noRegimes) 
      accum <- accum + 
        (c(1,newdata) %*% phi_1[i,]) * G(z, phi_2[i,2], phi_2[i,1])
    
    result <- accum
  }
  
#  z <- plogis(z, th, 1/gamma)

#  if(nrow(newdata)>1) {
#    xL <- cbind(1,newdata[,1:mL])
#    xH <- cbind(1,newdata[,1:mH])
#  } else {
#    xL <- c(1,newdata[,1:mL])
#    xH <- c(1,newdata[,1:mH])
#  }
#  (xL %*% phi1) * z + (xH %*% phi2) * (1-z)
}

