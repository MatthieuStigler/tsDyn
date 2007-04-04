## Copyright (C) 2005, 2006, 2007/2006  Antonio, Fabio Di Narzo
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

#Fitted values, given parameters
# phi1: vector of linear parameters
# phi2: vector of tr. functions' parameters
# noRegimes
# xx: matrix from an str object
F <- function(phi1, phi2, noRegimes, xx) {
  local <- array(0, c(noRegimes, T))
  int_xx <- cbind(1, xx)
  for (i in 1:noRegimes) 
    local[i,] <-
      (int_xx %*% phi1[i,]) * G(z, phi2[i,2], phi2[i,1])
  
  result <- apply(local, 2, sum)
  result
}

# Incremental STAR fitter
#
#   Builds a STAR model with as many regimes as needed, using the
#     bottom-up strategy proposed by Teräsvirta et al.
#
#   x
#   m
#   d
#   steps
#   series
#   rob
#   sig
star <- function(x, m=3, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay=1, thVar, sig=.05, trace=TRUE, control=list(), ...)
{

  # 1. Build the nlar object and associated variables.
  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  xx <- str$xx
  yy <- str$yy
  externThVar <- FALSE
  T <- NROW(xx)

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
  
  # 2. Linearity testing
  cat("Testing linearity...\n")
  testResults <- linearityTest(str, z, rob, sig)
  pValue <- testResults$pValue;

  cat("LM Linearity test: p-value = ", pValue,"\n")
  if(testResults$isLinear) {
    increase <- ! testResults$isLinear;
    cat("The series is linear. Use linear model instead.\n")
  }
  else {
    increase <- ! testResults$isLinear;
    cat("The series is nonlinear. Incremental building procedure:\n")
  }

  # 3. Add-regime loop
  m <- 0;
  while (increase) {

    m <- m + 1;

    estimateParams();
    addRegime();
    
  }
  
}

# Estimates the parameters of a given STAR model.
#
# object: a valid (estimated) STAR model.
#
# Estimates object$model.specific$gradient
#                   object$model.specific$phi_1
#                   object$model.specific$phi_2
estimateParams <- function(object, ...)
{

  # 1.- Find promising initial values
  # 2.- Estimate nonlinear parameters
  # 3.- Estimate linear parameters
  # 4.- Compute yhat and ehat
  # 5.- Compute the gradient

}

# Find promising initial values for next regime
#
# object: a valid (estimated) STAR model with m regimes
#
# returns good starting values for gamma and th of regime noRegime.
startingValues <- function(object, ...)
{

  s_t<- object$model.specific$thVar
  noRegimes <- object$model.specific$noRegimes
  xx <- object$str$xx
  th <- object$model.specific$phi2[,1]
  gamma <- object$model.specific$phi2[,2]
  phi1 <- object$model.specific$phi1
  
  bestCost <- 999999999999;
  
  # Maximum and minimum values for gamma
  maxGamma <- 40;
  minGamma <- 1;
  rateGamma <- 5;
  
  # Maximum and minimum values for c
  minTh <- quantile(s_t, .1) # percentil 10 of s_t
  maxTh <- quantile(s_t, .9) # percentil 90 of s_t
  rateTh <- (maxTh - minTh) / 100;
  
  gamma <- 0;
  c <- 0;
  for(newGamma in seq(minGamma, maxGamma, rateGamma)) {
    for(newTh in seq(minTh, maxTh, rateTh)) {

      # We fix the linear parameters.
      fX <- array(1, c(length(s_t), noRegimes - 1))
      int_xx <- cbind(1, xx)
      tmp <- int_xx;
      for(i in 2:(noRegimes - 1)) { # leave out the first and last regime
        fX[,i - 1] <- sigmoid(gamma[i - 1] * (s_t - th[i - 1]))
        tmp <- cbind(tmp, int_xx * fX[,i - 1])
      }
      fX[,noRegimes - 1] <- sigmoid(newGamma * (s_t - newTh))
      tmp <- cbind(tmp, int_xx * fX[,noRegimes - 1])
      
      newphi1 <- lm(yy ~ . - 1, data.frame(tmp))$coefficients;
      dim(newphi1) <- dim(phi1)
      
      local <- array(0, c(noRegimes, T))
      local[1,] <- int_xx %*% newphi1[1,]
      for (i in 2:noRegimes) 
        local[i,] <-
          (int_xx %*% new_phi1[i,]) * fX[,i - 1]
      
      y.hat <- apply(local, 2, sum)
      
      cost <- crossprod(yy - y.hat)
      
      if(cost <= bestCost) {
        bestCost <- cost;
        gamma <- newGamma;
        th <- newTh;
        phi1 <- new_phi[1:(mL+1)]
        phi2 <- new_phi[(mL+2):(mL+mH+2)]
      }
    }
  }
  
  if (trace)
    cat("Starting values fixed for regime ", noRegimes, ": th = ", th,
        ", gamma = ", gamma,"; SSE = ", bestCost, "\n");

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

  result
  
}


# Tests (within the LM framework), being the null hypothesis H_0 that the
#    model 'object' is complex enough and the alternative H_1 that an extra
#    regime should be considered.
#
# object: a STAR model already built with at least 2 regimes.
#
# returns the p-values of the F statistics (from 1st to 5th order)
addRegime.star <- function(object, ...)
{

  str <- object$str
  
  T <- NROW(str$xx);  # The number of lagged samples
  
  # Build the regressand vector
  y_t <- str$yy;
  
  # Build the regressors matrix
  if (object$model.specific$externThVar) {
    x_t <- cbind(1, str$xx)
  } else x_t <- str$xx

  # Get the transition variable
  s_t <- object$model.specific$thVar

  # Get the residuals
  e_t <- object$residuals
  
  # Build the Gradient Matrix (w/o derivative terms)
  gmatrix <- x_t
  for (i in 2:object$model.specific$noRegimes) {
    gmatrix <- cbind(gmatrix, x_t *
                     G(s_t, object$model$phi_2[i,2],
                       object$model$phi_2[i,1])
#                     , phi_1 * x_t * dGdc() deriv G / deriv c
#                     , phi_1 * x_t * dGds() dG / d sigma
                   )
  }
  
  # 1. Regress the residuals on the gradient matrix, get SSR0
  regression1 <- lm(e_t ~ ., data=data.frame(gmatrix));
  SSR0 <- sum(regression1$residuals^2);
  
  # 2. Regress the residuals of regression1 on the gmatrix and on 
  aux_data1 <- data.frame(gmatrix = gmatrix, a = x_t, b = x_t * s_t);
  aux_regression1 <- lm(regression1$residuals ~ ., data=aux_data1);
  SSR1 <- sum(aux_regression1$residuals^2);
  
  # 3. Compute the first order statistic
  n <- object$k # (number of parameters under the null)
  m <- dim(aux_data1)[2] - dim(gmatrix)[2];
  F_1 <- ((SSR0 - SSR1) / m) / (SSR1 / (T - n - m));
  cat("Degrees of freedom: n = ",n,", m = ",m,"\n")
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor1 <- pf(F_1, m, T - m - n, lower.tail = FALSE);
  
  # Regress y_t on the restrictions and compute the RSS
  aux_data3 <- data.frame(gmatrix = gmatrix, a = x_t, b = x_t * s_t,
                          c = x_t * s_t^2, d = x_t * s_t^3)
  aux_regression3 <- lm(y_t ~ ., data=aux_data3)
  SSR3 <- sum(aux_regression3$residuals^2);
  
  # Compute the third order statistic
  n <- object$k;
  m <- dim(aux_data1)[2] - dim(gmatrix)[2];
#  m <- dim(aux_data3)[2] - n;
  F_3 = ((SSR0 - SSR3) / m) / (SSR3 / (T - m - n));
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor3 <- pf(F_3, m, T - m - n, lower.tail = FALSE);
  
  # Regress y_t on the restrictions and compute the RSS
  aux_data5 <- data.frame(gmatrix = gmatrix, a = x_t, b = x_t * s_t,
                          c = x_t * s_t^2, d = x_t * s_t^3,
                          e = x_t * s_t^4, d = x_t * s_t^5)
  aux_regression5 <- lm(y_t ~ ., data=aux_data5)
  SSR5 <- sum(aux_regression5$residuals^2);
  
  # Compute the fifth order statistic
  n <- object$k;
  m <- dim(aux_data1)[2] - dim(gmatrix)[2];
#  m <- dim(aux_data5)[2] - n;
  F_5 = ((SSR0 - SSR5) / m) / (SSR5 / (T - m - n));
  
  # Look up the statistic in the table, get the p-value
  lmStatTaylor5 <- pf(F_5, m, T - m - n, lower.tail = FALSE);
  
  c(firstOrderTest = lmStatTaylor1,
    thirdOrderTest = lmStatTaylor3,
    fifthOrderTest = lmStatTaylor5)

}


# LM linearity testing against 2 regime STAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   object: a star object
#   rob
#   sig
linearityTest.star <- function(object, rob=FALSE, sig=0.05, ...)
{

  str <- object$str
  T <- NROW(str$xx);  # The number of lagged samples

  # Build the regressand vector
  y_t <- str$yy;
  
  # Regressors under the null
  xH0 <- cbind(1, str$xx)

  # Get the transition variable
  s_t <- object$model.specific$thVar

  # Linear Model (null hypothesis)
  linearModel <- lm(y_t ~ . , data=data.frame(xH0))
  u_t <- linearModel$residuals;
  SSE0 <- sum(u_t^2)

  # Regressors under the alternative
  if (object$model.specific$externThVar) {
    tmp <- rep(s_t, NCOL(str$xx) + 1)
    dim(tmp) <- c(length(s_t), NCOL(str$xx) + 1)
    xH1 <- cbind(cbind(1, str$xx) * tmp, cbind(1, str$xx) * (tmp^2),
                 cbind(1, str$xx) * (tmp^3))
  } else {
    tmp <- rep(s_t, NCOL(str$xx))
    dim(tmp) <- c(length(s_t), NCOL(str$xx))
    xH1 <- cbind(str$xx * tmp, str$xx * (tmp^2), str$xx * (tmp^3))
  }

  # Standarize the regressors
  Z <- cbind(xH0, xH1);
  nZ <- NCOL(Z);
  sdZ <- sd(Z)
  dim(sdZ) <- c(1, nZ)
  sdZ <- kronecker(matrix(1,T,1), sdZ) # repeat sdZ T rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Nonlinear model (alternative hypothesis)
  nonlinearModel <- lm(u_t ~ ., data=data.frame(Z));
  e_t <- nonlinearModel$residuals;
  SSE1 <- sum(e_t^2)

  # Compute the test statistic
  n <- dim(xH0)[2];
  m <- dim(xH1)[2];
  
  F = ((SSE0 - SSE1) / m) / (SSE1 / (T - m - n));
    
  # Look up the statistic in the table, get the p-value
  pValue <- pf(F, m, T - m - n, lower.tail = FALSE);

  if (pValue < sig) {
    return(list(isLinear = TRUE, pValue = pValue));
  }
  else {
    return(list(isLinear = FALSE, pValue = pValue));
  }

}

# Predefined STAR model
#
#   Builds a STAR with a fixed number of regimes ('noRegimes') and
#     fixed parameters phi1 (linear) and phi2 (nonlinear). If no
#     parameters are given they are set to random and evenly
#     distributed, respectively.
#                            TO-DO: proper initial values (grid serch)!!
#
#   mTh, thDelay, thVar: as in setar
#   maxRegressors[i]: maximum number of autoregressors in regime i
#   noRegimes: number of regimes of the model
#   phi1[i]: vector with the maxRegressors[i]+1 parameters of regime i
#   phi2[i]: vector with the parameters of tr. function i.
#   trace: should infos be printed?
#   control: 'control' options to be passed to optim
#
star.predefined <- function(x, m, noRegimes, d=1, steps=d, series,
                  maxRegressors, phi1, phi2,
                  mTh, thDelay=1, thVar, trace=TRUE, control=list())
{

  if(noRegimes == 1) 
   stop("A STAR with 1 regime is an AR model: use the linear model instead.")
  
  if(noRegimes == 2) {
    temp <- lstar(x, m, d, steps, series, mTh, mH=m, mL=m, thDelay, thVar, control=list())
    
    return(extend(nlar(str=temp$str, 
                       coefficients = temp$model.specific$coefficients,
                       fitted.values = temp$fitted,
                       residuals = temp$residuals,
                       k   = temp$k,
                       model.specific=temp$model.specific), "star"))
  } 

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
# TO DO: grid search over phi2.
  if(missing(phi2)) {         
    phi2 <- array(0, c(noRegimes, 2))
    range <- range(z)
    phi2[1:noRegimes,1] <- # phi2[,1] = th
      range[1] + ((range[2] - range[1]) / (noRegimes)) * (0:(noRegimes-1))
    phi2[,2] <- 4                   # phi2[,2] = gamma
    if(trace) {
      cat('Missing starting transition values. Using uniform distribution.\n')
    }
  }

  if(missing(phi1)) {
    phi1 <- array(rnorm(noRegimes * m+1), c(noRegimes, m+1))
    if(trace) {
      cat('Missing starting linear values. Using random values.\n')
    }
  }

#############################################

  #Sum of squares function
  #p: vector of parameters
  SS <- function(phi2) {
    dim(phi2) <- c(noRegimes, 2)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 1:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, phi2[i,2], phi2[i,1])

    new_phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(new_phi1) <- c(noRegimes, m + 1)

    # Return the sum of squares
    y.hat <- F(new_phi1, phi2)
    crossprod(yy - y.hat)
  }
  
  #Fitted values, given parameters
  # phi1: vector of linear parameters
  # phi2: vector of tr. functions' parameters
  F <- function(phi1, phi2) {
    local <- array(0, c(noRegimes, T))
    int_xx <- cbind(1, xx)
    for (i in 1:noRegimes) 
      local[i,] <-
        (int_xx %*% phi1[i,]) * G(z, phi2[i,2], phi2[i,1])

    result <- apply(local, 2, sum)
    result
  }
  
#Numerical optimization##########
  if(trace) 
    cat('Optimizing...')

  p <- as.vector(phi2)
  res <- optim(p, SS, hessian = TRUE, control = control)
  phi2 <- res$par
  dim(phi2) <- c(noRegimes,2)

  # A last linear optimization
  tmp <- rep(cbind(1,xx), noRegimes)
  dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
  for (i in 1:noRegimes) 
    tmp[,,i] <- tmp[,,i] * G(z, phi2[i,2], phi2[i,1])
  
  dim(tmp) <- c(NROW(xx), (NCOL(xx) + 1) * noRegimes)
  phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
  dim(phi1) <- c(noRegimes, m+1)

  if(trace) 
    cat(' Done.\n')

  if(trace)
    if(res$convergence!=0)
      cat("Convergence problem. Convergence code: ",res$convergence,"\n")
    else
      cat("Optimization algorithm converged\n")
################################
  
  #Results storing################
  res$phi1 <- phi1
  res$phi2 <- phi2
  res$externThVar <- externThVar

  if(!externThVar) {
    if(missing(mTh)) {
      mTh <- rep(0,m)
      mTh[thDelay+1] <- 1
    }
    res$mTh <- mTh
  }

  res$thVar <- z
  
  fitted <- F(phi1, phi2)
  residuals <- yy - fitted
  dim(residuals) <- NULL	#this should be a vector, not a matrix

  res$noRegimes <- noRegimes
  
  return(extend(nlar(str, 
                     coef= c(phi1, phi2),
                     fit = fitted,
                     res = residuals,
                     k   = length(as.vector(phi1)) +
                                 length(as.vector(phi2)),
                     model.specific=res), "star"))
}
