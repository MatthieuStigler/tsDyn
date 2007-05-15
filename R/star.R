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

## This program is a translation of Prof. Marcelo Medeiros's Matlab codes
##    and is indebted to him.

#Logistic transition function
# y: variable
# gamma: smoothing parameter
# th: threshold value
G <- function(z, gamma, th) {
  if((length(th) > 1) && (length(gamma) > 1))
    t( apply( as.matrix(z) , 1, plogis, th, 1/gamma) )
  else
    plogis(z, th, 1/gamma)
}

# Incremental STAR fitter
#
#   Builds a STAR model with as many regimes as needed, using the
#     bottom-up strategy proposed by Teräsvirta et al.
#
#   x: the time series 
#   m: the order of the autoregressive terms
#   d: 
#   steps
#   series
#   rob
#   sig
star <- function(x, m=3, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay=1, thVar, sig=0.05, trace=TRUE, control=list(), ...)
{

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  xx <- str$xx
  yy <- str$yy
  externThVar <- FALSE
  n.used <- NROW(xx)

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

      if(trace) cat("Using only first", nrow(xx), "elements of thVar\n")

    }
    else 

      z <- thVar

    externThVar <- TRUE

  } else {

    if(trace) cat("Using default threshold variable: thDelay=0\n")
    z <- xx[,1]

  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2. Linearity testing
  if (trace) cat("Testing linearity...   ")
  testResults <- linearityTest(str, z, rob=rob, sig=sig, trace = trace)
  pValue <- testResults$pValue;
  increase <- ! testResults$isLinear;

  if(trace) cat("p-Value = ", pValue,"\n")
  if(!increase) {
    if(trace) cat("The series is linear. Use linear model instead.\n")
    return(str);
  }
  else {
    if(trace) cat("The series is nonlinear. Incremental building procedure:\n")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Build the 2-regime star
    if(trace) cat("Building a 2 regime STAR.\n")
    object <- star.predefined(x, noRegimes=2, m=m, d=d, steps=steps,
                              series=series, mTh=mTh, thDelay=thDelay,
                              thVar=thVar, trace=trace, control=control)
  
    if(trace) cat("Testing for addition of regime 3.\n");
    
    if(trace) cat("  Estimating gradient matrix...\n");
    G <- computeGradient(object);
    if(trace) cat("  Done. Computing the test statistic...\n");
    testResults <- addRegime(object, G = G, rob = rob, sig = sig, trace=trace);
    increase <- testResults$remainingNonLinearity;
    if(increase && trace) { cat("  Done. Regime 3 is needed (p-Value = ",
                                testResults$pValue,").\n"); }
    else {if(trace) cat("  Done. Regime 3 is NOT accepted (p-Value = ",
                                testResults$pValue,").\n");}
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 4. Add-regime loop
    count <- 2; # Number of nonlinear terms in the model
    while (increase) {
      
      object <- estimateParams(object, control=control, trace=trace);
      count <- count + 1;

      if(trace) cat("Testing for addition of regime ", count + 1, "...  \n");
      
      if(trace) cat("  Estimating gradient matrix...\n");
      G <- computeGradient(object);
      if(trace) cat("  Done. Computing the test statistic...\n");
      testResults <- addRegime(object, G = G, rob = rob, sig = sig, trace=trace);
      increase <- testResults$remainingNonLinearity;
      if(increase && trace) {
        cat("  Done. Regime ", count + 1," is needed (p-Value = ",
                                  testResults$pValue,").\n"); 
      }
      else {
        cat("  Done. Regime ", count + 1," is NOT accepted (p-Value = ",
                testResults$pValue,").\n");
      }            
    }
    if(trace) cat("Finished building a MRSTAR with ",count, " regimes\n");
    return(object);
  }
  
}

# Tests (within the LM framework), being the null hypothesis H_0 that the
#    model 'object' is complex enough and the alternative H_1 that an extra
#    regime should be considered.
#
# object: a STAR model already built with at least 2 regimes.
# G: the gradient matrix obtained by using computeGradient()
# rob: boolean indicating if robust tests should be used or not.
# sig: significance level of the tests
# 
# returns a list containing the p-value of the F statistic and a boolean,
#      true if there is some remaining nonlinearity and false otherwise.
addRegime.star <- function(object, G, rob=FALSE, sig=0.05, trace = TRUE, ...)
{

  e <-  object$residuals;
  s_t <- object$model.specific$thVar;
  nG <- NCOL(G);
  normG <- norm(Matrix(G * e))
  n.used <- length(e);

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
  rG <- qr(G2, tol)$rank
#  rG <- sum(svd(G2)$d > tol);

  if (normG > 1e+6) {
    if(rG < nG) {
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
#      GPCA <- PCtmp$scores;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda> 0.99999));
      GPCA <- t(PC%*%t(G));
      GPCA <- GPCA[, 1:indmin];
      b <- solve(t(GPCA) %*% GPCA) %*% t(GPCA) %*% e;
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
      b <- solve(t(G) %*% G) %*% t(G) %*% e;
      u <- e - G %*% b;
      xH0 <- G;
    }
  } else {
    u <- e;
    if(rG < nG) {
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
#      GPCA <- PCtmp$scores;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda> 0.99999));
      GPCA <- t(PC%*%t(G));
      GPCA <- GPCA[, 1:indmin];
      xH0 <- GPCA;
    } else {
      xH0 = G;
    }
  }

  SSE0 <- sum(u^2);

  # Regressors under the alternative:
  xx <- object$str$xx;
  if (object$model.specific$externThVar) {
    x_t <- cbind(1, xx);
    tmp <- rep(s_t, NCOL(x_t))
    dim(tmp) <- c(length(s_t), NCOL(x_t))
    xH1 <- cbind(x_t * tmp, x_t * (tmp^2), x_t * (tmp^3))
  } else {
    tmp <- rep(s_t, NCOL(xx))
    dim(tmp) <- c(length(s_t), NCOL(xx))
    xH1 <- cbind(xx * tmp, xx * (tmp^2), xx * (tmp^3))
  }

  Z <- cbind(xH0, xH1)

  # Standarize the regressors
  nZ <- NCOL(Z);
  sdZ <- sd(Z)
  dim(sdZ) <- c(1, nZ)
  sdZ <- kronecker(matrix(1, n.used, 1), sdZ) # repeat sdZ n.used rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Nonlinear model (Alternative hypothesis)
  c <- solve(t(Z) %*% Z) %*% t(Z) %*% u;
  v <- u - Z %*% c;
  SSE <- sum(v^2);

  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  # Compute the third order statistic
  F = ((SSE0 - SSE) / nxH1) / (SSE / (n.used - nxH0 - nxH1));

  pValue <- 1 - pf(F, nxH1, n.used - nxH0 - nxH1, lower.tail = FALSE);

  if (pValue > 1- sig) {
    return(list(remainingNonLinearity = FALSE, pValue = pValue));
  }
  else {
    return(list(remainingNonLinearity = TRUE, pValue = pValue));
  }
 
}

# Computes the gradient under the null hypothesis
#
# object: a valid STAR model.
#
# Returns a list of the gradients with respect to  the linear and
#     nonlinear parameters
computeGradient <- function(object, ...)
{
 
  noRegimes <- object$model.specific$noRegimes;
  if(noRegimes == 2) {
    gamma <- object$model.specific$phi2[1];
    th <- object$model.specific$phi2[2];
  }
  else {
    gamma <- object$model.specific$phi2[,1];
    th <- object$model.specific$phi2[,2];
  }

  phi1 <- object$model.specific$phi1;

  n.used <- NROW(object$str$xx);
  s_t<- object$model.specific$thVar;
  x_t <- cbind(1, object$str$xx);

  fX <- array(0, c(noRegimes - 1, n.used));
  dfX <- array(0, c(noRegimes - 1, n.used));
  gPhi <- x_t;
  for (i in 1:(noRegimes - 1)) {
    fX[i,] <- sigmoid(gamma[i] * (s_t - th[i]));
    dfX[i,] <- dsigmoid(fX[i,]);
    gPhi <- cbind(gPhi, kronecker(matrix(1, 1, NCOL(x_t)), fX[i,]) * x_t)
  }
  
  if (noRegimes > 2) {
    gGamma <- array(0, c(n.used, noRegimes-1));
    gTh <- array(0, c(n.used, noRegimes-1))
    for (i in 1:(noRegimes - 1)) {
      gGamma[, i] <- as.vector(x_t %*% phi1[i + 1,]) *
                                                                             as.vector(dfX[i,] * (s_t - th[i]));
      gTh[,i] <- - as.vector(x_t %*% phi1[i+1,]) *
                                                                          as.vector(gamma[i] * dfX[i,]);
    }
  } else {
    gGamma <- as.vector(x_t %*% phi1[2,]) * as.vector(dfX * (s_t - th));
    gTh <- - as.vector(x_t %*% phi1[2,]) * as.vector(gamma * dfX);
  }
  
  return(cbind(gPhi, gGamma, gTh))

}

# Estimates the parameters of a given STAR model.
#
# object: a valid STAR model.
#
# Estimates object$model.specific$phi1
#                   object$model.specific$phi2
estimateParams <- function(object, trace=TRUE, control=list(), ...)
{

  if(trace) cat("Adding regime ",object$model.specific$noRegimes + 1,"...");
  
  object$model.specific$noRegimes <- object$model.specific$noRegimes + 1;
  noRegimes <- object$model.specific$noRegimes;
  
  gamma <- object$model.specific$phi2[,1];
  th <- object$model.specific$phi2[,2];
  gamma[noRegimes - 1] <- 0.0;
  th[noRegimes - 1] <- 0.0;
     
  object$model.specific$phi2 <- cbind(gamma, th);
    
  phi1 <- object$model.specific$phi1;
  phi1 <- rbind(phi1, rnorm(3));
  object$model.specific$phi1 <- phi1;
  
  object$model.specific$coefficients <- c(phi1, cbind(gamma, th));
  object$coefficients <- c(phi1, cbind(gamma, th));
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients);
  
  if(trace) cat(" OK.\n");
  
  s_t<- object$model.specific$thVar;
  xx <- object$str$xx
  x_t <- cbind(1, xx)
  yy <- object$str$yy
  n.used <- NROW(object$str$xx);

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1.- Find promising initial values
  object <- startingValues(object, trace=trace, control=control);    
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2.- Estimate nonlinear parameters

  #Fitted values, given parameters
  # phi1: vector of linear parameters
  # phi2: vector of tr. functions' parameters
  F <- function(phi1, phi2) {
    local <- array(0, c(noRegimes, n.used))
    local[1,] <- x_t %*% phi1[1,];

    for (i in 2:noRegimes) 
      local[i,] <-
        (x_t %*% phi1[i,]) * G(s_t, gamma= phi2[i - 1,1], th= phi2[i - 1,2]);
    
    result <- apply(local, 2, sum)
    result
  }
  
  # Sum of squares function
  # phi2: vector of parameters
  SS <- function(phi2) {
    dim(phi2) <- c(noRegimes - 1, 2)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 2:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(s_t, gamma = phi2[i - 1,1], th = phi2[i - 1,2])

    newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(newPhi1) <- c(noRegimes, NCOL(xx) + 1)

    # Return the sum of squares
    y.hat <- F(newPhi1, phi2)
    crossprod(yy - y.hat)
  }
    
  if(trace) cat('Optimizing...\n')
    
  p <- as.vector(object$model.specific$phi2)

  res <- optim(p, SS, control = control)

  newPhi2 <- res$par
  dim(newPhi2) <- c(noRegimes - 1, 2)
    
  for (i in 1:(object$model.specific$noRegimes - 1))
    if(newPhi2[i,1] <0)
      newPhi2[i,1] <- - newPhi2[i,1];
  
  gamma <- newPhi2[,1]
  th <- newPhi2[,2]
  if (trace) cat("Optimized values fixed for regime ", noRegimes,
                 ": gamma = ", gamma[noRegimes - 1],
                 ", th = ", th[noRegimes - 1],"\n");

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3.- Estimate linear parameters
  tmp <- rep(cbind(1,xx), noRegimes)
  dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
  for (i in 2:noRegimes) 
    tmp[,,i] <- tmp[,,i] * G(s_t, gamma = gamma[i - 1], th = th[i - 1])
  
  newPhi1 <- lm(yy ~ . - 1, data.frame(tmp))$coefficients;
  dim(newPhi1) <- c(noRegimes, NCOL(xx) + 1);
  
  if(trace)
    if(res$convergence!=0)
      cat("Convergence problem. Convergence code: ",res$convergence,"\n")
    else
      cat("Optimization algorithm converged\n")


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 4.- Compute yhat and ehat
  local <- array(0, c(noRegimes, n.used))
  local[1,] <- x_t %*% newPhi1[1,]
  for (i in 2:noRegimes) 
    local[i,] <-
      (x_t %*% newPhi1[i,]) * sigmoid(gamma[i - 1] * (s_t - th[i - 1]))
  
  object$fitted.values <- apply(local, 2, sum);
  object$residuals <- yy - object$fitted.values;
  object$coefficients <- c(newPhi1, newPhi2)
  object$model.specific$phi1 <- newPhi1
  object$model.specific$phi2 <- cbind(gamma, th);
  
  return(object)
  
}

# Find promising initial values for next regime
#
# object: a valid (estimated) STAR model with m regimes
#
# returns a modified copy of 'object' with good starting values for
#      gamma[noRegime-1] and th[noRegime-1]. It also modifies the linear
#      parameters phi1 (as they are estimated for the new gamma and th).
startingValues <- function(object, trace=TRUE, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  gamma <- object$model.specific$phi2[,1];
  th <- object$model.specific$phi2[,2];

  s_t<- object$model.specific$thVar;
  noRegimes <- object$model.specific$noRegimes;

  xx <- object$str$xx;
  yy <- object$str$yy;
  n.used <- NROW(object$str$xx);
  
  # Maximum and minimum values for gamma
  maxGamma <- 400;
  minGamma <- 1;
  rateGamma <- 40;
  
  # Maximum and minimum values for c
  minTh <- quantile(s_t, .1) # percentil 10 of s_t
  maxTh <- quantile(s_t, .9) # percentil 90 of s_t
  rateTh <- (maxTh - minTh) / 100;
  
  bestCost <- Inf;

  if(trace) cat("Performing grid search for starting values\n");
  
  for(newGamma in seq(minGamma, maxGamma, rateGamma)) {
    for(newTh in seq(minTh, maxTh, rateTh)) {

      # We fix the linear parameters.
      x_t <- cbind(1, xx)
      tmp <- x_t;
      if (noRegimes > 2) {
        for(i in 1:(noRegimes - 2)) { # leave out the first and last regime
          tmp <- cbind(tmp, x_t * sigmoid(gamma[i] * (s_t - th[i])))
        }
      }
      tmp <- cbind(tmp, x_t * sigmoid(newGamma * (s_t - newTh)))

      newPhi1 <- lm(yy ~ . - 1, data.frame(tmp))$coefficients;
      dim(newPhi1) <- c(noRegimes, NCOL(x_t));
      
      local <- array(0, c(noRegimes, n.used))
      local[1,] <- x_t %*% newPhi1[1,]
      for (i in 2:(noRegimes - 1)) 
        local[i,] <-
          (x_t %*% newPhi1[i,]) * sigmoid(gamma[i - 1] * (s_t - th[i - 1]))
      local[noRegimes,] <- (x_t %*% newPhi1[noRegimes,]) *
        sigmoid(newGamma * (s_t - newTh))
      
      y.hat <- apply(local, 2, sum)
      
      cost <- crossprod(yy - y.hat)
      
      if(cost <= bestCost) {
        bestCost <- cost;
        gamma[noRegimes - 1] <- newGamma;
        th[noRegimes - 1] <- newTh;
        phi1 <- newPhi1;
      }
    }
  }
  
  if (trace) cat("Starting values fixed for regime ", noRegimes,
                 ": gamma = ", gamma[noRegimes - 1],
                 ", th = ", th[noRegimes - 1], 
                 "; SSE = ", bestCost, "\n");

  object$model.specific$phi2 <- cbind(gamma, th);
    
  object$model.specific$phi1 <- phi1;
  object$model.specific$coefficients <- c(phi1, cbind(gamma, th));
  object$coefficients <- c(phi1, cbind(gamma, th));
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients)
  
  return(object);
  
}


# LM linearity testing against 2 regime STAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   object: a star object
#   rob
#   sig
linearityTest.star <- function(object, rob=FALSE, sig=0.05, trace=TRUE,...)
{

  str <- object$str
  n.used <- NROW(str$xx);  # The number of lagged samples

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
  sdZ <- kronecker(matrix(1, n.used, 1), sdZ) # repeat sdZ n.used rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Nonlinear model (alternative hypothesis)
  nonlinearModel <- lm(u_t ~ ., data=data.frame(Z));
  e_t <- nonlinearModel$residuals;
  SSE1 <- sum(e_t^2)

  # Compute the test statistic
  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  F = ((SSE0 - SSE1) / nxH1) / (SSE1 / (n.used - nxH1 - nxH0));
    
  # Look up the statistic in the table, get the p-value
  pValue <- pf(F, m, n.used - m - n, lower.tail = FALSE);

  if (pValue < sig) {
    return(list(isLinear = TRUE,  pValue = pValue));
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
  
  if((noRegimes == 2) && missing(phi1) && missing(phi2)) {
    temp <- lstar(x, m=m, d=d, steps=steps, series=series, mTh=mTh,
                  mH=m, mL=m, thDelay=thDelay, thVar=thVar,
                  trace=trace, control=list())

    # Cast the lstar into a valid star object
    temp$model.specific$noRegimes <- noRegimes;
    vecLength <- NCOL(temp$str$xx); # number of linear parameters on each rule

    temp$model.specific$phi1 <- # Linear parameters
      temp$model.specific$coefficients[1:((vecLength + 1) * 2)];
    dim(temp$model.specific$phi1) <- c(2, vecLength + 1)

    temp$model.specific$phi2 <- # Non linear parameters
      temp$model.specific$coefficients[(((vecLength + 1) * 2) + 1):
                                       (((vecLength+1) * 2) + 2)];
    dim(temp$model.specific$phi2) <- c(1, 2)

    if(trace) cat("Optimized values: gamma = ",temp$model.specific$phi2[1], ", th =", temp$model.specific$phi2[2], "\n")
    
    return(extend(nlar(str=temp$str, 
                       coefficients = temp$model.specific$coefficients,
                       fitted.values = temp$fitted,
                       residuals = temp$residuals,
                       noRegimes = 2,
                       k   = temp$k,
                       model.specific = temp$model.specific), "star"))
  } 

  if(missing(m))
    m <- max(maxRegressors, thDelay+1)

  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  xx <- str$xx *3
  yy <- str$yy
  externThVar <- FALSE
  n.used <- NROW(xx)

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
    phi2 <- array(0, c(noRegimes - 1, 2))
    range <- range(z)
    phi2[,1] <- 0.5                   # phi2[,1] = gamma
    phi2[1:(noRegimes - 1),2] <- # phi2[,2] = th
      range[1] + ((range[2] - range[1]) / (noRegimes - 1)) * (0:(noRegimes-2))
    if(trace) {
      cat('Missing starting transition values. Using uniform distribution.\n')
    }
    optimize <- TRUE
  }
  else {
    optimize <- FALSE
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
    dim(phi2) <- c(noRegimes - 1, 2)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 2:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, gamma = phi2[i - 1,1], th = phi2[i - 1,2])

    newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(newPhi1) <- c(noRegimes, m + 1)

    # Return the sum of squares
    y.hat <- F(newPhi1, phi2)
    crossprod(yy - y.hat)
  }
  
  #Fitted values, given parameters
  # phi1: vector of linear parameters
  # phi2: vector of tr. functions' parameters
  F <- function(phi1, phi2) {
    x_t <- cbind(1, xx)
    local <- array(0, c(noRegimes, n.used))
    local[1,] <- x_t %*% phi1[1,];

    for (i in 2:noRegimes) 
      local[i,] <-
        (x_t %*% phi1[i,]) * G(z, gamma= phi2[i - 1,1], th= phi2[i - 1,2]);
    
    result <- apply(local, 2, sum)
    result
  }
  
#Numerical optimization##########
  res <- list()
  res$convergence <- NA
  res$hessian <- NA
  res$message <- NA
  res$value <- NA

  if(optimize) {
    if(trace) cat('Optimizing...')
    
    p <- as.vector(phi2)
    
    res <- optim(p, SS, hessian = TRUE, control = control)
    
    if(trace) 
      cat(' Nonlinear optimisation done.')

    phi2 <- res$par
    dim(phi2) <- c(noRegimes - 1,2)
    
    for (i in 1:(noRegimes - 1))
      if(phi2[i,1] <0)
        phi2[i,1] <- - phi2[i,1];

#    x_t <- cbind(1,xx);
#    local <- array(0, c(noRegimes, n.used))
#    local[1,] <- x_t %*% phi1[1,]
#    for (i in 2:noRegimes) 
#      local[i,] <-
#        (x_t %*% phi1[i,]) * sigmoid(phi2[i - 1,1] * (z - phi2[i - 1, 2]))
    
#    phi1<- lm(yy ~ . - 1, as.data.frame(local))$coefficients
#    dim(phi1) <- c(noRegimes, m+1)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 2:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, gamma = phi2[i - 1,1], th = phi2[i - 1,2])

    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(phi1) <- c(noRegimes, m + 1)

    if(trace) 
      cat(' Linear optimisation done.\n')

    if(trace)
      if(res$convergence!=0)
        cat("Convergence problem. Convergence code: ",res$convergence,"\n")
      else
       cat("Optimization algorithm converged\n")
  }
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
  
  if(!optimize) {
    res$convergence=NA
    res$hessian=NA
    res$message=NA
    res$value=NA
    res$fitted <- fitted
    res$residuals <- residuals
    dim(res$residuals) <- NULL #this should be a vector, not a matrix
    res$k <- length(as.vector(phi1)) + length(as.vector(phi2))
  }
  
  return(extend(nlar(str, 
                     coef= c(phi1, phi2),
                     fit = fitted,
                     res = residuals,
                     k   = length(as.vector(phi1)) +
                                 length(as.vector(phi2)),
                     model.specific=res), "star"))
}

oneStep.star <- function(object, newdata, itime, thVar, ...){
  m <- object$model.specific$m;
  noRegimes <- object$model.specific$noRegimes;

  phi1 <- object$model.specific$phi1
  phi2 <- object$model.specific$phi2
  
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
    accum[1,] <- (cbind(1,newdata) %*% phi1[1,]);
    for (i in 2:noRegimes) 
      accum[i,] <-
        (cbind(1,newdata) %*% phi1[i,]) * G(z, phi2[i - 1,2], phi2[i - 1,1])
    
    result <- apply(accum, 2, sum)
  }
  else {
    accum <- c(1, newdata) %*% phi1[1,];
    for (i in 2:noRegimes) 
      accum <- accum + 
        (c(1, newdata) %*% phi1[i,]) * G(z, phi2[i - 1,2], phi2[i - 1,1])
    
    result <- accum
  }

  result
  
}
