## Copyright (C) 2005, 2006, 2007, 2008  Antonio, Fabio Di Narzo;
##                                                           José Luis Aznarte M.
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

# phi2 = c(gamma, th[1,], th[2,], ..., th[m,])

# Gaussian transition function
# z: threshold variable(s)
# gamma: smoothing parameter(s)
# th: threshold parameters, one vector for each regime
GG.gaussian <- function(z, gamma, th) {
  regimes <- length(gamma)
  n.used <- NROW(z)
  result <- array(NA, c(n.used, regimes))

  for(i in 1:regimes)
      result[,i] <- apply(exp( - gamma[i] *  t(t(z) - th[i,])^2), 1, prod)
  result
}

FF <- function(object, ...)
  UseMethod("FF")

#Fitted values, given parameters
# phi1: matrix of linear parameters
# phi2omega: vector of tr. functions' parameters
#                                     (gamma_{1...p}, th{1...p}, omega_{1...p})
FF.ncgstar <- function(xx, phi1, phi2) {
  noRegimes <- NROW(phi1)
  x_t <- cbind(1,xx)

  gamma <- phi2[1:noRegimes - 1]
  th <- phi2[noRegimes:length(phi2)]
  dim(th) <- c(noRegimes - 1, NCOL(xx))

  trfun <- cbind(1, GG.gaussian(xx, gamma, th))
  local <- (x_t %*% t(phi1)) * trfun;
  result <- apply(local, 1, sum)
  result
}

linearityTest <- function(object, ...)
  UseMethod("linearityTest")

# LM linearity testing against 2 regime NCSTAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   object: a star object
#   rob
#   sig
linearityTest.ncgstar <- function(str, rob=FALSE, sig=0.05, trace=TRUE,...)
{

  n.used <- NROW(str$xx)  # The number of lagged samples

  # Build the regressand vector
  y_t <- str$yy
  
  # Regressors under the null
  xH0 <- cbind(1, str$xx)

  # Get the transition variable
  xx <- str$xx

  # Linear Model (null hypothesis)
  linearModel <- lm(y_t ~ . , data=data.frame(xH0))
  u_t <- linearModel$residuals;
  SSE0 <- sum(u_t^2)

  # Regressors under the alternative
  xH1 <- cbind(xx^ 2, xx^3)#, xx^4)

  # Standarize the regressors
  Z <- cbind(xH0, xH1);
  nZ <- NCOL(Z);
  sdZ <- sd(Z)
  dim(sdZ) <- c(1, nZ)
  sdZ <- kronecker(matrix(1, n.used, 1), sdZ) # repeat sdZ n.used rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Compute the rank of Z
  s <- svd(Z);
  tol <- max(dim(Z)) * s[1]$d * 2.2204e-16
  rZ <- qr(Z, tol)$rank
  if(rZ < NCOL(Z)) stop("Multicollinearity problem.\n")
      
  # Nonlinear model (alternative hypothesis)
  nonlinearModel <- lm(u_t ~ ., data=data.frame(Z));
  e_t <- nonlinearModel$residuals;
  SSE1 <- sum(e_t^2)
 
  # Compute the test statistic
  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  F = ((SSE0 - SSE1) / nxH1) / (SSE1 / (n.used - nxH1 - nxH0));
    
  # Look up the statistic in the table, get the p-value
  pValue <- pf(F, nxH1, n.used - nxH0 - nxH1, lower.tail = FALSE);
  
  if (pValue >= sig) {
    return(list(isLinear = TRUE,  pValue = pValue));
  }
  else {
    return(list(isLinear = FALSE, pValue = pValue));
  }

}


computeGradient <- function(object, ...)
  UseMethod("computeGradient")

# Computes the gradient
#
# object: a valid STAR model.
#
# Returns a list of the gradients with respect to  the linear and
#     nonlinear parameters
computeGradient.ncgstar <- function(object, ...)
{
 
  noRegimes <- object$model.specific$noRegimes;
  n.used <- NROW(object$str$xx);
  x_t <- cbind(1, object$str$xx);
  xx <- object$str$xx;
  p <- NCOL(x_t)
  q <- NCOL(xx)

  phi2 <- object$model.specific$phi2
  gamma <- phi2[1:noRegimes - 1]
  th <- phi2[noRegimes:length(phi2)]
  dim(th) <- c(noRegimes - 1, q)

  # linear parameters
  phi1 <- object$model.specific$phi1;

  fX <- GG.gaussian(xx, gamma, th);
#  dfX <- GGd(xx, gamma, th);

  tsum <- array(NA, c(n.used, noRegimes - 1, q))
  tsum2 <- array(NA, c(n.used, noRegimes - 1))
  for (i in 1:(noRegimes-1)) 
    tsum[,i,] <- t(apply(xx, 1, "-", th[i,]))
  for (i in 1:(noRegimes-1)) 
    tsum2[,i] <- apply(tsum[,i,]^2, 1, sum)
  
  gPhi1 <- array(x_t, c(n.used, noRegimes * p));
  for (i in 2:noRegimes) 
    gPhi1[,(i * p - p + 1):(i * p)] <- x_t * fX[,i-1] ;

  ggamma <- array(NA, c(n.used, noRegimes-1));
  for (i in 1:(noRegimes-1))
    ggamma[,i] <- - fX[,i] * tsum2[,i];

  gth <- array(0, c(n.used, (noRegimes-1), q));
  for (i in 1:(noRegimes - 1))
    for (j in 1:q)
      gth[, i, j] <- 2 * gamma[i] * tsum[,i,j] * fX[,i]
  dim(gth) <- c(n.used, (noRegimes - 1) * q)
  
  return(cbind(gPhi1, ggamma, gth))

}

testRegime <- function(object, ...)
  UseMethod("testRegime")

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
testRegime.ncgstar <- function(object, G, rob=FALSE, sig=0.05, trace = TRUE, ...)
{

  e <-  object$residuals;
  nG <- NCOL(G);
  T <- length(e);
  xx <- object$str$xx

  normG <- norm(Matrix(t(G)) %*% e)
  norm2G <-  sum(abs(t(G) %*% e)^2)^(1/2) 
#  cat("norm2G = ", norm2G, "normG = ", normG, "\n");

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
#  rG <- sum(svd(G2)$d > tol); # Alternative methods
  rG <- qr(G2, tol)$rank
#  cat("rango G = ", rG, "\n");

  if (normG > 1e-6) {
    if(rG < nG) {
#      cat("A1\n")
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda > 0.99999));
      GPCA <- t(PC %*% t(G));
      GPCA <- GPCA[, 1:indmin];
#      b <- lm.fit(x = GPCA, y = e)$coefficients
      b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
#      cat("A2\n")
#      b <- lm.fit(x = G, y = e)$coefficients;
      b <- lm(e ~ . -1, data=data.frame(G))$coefficients;
      u <- e - G %*% b;
      xH0 <- G;
    }
  } else {
    u <- e;
    if(rG < nG) {
#     cat("B1\n")
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
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
  xH1 <- cbind(xx^2, xx^3)#, xx^4)
  
  Z <- cbind(xH0, xH1)

  # Standarize the regressors
  nZ <- NCOL(Z);
  sdZ <- sd(Z)
  dim(sdZ) <- c(1, nZ)

  sdZ <- kronecker(matrix(1, T, 1), sdZ) # repeat sdZ T rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Compute the rank of Z
  s <- svd(Z);
  tol <- max(dim(Z)) * s[1]$d * 2.2204e-16
  rZ <- qr(Z, tol)$rank
  if(rZ < NCOL(Z)) stop("Multicollinearity problem.\n")

 # Nonlinear model (Alternative hypothesis)
  c <- lm.fit(x = Z, y = u)$coefficients
#  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
  dim(c) <- c(NCOL(Z), 1);
  v <- u - Z %*% c;
  SSE <- sum(v^2);

  # Compute the third order statistic
  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  F = ((SSE0 - SSE) / nxH1) / (SSE / (T - nxH0 - nxH1));

  pValue <- pf(F, nxH1, T - nxH0 - nxH1, lower.tail = TRUE);

  if (pValue >= sig) {
    return(list(remainingNonLinearity = TRUE, pValue = pValue));
  }
  else {
    return(list(remainingNonLinearity = FALSE, pValue = pValue));
  }
}

addRegime <- function(object, ...)
  UseMethod("addRegime")

addRegime.ncgstar <- function(object)
{

  noRegimes <- object$model.specific$noRegimes
  xx <- object$str$xx
  phi2 <- object$model.specific$phi2

  if(noRegimes == 1) {

    gamma <- 0
    th <- rnorm(NCOL(xx))
    
    phi2 <- c(gamma, th);

    phi1 <-   object$model.specific$phi1
    phi1 <-   rbind(object$model.specific$phi1, rnorm(NCOL(xx) + 1));
    
  }
  else {
    # phi2 contains parameters gamma and th for each regime
    gamma <- phi2[1:noRegimes - 1]
    th <- phi2[noRegimes:length(phi2)]
    dim(th) <- c(noRegimes - 1, NCOL(xx))
      
    gamma[noRegimes] <- 0.0; # Add one gamma
    th <- rbind(th, rnorm(NCOL(xx)));        # Add one th

    phi2 <- c(gamma, th);

    # phi1 (one vector for each regime)
    phi1 <-   object$model.specific$phi1
    phi1 <-   rbind(object$model.specific$phi1, rnorm(NCOL(xx) + 1));
  }

  # Update the object
  object$model.specific$noRegimes <- object$model.specific$noRegimes + 1;
  object$noRegimes <- object$noRegimes + 1;
  
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2 <- phi2
  
  object$model.specific$coefficients <- c(phi1, phi2);
  object$coefficients <- c(phi1, phi2);
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients);
  
  return(object);
  
}

startingValues <- function(object, ...)
  UseMethod("startingValues")

startingValues.ncgstar <- function(object, trace=TRUE, svIter, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  xx <- object$str$xx;
  x_t <- cbind(1, xx)
  yy <- object$str$yy;
  n.used <- NROW(object$str$xx);
  p <- NCOL(x_t)
  q <- NCOL(xx)
  local1 <- array(x_t,  c(n.used * NCOL(x_t), noRegimes))
  local2 <- array(0, c(NROW(xx), noRegimes))
  
  phi2 <- object$model.specific$phi2
  gamma <- phi2[1:noRegimes - 1]
  th <- phi2[noRegimes:length(phi2)]
  dim(th) <- c(noRegimes - 1, NCOL(xx))
  
  # Fix 'svIter' random starting values for th
#  newTh <- rnorm(svIter * q, mean=mean(xx), sd=sd(xx))
#  dim(newTh) <- c(svIter, q);

  rango <- range(xx)
  div <- seq(rango[1], rango[2], length.out = svIter^(1/q))
  newTh <- as.matrix( do.call( expand.grid, rep( list(div), q ) ) )

  max.iter <- NROW(newTh)
  
  maxGamma <- 40; # abs(8 / ((max(xx %*% newOmega) - newTh)))
  minGamma <- 1; # abs(1 / ((min(xx %*% newOmega) - newTh)));
  rateGamma <- 2; # (maxGamma - minGamma) / 20;     
  bestCost <- Inf;

  for(i in 1:max.iter) { 

    if ((i %% 25 == 0) && trace) cat(".")

    th[noRegimes - 1,] <- newTh[i,]
        
    for(newGamma in seq(minGamma, maxGamma, rateGamma)) {

      gamma[noRegimes - 1] <- newGamma;

      trfun <- cbind(1, GG.gaussian(xx, gamma, th))
      
      # We fix the linear parameters here, before optimizing the nonlinear.
      local1 <- apply(trfun, 2, "*", x_t)
      dim(local1) <- c(n.used, noRegimes * p)
      newPhi1<- lm.fit(x = local1, y = yy)$coefficients
      dim(newPhi1) <- c(NCOL(xx) + 1, noRegimes)
#      newPhi1 <- t(newPhi1)

      local2 <- (x_t %*% newPhi1) * trfun;
 
      y.hat <- apply(local2, 1, sum)
      cost <- crossprod(yy - y.hat) / sqrt(n.used);
#      cat("cost for iter ",i,": ",cost,"\n")
#      cat("Phi1:\n")
#      print(t(newPhi1))
#      cat("gamma:\n")
#      print(gamma)
#      cat("th:\n")
#      print(th)

      if(! is.na(cost)) {
        if(cost < bestCost) {
#          cat("Updating bestCost, iter ",i,"\n")
          bestCost <- cost;

          bestGamma <- gamma
          bestTh <- th
          bestPhi1 <- newPhi1
        }
      }
    }
  }

  gamma <- bestGamma
  th <- bestTh
  phi1 <- bestPhi1
  
  if (trace) cat("\n  Starting values fixed for regime ", noRegimes,
                 ":\n\tgamma = ", gamma[noRegimes - 1],
                 "\n\tth = ", th[noRegimes - 1,], "\n");
  
  object$model.specific$phi2 <- c(gamma, th)
  object$model.specific$phi1 <- phi1;
  object$model.specific$coefficients <- c(phi1, gamma, th);
  object$coefficients <- c(phi1, gamma, th);
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients)
  
  return(object);
  
}

estimateParams <- function(object, ...)
  UseMethod("estimateParams")

# Estimates the parameters of a given STAR model.
#
# object: a valid STAR model.
#
# Estimates object$model.specific$phi1
#                   object$model.specific$phi2
estimateParams.ncgstar <- function(object, trace=TRUE,
                                  control=list(), alg="LM", cluster=NULL, ...)
{

  xx <- object$str$xx
  x_t <- cbind(1, xx)
  yy <- object$str$yy
  n.used <- NROW(object$str$xx);
  noRegimes <- object$model.specific$noRegimes;
  q <- NCOL(xx)
  p <- NCOL(x_t)

  # Function to compute the gradient 
  #
  # Returns the gradient with respect to the error
  gradEhat.ncgstar <- function(par) {
    gamma <- par[1:noRegimes - 1]
    th <- par[noRegimes:length(phi2)]
    dim(th) <- c(noRegimes - 1, q)
    
    # We fix the linear parameters because they're not known here
    trfun <- GG.gaussian(xx, gamma, th)
    tmp <- apply(cbind(1, trfun), 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
    phi1 <- t(phi1)

    # Compute the gradients
    fX <- GG.gaussian(xx, gamma, th);

    tsum <- array(NA, c(n.used, noRegimes - 1, q))
    for (i in 1:(noRegimes-1)) 
      tsum[,i,] <- t(apply(xx, 1, "-", th[i,]))

    tsum2 <- array(NA, c(n.used, noRegimes - 1))
    for (i in 1:(noRegimes-1)) 
      tsum2[,i] <- apply(tsum[,i,]^2, 1, sum)
    
    ggamma <- array(NA, c(n.used, noRegimes-1));
    for (i in 1:(noRegimes-1))
      ggamma[,i] <- - fX[,i] * tsum2[,i];
    
    gth <- array(NA, c(n.used, (noRegimes-1), q));
    for (i in 1:(noRegimes - 1))
      for (j in 1:q)
        gth[, i, j] <- 2 * gamma[i] * tsum[,i,j] * fX[,i]
    dim(gth) <- c(n.used, (noRegimes - 1) * q)
      
    J = - cbind(ggamma, gth) / sqrt(n.used)

    local <- (x_t %*% t(phi1)) * cbind(1, trfun);
    y.hat <- apply(local, 1, sum)
    e.hat <- yy - y.hat;
    return(2 * t(e.hat) %*% J)
    
  }

  #Sum of squares function
  #p: vector of parameters
  SS.ncgstar <- function(par) {
    gamma <- par[1:noRegimes - 1]
    th <- par[noRegimes:length(par)]
    dim(th) <- c(noRegimes - 1, q)
    
    # We fix the linear parameters because they're not known here
    trfun <- cbind(1, GG.gaussian(xx, gamma,th))
    tmp <- apply(trfun, 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
#    phi1 <- t(phi1)

    # Return the sum of squares
    local <- (x_t %*% phi1) * trfun;
    y.hat <- apply(local, 1, sum)
    
    return(crossprod(yy - y.hat))

  }

  # Function to compute the gradient 
  #
  # Returns the gradient with respect to the error
  gradEhat.lm <- function(par) {
    gamma <- par[1:noRegimes - 1]
    th <- par[noRegimes:length(phi2)]
    dim(th) <- c(noRegimes - 1, q)
    
    # We fix the linear parameters because they're not known here
    trfun <- GG.gaussian(xx, gamma, th)
    tmp <- apply(cbind(1, trfun), 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
    phi1 <- t(phi1)

    # Compute the gradients
    fX <- GG.gaussian(xx, gamma, th);
    tsum <- array(NA, c(n.used, noRegimes - 1, q))
    tsum2 <- array(NA, c(n.used, noRegimes - 1))
    for (i in 1:(noRegimes-1)) 
      tsum[,i,] <- t(apply(xx, 1, "-", th[i,]))
    for (i in 1:(noRegimes-1)) 
      tsum2[,i] <- apply(tsum[,i,]^2, 1, sum)
    
    ggamma <- array(0, c(n.used, noRegimes-1));
    for (i in 1:(noRegimes-1))
      ggamma[,i] <- - fX[,i] * tsum2[,i];
    
    gth <- array(0, c(n.used, (noRegimes-1), q));
    for (i in 1:(noRegimes - 1))
      for (j in 1:q)
        gth[, i, j] <- 2 * gamma[i] * tsum[,i,j] * fX[,i]
    dim(gth) <- c(n.used, (noRegimes - 1) * q)
      
    J = - cbind(ggamma, gth) / sqrt(n.used)

    return(J)
    
  }

  #Sum of squares function
  #p: vector of parameters
  SS.lm <- function(par) {
    gamma <- par[1:noRegimes - 1]
    th <- par[noRegimes:length(phi2)]
    dim(th) <- c(noRegimes - 1, q)
    
    # We fix the linear parameters because they're not known here
    trfun <- cbind(1, GG.gaussian(xx, gamma,th))
    tmp <- apply(trfun, 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
#    phi1 <- t(phi1)

    # Return the sum of squares
    local <- (x_t %*% phi1) * trfun;
    y.hat <- apply(local, 1, sum)
    
    return((yy - y.hat) / sqrt(n.used))

  }

  # Reorder the regimes according to the values of th
  if (noRegimes > 2) {
    if(trace) cat("  Reordering regimes...\n")

    ordering <-  apply(th, 2, order)[,1]
    
    th <- th[ordering,]
    gamma <- gamma[ordering]
    phi1 <- phi1[ordering,]
  }

  phi2 <- object$model.specific$phi2;

  domainsGamma <- t(matrix(rep(c(0,100), noRegimes - 1), c(2,noRegimes - 1)))
  domainsTh <- t(matrix(rep(c(min(xx) - sd(xx[,1]), max(xx) + sd(xx[,1])),
                            q * (noRegimes - 1)), c(2, q * (noRegimes - 1))))
  domainsGA <- rbind(domainsGamma, domainsTh)

  if (alg == "LM") {
    
    res <- nls.lm(phi2, fn=SS.lm, jac = gradEhat.lm,
                  control=nls.lm.control(nprint=100, maxfev=1e+9, maxiter=4000));

  } else if (alg == "BFGS") {

    res <- optim(phi2, fn=SS.ncgstar, gr = gradEhat.ncgstar,
                 method="BFGS",
                  control = list(trace=0, maxit=1000));

  } else if (alg == "SANN") {

    res <- optim(phi2, fn=SS.ncgstar, method="SANN",
                  control = list(trace=10, maxit=20000, temp=200));

  } else if (alg == "GAD") {

    res <- genoud(fn = SS.ncgstar, starting.values=phi2, nvars=length(phi2),
                  gr = gradEhat.ncgstar, cluster=cluster,
                  control=list(maxit=1000), print.level = 1,
                  pop.size=1000, max.generations=200,
                  wait.generations=50, Domains = domainsGA)

  } else if (alg == "GA") {

    res <- genoud(fn = SS.ncgstar, starting.values=phi2, nvars=length(phi2),
                  BFGS=FALSE, cluster=cluster, print.level = 1,
                  pop.size=1000, max.generations=200,
                  wait.generations=50, Domains = domainsGA)

  }
  
  gamma <- res$par[1:noRegimes - 1]
  th <- res$par[noRegimes:length(phi2)]
  dim(th) <- c(noRegimes - 1, q)

  gamma <- abs(gamma) # We don't want negative values for gamma.
  
  if (trace) cat("  Optimized values fixed for regime ", noRegimes,
                 ": \n    gamma = ", gamma,
                 "\n    th = ", th, "\n");

#  if(trace) 
#    if(res$convergence != 0)
#      cat("  *** Convergence problem. Code: ",res$convergence,"\n")
#    else {
#      cat("  Optimization algorithm converged\n")
#    }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3.- Estimate linear parameters again
  trfun <- GG.gaussian(xx, gamma, th)
  tmp <- apply(cbind(1, trfun), 2, "*", x_t)
  dim(tmp) <- c(n.used, noRegimes * p)
  phi1 <- lm.fit(x = tmp, y = yy)$coefficients
  dim(phi1) <- c(p, noRegimes)
  phi1 <- t(phi1)

  if(trace) {
    cat("  Optimized linear values: \n");
    print(phi1);
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 4.- Compute yhat and ehat
  
  object$fitted.values <- FF.ncgstar(xx, phi1, c(gamma, th)); # y.hat
  object$residuals <- yy - object$fitted.values; # e.hat
  object$coefficients <- c(phi1, gamma, th)
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2 <- c(gamma, th);
  
  return(object)
  
}

# Incremental STAR fitter
#
#   Builds a STAR model with as many regimes as needed, using the
#     bottom-up strategy proposed by Terï¿½svirta et al.
#
#   x: the time series 
#   m: the order of the autoregressive terms
#   d: 
#   steps
#   series
#   rob
#   sig
ncgstar <- function(x, m=2, noRegimes, d = 1, steps = d, series, 
                   mTh, thDelay, thVar, sig=0.05, trace=TRUE, svIter = 1000,
                   cluster= NULL, alg="LM", control=list(), ...)
{

#  debug(startingValues.ncgstar)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))   series <- deparse(substitute(x))

  # Normalize the series.
#  if ((mean(x) >= 0.05) && (sd(x) >= 0.05) {
#    if (trace) cat("Normalizing the series.\n")
#    x <- (x - mean(x)) / sd(x);
#  }
  
  if(!is.null(cluster)) {
    setDefaultClusterOptions(master="localhost", port=10187)
    cl <- makeCluster(cluster, type="SOCK")
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(Matrix))
  }

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  cl <- makeCluster(cluster, "SOCK")
  
  xx <- str$xx
  x_t <- cbind(1, xx)
  yy <- str$yy

  externThVar <- FALSE
  n.used <- NROW(xx)

  if(missing(noRegimes))
     noRegimes <- Inf;     
    
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2. Linearity testing
  if (trace) cat("Testing linearity...   ")
  testResults <- linearityTest.ncgstar(str, rob=rob, sig=sig, trace = trace)
  pValue <- testResults$pValue;
  increase <- ! testResults$isLinear;

  if(trace) cat("p-Value = ", pValue,"\n")
  if(!increase) {
    if(trace) cat("The series is linear.\n")
    # Build the linear model
    linearModel <- lm.fit(x = x_t, y = yy);
#    linearModel <- lm(yy ~ .  - 1, data = data.frame(x_t));

    phi1= linearModel$coefficients
    dim(phi1) = c(1, NCOL(x_t)) # one row, one regime
    
    # Create the ncstar object
    list = list(noRegimes=1, m = m, fitted = linearModel$fitted.values,  
                     residuals = linearModel$residuals,
                     coefficients = linearModel$coefficients,
                     k = length(linearModel$coefficients),
                     convergence=NA, counts=NA, hessian=NA, message=NA,
                     value=NA, par=NA,
                     phi2 = NA, phi1= phi1) 
    object <- extend(nlar(str, coefficients = linearModel$coefficients,
                          fitted = linearModel$fitted.values,
                          residuals = linearModel$residuals,
                          k=NCOL(x_t), noRegimes=1,
                          model.specific=list),
                     "ncgstar")
    return(object);
  }
  else {
    if(trace) cat("The series is nonlinear. Incremental building procedure:\n")
    
    # Build the linear model
    linearModel <- lm.fit(x=x_t, y = yy);

    phi1= linearModel$coefficients
    dim(phi1) = c(1, NCOL(x_t)) # one row, one regime
    
    # Create the ncstar object
    list = list(noRegimes=1, m = m, fitted = linearModel$fitted.values,  
                     residuals = linearModel$residuals,
                     coefficients = linearModel$coefficients,
                     k = length(linearModel$coefficients),
                     convergence=NA, counts=NA, hessian=NA, message=NA,
                     value=NA, par=NA,
                     phi2 = NA, phi1= phi1) 
    object <- extend(nlar(str, coefficients = linearModel$coefficients,
                          fitted = linearModel$fitted.values,
                          residuals = linearModel$residuals,
                          k=NCOL(x_t), noRegimes=1,
                          model.specific=list),
                     "ncgstar")
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Add-regime loop
    while (increase && (object$model.specific$noRegimes < noRegimes)) {
      
      nR <- object$model.specific$noRegimes + 1
      
      if(trace) cat("- Adding regime ", nR, ".\n");
      object <- addRegime.ncgstar(object);
      
      if(trace) cat("- Fixing good starting values for regime ", nR);
      if(is.null(cluster)) {
        object <- startingValues.ncgstar(object, trace=trace, svIter = svIter);
      } else {
        if(trace) cat("\n   + Doing distributed computations... ")
        solutions <- clusterCall(cl, startingValues.ncgstar, object, trace=trace,
                                 svIter = svIter %/% length(cluster))
        cost <- rep(Inf, length(solutions))
        if(trace) cat("\n   + Gathering results...\n")
        for (i in 1:length(solutions)) {
          cost[i] <- sum(solutions[[i]]$residuals^2) / sqrt(n.used)
        }
        
        object <- solutions[[which.min(cost)]]

        if (trace) cat("\n  Starting values fixed for regime ", nR,
                       ":\n", object$model.specific$phi2, "\n");
#                       ":\n\tgamma = ", object$model.specific$phi2omega[noRegimes - 1],
#                       "\n\tth = ", object$model.specific$phi2omega[(noRegimes - 1) * 2],
#                       "\n\tomega = ",
#                       object$model.specific$phi2omega[((noRegimes - 1) * 2 + 1):((noRegimes - 1) * 2 + NCOL(xx))], "\n");
      }

      if(trace) cat('- Estimating parameters of regime', nR, ' (alg: ', alg, ')...\n')
      object <- estimateParams(object, alg=alg, cluster=cluster,
                               control=control, trace=trace);
      
      if(trace) cat("\n- Testing for addition of regime ", nR + 1, ".\n");
      if(trace) cat("  Estimating gradient matrix...\n");
      G <- computeGradient.ncgstar(object);

#      sig <- sig / 2 # We halve the significance level
      if(trace) cat("  Computing the test statistic (sig = ", 1 - 0.05 / nR, ")...\n");
      testResults <- testRegime(object, G = G,
                                       rob = rob, sig = 1 - 0.05 / nR, trace=trace);

      increase <- testResults$remainingNonLinearity;
      if(increase) {
        if(trace) cat("  Regime ", nR + 1, " is needed (p-Value = ",
                      testResults$pValue,").\n"); 
      }
      else {
        if(trace) cat("  Regime ", nR + 1, " is NOT accepted (p-Value = ",
                      testResults$pValue,").\n");
      }            
    }
    
    if(trace) cat("\n- Finished building an NCGSTAR with ",
                  object$model.specific$noRegimes, " regimes\n");
    return(object);
  }
  
}

# Predefined NCGSTAR model
#
#   Builds a NCGSTAR with a fixed number of regimes ('noRegimes') and
#     fixed parameters phi1 (linear) and phi2 (nonlinear). If no
#     parameters are given they are set to random and evenly
#     distributed, respectively.
#
#   mTh, thDelay, thVar: as in setar
#   maxRegressors[i]: maximum number of autoregressors in regime i
#   noRegimes: number of regimes of the model
#   phi1[i,]: vector with the maxRegressors[i]+1 linear parameters of regime i
#   phi2: vector with the nonlinear parameters of transition 
#          function i, c(gamma, th)
#   trace: should infos be pinted?
#   control: 'control' options to be passed to optim
#
ncgstar.predefined <- function(x, m=2, noRegimes, d = 1, steps = d, series,
                              mTh, thDelay, thVar, phi1, phi2,
                              alg="BFGS", cluster= NULL, tests=F,
                              sig=0.05, trace=TRUE, svIter = 1000, control=list(), ...)
{

  if(!is.null(cluster)) {
    setDefaultClusterOptions(master="localhost", port=10187)
    cl <- makeCluster(cluster, type="SOCK")
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(Matrix))
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))   series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  xx <- str$xx
  x_t <- cbind(1, xx)
  yy <- str$yy

  externThVar <- FALSE
  n.used <- NROW(xx)

  if(missing(noRegimes)) {
    if(!missing(phi1)) {
      noRegimes <- NROW(phi1)
    }
    else{
     cat("Number of regimes not provided.");
     exit(1);
   }
  }
    
   # Create the ncgstar object
  list = list(noRegimes=noRegimes, m = m, fitted = NA,  
    residuals = NA,
    coefficients = c(phi1, phi2),
    k = length(c(phi1, phi2)),
    convergence=NA, counts=NA, hessian=NA, message=NA,
    value=NA, par=NA,
    phi2 = phi2, phi1= phi1) 
  object <- extend(nlar(str, coefficients = c(phi1, phi2),
                        fitted = NA,
                        residuals = NA,
                        k =length(c(phi1, phi2)),  noRegimes=noRegimes,
                        model.specific=list),
                   "ncgstar")

#      cat("\nStarting values for linear parameters:\n    ")
#      print(phi1_median)
      
#      cat("\nStarting values for nonlinear parameters:\n")
#      cat("gamma = ", phi2_median[1:(noRegimes - 1)])
#      cat("; th = ", phi2_median[noRegimes:(2*(noRegimes - 1))])
#      cat("\nomega=", phi2_median[(2 * (noRegimes - 1) + 1):length(phi2_median)])

  if(trace) cat('- Estimating parameters for all ', noRegimes, ' regimes, (alg: ', alg, ')...\n')
  object <- estimateParams.ncgstar(object, control=control, trace=trace,
                                  alg=alg, cluster=cluster);    
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 4. Diagnostic checking
    
  if(tests) {
    if(trace) cat("== Testing for linear independence of the residuals:\n")

    isIID <- array(NA, 12)
    pValue <- array(NA, 12)
    for(r in 1:12) {
      test1 <- testIID.ncgstar(object, G, r, rob=rob, sig=sig, trace = trace)
      isIID[r] <- test1$isIID
      pValue[r] <- test1$pValue
    }
    object$model.specific$testIID <- data.frame(isIID, pValue)
    if(trace) print(data.frame(isIID, pValue))
    
    if(trace) cat("== Testing for constant variance of the residuals:\n")
    test2 <- testConstVar.ncgstar(object, G, rob=rob, sig=sig, trace = trace)
    object$model.specific$testConstVar <- test2; 
    
    if(test2$isConstVar) 
      cat("      Constant variance detected, pValue = ", test2$pValue, "\n")
    else cat("     Smoothly changing variance detected, pvalue = ", test2$pValue, "\n")
    
    if(trace) cat("== Testing for parameter constancy:\n")
    test3 <- testParConst.ncgstar(object, G, rob=rob, sig=sig, trace = trace)
    object$model.specific$testParConst <- test3;
    
    if(test3$isConstVar) 
      cat("      Constant parameters detected, pValue = ", test3$pValue, "\n")
    else cat("     Smoothly changing parameters detected, pvalue = ", test3$pValue, "\n")
  }
  
  if(trace) cat("\n- Finished building an NCGSTAR with ",
                object$model.specific$noRegimes, " regimes\n");
  
  if(!is.null(cluster)) stopCluster(cl)
  return(object);

}

oneStep.ncgstar <- function(object, newdata, itime, thVar, ...)
{

  phi1 <- object$model.specific$phi1;
  phi2 <- object$model.specific$phi2;

  noRegimes <- NROW(phi1)
  x_t <- cbind(1,xx)

  gamma <- phi2[1:noRegimes - 1]
  th <- phi2[noRegimes:length(phi2)]
  dim(th) <- c(noRegimes - 1, NCOL(newdata))

  trfun <- cbind(1, GG.gaussian(newdata, gamma, th))
  local <- (cbind(1,newdata) %*% t(phi1)) * trfun;
  result <- apply(local, 1, sum)

  result
  
}


print.star <- function(x, ...) {
  NextMethod(...)
  cat("\nMultiple regime STAR model\n\n")
  x <- x$model.specific
  
  for (i in 1:x$noRegimes) {
    cat("Regime ", i, ":\n")

    cat("    Linear parameters: ")
    print(x$phi1[i,])
    cat("\n")

    if(i > 1) {
      cat("    Non-linear parameters:\n")
      print(x$phi2[i - 1,])
    }  

    cat("\n")

  }
  
  invisible(x)
}
