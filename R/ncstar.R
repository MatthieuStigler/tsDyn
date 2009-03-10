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


#Logistic transition function
# z: threshold variable(s)
# gamma: smoothing parameter(s)
# th: threshold parameter(s)
GG <- function(z, gamma, th) {
  regimes <- length(gamma)
  m <- NROW(z)
  result <- array(0, c(m, regimes))

  for(i in 1:regimes) 
#    result[,i] <- plogis(z[,i], th[i], 1/(gamma[i]))
    result[,i] <- sigmoid(gamma[i] * (z[,i] - th[i]))
    
  result
}

#Logistic transition function
# z: threshold variable(s)
# gamma: smoothing parameter(s)
# th: threshold parameter(s)
GG1 <- function(z, th) {
  regimes <- length(th)
  m <- NROW(z)
  result <- array(0, c(m, regimes))

  for(i in 1:regimes) 
#    result[,i] <- plogis(z[,i], th[i], 1/(gamma[i]))
    result[,i] <- sigmoid(z[,i] - th[i])
    
  result
}

FF <- function(object, ...)
  UseMethod("FF")

#Fitted values, given parameters
# phi1: matrix of linear parameters
# phi2omega: vector of tr. functions' parameters
#                                     (gamma_{1...p}, th{1...p}, omega_{1...p})
FF.ncstar <- function(phi1, phi2omega, xx) {
  noRegimes <- NROW(phi1)
  x_t <- cbind(1,xx)
  
  phi2 <- phi2omega[1:((noRegimes - 1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  gamma <- phi2[,1]
  th <- phi2[,2]
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  trfun <- cbind(1, GG(xx %*% omega, gamma, th))
  local <- (x_t %*% t(phi1)) * trfun;

# Variante original (ídem, en bucle)
#  local <- array(0, c(noRegimes, NROW(xx)))
#  for (i in 1:noRegimes) 
#    local[i,] <- (x_t %*% phi1[i,]) * trfun[,i]
  
# Variante 1 (Medeiros)
#  local <- t(t(phi1) %*% t(trfun)) * x_t

  result <- apply(local, 1, sum)
  result
}

# LM linearity testing against 2 regime NCSTAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   object: a star object
#   rob
#   sig
linearityTest.ncstar <- function(str, rob=FALSE, sig=0.95, trace=TRUE,...)
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
  xH1 <- cbind(xx^ 2, xx^3, xx^4)

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
testRegime.ncstar <- function(object, G, rob=FALSE, sig=0.95, trace = TRUE, ...)
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
  xH1 <- cbind(xx^2, xx^3, xx^4)
  
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

testIID.ncstar <- function(object, G, r, rob=FALSE, sig=0.95, trace = TRUE, ...)
{

  e <-  object$residuals;
  n <- NCOL(G);
  T <- length(e);
  t0 <- r+1;
  G <- G[t0:T,]
  selectr <- 1:r

  v <- array(NA, c(T - t0 + 1, r))
  for (t in t0:T) 
    v[t-t0+1,] = e[t - selectr];

  e <- e[t0:T]

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
  rG <- qr(G2, tol)$rank

  if(rG < n) {
    PCtmp <- princomp(G);
    PC <- PCtmp$loadings;
    lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
    indmin <- min(which(lambda > 0.99999));
    GPCA <- t(PC %*% t(G));
    GPCA <- GPCA[, 1:indmin];
    b <- lm.fit(x = GPCA, y = e)$coefficients
#    b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
    u <-  e - GPCA %*% b;
    xH0 <- GPCA;
  } else {
#    b <- solve(t(G) %*% G) %*% t(G) * e;
#    lm.fit(x = G, y = e)$coefficients;
    b <- lm(e ~ . -1, data=data.frame(G))$coefficients;
    u <- e - G %*% b;
  }

  SSE0 <- sum(u^2) / T;

  if (rG < n) {
    xH1 = cbind(v, GPCA);
    n = NCOL(GPCA);
  } else {
    xH1 = cbind(v, G);
  }

 # Nonlinear model (Alternative hypothesis)
  c <- lm.fit(x = xH1, y = u)$coefficients
#  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
#  dim(c) <- c(NCOL(xH1), 1);
  e <- u - xH1 %*% c #e=u-X*c;
  SSE <- sum(e^2) / T;

  # Compute the statistic
  F = ((SSE0 - SSE) / r) / (SSE / (T - n - r));

  pValue <- pf(F, r, T - n - r, lower.tail = FALSE);

  if (pValue >= sig) {
    return(list(isIID = FALSE, pValue = pValue));
  }
  else {
    return(list(isIID = TRUE, pValue = pValue));
  }
}

testConstVar.ncstar <- function(object, G, rob=FALSE, sig=0.95, trace = TRUE, ...)
{

  e <-  object$residuals;
  n <- NCOL(G);
  T <- length(e);
  q <- NCOL(object$str$xx)
  normG <- norm(t(G) %*% e);

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
  rG <- qr(G2, tol)$rank

  if(normG > 1e-6) {
    if(rG < n) {
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda > 0.99999));
      GPCA <- t(PC %*% t(G));
      GPCA <- GPCA[, 1:indmin];
      b <- lm.fit(x = GPCA, y = e)$coefficients
#    b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
#    b <- solve(t(G) %*% G) %*% t(G) * e;
#    lm.fit(x = G, y = e)$coefficients;
      b <- lm(e ~ . -1, data=data.frame(G))$coefficients;
      u <- e - G %*% b;
    }
  } else {
    u <-  e;
  }

  sigma_u <- var(u)

  u2 <- (u^2 / as.numeric(sigma_u)) - 1;

  SSE0 <- sum(u2^2) / T;

  x_t <- cbind(1, object$str$xx / sd(object$str$xx));

  nXH0 <- 1;
  nXH1 <- NCOL(object$str$xx);
           
 # Nonlinear model (Alternative hypothesis)
  c <- lm.fit(x = x_t, y = u2)$coefficients
#  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
#  dim(c) <- c(NCOL(xH1), 1);
  e <- u2 - x_t %*% c;
  SSE <- sum(e^2) / T;

  # Compute the third order statistic
  F = ((SSE0 - SSE) / nXH0) / (SSE / (T - nXH0 - nXH1));

  pValue <- pf(F, nXH0, T - nXH0 - nXH1, lower.tail = FALSE);

  if (pValue >= sig) {
    return(list(isConstVar = FALSE, pValue = pValue));
  }
  else {
    return(list(isConstVar = TRUE, pValue = pValue));
  }
}

testParConst.ncstar <- function(object, G, rob=FALSE, sig=0.95, trace = TRUE, ...)
{

  e <-  object$residuals;
  n <- NCOL(G);
  T <- length(e);
  t <- 1:T;
  q <- NCOL(object$str$xx)
  normG <- norm(t(G) %*% e);
  noRegimes <- object$noRegimes
  xx <- object$str$xx
  x_t <- cbind(1, xx)

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
  rG <- qr(G2, tol)$rank

  if(normG > 1e-6) {
    if(rG < n) {
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda > 0.99999));
      GPCA <- t(PC %*% t(G));
      GPCA <- GPCA[, 1:indmin];
      b <- lm.fit(x = GPCA, y = e)$coefficients
#    b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
#    b <- solve(t(G) %*% G) %*% t(G) * e;
#    lm.fit(x = G, y = e)$coefficients;
      b <- lm(e ~ . -1, data=data.frame(G))$coefficients;
      u <- e - G %*% b;
    }
  } else {
    u <-  e;
  }

  SSE0 <- sum(u^2) / T;

  phi2omega <- object$model.specific$phi2omega
  phi2 <- phi2omega[1:((noRegimes - 1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  gamma <- phi2[,1]
  th <- phi2[,2]
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  trfun <-  GG(xx %*% omega, gamma, th)

  XH0 <- xx;
  for (i in 2:noRegimes)
    XH0 <- cbind(XH0, repmat(trfun[,i-1], 1 , q) * xx);   

  XH1 <- repmat(t, 1 , q) * xx;
  for (i in 2:noRegimes)
    XH1 <- cbind(XH1, repmat(trfun[,i-1], 1, q) * repmat(t, 1, q)  * xx);   

  X <- cbind(XH0, XH1);
  X <- X / sd(X)

  nXH0 <- NCOL(XH0);
  nXH1 <- NCOL(XH1);
           
 # Nonlinear model (Alternative hypothesis)
  c <- lm.fit(x = X, y = u)$coefficients
#  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
#  dim(c) <- c(NCOL(xH1), 1);
  e <- u - X %*% c;
  SSE <- sum(e^2) / T;

  # Compute the third order statistic
  F = ((SSE0 - SSE) / nXH0) / (SSE / (T - nXH0 - nXH1));

  pValue <- pf(F, nXH0, T - nXH0 - nXH1, lower.tail = FALSE);

  if (pValue >= sig) {
    return(list(isConstVar = FALSE, pValue = pValue));
  }
  else {
    return(list(isConstVar = TRUE, pValue = pValue));
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
computeGradient.ncstar <- function(object, ...)
{
 
  noRegimes <- object$model.specific$noRegimes;
  n.used <- NROW(object$str$xx);
  x_t <- cbind(1, object$str$xx);
  xx <- object$str$xx;
  p <- NCOL(x_t)
  q <- NCOL(xx)

  phi2omega <- object$model.specific$phi2omega
  
  # phi2 contains parameters gamma and th for each regime
  phi2 <- phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  
  gamma <- phi2[,1]
  th <- phi2[,2]

  # omega (one vector for each regime)
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  b1 <- gamma * th
  w1 <- array(NA, c(NCOL(xx), noRegimes - 1))
  for (i in 1:NCOL(omega))
    w1[,i] <- gamma[i] * omega[,i]
  
  # linear parameters
  phi1 <- object$model.specific$phi1;

  fX <- GG1(xx %*% w1, b1);
  dfX <- dsigmoid(fX);
  
  gPhi1 <- array(x_t, c(n.used, noRegimes * p));
  for (i in 2:noRegimes) 
    gPhi1[,(i * p - p + 1):(i * p)] <- x_t * fX[,i-1] ;
  
  gb1 <- array(0, c(n.used, noRegimes-1));
  for (i in 1:(noRegimes-1)) 
    gb1[,i] <- - (x_t %*% phi1[i+1,]) * dfX[,i]

  gw1 <- array(0, c(n.used, (noRegimes-1) * q));
  for (i in 1:(noRegimes - 1)) 
    gw1[, (i*q - q+1):(i*q)] <- repmat(dfX[,i] * (x_t %*% phi1[i + 1,]), 1, q) * xx

#  gGamma <- array(0, c(n.used, noRegimes-1));
#  for (i in 1:(noRegimes - 1)) 
#    gGamma[, i] <- (x_t %*% phi1[i + 1,]) * (dfX[,i] * (xx %*% omega[,i] - th[i]));

#  gTh <- array(0, c(n.used, noRegimes-1))
#  for (i in 1:(noRegimes - 1)) 
#    gTh[,i] <- - (x_t %*% phi1[i + 1,]) * (gamma[i] * dfX[,i]);
                 
#  gOmega <- array(0, c(n.used, (noRegimes - 1) * q))
#  for (i in 1:(noRegimes - 1)) {
##    for (j in 1:NCOL(xx)) 
#      gOmega[, (i*q - q+1):(i*q)] <-
#        repmat(dfX[,i] * gamma[i] * (x_t %*% phi1[i + 1,]), 1, q) * xx;
#  }
##  dim(gOmega) <- c(n.used, (noRegimes - 1) * NCOL(xx))

#  return(cbind(gPhi1, gGamma, gTh, gOmega))
  return(cbind(gPhi1, gw1, gb1))

}

addRegime <- function(object, ...)
  UseMethod("addRegime")

addRegime.ncstar <- function(object)
{

  noRegimes <- object$model.specific$noRegimes
  xx <- object$str$xx
  phi2omega <- object$model.specific$phi2omega

  if(noRegimes == 1) {

    # omega is random but has norm equal to 1
    omega <- c(runif(1), runif(NCOL(xx) - 1, min=-1,max=1))
#    omega <- omega * (1 / norm(Matrix(omega), "f"))
    omega <- omega * (1 / norma(omega))
    dim(omega) <- c(NCOL(xx), 1)

    gamma <- 0
    th <- rnorm(1, mean(xx %*% omega), sd(xx %*% omega))
    
    phi2 <- cbind(gamma, th);

    phi1 <-   object$model.specific$phi1
    phi1 <-   rbind(object$model.specific$phi1, rnorm(NCOL(xx) + 1));
    
  }
  else {
    # phi2 contains parameters gamma and th for each regime
    phi2 <- phi2omega[1:((noRegimes-1) * 2)];
    dim(phi2) <- c(noRegimes - 1, 2)
    
    gamma <- phi2[,1];
    th <- phi2[,2];
    
    gamma[noRegimes] <- 0.0; # Add one gamma
    th[noRegimes] <- 0.0;        # Add one th

    phi2 <- cbind(gamma, th);

    # omega (one column vector for each regime)
    omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
    dim(omega) <- c(NCOL(xx), noRegimes - 1)

#    newOmega <- rep(0.0, NCOL(xx)) #c(runif(1), runif(NCOL(xx) - 1, min=-1,max=1))
#    newOmega <- newOmega * (1 / norm(Matrix(newOmega), "f"))
#    dim(newOmega) <- c(NCOL(xx), 1)
    
    omega <- cbind(omega, rnorm(NCOL(xx))); # Add one (random) vector

    # phi1 (one vector for each regime)
    phi1 <-   object$model.specific$phi1
    phi1 <-   rbind(object$model.specific$phi1, rnorm(NCOL(xx) + 1));
  }

  # Update the object
  object$model.specific$noRegimes <- noRegimes + 1;
  object$noRegimes <- noRegimes + 1;
  
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2omega <- c(phi2, omega)
  
  object$model.specific$coefficients <- c(phi1, phi2, omega);
  object$coefficients <- c(phi1, phi2, omega);
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients);
  
  return(object);
  
}

startingValues <- function(object, ...)
  UseMethod("startingValues")

startingValues.ncstar <- function(object, trace=TRUE, svIter, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  xx <- object$str$xx;
  x_t <- cbind(1, xx)
  yy <- object$str$yy;
  n.used <- NROW(object$str$xx);
  local1 <- array(x_t,  c(n.used * NCOL(x_t), noRegimes))
  local2 <- array(0, c(NROW(xx), noRegimes))
  
  phi2omega <- object$model.specific$phi2omega
  
  # phi2 contains parameters gamma and th for each regime
  phi2 <- phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)

  gamma <- phi2[,1];
  th <- phi2[,2];

  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  # Fix 'svIter' random starting values for omega and th
  newOmega <- c(runif(svIter), runif(svIter * (NCOL(xx) - 1), min=-1, max=1))
  dim(newOmega) <- c(svIter, NCOL(xx));
  newOmega <- t(newOmega) # for coherence with omega's structure
 
  newTh <- array(NA, svIter)
  for(i in 1:svIter) {
#    newOmega[,i] <- newOmega[,i] / norm(Matrix(newOmega[,i]), "f")
    newOmega[,i] <- newOmega[,i] / norma(newOmega[,i])
#    newTh[i] <- rnorm(1, mean(xx %*% newOmega[,i]), sd(xx %*% newOmega[,i]))
    newTh[i] <- median(xx %*% newOmega[,i])
  }
  
  maxGamma <- 40; # abs(8 / ((max(xx %*% newOmega) - newTh)))
  minGamma <- 1; # abs(1 / ((min(xx %*% newOmega) - newTh)));
  rateGamma <- 2; # (maxGamma - minGamma) / 20;     
  bestCost <- Inf;

  for(i in 1:svIter) { 

    if ((i %% 25 == 0) && trace) cat(".")

#    newOmega <- c(runif(1, min=0, max=1),
#                                runif(NCOL(xx) - 1, min=-1, max=1))
#    newOmega <- newOmega * (1 / norma(newOmega))
#    dim(newOmega) <- c(NCOL(xx), 1)
#    newTh <- rnorm(1, mean(xx %*% newOmega), sd(xx %*% newOmega))

    omega[, noRegimes - 1] <- newOmega[,i]
    th[noRegimes - 1] <- newTh[i]
    z <- xx %*% omega
        
    for(newGamma in seq(minGamma, maxGamma, rateGamma)) {

      gamma[noRegimes - 1] <- newGamma;

      trfun <- cbind(1, GG(z, gamma, th))
      
      # We fix the linear parameters here, before optimizing the nonlinear.
      local1 <- apply(trfun, 2, "*", x_t)
      dim(local1) <- c(n.used, noRegimes * NCOL(x_t))
      newPhi1<- lm.fit(x = local1, y = yy)$coefficients
      dim(newPhi1) <- c(NCOL(xx) + 1, noRegimes)
#      newPhi1 <- t(newPhi1)

      local2 <- (x_t %*% newPhi1) * trfun;
      
      y.hat <- apply(local2, 1, sum)
      cost <- crossprod(yy - y.hat) / sqrt(n.used);

      if(! is.na(cost)) {
        if(cost < bestCost) {
          bestCost <- cost;

          bestGamma <- gamma
          bestTh <- th
          bestOmega <- omega
          bestPhi1 <- newPhi1
        }
      }
    }
  }

  gamma <- bestGamma
  th <- bestTh
  omega <- bestOmega
  phi1 <- t(bestPhi1)
  
   if (trace) cat("\n  Starting values fixed for regime ", noRegimes,
                 ":\n\tgamma = ", gamma[noRegimes - 1],
                 "\n\tth = ", th[noRegimes - 1],
                 "\n\tomega = ", omega[, noRegimes - 1], "\n");

  object$fitted.values <- FF.ncstar(phi1, c(gamma, th, omega), xx); # y.hat
  object$residuals <- yy - object$fitted.values; # e.hat
  object$model.specific$phi2omega <- c(gamma, th, omega)
  object$model.specific$phi1 <- phi1;
  object$model.specific$coefficients <- c(phi1, gamma, th, omega);
  object$coefficients <- c(phi1, c(gamma, th, omega));
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
estimateParams.ncstar <- function(object, trace=TRUE,
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
  gradEhat <- function(par) {
    b1 <- par[1:(noRegimes-1)];                            # th * gamma
    w1 <- par[noRegimes:length(par)]  # omega * gamma
    dim(w1) <- c(q, noRegimes - 1)
    
    # We fix the linear parameters because they're not known here
    trfun <- GG1(xx %*% w1, b1)
    tmp <- apply(cbind(1, trfun), 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
    phi1 <- t(phi1)

    # Compute the gradients
    dfX <- dsigmoid(trfun);
    gb1 <- array(NA, c(n.used, noRegimes - 1)) 
    for (i in 1:(noRegimes - 1)) 
      gb1[,i] <- - (x_t %*% phi1[i+1,]) * dfX[,i];
    gw1 <- array(NA, c(n.used, (noRegimes - 1) * q)) 
    for (i in 1:(noRegimes - 1)) 
      gw1[, (i*q - q+1):(i*q)] <- repmat((x_t %*% phi1[i+1,]) * dfX[,i], 1, q) * xx;
    
    J = - cbind(gb1, gw1) / sqrt(n.used)

    local <- (x_t %*% t(phi1)) * cbind(1, trfun);
    y.hat <- apply(local, 1, sum)
    e.hat <- yy - y.hat;
    return(2 * t(e.hat) %*% J)
    
  }

  #Sum of squares function
  #p: vector of parameters
  SS <- function(par) {
    b1 <- par[1:(noRegimes-1)];
    w1 <- par[noRegimes:length(par)]
    dim(w1) <- c(q, noRegimes - 1)
    
    # We fix the linear parameters because they're not known here
    trfun <- cbind(1, GG1(xx %*% w1, b1))
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

  #Sum of squares function
  #p: vector of parameters
  SS.ga <- function(par) {
    b1 <- par[1:(noRegimes-1)];
    w1 <- par[noRegimes:length(par)]
    dim(w1) <- c(q, noRegimes - 1)

    # We check the restrictions
    for (i in 1:(noRegimes - 1)) 
      if( abs(norma(w1[,i]) - 1) < 0.05 )
#      if( abs(norm(Matrix(w1[,i]), "f") - 1) < 0.05 )
        return(Inf)
    
    if( noRegimes > 2) {
      sorting <- sort(b1, index.return = TRUE)
      if (!prod((1:(noRegimes - 1)) == sorting$ix))
        return(Inf)
    }

    # We fix the linear parameters because they're not known here
    trfun <- cbind(1, GG1(xx %*% w1, b1))
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
    b1 <- par[1:(noRegimes-1)];                            # th * gamma
    w1 <- par[noRegimes:length(par)]  # omega * gamma
    dim(w1) <- c(q, noRegimes - 1)
    
    # We fix the linear parameters because they're not known here
    trfun <- GG1(xx %*% w1, b1)
    tmp <- apply(cbind(1, trfun), 2, "*", x_t)
    dim(tmp) <- c(n.used, noRegimes * p)
    phi1 <- lm.fit(x = tmp, y = yy)$coefficients
    dim(phi1) <- c(p, noRegimes)
    phi1 <- t(phi1)

    # Compute the gradients
    dfX <- dsigmoid(trfun);
    gb1 <- array(NA, c(n.used, noRegimes - 1)) 
    for (i in 1:(noRegimes - 1)) 
      gb1[,i] <- - (x_t %*% phi1[i+1,]) * dfX[,i];
    gw1 <- array(NA, c(n.used, (noRegimes - 1) * q)) 
    for (i in 1:(noRegimes - 1)) 
      gw1[, (i*q - q+1):(i*q)] <- repmat((x_t %*% phi1[i+1,]) * dfX[,i], 1, q) * xx;
    
    J = - cbind(gb1, gw1) / sqrt(n.used)

    return(J)
    
  }

  #Sum of squares function
  #p: vector of parameters
  SS.lm <- function(par) {
    b1 <- par[1:(noRegimes-1)];
    w1 <- par[noRegimes:length(par)]
    dim(w1) <- c(q, noRegimes - 1)
    
    # We fix the linear parameters because they're not known here
    trfun <- cbind(1, GG1(xx %*% w1, b1))
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

  # Buid 'par', without gamma
  phi2 <- object$model.specific$phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)

  gamma <- phi2[,1]
  th <- phi2[,2]

    # omega (one column vector for each regime)
  omega <- object$model.specific$phi2omega[(((noRegimes - 1) * 2) + 1):
                                           length(object$model.specific$phi2omega)]
  dim(omega) <- c(q, noRegimes - 1)

  # Reorder the regimes according to the values of th
  sorting <- sort(th, index.return=TRUE)
  if ((noRegimes > 2) &&
      ! prod((1:(noRegimes - 1)) == sorting$ix)) {
    if(trace) cat("  Reordering regimes... ")
    
    ordering <-  sorting$ix
    
    th <- sorting$x
    gamma <- gamma[ordering]
    omega <- omega[, ordering]
#    phi1 <- phi1[ordering,]

    if(trace) cat("OK.\n")
    
    z <- xx %*% omega
    trfun <- cbind(1, GG(z, gamma, th))
      
    # We fix the linear parameters again.
    local1 <- apply(trfun, 2, "*", x_t)
    dim(local1) <- c(n.used, noRegimes * NCOL(x_t))
    phi1<- lm.fit(x = local1, y = yy)$coefficients
    dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
    phi1 <- t(phi1)
  }

  w1 <- array(NA, c(q, noRegimes - 1))
  for (i in 1:NCOL(omega))
    w1[,i] <- omega[,i] * gamma[i]
  
  par <- c(th * gamma, w1)

#  domainsGamma <- t(matrix(rep(c(0,100), noRegimes - 1), c(2,noRegimes - 1)))
#  domainsTh <- t(matrix(rep(c(min(xx) - sd(xx[,1]), max(xx) + sd(xx[,1])),
#                            noRegimes - 1), c(2, noRegimes - 1)))
#  domainsOmega <- t(matrix(rep(c(-1,1), noRegimes - 1), c(2, noRegimes - 1)))

  domainsb1 <- t(matrix(rep(c(100 * (min(xx) - sd(xx[,1])),
                              100 * (max(xx) + sd(xx[,1]))), noRegimes - 1),
                        c(2, noRegimes - 1)))
  domainsw1 <- t(matrix(rep(c(-100, 100), q * (noRegimes - 1)),
                        c(2,q * (noRegimes - 1))))
  domainsGA <- rbind(domainsb1, domainsw1)

  if (alg == "LM") {
    res <- try(nls.lm(par=par, fn=SS.lm, jac = gradEhat.lm,
                  control=nls.lm.control(nprint=100, maxfev=1e+9, maxiter=1024)));
  } else if (alg == "BFGS") {
    res <- optim(par, fn=SS, gr = gradEhat, method = alg,
                  control = list(trace=0, hessian = FALSE, maxit=1000));
  } else if (alg == "SANN") {
    res <- optim(par, fn=SS, method = alg,
                  control = list(trace=0, maxit=20000, temp=200,tmax=50));
  } else if (is.null(cluster)) {
    if (alg == "GAD") {
      res <- genoud(fn = SS.ga, starting.values=par, nvars=length(par),
                    gr = gradEhat, cluster=FALSE, control=list(maxit=1000),
                    print.level = 0, pop.size=1000, max.generations=200,
                    wait.generations=50, Domains = domainsGA)
    } else if (alg == "GA") {
      res <- genoud(fn = SS.ga, starting.values=par, nvars=length(par),
                    BFGS=FALSE, cluster= FALSE, control=list(maxit=1000),
                    print.level = 0, pop.size=1000, max.generations=200,
                    wait.generations=50, Domains = domainsGA)
    }
  } else {
    if (alg == "GAD") {
      res <- genoud(fn = SS.ga, starting.values=par, nvars=length(par),
                    gr = gradEhat, cluster=cluster, control=list(maxit=1000),
                    print.level = 0, pop.size=1000, max.generations=200,
                    wait.generations=50, Domains = domainsGA)
    } else if (alg == "GA") {
      res <- genoud(fn = SS.ga, starting.values=par, nvars=length(par),
                    BFGS=FALSE, cluster=cluster, control=list(maxit=1000),
                    print.level = 0, pop.size=1000, max.generations=200,
                    wait.generations=50, Domains = domainsGA)
    }
  }


  b1 <- res$par[1:(noRegimes - 1)]
  w1 <- res$par[noRegimes:length(res$par)]
  dim(w1) <- c(q, noRegimes - 1)

  for (i in 1:(noRegimes - 1)) {
    if(w1[1,i] < 0) {
      w1[, i] <- -1 * w1[,i] 
      b1[i] <- -1 * b1[i];
    }

#    gamma[i] <- norm(Matrix(w1[,i]), "f")
    gamma[i] <- norma(w1[,i])
    omega[, i] <- w1[,i] / gamma[i]; 
    th[i] <- b1[i] / gamma[i];

  }
  
  if (trace) cat("  Optimized values fixed for regime ", noRegimes,
                 ": \n    gamma = ", gamma,
                 "\n    th = ", th,
                 "\n    omega = ", omega, "\n");

#  if(trace) 
#    if(res$convergence != 0)
#      cat("  *** Convergence problem. Code: ",res$convergence,"\n")
#    else {
#      cat("  Optimization algorithm converged\n")
#    }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3.- Estimate linear parameters again
  trfun <- GG(xx %*% omega, gamma, th)
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
  
  object$fitted.values <- FF.ncstar(phi1, c(gamma, th, omega), xx); # y.hat
  object$residuals <- yy - object$fitted.values; # e.hat
  object$coefficients <- c(phi1, c(gamma, th, omega))
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2omega <- c(gamma, th, omega);

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
ncstar <- function(x, m=2, noRegimes, d = 1, steps = d, series, tests = FALSE,
                   mTh, thDelay, thVar, sig=0.95, trace=TRUE, svIter = 1000,
                   cluster= NULL, control=list(), alg="LM", ...)
{

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))   series <- deparse(substitute(x))

  # Normalize the series.
#  if ((mean(x) >= 0.05) && (sd(x) >= 0.05) {
#    if (trace) cat("Normalizing the series.\n")
#    x <- (x - mean(x)) / sd(x);
#  }
  
  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  if(!is.null(cluster)) {
    setDefaultClusterOptions(master="localhost", port=10187)
    cl <- makeCluster(cluster, type="SOCK")
    clusterSetupRNG(cl)
#    clusterEvalQ(cl, library(Matrix))
  }

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
  testResults <- linearityTest.ncstar(str, rob=rob, sig=sig, trace = trace)
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
                     phi2omega = NA, phi1= phi1) 
    object <- extend(nlar(str, coefficients = linearModel$coefficients,
                          fitted = linearModel$fitted.values,
                          residuals = linearModel$residuals,
                          k=NCOL(x_t), noRegimes=1,
                          model=data.frame(yy,xx),
                          model.specific=list),
                     "ncstar")

    if(!is.null(cluster)) stopCluster(cl)
    return(object);
  }
  else {
    if(trace) cat("The series is nonlinear. Incremental building procedure:\n")
    
    # Build the linear model
    linearModel <- lm.fit(x = x_t, y = yy);

    phi1= linearModel$coefficients
    dim(phi1) = c(1, NCOL(x_t)) # one row, one regime
    
    # Create the ncstar object
    list = list(noRegimes=1, m = m, fitted = linearModel$fitted.values,  
                     residuals = linearModel$residuals,
                     coefficients = linearModel$coefficients,
                     k = length(linearModel$coefficients),
                     convergence=NA, counts=NA, hessian=NA, message=NA,
                     value=NA, par=NA,
                     phi2omega = NA, phi1= phi1) 
    object <- extend(nlar(str, coefficients = linearModel$coefficients,
                          fitted = linearModel$fitted.values,
                          residuals = linearModel$residuals,
                          k=NCOL(x_t), noRegimes=1,
                          model=data.frame(yy,xx),
                          model.specific=list),
                     "ncstar")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Add-regime loop
    while (increase && (object$model.specific$noRegimes < noRegimes)) {
      
      nR <- object$model.specific$noRegimes + 1
      
      if(trace) cat("- Adding regime ", nR, ".\n");
      object <- addRegime(object);

#      if(nR==4) debug(startingValues.ncstar)
      
      if(trace) cat("- Fixing good starting values for regime ", nR);
      if(is.null(cluster)) {
        object <- startingValues.ncstar(object, trace=trace, svIter = svIter);
      } else {
        if(trace) cat("\n   + Doing distributed computations... ")
        solutions <- clusterCall(cl, startingValues.ncstar, object, trace=trace,
                                 svIter = svIter %/% length(cluster))
        cost <- rep(Inf, length(solutions))
        if(trace) cat("\n   + Gathering results...\n")
        for (i in 1:length(solutions)) {
          cost[i] <- sum(solutions[[i]]$residuals^2) / sqrt(n.used)
        }
        
        object <- solutions[[which.min(cost)]]

        if (trace) cat("\n  Starting values fixed for regime ", nR,
#                       ":\n", object$model.specific$phi2omega, "\n");
             ":\n\tgamma = ", object$model.specific$phi2omega[1:(nR - 1)],
             "\n\tth = ", object$model.specific$phi2omega[nR:(2*(nR - 1))],
             "\n\tomega = ",
                object$model.specific$phi2omega
                       [(2 * (nR - 1) + 1):length(object$model.specific$phi2omega)], "\n");
      }

      if(trace) cat('- Estimating parameters of regime', nR, ', (alg: ', alg, ')...\n')
      object <- estimateParams.ncstar(object, control=control, trace=trace,
                               alg=alg, cluster=cluster);
      
      if(trace) cat("\n- Testing for addition of regime ", nR + 1, ".\n");
      if(trace) cat("  Estimating gradient matrix...\n");
      G <- computeGradient.ncstar(object);

      if(trace) cat("  Computing the test statistic (sig = ", 1 - 0.05 / nR,
                    ")...\n");
      testResults <- testRegime.ncstar(object, G = G,
                                       rob = rob, sig = 1 - 0.05 / nR,
                                       trace=trace);

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
    
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 4. Diagnostic checking
    
    if(tests) {
      if(trace) cat("== Testing for linear independence of the residuals:\n")

      isIID <- array(NA, 12)
      pValue <- array(NA, 12)
      for(r in 1:12) {
        test1 <- testIID.ncstar(object, G, r, rob=rob, sig=sig, trace = trace)
        isIID[r] <- test1$isIID
        pValue[r] <- test1$pValue
      }
      object$model.specific$testIID <- data.frame(isIID, pValue)
      if(trace) print(data.frame(isIID, pValue))
      
      if(trace) cat("== Testing for constant variance of the residuals:\n")
      test2 <- testConstVar.ncstar(object, G, rob=rob, sig=sig, trace = trace)
      object$model.specific$testConstVar <- test2; 
      
      if(test2$isConstVar) 
        cat("      Constant variance detected, pValue = ", test2$pValue, "\n")
      else cat("     Smoothly changing variance detected, pvalue = ", test2$pValue, "\n")
      
      if(trace) cat("== Testing for parameter constancy:\n")
      test3 <- testParConst.ncstar(object, G, rob=rob, sig=sig, trace = trace)
      object$model.specific$testParConst <- test3;
      
      if(test3$isConstVar) 
        cat("      Constant parameters detected, pValue = ", test3$pValue, "\n")
      else cat("     Smoothly changing parameters detected, pvalue = ", test3$pValue, "\n")
    }
    
    if(trace) cat("\n- Finished building an NCSTAR with ",
                  object$model.specific$noRegimes, " regimes\n");

    if(!is.null(cluster)) stopCluster(cl)
    return(object);
  }
  
}

ncstar.predefined <- function(x, m=2, d = 1, steps = d,
                              noRegimes, phi1, phi2omega, series, tests = FALSE,
                   mTh, thDelay, thVar, sig=0.95, trace=TRUE, svIter = 1000,
                   cluster= NULL, control=list(), alg="LM", ...)
{

  if(!is.null(cluster)) {
    setDefaultClusterOptions(master="localhost", port=10187)
    cl <- makeCluster(cluster, type="SOCK")
    clusterSetupRNG(cl)
#    clusterEvalQ(cl, library(Matrix))
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
    
  # Create the ncstar object
  list = list(noRegimes=noRegimes, m = m, fitted = NA,  
    residuals = NA,
    coefficients = c(phi1, phi2omega),
    k = length(c(phi1, phi2omega)),
    convergence=NA, counts=NA, hessian=NA, message=NA,
    value=NA, par=NA,
    phi2omega = phi2omega, phi1= phi1) 
  object <- extend(nlar(str, coefficients = c(phi1, phi2omega),
                        fitted = NA,
                        residuals = NA,
                        k =length(c(phi1, phi2omega)),  noRegimes=noRegimes,
                        model=NULL,
                        model.specific=list),
                   "ncstar")

#      cat("\nStarting values for linear parameters:\n    ")
#      print(phi1_median)
      
#      cat("\nStarting values for nonlinear parameters:\n")
#      cat("gamma = ", phi2_median[1:(noRegimes - 1)])
#      cat("; th = ", phi2_median[noRegimes:(2*(noRegimes - 1))])
#      cat("\nomega=", phi2_median[(2 * (noRegimes - 1) + 1):length(phi2_median)])

  if(trace) cat('- Estimating parameters for all ', noRegimes, ' regimes, (alg: ', alg, ')...\n')
  object <- estimateParams.ncstar(object, control=control, trace=trace,
                                  alg=alg, cluster=cluster);    
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 4. Diagnostic checking
    
  if(tests) {
    if(trace) cat("== Testing for linear independence of the residuals:\n")

    isIID <- array(NA, 12)
    pValue <- array(NA, 12)
    for(r in 1:12) {
      test1 <- testIID.ncstar(object, G, r, rob=rob, sig=sig, trace = trace)
      isIID[r] <- test1$isIID
      pValue[r] <- test1$pValue
    }
    object$model.specific$testIID <- data.frame(isIID, pValue)
    if(trace) print(data.frame(isIID, pValue))
    
    if(trace) cat("== Testing for constant variance of the residuals:\n")
    test2 <- testConstVar.ncstar(object, G, rob=rob, sig=sig, trace = trace)
    object$model.specific$testConstVar <- test2; 
    
    if(test2$isConstVar) 
      cat("      Constant variance detected, pValue = ", test2$pValue, "\n")
    else cat("     Smoothly changing variance detected, pvalue = ", test2$pValue, "\n")
    
    if(trace) cat("== Testing for parameter constancy:\n")
    test3 <- testParConst.ncstar(object, G, rob=rob, sig=sig, trace = trace)
    object$model.specific$testParConst <- test3;
    
    if(test3$isConstVar) 
      cat("      Constant parameters detected, pValue = ", test3$pValue, "\n")
    else cat("     Smoothly changing parameters detected, pvalue = ", test3$pValue, "\n")
  }
  
  if(trace) cat("\n- Finished building an NCSTAR with ",
                object$model.specific$noRegimes, " regimes\n");
  
  if(!is.null(cluster)) stopCluster(cl)
  return(object);
}


oneStep.ncstar <- function(object, newdata, itime, thVar, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  xx <- object$str$xx

  phi1 <- object$model.specific$phi1;
  phi2omega <- object$model.specific$phi2;

  phi2 <- phi2omega[1:((noRegimes - 1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  gamma <- phi2[,1]
  th <- phi2[,2]
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

#  if(nrow(newdata) > 1) {
    trfun <- cbind(1, GG(newdata %*% omega, gamma, th))
    local <- (cbind(1, newdata) %*% t(phi1)) * trfun;
    result <- apply(local, 1, sum)
#  }
#  else {
#    accum <- c(1, newdata) %*% phi1[1,];
#    for (i in 2:noRegimes) 
#      accum <- accum + 
#        (c(1, newdata) %*% phi1[i,]) * G(z, phi2[i - 1,1], phi2[i - 1,2])
    
#    result <- accum
#  }

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
