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
# z: variable
# gamma: smoothing parameter
# th: threshold value
GG <- function(z, gamma, th) {
  p <- length(gamma)
  m <- NROW(z)
  result <- array(0, c(m, p))

  for(i in 1:p) 
#    result[,i] <- plogis(z[,i], th[i], 1/(gamma[i]))
    result[,i] <- sigmoid(gamma[i] * (z[,i] - th[i]))
    
  result
}

#Fitted values, given parameters
# phi1: matrix of linear parameters
# phi2omega: vector of tr. functions' parameters
#                                     (gamma_{1...p}, th{1...p}, omega_{1...p})
FF <- function(phi1, phi2omega, xx) {
  noRegimes <- NROW(phi1)
  
  phi2 <- phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  gamma <- phi2[,1]
  th <- phi2[,2]
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  trfun <- GG(xx %*% omega, gamma, th)

  local <- array(0, c(noRegimes, NROW(xx)))
  local[1,] <- cbind(1,xx) %*% phi1[1,];
  for (i in 2:noRegimes) 
    local[i,] <- (cbind(1,xx) %*% phi1[i,]) * trfun[,i - 1]
#        sigmoid(gamma[i-1] * (xx %*% omega[,i-1] - th[i-1]));
  
  result <- apply(local, 2, sum)
  result
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

  phi2omega <- object$model.specific$phi2omega
  
  # phi2 contains parameters gamma and th for each regime
  phi2 <- phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)
  
  gamma <- phi2[,1]
  th <- phi2[,2]
  
  # omega (one vector for each regime)
  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  # linear parameters
  phi1 <- object$model.specific$phi1;

  fX <- array(0, c(noRegimes - 1, n.used));
  dfX <- array(0, c(noRegimes - 1, n.used));
  gPhi1 <- x_t;
  for (i in 1:(noRegimes - 1)) {
    fX[i,] <- plogis(xx %*% omega[,i], th[i], 1/gamma[i])
    dfX[i,] <- plogis(xx %*% omega[,i], th[i], 1/gamma[i])
                                        (1 - plogis(xx %*% omega[,i], th[i], 1/gamma[i]));
    gPhi1 <- cbind(gPhi1, kronecker(matrix(1, 1, NCOL(x_t)), fX[i,]) * x_t)
  }
  
  gTh <- array(0, c(n.used, noRegimes-1))
  gGamma <- array(0, c(n.used, noRegimes-1));
  gOmega <- array(0, c(n.used, (noRegimes - 1), NCOL(xx)))
  for (i in 1:(noRegimes - 1)) {
    gTh[,i] <-           - (x_t %*% phi1[i + 1,]) * (gamma[i] * dfX[i,]);
    gGamma[, i] <- (x_t %*% phi1[i + 1,]) * (dfX[i,] * (xx %*% omega[,i] - th[i]));
    for (j in 1:NCOL(xx)) 
      gOmega[, i, j] <-  (x_t %*% phi1[i + 1,]) * (dfX[i,] * gamma[i] * xx[i,j])
  }
  dim(gOmega) <- c(n.used, (noRegimes - 1) * NCOL(xx))

  return(cbind(gPhi1, gGamma, gTh, gOmega))

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
testRegime.ncstar <- function(object, G, rob=FALSE, sig=0.05, trace = TRUE, ...)
{

  e <-  object$residuals;
 # x_t <- object$model.specific$thVar;
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
      GPCA <- t(PC%*%t(G));
      GPCA <- GPCA[, 1:indmin];
      b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
#      cat("A2\n")
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
  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
  dim(c) <- c(NCOL(Z), 1);
  v <- u - Z %*% c;
  SSE <- sum(v^2);

  # Compute the third order statistic
  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  F = ((SSE0 - SSE) / nxH1) / (SSE / (T - nxH0 - nxH1));

  pValue <- pf(F, nxH1, T - nxH0 - nxH1, lower.tail = FALSE);

  if (pValue >= sig) {
    return(list(remainingNonLinearity = FALSE, pValue = pValue));
  }
  else {
    return(list(remainingNonLinearity = TRUE, pValue = pValue));
  }
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
    omega <- c(runif(1, 0.000000001, 1), runif(NCOL(xx) - 1, -1,1))
    omega <- omega * (1 / norm(Matrix(omega), "f"))
    dim(omega) <- c(NCOL(xx), 1)

    gamma <- 0
    th <- median(xx %*% omega)

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

    newOmega <- c(runif(1, 0.000000001, 1), runif(NCOL(xx) - 1, -1,1))
    newOmega <- newOmega * (1 / norm(Matrix(newOmega), "f"))
    dim(newOmega) <- c(NCOL(xx), 1)
    
    omega <- cbind(omega, newOmega); # Add one (random) vector

    # phi1 (one vector for each regime)
    phi1 <-   object$model.specific$phi1
    phi1 <-   rbind(object$model.specific$phi1, rnorm(NCOL(xx) + 1));
  }

  # Update the object
  object$model.specific$noRegimes <- object$model.specific$noRegimes + 1;
  object$noRegimes <- object$noRegimes + 1;
  
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2omega <- c(phi2, omega)
  
  object$model.specific$coefficients <- c(phi1, c(phi2, omega));
  object$coefficients <- c(phi1, c(phi2, omega));
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients);
  
  return(object);
  
}

startingValues <- function(object, ...)
  UseMethod("startingValues")

startingValues.ncstar <- function(object, trace=TRUE, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  xx <- object$str$xx;
  x_t <- cbind(1, xx)
  yy <- object$str$yy;
  n.used <- NROW(object$str$xx);
  
  phi2omega <- object$model.specific$phi2omega
  
  # phi2 contains parameters gamma and th for each regime
  phi2 <- phi2omega[1:((noRegimes-1) * 2)];
  dim(phi2) <- c(noRegimes - 1, 2)

  gamma <- phi2[,1];
  th <- phi2[,2];

  omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)
         
  bestCost <- Inf;

<<<<<<< .mine
  for(i in 1:1000) {  ######################    1 0 0 0
=======
  # Maximum and minimum values for gamma
  maxGamma <- 40;
  minGamma <- 1;
  rateGamma <- 2; ############################################!!!
  
  for(i in 1:1000) {
>>>>>>> .r139

<<<<<<< .mine
    if ((i %% 25 == 0) && trace) cat(".")
=======
    if ((i %% 20 == 0) && trace) cat(".")
>>>>>>> .r139

    newOmega <- c(runif(1, min=0, max=1),
                                runif(NCOL(xx) - 1, min=-1, max=1))
    newOmega <- newOmega * (1 / norm(Matrix(newOmega), "f"))
    dim(newOmega) <- c(NCOL(xx), 1)
<<<<<<< .mine
    omega[, noRegimes - 1] <- newOmega
=======
    omega[, noRegimes - 1] <- newOmega;
>>>>>>> .r139

#    if(i %% 2 == 0)  newTh <- -1
#    else newTh <- 1
    
#    newTh <- median(xx %*% newOmega)

    newTh <- rnorm(1, mean(xx %*% newOmega), sd(xx %*% newOmega))
    th[noRegimes - 1] <- newTh
    
    # Maximum and minimum values for gamma
    maxGamma <- 40; #abs(8 / ((max(xx %*% newOmega) - newTh)))
    minGamma <- 1; #abs(1 / ((min(xx %*% newOmega) - newTh)));
    rateGamma <- 2#(maxGamma-minGamma) / 20;     
    
    for(newGamma in seq(minGamma, maxGamma, rateGamma)) {

      gamma[noRegimes - 1] <- newGamma;

      # We fix the linear parameters here, before optimizing the nonlinear.
      tmp <- cbind(x_t, array(apply(GG(xx %*% omega, gamma, th),
                                    2, "*", x_t), c(n.used, (noRegimes - 1) * NCOL(x_t))))
      newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
      dim(newPhi1) <- c(NCOL(xx) + 1, noRegimes)
      newPhi1 <- t(newPhi1)

      y.hat <- FF(newPhi1, c(gamma, th, omega), xx);
      cost <- crossprod(yy - y.hat)

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
  phi1 <- bestPhi1
  
  if (trace) cat("\n  Starting values fixed for regime ", noRegimes,
                 ": gamma = ", gamma[noRegimes - 1],
                 ", th = ", th[noRegimes - 1],"\n");
  
  # Reorder the regimes according to the values of th
<<<<<<< .mine
  if ((noRegimes > 2) &&
                   prod((1:(noRegimes - 1)) != sort(th, index.return=TRUE)$ix)) {
    if(trace) cat("  Reordering regimes...\n")

=======
  if ((noRegimes > 2)  &&
                    prod((1:(noRegimes - 1)) != sort(th, index.return=TRUE)$ix))  {

    if(trace) cat("  Reordering regimes...\n")
>>>>>>> .r139
    ordering <-  sort(th, index.return=TRUE)$ix

    th <- sort(th, index.return=TRUE)$x
    gamma <- gamma[ordering]
    omega <- omega[, ordering]
    phi1 <- phi1[ordering,]

    # reestimate phi's
#    tmp <- cbind(x_t, array(apply(GG(xx %*% omega, gamma, th),
#                                  2, "*", x_t), c(n.used, (noRegimes - 1) * NCOL(x_t))))
#    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
#    dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
#    phi1 <- t(phi1)
  }
  
  object$model.specific$phi2omega <- c(cbind(gamma, th), omega)
  object$model.specific$phi1 <- phi1;
  object$model.specific$coefficients <- c(phi1, c(cbind(gamma, th), omega));
  object$coefficients <- c(phi1, c(cbind(gamma, th), omega));
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
estimateParams.ncstar <- function(object, trace=TRUE, control=list(), ...)
{

  xx <- object$str$xx
  x_t <- cbind(1, xx)
  yy <- object$str$yy
  n.used <- NROW(object$str$xx);
  noRegimes <- object$model.specific$noRegimes;

  # Function to compute the gradient 
  #
  # Returns the gradient with respect to the error
  gradEhat <- function(phi2omega) {
    # phi2 contains parameters gamma and th for each regime
    phi2 <- phi2omega[1:((noRegimes-1) * 2)];
    dim(phi2) <- c(noRegimes - 1, 2)

    gamma <- phi2[,1]
    th <- phi2[,2]

    # omega (one vector for each regime)
    omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
    dim(omega) <- c(NCOL(xx), noRegimes - 1)

    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- cbind(x_t, array(apply(GG(xx %*% omega, gamma, th),
                                  2, "*", x_t), c(n.used, (noRegimes - 1) * NCOL(x_t))))
    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
    phi1 <- t(phi1)

    y.hat <- FF(phi1, phi2omega, xx)
    e.hat <- yy - y.hat;

    fX <- array(0, c(noRegimes - 1, n.used));
    dfX <- array(0, c(noRegimes - 1, n.used));
#    gPhi1 <- x_t;
    gOmega <- array(0, c(n.used, (noRegimes - 1), NCOL(xx)))
    for (i in 1:(noRegimes - 1)) {
      fX[i,] <- sigmoid(gamma[i] * (xx %*% omega[,i] - th[i]));
      dfX[i,] <- sigmoid(gamma[i] * (xx %*% omega[,i] - th[i])) *
                                              (1 - sigmoid(gamma[i] * (xx %*% omega[,i] - th[i])));
#      gPhi1 <- cbind(gPhi1, kronecker(matrix(1, 1, NCOL(x_t)), fX[i,]) * x_t)
    }
    
    gGamma <- array(0, c(n.used, noRegimes-1));
    gTh <- array(0, c(n.used, noRegimes-1))
    gOmega <- array(0, c(n.used, (noRegimes - 1), NCOL(xx)))
    for (i in 1:(noRegimes - 1)) {
      gTh[,i] <-           - (x_t %*% phi1[i + 1,]) * (gamma[i] * dfX[i,]);
      gGamma[, i] <- (x_t %*% phi1[i + 1,]) * (dfX[i,] * (xx %*% omega[,i] - th[i]));
      for (j in 1:NCOL(xx)) 
        gOmega[, i, j] <-  (x_t %*% phi1[i + 1,]) * (dfX[i,] * gamma[i] * xx[i,j])
    }
    dim(gOmega) <- c(n.used, (noRegimes - 1) * NCOL(xx))
    
    J = - cbind(gGamma, gTh, gOmega) / sqrt(n.used)
      
    return(2 * t(e.hat) %*% J)
    
  }

  #Sum of squares function
  #p: vector of parameters
  SS <- function(phi2omega) {
    # phi2 contains parameters gamma and th for each regime
    phi2 <- phi2omega[1:((noRegimes-1) * 2)];
    dim(phi2) <- c(noRegimes - 1, 2)

    gamma <- phi2[,1]
    th <- phi2[,2]

    # omega (one column vector for each regime)
    omega <- phi2omega[(((noRegimes - 1) * 2) + 1):length(phi2omega)]
    dim(omega) <- c(NCOL(xx), noRegimes - 1)

    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- cbind(x_t, array(apply(GG(xx %*% omega, gamma, th),
                                  2, "*", x_t), c(n.used, (noRegimes - 1) * NCOL(x_t))))
    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
    phi1 <- t(phi1)

    # Return the sum of squares
    y.hat <- FF(phi1, phi2omega, xx)
    crossprod(yy - y.hat)
  }

  res <- optim(object$model.specific$phi2omega, SS, gr = gradEhat,
<<<<<<< .mine
               method="CG", control = control)
=======
               method="BFGS", control = control)
>>>>>>> .r139

  newPhi2 <- res$par[1:((noRegimes-1) * 2)];
  dim(newPhi2) <- c(noRegimes - 1, 2)
  
  gamma <- newPhi2[,1]
  th <- newPhi2[,2]
  
  omega <- res$par[(((noRegimes - 1) * 2) + 1):length(res$par)]
  dim(omega) <- c(NCOL(xx), noRegimes - 1)

  for (i in 1:(noRegimes - 1))
    if(gamma[i] <0)
      gamma[i] <- - gamma[i];
  
<<<<<<< .mine
  if (trace) cat("  Optimized values fixed for regime ", noRegimes,
                 ": \n    gamma = ", gamma,
                 "\n    th = ", th,
                 "\n    omega = ", omega, "\n");
=======
  if (trace) cat("  Optimized values fixed for regime ", noRegimes,
                 ": \n    gamma = ", gamma,
                 "\n    th = ", th,"\n");
>>>>>>> .r139

  if(trace) 
    if(res$convergence != 0)
      cat("  *** Convergence problem. Code: ",res$convergence,"\n")
    else {
      cat("  Optimization algorithm converged\n")
    }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3.- Estimate linear parameters again
  for (i in 1:(noRegimes - 1))
    if(gamma[i] <0)
      gamma[i] <- - gamma[i];
  
  tmp <- cbind(x_t, array(apply(GG(xx %*% omega, gamma, th),
                                2, "*", x_t), c(n.used, (noRegimes - 1) * NCOL(x_t))))
  phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
  dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
  phi1 <- t(phi1)

  if(trace) {
    cat("  Optimized linear values: \n");
    for(i in 1:noRegimes)
      cat("    ", phi1[i,], "\n")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 4.- Compute yhat and ehat
  
  object$fitted.values <- FF(phi1, c(cbind(gamma, th), omega), xx); # y.hat
  object$residuals <- yy - object$fitted.values; # e.hat
  object$coefficients <- c(phi1, c(cbind(gamma, th), omega))
  object$model.specific$phi1 <- phi1
  object$model.specific$phi2omega <- c(cbind(gamma, th), omega);
  
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
ncstar <- function(x, m=2, noRegimes, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay, thVar, sig=0.05, trace=TRUE, control=list(), ...)
{

#  debug(computeGradient.ncstar)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))   series <- deparse(substitute(x))

  # Normalize the series.
<<<<<<< .mine
#  if ((mean(x) >= 0.05) && (sd(x) >= 0.05) {
#    if (trace) cat("Normalizing the series.\n")
#    x <- (x - mean(x)) / sd(x);
#  }
=======
#  x <- (x - mean(x)) / sd(x);
>>>>>>> .r139
  
  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

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
    if(trace) cat("The series is linear. Use linear model instead.\n")
    return(str);
  }
  else {
    if(trace) cat("The series is nonlinear. Incremental building procedure:\n")
    
    # Build the linear model
    linearModel <- lm(yy ~ .  - 1, data = data.frame(x_t));

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
                          model.specific=list),
                     "ncstar")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Add-regime loop
    while (increase && (object$model.specific$noRegimes < noRegimes)) {
      
      nR <- object$model.specific$noRegimes + 1
      
      if(trace) cat("- Adding regime ", nR, ".\n");
      object <- addRegime(object);

      if(trace) cat("- Fixing good starting values for regime ", nR);
      object <- startingValues.ncstar(object, control=control, trace=trace);

      if(trace) cat('- Estimating parameters of regime', nR, '...\n')
      object <- estimateParams(object, control=control, trace=trace);
      
      if(trace) cat("\n- Testing for addition of regime ", nR + 1, ".\n");
      if(trace) cat("  Estimating gradient matrix...\n");
      G <- computeGradient(object);

      sig <- sig / 2 # We halve the significance level
      if(trace) cat("  Computing the test statistic (sig = ", sig, ")...\n");
      testResults <- testRegime(object, G = G,
                                       rob = rob, sig = sig, trace=trace);

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
    
<<<<<<< .mine
    if(trace) cat("\n- Finished building an NCSTAR with ",
=======
    if(trace) cat("\n- Finished building a MRSTAR with ",
>>>>>>> .r139
                  object$model.specific$noRegimes, " regimes\n");
    return(object);
  }
  
}

# Predefined NCSTAR model
#
#   Builds a NCSTAR with a fixed number of regimes ('noRegimes') and
#     fixed parameters phi1 (linear) and phi2 (nonlinear). If no
#     parameters are given they are set to random and evenly
#     distributed, respectively.
#
#   mTh, thDelay, thVar: as in setar
#   maxRegressors[i]: maximum number of autoregressors in regime i
#   noRegimes: number of regimes of the model
#   phi1[i]: vector with the maxRegressors[i]+1 parameters of regime i
#   phi2omega[i]: vector with the parameters of tr. function i.
#   trace: should infos be printed?
#   control: 'control' options to be passed to optim
#
ncstar.predefined <- function(x, m, noRegimes, d=1, steps=d, series,
                  maxRegressors, phi1, phi2omega,
                  mTh, thDelay, thVar, trace=TRUE, control=list())
{

  if(noRegimes == 1) 
   stop("An NCSTAR with 1 regime is an AR model: use the linear model instead.")
  
  if(missing(m))
    m <- max(maxRegressors, thDelay+1)

  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  xx <- str$xx 
  x_t <- cbind(1, xx)
  yy <- str$yy
  externThVar <- FALSE
  n.used <- NROW(xx)

  if (missing(maxRegressors))
  {
    maxRegressors <- rep(m, times = noRegimes)
    if (trace) 
      cat("Using maximum autoregressive order for all regimes: ", m,"\n")
  }

#  if(!missing(thDelay))
#  {
#    if(thDelay>=m)  stop(paste("thDelay too high: should be < m (=",m,")"))
#    z <- xx[,thDelay+1]
#  }
#  else
#  {
#    if(!missing(mTh))
#    {
#      if(length(mTh) != m) stop("length of 'mTh' should be equal to 'm'")
#      z <- xx %*% mTh #threshold variable
#      dim(x) <- NULL
# }
#    else
#    {
#      if(!missing(thVar))
#     {
#        if(length(thVar)>nrow(xx))
#        {
#          z <- thVar[1:nrow(xx)]
#          if(trace) cat("Using only first", nrow(xx), "elements of thVar\n")
#        }
#        else
#        {
#          z <- thVar
#        }
#        externThVar <- TRUE
#      }
#      else
#      {
#        if(trace) cat("Using default threshold variable: thDelay=0\n")
#        z <- xx[,1]
#      }
#    }
#  }
  
#Automatic starting values####################
# TO DO: grid search over phi2.

  if(missing(phi1)) {
    phi1 <- array(rnorm(noRegimes * m+1), c(noRegimes, m+1))
    if(trace) {
      cat('Missing starting linear values. Using random values.\n')
    }
  }
  
  if(missing(phi2omega)) {
    omega <- phi1[2:noRegimes,1:m] # for example
                                                       # (is random now, TO DO grid search)
    dim(omega) <- c(m, noRegimes - 1)
    
    phi2 <- array(0, c(noRegimes - 1, 2))
    range <- range(x_t[,2])

    phi2[,1] <- 5                   # phi2[,1] = gamma
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

#############################################

  #Sum of squares function
  #p: vector of parameters
  SS <- function(phi2omega) {
    # phi2 contains parameters gamma and th for each regime
    phi2 <- phi2omega[1:((noRegimes-1) * 2)];
    dim(phi2) <- c(noRegimes - 1, 2)

    gamma <- phi2[,1]
    th <- phi2[,2]

    # omega (one column vector for each regime)
    omega <- phi2omega[(((noRegimes-1) * 2) + 1):length(phi2omega)]
    dim(omega) <- c(NCOL(xx), noRegimes - 1)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- cbind(x_t, matrix(apply(G(xx %*% omega, gamma, th),
                                   2, "*", x_t), nrow = n.used,
                             ncol = (noRegimes - 1) * NCOL(x_t)))
    newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(newPhi1) <- c(m + 1, noRegimes)
    newPhi1 <- t(newPhi1)
    
    # Return the sum of squares
    y.hat <- FF(newPhi1, phi2omega, xx)
    crossprod(yy - y.hat)
  }
  
#Numerical optimization##########
  res <- list()
  res$convergence <- NA
  res$hessian <- NA
  res$message <- NA
  res$value <- NA

  if(optimize) {
    if(trace) cat('Optimizing...')
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- cbind(x_t, matrix(apply(G(xx %*% omega, gamma, th),
                                   2, "*", x_t), nrow = n.used,
                             ncol = (noRegimes - 1) * NCOL(x_t)))
    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(phi1) <- c(m + 1, noRegimes)
    phi1 <- t(phi1)

    # Optimize
    p <- c(phi2, omega)
    res <- optim(p, SS, hessian = TRUE, control = control)
    
    if(trace) cat(' Nonlinear optimisation done.')

    phi2omega <- res$par

    # phi2 contains parameters gamma and th for each regime
    phi2 <- phi2omega[1:((noRegimes-1) * 2)];
    dim(phi2) <- c(noRegimes - 1, 2)

    gamma <- phi2[,1]
    th <- phi2[,2]
    
    # omega (one vector for each regime)
    omega <- phi2omega[(((noRegimes-1) * 2) + 1):length(phi2omega)]
    dim(omega) <- c(NCOL(xx), noRegimes - 1)
    
    for (i in 1:(noRegimes - 1))
      if(gamma[i] <0)
        gamma[i] <- - gamma[i];
    
    # We fix the linear parameters again.
    tmp <- cbind(x_t, matrix(apply(G(xx %*% omega, gamma, th),
                                   2, "*", x_t), nrow = n.used,
                             ncol = (noRegimes - 1) * NCOL(x_t)))
    phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    dim(phi1) <- c(m+1, noRegimes)
    phi1 <- t(phi1)

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
  res$phi2omega <- c(phi2, omega)
  
  res$thVar <- x_t
  
  fitted <- FF(phi1, c(phi2, omega), xx)
  
  residuals <- yy - fitted
  dim(residuals) <- NULL	#this should be a vector, not a matrix

  res$noRegimes <- noRegimes
  res$m <- m;
  
  if(!optimize) {
    res$convergence=NA
    res$hessian=NA
    res$message=NA
    res$value=NA
    res$fitted <- fitted
    res$residuals <- residuals
    dim(res$residuals) <- NULL #this should be a vector, not a matrix
    res$k <- length(as.vector(phi1)) + length(as.vector(phi2omega))
  }
  
  return(extend(nlar(str, 
                     coef= c(phi1, phi2omega),
                     fit = fitted,
                     res = residuals,
                     k   = length(as.vector(phi1)) +
                                 length(as.vector(phi2omega)),
                     model.specific=res), "ncstar"))
}

oneStep.star <- function(object, newdata, itime, thVar, ...)
{

  noRegimes <- object$model.specific$noRegimes;

  phi1 <- object$model.specific$phi1;
  phi2 <- object$model.specific$phi2;

  if(object$model.specific$externThVar) {
    z <- thVar[itime]
  } else {
    z <- newdata %*% object$model$mTh;
    dim(z) <- NULL;
  }

  if(nrow(newdata) > 1) {
    accum <- array(0, c(noRegimes, nrow(newdata)))
    accum[1,] <- (cbind(1,newdata) %*% phi1[1,]);
    for (i in 2:noRegimes) 
      accum[i,] <-
        (cbind(1,newdata) %*% phi1[i,]) * G(z, phi2[i - 1,1], phi2[i - 1,2])
    
    result <- apply(accum, 2, sum)
  }
  else {
    accum <- c(1, newdata) %*% phi1[1,];
    for (i in 2:noRegimes) 
      accum <- accum + 
        (c(1, newdata) %*% phi1[i,]) * G(z, phi2[i - 1,1], phi2[i - 1,2])
    
    result <- accum
  }

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
