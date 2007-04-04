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

nlar.struct <- function(x, m, d=1, steps=d, series) {
  
  if(missing(series))
    series <- deparse(substitute(x))
  
  x <- as.ts(x)
  n <- length(x)

  if(NCOL(x)>1)
    stop("only univariate time series are allowed")

  if(any(is.na(x)))
    stop("missing values not allowed")

  if( missing(m) )
    stop("missing argument: 'm' is necessary")

  if((m*d + steps) > n)
    stop("time series too small to handle these embedding parameters")

  xxyy <- embedd(x, lags=c((0:(m-1))*(-d), steps) )

  extend(list(), "nlar.struct",
         x=x, m=m, d=d, steps=steps, series=series,
         xx=xxyy[,1:m,drop=FALSE], yy=xxyy[,m+1], n.used=length(x))
}

getXXYY <- function(obj, ...) UseMethod("getXXYY")

getXXYY.nlar.struct <- function(obj, ...) {
	x <- obj$x
	m <- obj$m
	d <- obj$d
	steps <- obj$steps
	embedd(x, lags=c((0:(m-1))*(-d), steps) )
}

getXX <- function(obj, ...)
	getXXYY(obj,...)[ , 1:obj$m , drop=FALSE]

getYY <- function(obj, ...)
	getXXYY(obj, ...)[ , obj$m+1]

getNUsed <- function(obj, ...)
	UseMethod("getNUsed")

getNUsed.nlar.struct <- function(obj, ...)
	length(obj$x)

#non-linear autoregressive model fitting
#str: result of a call to nlar.struct
nlar <- function(str, coefficients, fitted.values, residuals, k,
                 model.specific=NULL, ...) {
  return(extend(list(), "nlar",
                str=str,
                coefficients = coefficients,
                fitted.values= fitted.values,
                residuals = residuals,
                k=k,
                model.specific=model.specific,
                ...
                ))
}

#Print nlar object
print.nlar <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nNon linear autoregressive model\n")
  invisible(x)
}

#Coefficients of a nlar.fit object
coef.nlar <- function(object, ...)
  object$coefficients

#Fitted values for the fitted nlar object
fitted.nlar <- function(object, ...) {
  ans <- c(rep(NA, object$str$n.used - length(object$fitted.values)), object$fitted.values)
  tsp(ans) <- tsp( object$str$x )
  ans <- as.ts(ans)
  ans
}

#Observed residuals for the fitted nlar object
residuals.nlar <- function(object, ...) {
  str <- object$str
  data <- str$x
  ans <- c(rep(NA, str$n.used - length(object$residuals) ), object$residuals)
  tsp(ans) <- tsp(data)
  ans <- as.ts(ans)
  ans
}

#Mean Square Error for the specified object
mse <- function (object, ...)  
  UseMethod("mse")

mse.default <- function(object, ...)
  NULL

mse.nlar <- function(object, ...)
  sum(object$residuals^2)/object$str$n.used

#AIC for the fitted nlar model
AIC.nlar <- function(object, ...){
  n <- object$str$n.used
  k <- object$k
  n * log( mse(object) ) + 2 * k
}

#Mean Absolute Percent Error
MAPE <- function(object, ...)
  UseMethod("MAPE")

MAPE.default <- function(object, ...)
  NULL

MAPE.nlar <- function(object, ...) {
  e <- abs(object$residuals/object$str$yy)
  mean( e[is.finite(e)] )
}

#Computes summary infos for the fitted nlar model
summary.nlar <- function(object, ...) {
  ans <- list()
  ans$object <- object
  ans$mse <- mse(object)
  ans$AIC <- AIC(object)
  ans$MAPE <- MAPE(object)
  ans$df <- object$k
  ans$residuals <- residuals(object)
  return(extend(list(), "summary.nlar", listV=ans))
}

#Prints summary infos for the fitted nlar model
print.summary.nlar <- function(x, ...) {
  print(x$object)
  cat("\nResiduals:\n")
  rq <- structure(quantile(x$residuals, na.rm=TRUE), names = c("Min","1Q","Median","3Q","Max"))
  print(rq, ...)
  cat("\nFit:\n")
  cat("residuals variance = ", format(x$mse, digits = 4), 
      ",  AIC = ", format(x$AIC, digits=2), ", MAPE = ",format(100*x$MAPE, digits=4),"%\n", sep="")
  invisible(x)
}

plot.nlar <- function(x, ask = interactive(), ...) {
  str <- x$str
  op <- par(no.readonly=TRUE)
  par(ask = ask, mfrow = c(2, 1), no.readonly=TRUE)
  series <- str$series
  data <- str$x
  plot(data, main = series, ylab = "Series")
  plot(residuals(x), main = "Residuals", ylab = "Series")
  tmp <- c(acf(data, na.action=na.remove, plot=FALSE)$acf[-1,,],
           acf(residuals(x), na.action=na.remove, plot=FALSE)$acf[-1,,])
  ylim <- range(tmp)
  acf.custom(data, main = paste("ACF of",series), na.action = na.remove, ylim=ylim)
  acf.custom(residuals(x), main = "ACF of Residuals", na.action = na.remove, ylim=ylim)
  tmp <- c(pacf(data, na.action=na.remove, plot=FALSE)$acf[-1,,],
           pacf(residuals(x), na.action=na.remove, plot=FALSE)$acf[-1,,])
  ylim<-range(tmp)
  pacf.custom(data, main = paste("PACF of",series), ylim=ylim)
  pacf.custom(residuals(x), main = "PACF of Residuals", na.action = na.remove, ylim=ylim)
  partitions <- length(na.remove(data))^(1/3)
  partitions <- max(2, partitions)
  tmp <- c(mutual(na.remove(data), partitions=partitions, plot=FALSE)[-1],
           mutual(na.remove(residuals(x)), partitions=partitions, plot=FALSE)[-1])
  ylim <- c(0,max(tmp))
  mutual.custom(na.remove(data), partitions=partitions, 
		main=paste("Average Mutual Information of",series), ylim=ylim)
  mutual.custom(na.remove(residuals(x)), partitions=partitions, 
		main="Average Mutual Information of residuals", ylim=ylim)
  par(op)
  invisible(x)
}

predict.nlar <- function(object, newdata, n.ahead=1, simulate=FALSE, ...)
{

  if(missing(newdata)) 
    newdata <- object$str$x

  res <- newdata
  n.used <- length(res)
  m <- object$str$m
  d <- object$str$d
  sd <- sqrt( mse(object) )
  steps <- object$str$steps
  
  tsp(res) <- NULL
  class(res) <- NULL

  res <- c(res, rep(0, n.ahead))
  xrange <- (m-1)*d + steps - ((m-1):0)*d

  for(i in (n.used + 1:n.ahead)) {
    res[i] <- oneStep(object, newdata = t(as.matrix(res[i - xrange])),
                      itime=(i - n.used), ...)
    if(simulate)
      res[i] <- res[i] + rnorm(1, mean=0, sd=sd)
  }

  pred <- res[n.used + 1:n.ahead]
  pred <- ts(pred, start = tsp(newdata)[2] + deltat(newdata),
             frequency=frequency(newdata))
  return(pred)
  
}

oneStep <- function(object, newdata, ...)
  UseMethod("oneStep")

toLatex.nlar <- function(object, ...) {
  obj <- object
  str <- obj$str
  m <- str$m
  d <- str$d
  steps <- str$steps
  res <- character()
  res[1] <- "\\["
  res[2] <- paste("X_{t+",steps,"} = F( X_{t}",sep="")
  if(m>1) for(j in 2:m)
    res[2] <- paste(res[2], ", X_{t-",(j-1)*d,"}",sep="")
  res[2] <- paste(res[2], " )")
  res[3] <- "\\]"
  res[4] <- ""
  return(structure(res, class="Latex"))
}


addRegime <- function(object, ...)
  UseMethod("addRegime")

linearityTest <- function(object, ...)
	UseMethod("linearityTest")

# LM linearity testing against 2 regime STAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   str: an nlar.struct object
#   rob
#   sig
linearityTest.nlar.struct <- function(str, thVar, externThVar=FALSE,
                                      rob=FALSE, sig=0.05, ...)
{

  T <- NROW(str$xx);  # The number of lagged samples

  # Build the regressand vector
  y_t <- str$yy;
  
  # Regressors under the null
  xH0 <- cbind(1, str$xx)

  # Get the transition variable
  s_t <- thVar

  # Linear Model (null hypothesis)
  linearModel <- lm(y_t ~ . , data=data.frame(xH0))
  u_t <- linearModel$residuals;
  SSE0 <- sum(u_t^2)

  # Regressors under the alternative
  if (externThVar) {
#    if (rob) {} else {
    tmp <- rep(s_t, NCOL(str$xx) + 1)
    dim(tmp) <- c(length(s_t), NCOL(str$xx) + 1)
    xH1 <- cbind(cbind(1, str$xx) * tmp, cbind(1, str$xx) * (tmp^2),
                 cbind(1, str$xx) * (tmp^3))
  } else {
#    if (rob) {} else {
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
  pValue <- 1 - pf(F, m, T - m - n, lower.tail = FALSE);

  if (pValue < sig) {
    return(list(isLinear = TRUE, pValue = pValue));
  }
  else {
    return(list(isLinear = FALSE, pValue = pValue));
  }

}

