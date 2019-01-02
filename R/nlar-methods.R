#' NLAR methods
#' 
#' Generic \sQuote{nlar} methods. Method \sQuote{nlar} is described in a
#' separate page: \code{\link{nlar}}
#' 
#' \describe{ 
#'   \item{MAPE}{ Mean Absolute Percent Error } 
#'   \item{mse}{ Mean Square Error } 
#'   \item{plot}{ Diagnostic plots } }
#' 
#' @param x,object fitted \sQuote{nlar} object
#' @param ask graphical option. See \code{\link{par}}
#' @param digits For print method, see \code{\link{printCoefmat}}.
#' @param label LaTeX label passed to the equation
#' @param \dots further arguments to be passed to and from other methods
#' @author Antonio, Fabio Di Narzo
#' @seealso \code{\link{availableModels}} for listing all currently available
#' models.
#' @keywords ts
#' @examples
#' 
#' x <- log10(lynx)
#' mod.setar <- setar(x, m=2, thDelay=1, th=3.25)
#' mod.setar
#' AIC(mod.setar)
#' mse(mod.setar)
#' MAPE(mod.setar)
#' coef(mod.setar)
#' summary(mod.setar)
#' 
#' e <- residuals(mod.setar)
#' e <- e[!is.na(e)]
#' plot(e)
#' acf(e)
#' 
#' plot(x)
#' lines(fitted(mod.setar), lty=2)
#' legend(x=1910, y=3.9,lty=c(1,2), legend=c("observed","fitted"))
#' 
#' plot(mod.setar)
#'
#' @name nlar-methods
NULL

## Copyright (C) 2005, 2006, 2007/2006, 2008, 2018  Antonio, Fabio Di Narzo
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

getdXXYY <- function(obj, ...) UseMethod("getdXXYY")

getdXXYY.nlar.struct <- function(obj,same.dim=FALSE, ...) {
	x <- obj$x
	m<-if(same.dim) obj$m-1 else obj$m
	d <- obj$d
	steps <- obj$steps
	embedd(x, lags=c((0:(m-1))*(-d), steps) )
}

getdXX <- function(obj, ...)
	diff(getdXXYY(obj,...))[ , 1:obj$m , drop=FALSE]

getdYY <- function(obj, ...)
	diff(getdXXYY(obj, ...))[ , obj$m+1]

getdX1 <- function(obj, ...)
	getdXXYY(obj,...)[ -1, 1, drop=FALSE]

getNUsed <- function(obj, ...)
	UseMethod("getNUsed")

getNUsed.nlar.struct <- function(obj, ...)
	length(obj$x)

getNUsed.nlar <- function(obj, ...)
	length(obj$str$x)
	
#non-linear autoregressive model fitting
#str: result of a call to nlar.struct


#'Non-linear time series model, base class definition
#'
#'Generic non-linear autogregressive model class constructor.
#'
#'Constructor for the generic \code{nlar} model class. On a fitted object you
#'can call some generic methods. For a list of them, see
#'\code{\link{nlar-methods}}.
#'
#'An object of the \code{nlar} class is a list of (at least) components:
#'\describe{ \item{str}{ \code{\link{nlar.struct}} object, encapsulating
#'general infos such as time series length, embedding parameters, forecasting
#'steps, model design matrix } \item{coefficients}{ a named vector of model
#'estimated/fixed coefficients } \item{k}{ total number of estimated
#'coefficients } \item{fitted.values}{ model fitted values } \item{residuals}{
#'model residuals } \item{model}{ data frame containing the variables used }
#'\item{model.specific}{ (optional) model specific additional infos} }
#'
#'A \code{\link{nlar}} object normally should also have a model-specific
#'subclass (i.e., \code{nlar} is a virtual class).
#'
#'Each subclass should define at least a \code{print} and, hopefully, a
#'\code{oneStep} method, which is used by \code{\link{predict.nlar}} to
#'iteratively extend ahead the time series.
#'
#'@param str a \code{nlar.struct} object, i.e. the result of a call to
#'\code{\link{nlar.struct}}
#'@param coefficients,fitted.values,residuals,k,model,model.specific internal
#'structure
#'@param \dots further model specific fields
#'@return An object of class \code{nlar}. \link{nlar-methods} for a list of
#'available methods.
#'@author Antonio, Fabio Di Narzo
#'@seealso \code{\link{availableModels}} for currently available built-in
#'models.  \link{nlar-methods} for available \code{nlar} methods.
#'@references Non-linear time series models in empirical finance, Philip Hans
#'Franses and Dick van Dijk, Cambridge: Cambridge University Press (2000).
#'
#'Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford:
#'Oxford University Press (1990).
#'@keywords ts internal
#'@export
nlar <- function(str, coefficients, fitted.values, residuals, k, model,
                 model.specific=NULL, ...) {
  return(extend(list(), "nlar",
                str=str,
                coefficients = coefficients,
                fitted.values= fitted.values,
                residuals = residuals,
                k=k,
		model=model,
                model.specific=model.specific,
                ...
                ))
}


#' @export
#Print nlar object
print.nlar <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nNon linear autoregressive model\n")
  invisible(x)
}

#'@rdname nlar-methods
#' @export
#Coefficients of a nlar.fit object (actually only for aar, as setar and lstar have dedicated one)
coef.nlar <- function(object,  ...){
  co <- object$coefficients
  co
}
  

#'@rdname nlar-methods
#' @export
#Fitted values for the fitted nlar object
fitted.nlar <- function(object, ...) {
  ans <- c(rep(NA, object$str$n.used - length(object$fitted.values)), object$fitted.values)
  tsp(ans) <- tsp( object$str$x )
  ans <- as.ts(ans)
  ans
}

#'@rdname nlar-methods
#' @template param_timeAttr
#' @template param_initVal
#' @export
#Observed residuals for the fitted nlar object
residuals.nlar <- function(object, initVal=TRUE, timeAttr = TRUE, ...) {
  str <- object$str
  data <- str$x
  resids <- object$residuals
  n_init <- str$n.used - length(resids)
  
  if(initVal) {
    res <- c(rep(NA,  n_init), resids)  
  } else {
    res <- resids
  }
  
  if(timeAttr) {
    if(!initVal) {
      data <- tail(data, -n_init)
    }
    tsp(res) <- tsp(data)
    res <- as.ts(res)
  }
  
  res
}


if(FALSE) {
  library(tsDyn)
  linear_l2_none <- linear(lh, m = 2, include = "none")
  residuals(linear_l2_none)
  residuals(linear_l2_none, timeAttr = FALSE)
  residuals(linear_l2_none, initVal = FALSE)
  residuals(linear_l2_none, initVal = FALSE, timeAttr = FALSE)
  
  linear_l2_none_num <- linear(as.numeric(lh), m = 2, include = "none")
  residuals(linear_l2_none_num)
  residuals(linear_l2_none_num, timeAttr = FALSE)
  residuals(linear_l2_none_num, initVal = FALSE)
  residuals(linear_l2_none_num, initVal = FALSE, timeAttr = FALSE)
}

#get the threshold for setar and nlVar


#'Extract threshold(s) coefficient
#'
#'Extract threshold coefficient(s)
#'
#'
#'@aliases getTh
#'@param object object of class \code{setar}, \code{summary.setar},
#'\code{nlVar}
#'@param \dots additional arguments to \code{getTh}
#'@return Threshold value.
#'@author Matthieu Stigler
#'@keywords ts
#'@export
#'@examples
#'
#'set<-setar(lynx, m=3)
#'getTh(set)
#'getTh(summary(set))
#'
getTh<- function (object, ...)  
  UseMethod("getTh")

#' @rdname getTh
#' @export
getTh.default <- function(object, ...){
  # look first just in object
  allth<-object[grep("th",names(object))]
  # look then in coef(object)
  if(length(allth)==0){
    allth<-coef(object)[grep("th",names(coef(object)))]
  }
  if(length(allth)==0) allth <- NULL
  # remove thD if there
  if(any(grepl("thD",names(allth))))
    allth<-allth[-grep("thD",names(allth))]
  return(allth)
}

#' @export
getTh.setar<-function(object,...){
  object<-object$coefficients
  getTh.default(object)
}

#' @export
getTh.lstar<-function(object,...){
  object<-object$coefficients
  getTh.default(object)
}

getTh.lstar<-function(object,...){
  object <-object$coef
  object["th"]
}

#' @export
getTh.summary.setar<-function(object,...){
  object$th
}

#' @export
getTh.nlVar<-function(object,...){
  nth <- object$model.specific$nthresh

  if(nth>0){
    th<-object$model.specific$Thresh
    if(length(th)==1){
      names(th)<-"th"
    } else{
      names(th)<-c("th1", "th2")
    }
  } else {
    th <- NULL
  }

  return(th)
}


#'@rdname nlar-methods
#' @export
#Observed residuals for the fitted nlar object
deviance.nlar<-function(object,...) crossprod(object$residuals)

#Mean Square Error for the specified object


#'Mean Square Error
#'
#'Generic function to compute the Mean Squared Error of a fitted model.
#'
#'
#'@aliases mse
#'@param object object of class \code{nlar.fit}
#'@param \dots additional arguments to \code{mse}
#'@return Computed MSE for the fitted model.
#'@author Antonio, Fabio Di Narzo
#'@keywords ts
#'@export
mse <- function (object, ...)  
  UseMethod("mse")


#' @rdname mse
#' @export
mse.default <- function(object, ...)
  NULL

#' @rdname nlar-methods
#' @export
mse.nlar <- function(object, ...)
  sum(object$residuals^2)/object$str$n.used

#' @rdname nlar-methods
#' @param k numeric, the penalty per parameter to be used for AIC/BIC; the default k = 2 is
#' the classical AIC
#' @export
#AIC for the fitted nlar model
AIC.nlar <- function(object,k=2, ...){
  n <- object$str$n.used
  npar <- object$k
  n * log( mse(object) ) + k * npar
}

#' @rdname nlar-methods
#' @export
#BIC for the fitted nlar model
BIC.nlar <- function(object, ...)
	AIC.nlar(object, k=log(getNUsed(object)))



#Mean Absolute Percent Error


#'Mean Absolute Percent Error
#'
#'Generic function to compute the Mean Absolute Percent Error of a fitted
#'model.
#'
#'
#'@aliases MAPE
#'@param object object of class \code{nlar.fit}
#'@param \dots additional arguments to \code{MAPE}
#'@return Computed Mean Absolute Percent Error for the fitted model.
#'@author Antonio, Fabio Di Narzo
#'@keywords ts
#' @export
MAPE <- function(object, ...)
  UseMethod("MAPE")

#' @rdname MAPE
#' @export
MAPE.default <- function(object, ...)
  NULL

#' @rdname nlar-methods
#' @export
MAPE.nlar <- function(object, ...) {
  e <- abs(object$residuals/object$str$yy)
  mean( e[is.finite(e)] )
}

#'@rdname nlar-methods
#' @export
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

#' @export
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

#'@rdname nlar-methods
#' @export
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



#'oneStep
#'
#'Doing one step forward within a \code{NLAR} model
#'
#'If in \code{object} is encapsulated the \code{NLAR} map, say, \code{F(X[t],
#'X[t-d], ..., X[t-(m-1)d])}, this function should return the value of \code{F}
#'(with already fitted parameters) applied to given new data, which can be a
#'single vector of length \code{m} or a matrix with \code{m} columns.
#'
#'@param object fitted \sQuote{nlar} object
#'@param newdata data from which to step forward
#'@param \dots further arguments to be passed to and from other methods
#'@return Computed value(s)
#'@note This is an internal function, and should not be called by the user
#'@author Antonio, Fabio Di Narzo
#'@keywords internal ts
#'
#'tsDyn:::oneStep.linear
#'
oneStep <- function(object, newdata, ...)
  UseMethod("oneStep")


#'@rdname nlar-methods
#'@export
toLatex.nlar <- function(object, digits, label,...) {
  obj <- object
  str <- obj$str
  m <- str$m
  d <- str$d
  steps <- str$steps
  res <- character()
  res[1] <- "\\["
  if(!missing(label)) res[1]<- paste(res[1], "\\label{", label, "}", sep="")
  res[2] <- paste("X_{t+",steps,"} = F( X_{t}",sep="")
  if(m>1) for(j in 2:m)
    res[2] <- paste(res[2], ", X_{t-",(j-1)*d,"}",sep="")
  res[2] <- paste(res[2], " )")
  res[3] <- "\\]"
  res[4] <- ""
  return(structure(res, class="Latex"))
}


# LM linearity testing against 2 regime STAR
#
#   Performs an 3rd order Taylor expansion LM test
#
#   str: an nlar.struct object
#   rob
#   sig
linearityTest.star <- function(str, thVar, externThVar=FALSE,
                                      rob=FALSE, sig=0.05, trace=TRUE, ...)
{
  n.used <- NROW(str$xx);  # The number of lagged samples

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
  sdZ <- apply(Z,2,sd)
  dim(sdZ) <- c(1, nZ)
  sdZ <- kronecker(matrix(1,n.used,1), sdZ) # repeat sdZ n.used rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

  # Nonlinear model (alternative hypothesis)
  nonlinearModel <- lm(u_t ~ ., data=data.frame(Z));
  e_t <- nonlinearModel$residuals;
  SSE1 <- sum(e_t^2)

  # Compute the test statistic
  n <- NCOL(xH0);
  m <- NCOL(xH1);
  
  F = ((SSE0 - SSE1) / m) / (SSE1 / (n.used - m - n));
  
  # Look up the statistic in the table, get the p-value
  pValue <- pf(F, m, n.used - m - n, lower.tail = FALSE);
 
  if (pValue >= sig) {
    return(list(isLinear = TRUE, pValue = pValue));
  }
  else {
    return(list(isLinear = FALSE, pValue = pValue));
  }

}

