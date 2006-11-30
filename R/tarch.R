## Copyright (C) 2006  Antonio, Fabio Di Narzo
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

#Treshold-ARCH model fitting, of order p
tarch <- function(x, m, d=1, steps=d, series, coef, thDelay=0, control=list(), ...) {
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- str$xx
	yy <- str$yy
  z <- xx[,thDelay+1]
  k <- ncol(xx)*2+2
  if(missing(coef))
    coef <- exp(rnorm(k))
  yy2 <- yy^2
  xx <- xx^2
  xx <- cbind(1,xx)
  isL <- z <= 0
  xxL <- xx[isL,]
  xxH <- xx[!isL,]
  yy2L <- yy2[isL]
  yy2H <- yy2[!isL]

  LInd <- 1:ncol(xx)
  HInd <- ncol(xx)+LInd
  lik <- function(coef) {
    sigma2L <- (xxL %*% coef[LInd])
    sigma2H <- (xxH %*% coef[HInd])
    rsL <- sum(log(sigma2L) + yy2L/sigma2L)
    rsH <- sum(log(sigma2H) + yy2H/sigma2H)
    rs <- rsL + rsH
    return(rs)
  }
  lbound <- rep(0,k)
  lbound[c(1,k/2+1)] <- .Machine$double.neg.eps
  res <- optim(coef, lik, control=control, method="L-BFGS-B", lower=lbound, hessian=TRUE, ...)
  coef <- res$par
  k <- length(coef)
  p <- k/2-1
  names(coef) <- c(paste("b0",0:p,sep="."), paste("b1",0:p,sep="."))
  sigma2 <- rep(NA, NROW(xx))
  sigma2[isL] <- xxL %*% coef[LInd]
  sigma2[!isL] <- xxH %*% coef[HInd]
  res$thDelay <- thDelay
  res$lowRegProp <- mean(isL)
  res$coefficients <- coef
  res$fitted.values <- sqrt(sigma2)
  res$residuals <- yy/res$fitted.values
  res$k <- k
  return(structure(list(str=str, model.specific=res), class="tarch"))
}

print.tarch <- function(x, digits = max(3, getOption("digits") - 3), ...) {
	NextMethod(digits=digits, ...)
  p <- x$k/2 - 1
  LInd <- 1:(p+1)
  HInd <- p + 1 + LInd
  coef <- x$coef
  cat("Treshold-ARCH(",p,") Model\n",sep="")
  cat("\nModel coefficients:\n")
  cat("Low regime:\n")
  print(coef[LInd], ...)
  cat("\nHigh regime:\n")
  print(coef[HInd], ...)
  cat("\nThreshold")
  cat("\nVariable: ")
  cat('Z(t) = X(t-', x$thDelay*x$d,')\n',sep="")
  cat("Proportion of points in low regime: ", format(x$lowRegProp*100, digits=3), "%\n", sep="")  
  cat("\n")
  invisible(x)
}

summary.tarch <- function(object, ...) {
	ans <- list()
  e <- object$model$residuals
  H <- object$model$hessian
  se <- sqrt(diag(solve(H)))
  ans$residuals <- na.remove(e)
  tval <- object$model$coef/se
  ans$coef <- cbind(object$model$coef, se, tval, 2 * (1 - pnorm(abs(tval))))
  dimnames(ans$coef) <- list(names(object$model$coef), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
#  return(extend(summary.nlar(object,...), "summary.tarch", listV=ans))
	return(structure(ans, class="summary.tarch"))
}

print.summary.tarch <- function(x, digits = max(3, getOption("digits") - 3), ...) {
	#NextMethod(digits=digits, ...)
  p <- nrow(x$coef)/2-1
  cat("\nTreshold-ARCH(",p,") Model\n",sep="")
  cat("\nEstimated coefficients:\n")
  printCoefmat(x$coef, digits = digits, ...)
  invisible(x)
}

oneStep.tarch <- function(object, newdata, ...) {
  thVal <- newdata[,object$thDelay+1]
  p <- object$k/2
  LInd <- 1:p
  HInd <- p+1:p
  xx <- cbind(1, newdata^2)
  isL <- thVal<=0
  if(length(isL) == 1)
    return(
           sqrt(xx %*% coef[ ifelse(isL,LInd,Hind) ] )
           )
  res <- rep(NA, length(isL))
  res[isL] <- xx[isL,] %*% coef[LInd]
  res[!isL] <- xx[!isL,] %*% coef[HInd]
  return(sqrt(res))
}