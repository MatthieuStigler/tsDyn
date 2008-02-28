## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
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

#linear model fitter (via OLS)
#str: call to nlar.struct
linear <- function(x, m, d=1, steps=d, series) {
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- getXX(str)
	yy <- getYY(str)
	xx <- cbind(1,xx)
	colnames(xx) <- c("(intercept)", paste("phi",1:(ncol(xx)-1), sep="."))
	res <- lm.fit(xx, yy)
	return(extend(nlar(str, coefficients=res$coefficients, fitted.values=res$fitted.values,
		residuals=res$residuals, k=res$rank, model.specific=res), 
		"linear"))
}

print.linear <- function(x, ...) {
	NextMethod(...)
	cat("\nAR model\n")
	cat("Coefficients:\n")
	print(x$coef, ...)
	invisible(x)
}

summary.linear <- function(object, ...) {
	ans <- list()
	obj <- c(object, object$model.specific)
	Qr <- obj$qr
	n <- nrow(Qr$qr)
	p <- obj$rank
	resvar <- mse(object)*n/(n-p)
	p1 <- 1:p
	R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
	se <- sqrt(diag(R) * resvar)
	est <- obj$coefficients[Qr$pivot[p1]]
	tval <- est/se
	coef <- cbind(est, se, tval, 2*pt(abs(tval), n-p, lower.tail = FALSE))
	dimnames(coef) <- list(names(est), c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
	ans$coef <- coef
	return( extend(summary.nlar(object, ...), "summary.linear", listV=ans) )
}

print.summary.linear <- function(x, digits=max(3, getOption("digits") - 2),
	signif.stars = getOption("show.signif.stars"), ...) {
	NextMethod(...)
	cat("\nCoefficient(s):\n")
	printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
	invisible(x)
}

oneStep.linear <- function(object, newdata, ...) {
	cbind(1,newdata) %*% object$coefficients
}

showDialog.linear <- function(x, ...) {
	frRoot <- Frame()
	vM <- tclVar(1)
	vD <- tclVar(1)
	vSteps <- tclVar(1)
	onFinish <- function() {
		res <- linear(x, m=as.numeric(tclObj(vM)), d=as.numeric(tclObj(vD)), steps=as.numeric(tclObj(vSteps)) )
		tkdestroy(frRoot$tkvar)
		assign("nlarModel", res, .GlobalEnv)
	}
	onCancel <- function()
		tkdestroy(frRoot$tkvar)
	frMain <- nlar.struct.Frame(vM, vD, vSteps)
	add(frRoot,
		frMain,
		makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
	)
	buildDialog(title="linear model", frRoot)
}
