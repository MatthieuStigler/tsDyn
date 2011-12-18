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

#Neural Network fitter
nnetTs <- function(x, m, d=1, steps=d, series, size, control=list(trace=FALSE)) {
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	args <- c(list(getXX(str), getYY(str), size=size, linout=TRUE), control)
	res <- do.call(nnet::nnet, args)
	res$k <- length(res$wts)
	res$fitted <- res$fitted.values	
	return(extend(nlar(str,
		coefficients = res$wts,
		residuals=res$residuals,
		fitted.values=res$fitted,
		k=res$k,
		model=NULL,
		model.specific=res), "nnetTs"))
}

print.nnetTs <- function(x, ...) {
	NextMethod(...)
	cat("\nNNET time series model\n")
	nnet:::print.nnet(x$model.specific, ...)
	invisible(x)
}

oneStep.nnetTs <- function(object, newdata, ...)
	nnet:::predict.nnet(object$model.specific, newdata)

selectNNET <- function(x, m, d=1, steps=d, size=1:(m+1), maxit=1e3) {
	IDS <- as.matrix( size )
	colnames(IDS) <- c("size")
	computedAIC <- function(j) 
		AIC( nnetTs(x=x, m=m, d=d, steps=steps, size=j, control=list(maxit=maxit)) )
	computedAIC <- apply(IDS, 1, computedAIC)
	res <- cbind(IDS, AIC = computedAIC)
	idSel <- sort(computedAIC, index=TRUE)$ix
	idSel <- idSel[1:min(10, length(idSel))]
	res <- data.frame(res[idSel,], row.names=NULL)
	return(res)
}

showDialog.nnetTs <- function(x, ...) {
	frRoot <- Frame()
	vM <- tclVar(1)
	vD <- tclVar(1)
	vSteps <- tclVar(1)
	vSize <- tclVar(1)
	onFinish <- function() {
		res <- nnetTs(x, m=as.numeric(tclObj(vM)), d=as.numeric(tclObj(vD)), steps=as.numeric(tclObj(vSteps)) , size=as.numeric(tclObj(vSize)))
		tkdestroy(frRoot$tkvar)
		assign("nlarModel", res, .GlobalEnv)
	}
	onCancel <- function()
		tkdestroy(frRoot$tkvar)
	frMain <- nlar.struct.Frame(vM, vD, vSteps)
	add(frMain,
		Widget(opts=list(type="label", text="Hidden units")),
		Widget(opts=list(type="spinbox", from=1, to=100, increment=1, textvariable=vSize, width=4))
	)
	add(frRoot,
		frMain,
		makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
	)
	buildDialog(title="linear model", frRoot)
}
