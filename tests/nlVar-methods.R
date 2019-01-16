library(tsDyn)

data(zeroyld)

###TVECM
tvecm <- TVECM(zeroyld, nthresh=2,lag=1, ngridBeta=20, ngridTh=20, plot=FALSE,trim=0.05, common="All")
tvar <- TVAR(zeroyld, nthresh=2,lag=1,  plot=FALSE, trim=0.05)
linVECM<-lineVar(zeroyld, lag=2, model="VECM")
lin<-lineVar(zeroyld,lag=2)


mods <- list(tvecm=tvecm, linVECM=linVECM, lin=lin, tvar=tvar)
mods_nonLIn <- mods[c("tvecm", "tvar")]



## Standard functions
sapply(mods, class)
sapply(mods, print)
sapply(mods, summary)

sapply(mods, coef)
sapply(mods, tsDyn:::coefMat.nlVar)
sapply(mods, tsDyn:::coefVec.nlVar)
sapply(mods_nonLIn, coef, regime = "L")
sapply(mods_nonLIn, coef, regime = "H")


sapply(mods, AIC)
sapply(mods, BIC)
sapply(mods, logLik)
sapply(mods, deviance)
sapply(mods, function(x) dim(residuals(x, initVal=FALSE)))
sapply(mods, function(x) dim(residuals(x, initVal=TRUE)))
sapply(mods, function(x) head(residuals(x), 3))
sapply(mods, function(x) head(residuals(x, initVal=TRUE), 3))
sapply(mods, function(x) tail(residuals(x), 3))

sapply(mods, function(x) head(fitted(x), 3))
sapply(mods, function(x) tail(fitted(x), 3))


## NOn linear functions
sapply(mods_nonLIn, function(x) head(regime(x), 3))
sapply(mods_nonLIn, function(x) tail(regime(x), 3))

sapply(mods_nonLIn, function(x) head(regime(x, initVal=FALSE), 3))
sapply(mods_nonLIn, function(x) tail(regime(x, initVal=FALSE), 3))

sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE), 3))
sapply(mods_nonLIn, function(x) head(regime(x, time=FALSE, initVal=FALSE), 3))


## toLatex
sapply(mods, toLatex)
options(show.signif.stars=FALSE)
sapply(mods, function(x) toLatex(summary(x)))

