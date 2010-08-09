library(tsDyn)

data(zeroyld)
dat<-zeroyld

###SETAR

mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)
mod.setar
summary(mod.setar)

setar(lynx,m=2, th=1426,trace=TRUE,thDelay=1)
setar(lynx,m=2, th=c(600,1426),trace=TRUE)
setar(lynx,m=2, trace=TRUE,nthresh=1)
setar(lynx,m=2, trace=TRUE,nthresh=2, trim=0.05)

mod.set<-setar(lynx,m=2)
deviance(mod.set)
AIC(mod.set)
residuals(mod.set)
fitted(mod.set)
coef(mod.set)
regime(mod.set)
regime(mod.set, initVal=FALSE)
regime(mod.set, time=FALSE)
regime(mod.set, time=FALSE, initVal=FALSE)


###SelectSetar
selectSETAR(lynx, m=2, d=1,  trace=TRUE, include = "const", common="none", model="TAR", nthresh=1,trim=0.15,criterion = "SSR",thSteps = 7,  plot=FALSE,max.iter=2)
selectSETAR(lynx, m=2, d=1,  trace=TRUE, include = "const", common="none", model="TAR", nthresh=2,trim=0.15,criterion = "SSR",thSteps = 7,  plot=FALSE,max.iter=3)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(ngrid="Half"), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(ngrid=20), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(exact=1388), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(around=1388, ngrid=20), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824)), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824), ngrid=5), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824), ngrid="Half"), common="none",criterion = "SSR", plot=FALSE)


###TVAR
tvar<-TVAR(dat[1:100,], lag=2, nthresh=2,thDelay=1,trim=0.1, plot=FALSE, include="const")
class(tvar)
tvar
print(tvar)
coefficients(tvar)
##FIXME
summary(tvar)$VAR
##FIXME
tvar$VAR

coefficients(summary(tvar))
logLik(tvar)
AIC(tvar)
BIC(tvar)
coef(tvar)
deviance(tvar)
head(residuals(tvar))
tail(residuals(tvar))
head(fitted(tvar))
tail(fitted(tvar))

regime(tvar)
regime(tvar, initVal=FALSE)
regime(tvar, time=FALSE)
regime(tvar, time=FALSE, initVal=FALSE)


##FIXME
options(show.signif.stars=TRUE)
summary(tvar)

options(show.signif.stars=FALSE)
summary(tvar)

print(summary(tvar), digits=3)

toLatex(tvar)

if(0) {##FIXME
toLatex(summary(tvar), digits=2)
tvar$coefficients
tvar$StDev
options(show.signif.stars=FALSE)
toLatex(summary(tvar), digits=2)
}

###TVECM
tvecm<-TVECM(dat, nthresh=2,lag=1, ngridBeta=20, ngridTh=20, plot=FALSE,trim=0.05, common="All")
class(tvecm)
tvecm
print(tvecm)
coef(tvecm)
logLik(tvecm)
AIC(tvecm)
BIC(tvecm)
deviance(tvecm)
head(residuals(tvecm))
tail(residuals(tvecm))
head(fitted(tvecm))
tail(fitted(tvecm))
summary(tvecm)

regime(tvecm)
regime(tvecm, initVal=FALSE)
regime(tvecm, time=FALSE)
regime(tvecm, time=FALSE, initVal=FALSE)



toLatex(tvecm)
options(show.signif.stars=FALSE)
toLatex(summary(tvecm))

###Linear
lin<-lineVar(dat,lag=2)
class(lin)
lin
print(lin)
logLik(lin)
AIC(lin)
BIC(lin)
deviance(lin)
coef(lin)
summary(lin)
toLatex(lin)
toLatex(summary(lin))

linVECM<-lineVar(dat,lag=2, model="VECM")
class(linVECM)
linVECM
print(linVECM)
logLik(linVECM)
AIC(linVECM)
BIC(linVECM)
deviance(linVECM)
coef(linVECM)
summary(linVECM)
toLatex(linVECM)
toLatex(summary(linVECM))
