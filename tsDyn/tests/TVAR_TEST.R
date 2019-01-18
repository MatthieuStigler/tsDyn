library(tsDyn)

data(zeroyld)

###TVAR
tvar<-TVAR(zeroyld[1:100,], lag=2, nthresh=2,thDelay=1,trim=0.1, plot=FALSE, include="const")
class(tvar)
tvar
print(tvar)
coefficients(tvar)

##FIXME ?
summary(tvar)$VAR
tvar$VAR



coefficients(summary(tvar))
logLik(tvar)
AIC(tvar)
BIC(tvar)
coef(tvar)
deviance(tvar)
vcov(tvar)
df.residual(tvar)
all.equal(summary(tvar)$VarCov[[1]], vcov(tvar))
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

if(FALSE) {##FIXME
  toLatex(summary(tvar), digits=2)
  tvar$coefficients
  tvar$StDev
  options(show.signif.stars=FALSE)
  toLatex(summary(tvar), digits=2)
}


#################################
### Many more TVARs
#################################
gr_args <- expand.grid(include= c("const", "trend", "none", "both"), 
                       nthresh=1:2,
                       lag=1:3,
                       commonInter=c(FALSE, TRUE),
                       thDelay=c(1,2),
                       mTh=c(1,2),
                       stringsAsFactors=FALSE)

## remove if commonInter and const
gr_args_ok <- subset(gr_args, !(commonInter & include!="const") & lag > thDelay) 


## run
TVAR_outs <- mapply(TVAR,
                    include=gr_args_ok$include,
                    lag=gr_args_ok$lag,
                    commonInter=gr_args_ok$commonInter,
                    thDelay=gr_args_ok$thDelay,
                    nthresh=gr_args_ok$nthresh,
                    mTh=gr_args_ok$mTh,
                    MoreArgs = list(data=zeroyld[1:100,], plot=FALSE, trace=FALSE),
                    SIMPLIFY=FALSE)

names(TVAR_outs) <- paste("mod", 1:length(TVAR_outs), sep="_")

## check outputs
lapply(TVAR_outs, print)
lapply(TVAR_outs, getTh)
lapply(TVAR_outs, summary)
sapply(TVAR_outs, AIC)
sapply(TVAR_outs, BIC)
sapply(TVAR_outs, logLik)
sapply(TVAR_outs, deviance)
sapply(TVAR_outs, df.residual)


## Boot
gr_args_ok_df <- gr_args_ok
gr_args_ok_df$models <- TVAR_outs 

gr_args_boot <- subset(gr_args_ok_df, !commonInter)

# sapply(gr_args_boot$models, function(x) head(TVAR.boot(x, seed=123)))

boots_runs <- lapply(gr_args_boot$models, function(x) TVAR.boot(x, boot.scheme="check", round_digits=5))
boots_checks <- sapply(boots_runs, function(x) all.equal(as.matrix(zeroyld[1:100,]), x, check.attributes = FALSE))

gr_args_boot$boots_checks <- boots_checks
gr_args_boot[,-7]


## predict
predict_runs <- lapply(gr_args_boot$models, function(x) predict(x, n.ahead=3))
