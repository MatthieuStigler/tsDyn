library(tsDyn)

mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar
summary(mod.lstar)
deviance(mod.lstar)
c(AIC(mod.lstar),BIC(mod.lstar))

mod.lstar2 <- lstar(log10(lynx), m=1, control=list(maxit=3000))
mod.lstar2
summary(mod.lstar2)
deviance(mod.lstar2)
c(AIC(mod.lstar2),BIC(mod.lstar2))
