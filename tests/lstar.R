library(tsDyn)

mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar
deviance(mod.lstar)

mod.lstar2 <- lstar(log10(lynx), m=1, control=list(maxit=3000))
mod.lstar2
deviance(mod.lstar2)
