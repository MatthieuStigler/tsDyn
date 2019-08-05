library(tsDyn)
suppressMessages(library(tidyverse))

###SETAR

mod.setar <- setar(log10(lynx), m=2, thDelay=1, th=3.25)
mod.setar
summary(mod.setar)

set0 <- setar(lynx,m=2)
set1 <- setar(lynx,m=2, th=1426,trace=TRUE,thDelay=1)
set2 <- setar(lynx,m=2, th=c(600,1426),trace=TRUE)
set3 <- setar(lynx,m=2, trace=TRUE,nthresh=1)
set4 <- setar(lynx,m=2, trace=TRUE,nthresh=2, trim=0.05)
set5 <- setar(lynx,m=2, th=c(600,1426),include="trend")
set6 <- setar(lynx,m=2, th=c(600,1426),include="both")
set7 <- setar(lynx,m=2, th=c(600,1426),include="none")

set_all <- list(set0=set0, set1=set1, set2=set2, set3=set3, set4=set4,
                set5=set5, set6=set6, set7=set7)

set_all
sapply(set_all, function(x) c(AIC=AIC(x), BIC=BIC(x), deviance=deviance(x)))
lapply(set_all, summary)
sapply(set_all, function(x) tail(coef(x),4))
sapply(set_all, function(x) tail(coef(x, hyperCo=FALSE),4))
lapply(set_all, function(x) round(vcov(x),9))

lapply(set_all, function(x) window(residuals(x), end=c(1830)))
lapply(set_all, function(x) window(fitted(x), end=c(1830)))

lapply(set_all, function(x) window(regime(x), end=c(1830)))
lapply(set_all, function(x) window(regime(x, initVal=FALSE), end=c(1830)))
sapply(set_all, function(x) head(regime(x, time=FALSE),10))
sapply(set_all, function(x) head(regime(x, initVal=FALSE, time=FALSE),10))

sapply(set_all, function(x) head(x$model,2))

lapply(set_all[c(1,2,4)], toLatex)

## coefs
sapply(set_all, function(x) tail(coef(x, regime = "L"),4))
sapply(set_all, function(x) tail(coef(x, regime = "H"),4))


###SelectSetar
selectSETAR(lynx, m=2, d=1,  trace=TRUE, include = "const", common="none", model="TAR", nthresh=1,trim=0.15,criterion = "SSR", plot=FALSE,max.iter=2)
selectSETAR(lynx, m=2, d=1,  trace=TRUE, include = "const", common="none", model="TAR", nthresh=2,trim=0.15,criterion = "SSR", plot=FALSE,max.iter=3)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(ngrid="Half"), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(ngrid=20), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(exact=1388), common="none",trim=0.15,criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(around=1388, ngrid=20), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824)), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824), ngrid=5), common="none",criterion = "SSR", plot=FALSE)
selectSETAR(lynx, m=2, d=1, th=MakeThSpec(int=c(957, 1824), ngrid="Half"), common="none",criterion = "SSR", plot=FALSE)


## setar.sim
set_all_boot <- lapply(set_all, function(x) setar.boot(x, boot.scheme = "check", returnStarting = TRUE, round_digits=10))
set_all_boot_check <-  sapply(set_all_boot, function(x) all.equal(x, as.numeric(lynx)))

vals_df <- bind_cols(set_all_boot) %>%  
  mutate(n_row = 1:n()) %>% 
  gather(model, value, -n_row) %>% 
  mutate(original = rep(as.numeric(lynx), length(set_all_boot)),
         diff = abs(value-original) < 1e-010)

vals_df %>% 
  filter(value!=original) %>% 
  group_by(model) %>%
  mutate(rmse = sqrt(mean(value-original)^2)) %>% 
  filter(n_row ==min(n_row)) %>% 
  dplyr::select(-rmse)
