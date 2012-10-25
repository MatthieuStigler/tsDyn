library(tsDyn)


data(zeroyld)
data(barry)




## Test a few VECM models
vecm_OLS_l1_co <-VECM(barry, lag=1)
vecm_OLS_l3_co <-VECM(barry, lag=3, include="const")
vecm_OLS_l3_co_betaGiven<-VECM(barry, lag=3, include="const", beta=c(0.1, -0.05))
vecm_OLS_l1_tr <-VECM(barry, lag=1, include="trend")
vecm_OLS_l1_bo <-VECM(barry, lag=1, include="both")
vecm_OLS_l1_no <-VECM(barry, lag=1, include="none")
vecm_OLS_l1_LRco <-VECM(barry, lag=1, LRinclude="const")
vecm_OLS_l1_LRtr <-VECM(barry, lag=1, LRinclude="trend")
vecm_OLS_l1_LRtr_noCo <-VECM(barry, lag=1, LRinclude="trend", include="none")
vecm_OLS_l1_LRbo <-VECM(barry, lag=1, LRinclude="both")

vecm_ML_l1_co <-VECM(barry, lag=1, estim="ML")
vecm_ML_l3_co <-VECM(barry, lag=3, include="const", estim="ML")
# vecm_ML_l3_co_betaGiven<-VECM(barry, lag=3, include="const", beta=-1, estim="ML")
vecm_ML_l1_tr <-VECM(barry, lag=1, include="trend", estim="ML")
vecm_ML_l1_bo <-VECM(barry, lag=1, include="both", estim="ML")
vecm_ML_l1_no <-VECM(barry, lag=1, include="none", estim="ML")
vecm_ML_l1_LRco <-VECM(barry, lag=1, LRinclude="const", estim="ML")
vecm_ML_l1_LRtr <-VECM(barry, lag=1, LRinclude="trend", estim="ML")
vecm_ML_l1_LRtr_noCo <-VECM(barry, lag=1, LRinclude="trend", include="none", estim="ML")
vecm_ML_l1_LRbo <-VECM(barry, lag=1, LRinclude="both", estim="ML")

vecm_all <- list(
		vecm_OLS_l1_co, vecm_OLS_l3_co, vecm_OLS_l3_co_betaGiven, vecm_OLS_l1_tr, 
		vecm_OLS_l1_bo, vecm_OLS_l1_no, vecm_OLS_l1_LRco, vecm_OLS_l1_LRtr, 
		vecm_OLS_l1_LRtr_noCo, vecm_OLS_l1_LRbo, 
		vecm_ML_l1_co, vecm_ML_l3_co,  vecm_ML_l1_tr, 
		vecm_ML_l1_bo, vecm_ML_l1_no, vecm_ML_l1_LRco, vecm_ML_l1_LRtr, 
		vecm_ML_l1_LRtr_noCo, vecm_ML_l1_LRbo)

names(vecm_all) <-c("vecm_OLS_l1_co", "vecm_OLS_l3_co", "vecm_OLS_l3_co_betaGiven", "vecm_OLS_l1_tr", 
		"vecm_OLS_l1_bo", "vecm_OLS_l1_no", "vecm_OLS_l1_LRco", "vecm_OLS_l1_LRtr", 
		"vecm_OLS_l1_LRtr_noCo", "vecm_OLS_l1_LRbo", 
		"vecm_ML_l1_co", "vecm_ML_l3_co", " vecm_ML_l1_tr", 
		"vecm_ML_l1_bo", "vecm_ML_l1_no", "vecm_ML_l1_LRco", "vecm_ML_l1_LRtr", 
		"vecm_ML_l1_LRtr_noCo", "vecm_ML_l1_LRbo")

vecm_ML <- vecm_all[grep("ML", names(vecm_all))]

lapply(vecm_all, print)
lapply(vecm_all, summary)

lapply(vecm_all, function(x) head(residuals(x), 3))
lapply(vecm_all, function(x) head(fitted(x), 3))
lapply(vecm_all, function(x) head(fitted(x, level="original"), 3))
sapply(vecm_all, deviance)



## logLik
sapply(vecm_all, logLik)
sapply(vecm_ML, logLik, r=0)
sapply(vecm_ML, logLik, r=1)
sapply(vecm_ML, logLik, r=2)

## AIC/BIC
sapply(vecm_all, AIC)
sapply(vecm_ML, AIC, r=0)
sapply(vecm_ML, AIC, r=1)
sapply(vecm_ML, AIC, r=2)
sapply(vecm_ML, AIC, r=0, fitMeasure="LL")
sapply(vecm_ML, AIC, r=1, fitMeasure="LL")
sapply(vecm_ML, AIC, r=2, fitMeasure="LL")

sapply(vecm_all, BIC)
sapply(vecm_ML, BIC, r=0)
sapply(vecm_ML, BIC, r=0, fitMeasure="LL")


## coint
sapply(vecm_all, function(x) x$model.specific$coint )
sapply(vecm_all, function(x) x$model.specific$beta)

### VARrep
lapply(vecm_all, VARrep)

### fevd
lapply(vecm_all, function(x) sapply(fevd(x, n.ahead=2), head))

### irf
vecm_irf <- vecm_all[-grep("l1_no|bo", names(vecm_all))] ## does not work for these models
lapply(vecm_irf, function(x) sapply(irf(x, runs=1)$irf,head,2))

## predict
lapply(vecm_all[-c(5, 6, 9, 10, 14, 15, 18, 19)], predict,  n.ahead=2)
lapply(vecm_all[-c(5, 6, 9, 10, 14, 15, 18, 19)], function(x) sapply(tsDyn:::predictOld.VECM(x, n.ahead=2)$fcst, function(y) y[,"fcst"]))

### rank test
vecm_ML_rtest <- vecm_ML[-grep("vecm_ML_l1_LRtr_noCo|vecm_ML_l1_LRbo", names(vecm_ML))] ## does not work for these models

rank.tests <- lapply(vecm_ML_rtest , rank.test)
rank.tests_rnull1 <- lapply(vecm_ML_rtest , rank.test, r_null=1)
rank.tests_tr <- lapply(vecm_ML_rtest , rank.test, type="trace")
rank.tests_tr_rnull1 <- lapply(vecm_ML_rtest , rank.test, r_null=1, type="trace")

rank.tests.all <- c(rank.tests , rank.tests_rnull1, rank.tests_tr,rank.tests_tr_rnull1 )

lapply(rank.tests.all, print)
lapply(rank.tests.all, summary)


### rank select
data(barry)
r_sel <- rank.select(barry)
r_sel_tre <- rank.select(barry, include="trend")
r_sel_none <- rank.select(barry, include="none")
r_sel_both <- rank.select(barry, include="both")

r_sel$LLs
r_sel$AICs

r_sel_tre$LLs
r_sel_tre$AICs

r_sel_none$LLs
r_sel_none$AICs

r_sel_both$LLs
r_sel_both$AICs
