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




###############################################################
### Check Johansen MLE: comparing with vars package
###############################################################

if(require(vars)){
data(Canada)

print(summary(VECM(Canada, lag=2, include="const", estim="ML")))
print(summary(VECM(Canada, lag=2, include="trend", estim="ML")))
print(summary(VECM(Canada, lag=2, include="both", estim="ML")))
print(summary(VECM(Canada, lag=2, include="none", estim="ML")))
print(summary(VECM(Canada, lag=2, LRinclude="const", estim="ML")))
print(summary(VECM(Canada, lag=2, LRinclude="trend", estim="ML")))
print(summary(VECM(Canada, lag=2, LRinclude="both", estim="ML")))

## VECM, l=1, inc=const
myVECM_l1<-VECM(Canada, lag=1, include="const", estim="ML")
ca_tes_l1 <- ca.jo(Canada,K=2, spec="trans")
ca_tes_l1_tr <- ca.jo(Canada,K=2, spec="trans", type="trace")
VECM_vars_l1<-cajorls(ca_tes_l1)
print(all.equal(ca_tes_l1@teststat, rev(rank.test(myVECM_l1)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(ca_tes_l1_tr@teststat, rev(rank.test(myVECM_l1)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l1$beta, myVECM_l1$model.specific$coint, check.attributes=FALSE))
print(all.equal(t(coefficients(myVECM_l1)), coefficients(VECM_vars_l1$rlm), check.attr=FALSE))

## VECM, l=2, inc=const
myVECM_l2<-VECM(Canada, lag=2, include="const", estim="ML")
ca_tes_l2 <- ca.jo(Canada,K=3, spec="trans")
ca_tes_l2_tr <- ca.jo(Canada,K=3, spec="trans", type="trace")
VECM_vars_l2<-cajorls(ca_tes_l2)
print(all.equal(ca_tes_l2_tr@teststat, rev(rank.test(myVECM_l2)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(ca_tes_l2@teststat, rev(rank.test(myVECM_l2)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l2$beta, myVECM_l2$model.specific$coint, check.attributes=FALSE))
print(all.equal(t(coefficients(myVECM_l2)), coefficients(VECM_vars_l2$rlm), check.attr=FALSE))

## VECM, l=1, inc const in coint
myVECM_l1_conCoint<-VECM(Canada, lag=1, LRinclude="const", estim="ML")
ca_tes_l1_conCoint <- ca.jo(Canada,K=2, spec="trans", ecdet="const")
ca_tes_l1_conCoint_tr <- ca.jo(Canada,K=2, spec="trans", ecdet="const", type="trace")
VECM_vars_l1_conCoint <-cajorls(ca_tes_l1_conCoint)
print(all.equal(ca_tes_l1_conCoint@teststat, rev(rank.test(myVECM_l1_conCoint)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(ca_tes_l1_conCoint_tr@teststat, rev(rank.test(myVECM_l1_conCoint)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l1_conCoint$beta, myVECM_l1_conCoint$model.specific$coint, check.attributes=FALSE))
print(all.equal(t(coefficients(myVECM_l1_conCoint)), coefficients(VECM_vars_l1_conCoint$rlm), check.attr=FALSE))

## VECM, l=2, inc const in coint
myVECM_l2_conCoint<-VECM(Canada, lag=2, LRinclude="const", estim="ML")
ca_tes_l2_conCoint <- ca.jo(Canada,K=3, spec="trans", ecdet="const")
ca_tes_l2_conCoint_tr <- ca.jo(Canada,K=3, spec="trans", ecdet="const", type="trace")
VECM_vars_l2_conCoint <-cajorls(ca_tes_l2_conCoint)
print(all.equal(ca_tes_l2_conCoint@teststat, rev(rank.test(myVECM_l2_conCoint)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(ca_tes_l2_conCoint_tr@teststat, rev(rank.test(myVECM_l2_conCoint)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l2_conCoint$beta, myVECM_l2_conCoint$model.specific$coint, check.attributes=FALSE))
print(all.equal(t(coefficients(myVECM_l2_conCoint)), coefficients(VECM_vars_l2_conCoint$rlm), check.attr=FALSE))

## VECM, l=1, inc trend in coint
myVECM_l2_trCoint<-VECM(Canada, lag=2, LRinclude="trend", estim="ML")
ca_tes_l2_trCoint <- ca.jo(Canada,K=3, spec="trans", ecdet="trend")
ca_tes_l2_trCoint_tr <- ca.jo(Canada,K=3, spec="trans", ecdet="trend", type="trace")
VECM_vars_l2_trCoint <-cajorls(ca_tes_l2_trCoint)
print(all.equal(ca_tes_l2_trCoint@teststat, rev(rank.test(myVECM_l2_trCoint)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(ca_tes_l2_trCoint_tr@teststat, rev(rank.test(myVECM_l2_trCoint)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l2_trCoint$beta, myVECM_l2_trCoint$model.specific$coint, check.attributes=FALSE))
print(all.equal(VECM_vars_l2_trCoint$beta, myVECM_l2_trCoint$model.specific$coint, check.attributes=FALSE, tol=1e-06))
print(all.equal(t(coefficients(myVECM_l2_trCoint)), coefficients(VECM_vars_l2_trCoint$rlm), check.attr=FALSE, tol=1e-02))

## VECM, l=1, inc both in coint
myVECM_l2_boCoint<-VECM(Canada, lag=2, LRinclude="trend", estim="ML")
ca_tes_l2_boCoint <- ca.jo(Canada,K=3, spec="trans", ecdet="trend")
ca_tes_l2_boCoint_tr <- ca.jo(Canada,K=3, spec="trans", ecdet="trend", type="trace")
VECM_vars_l2_boCoint <-cajorls(ca_tes_l2_boCoint)
print(all.equal(ca_tes_l2_boCoint@teststat, rev(rank.test(myVECM_l2_boCoint)$res_df[,"eigen"]), check.attributes=FALSE))
print(all.equal(ca_tes_l2_boCoint_tr@teststat, rev(rank.test(myVECM_l2_boCoint)$res_df[,"trace"]), check.attributes=FALSE))
print(all.equal(VECM_vars_l2_boCoint$beta, myVECM_l2_boCoint$model.specific$coint, check.attributes=FALSE))
print(all.equal(VECM_vars_l2_boCoint$beta, myVECM_l2_boCoint$model.specific$coint, check.attributes=FALSE, tol=1e-06))
print(all.equal(t(coefficients(myVECM_l2_boCoint)), coefficients(VECM_vars_l2_boCoint$rlm), check.attr=FALSE, tol=1e-02))


## Check LL
l1<-2*(logLik(myVECM_l1,r=4)-logLik(myVECM_l1,r=3))
l2<-2*(logLik(myVECM_l1,r=3)-logLik(myVECM_l1,r=2))
l3<-2*(logLik(myVECM_l1,r=2)-logLik(myVECM_l1,r=1))
l4<-2*(logLik(myVECM_l1,r=1)-logLik(myVECM_l1,r=0))
print(c(l1,l2,l3,l4))
print(ca.jo(Canada, spec="trans"))
print(all.equal(c(l1, l2, l3, l4),ca.jo(Canada, spec="trans")@teststat))

print(AIC(myVECM_l1,r=0, k=2*log(log(myVECM_l1$t))))
print(AIC(myVECM_l1,r=1, k=2*log(log(myVECM_l1$t))))
print(AIC(myVECM_l1,r=2, k=2*log(log(myVECM_l1$t))))
print(AIC(myVECM_l1,r=3, k=2*log(log(myVECM_l1$t))))
print(AIC(myVECM_l1,r=4, k=2*log(log(myVECM_l1$t))))


}
