library(tsDyn)
data(zeroyld)




## Test a few VECM models
myVECM1 <-VECM(zeroyld, lag=1)
myVECM2 <-VECM(zeroyld, lag=3, include="const")
myVECM2a<-VECM(zeroyld, lag=3, include="const", beta=-1)
myVECM3 <-VECM(zeroyld, lag=1, estim="ML")
myVECM4 <-VECM(zeroyld, lag=3, estim="ML")
myVECM5 <-VECM(zeroyld, lag=2, estim="ML", LRinclude="const")
myVECM6 <-VECM(zeroyld, lag=2, estim="ML", LRinclude="trend")
myVECM7 <-VECM(zeroyld, lag=2, estim="ML", LRinclude="both")


summary(myVECM1)
summary(myVECM2)
summary(myVECM2a)
summary(myVECM3)
summary(myVECM4)
summary(myVECM5)
summary(myVECM6)
summary(myVECM7)

logLik(myVECM1)
logLik(myVECM2)
logLik(myVECM2a)
logLik(myVECM3)
logLik(myVECM4)
logLik(myVECM5)
logLik(myVECM6)
logLik(myVECM7)

logLik(myVECM3, r=2)
logLik(myVECM4, r=2)

AIC(myVECM3)
AIC(myVECM4)
AIC(myVECM5)
AIC(myVECM6)
AIC(myVECM7)
AIC(myVECM3, r=2)
AIC(myVECM4, r=2)
AIC(myVECM3, fitMeasure="LL")
AIC(myVECM4, fitMeasure="LL")
AIC(myVECM3, r=2, fitMeasure="LL")
AIC(myVECM4, r=2, fitMeasure="LL")

BIC(myVECM3)
BIC(myVECM4)
BIC(myVECM5)
BIC(myVECM6)
BIC(myVECM7)
BIC(myVECM3, r=2)
BIC(myVECM4, r=2)
BIC(myVECM3, fitMeasure="LL")
BIC(myVECM4, fitMeasure="LL")
BIC(myVECM3, r=2, fitMeasure="LL")
BIC(myVECM4, r=2, fitMeasure="LL")

myVECM1$model.specific$coint
myVECM1$model.specific$beta

myVECM2a$model.specific$coint
myVECM2a$model.specific$beta

myVECM3$model.specific$coint
myVECM3$model.specific$beta

myVECM4$model.specific$coint
myVECM5$model.specific$coint
myVECM6$model.specific$coint
myVECM7$model.specific$coint


### rank test
r.test_VECM4 <- rank.test(myVECM4)
r.test_VECM4_tr <- rank.test(myVECM4, type="trace")
r.test_VECM4_spec <- rank.test(myVECM4, r_null=1)
r.test_VECM4_spec_tr <- rank.test(myVECM4,, r_null=1, type="trace")

r.test_VECM5 <- rank.test(myVECM5)
r.test_VECM6 <- rank.test(myVECM6)


r.test_VECM4
r.test_VECM4_tr
r.test_VECM4_spec
r.test_VECM4_spec_tr
r.test_VECM5
r.test_VECM6

summary(r.test_VECM4)
summary(r.test_VECM5)
summary(r.test_VECM6)

### rank select
data(barry)
r_sel <- rank.select(barry)
r_sel_tre <- rank.select(barry, include="trend")

r_sel$LLs
r_sel$AICs

r_sel_tre$LLs
r_sel_tre$AICs

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
