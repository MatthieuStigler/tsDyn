library(tsDyn)
data(zeroyld)




## Test a few VECM models
myVECM1<-VECM(zeroyld, lag=1)
myVECM2<-VECM(zeroyld, lag=3, include="const")
myVECM2a<-VECM(zeroyld, lag=3, include="const", beta=-1)
myVECM3<-VECM(zeroyld, lag=1, estim="ML")
myVECM4<-VECM(zeroyld, lag=3, estim="ML")


summary(myVECM1)
summary(myVECM2)
summary(myVECM2a)
summary(myVECM3)
summary(myVECM4)

myVECM1$model.specific$coint
myVECM1$model.specific$beta

myVECM2a$model.specific$coint
myVECM2a$model.specific$beta

myVECM3$model.specific$coint
myVECM3$model.specific$beta


###Check Johansen MLE: comparing with vars package
if(require(vars)){
data(Canada)

myVECM<-VECM(Canada, lag=1, include="const", estim="ML")
VECM_vars<-cajorls(ca.jo(Canada, spec="trans"))
all.equal(VECM_vars$beta, myVECM$model.specific$coint, check.attributes=FALSE)


## Check LL
l1<-2*(logLik(myVECM,r=4)-logLik(myVECM,r=3))
l2<-2*(logLik(myVECM,r=3)-logLik(myVECM,r=2))
l3<-2*(logLik(myVECM,r=2)-logLik(myVECM,r=1))
l4<-2*(logLik(myVECM,r=1)-logLik(myVECM,r=0))
l1;l2;l3;l4
ca.jo(Canada, spec="trans")
all.equal(c(l1, l2, l3, l4),ca.jo(Canada, spec="trans")@teststat)
logLik(myVECM,r=5)

AIC(myVECM,r=0, k=2*log(log(myVECM$t)))
AIC(myVECM,r=1, k=2*log(log(myVECM$t)))
AIC(myVECM,r=2, k=2*log(log(myVECM$t)))
AIC(myVECM,r=3, k=2*log(log(myVECM$t)))
AIC(myVECM,r=4, k=2*log(log(myVECM$t)))


}
