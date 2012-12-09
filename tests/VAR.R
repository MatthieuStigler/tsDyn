library(tsDyn)


data(zeroyld)
data(barry)




## Test a few VAR models
var_l1_co <-lineVar(barry, lag=1, include="const")
var_l1_tr <-lineVar(barry, lag=1, include="trend")
var_l1_bo <-lineVar(barry, lag=1, include="both")
var_l1_no <-lineVar(barry, lag=1, include="none")

var_l3_co <-lineVar(barry, lag=3, include="const")
var_l3_tr <-lineVar(barry, lag=3, include="trend")
var_l3_bo <-lineVar(barry, lag=3, include="both")
var_l3_no <-lineVar(barry, lag=3, include="none")

var_l2_diff_co <-lineVar(barry, lag=2, include="const", I="diff")
var_l2_diff_tr <-lineVar(barry, lag=2, include="trend", I="diff")
var_l2_diff_bo <-lineVar(barry, lag=2, include="both", I="diff")
var_l2_diff_no <-lineVar(barry, lag=2, include="none", I="diff")

var_l2_adf_co <-lineVar(barry, lag=2, include="const", I="ADF")
var_l2_adf_tr <-lineVar(barry, lag=2, include="trend", I="ADF")
var_l2_adf_bo <-lineVar(barry, lag=2, include="both", I="ADF")
var_l2_adf_no <-lineVar(barry, lag=2, include="none", I="ADF")

var_all <- list(
		var_l1_co, var_l1_tr, var_l1_bo, var_l1_no, 
		var_l3_co, var_l3_tr, var_l3_bo, var_l3_no,
		var_l2_diff_co, var_l2_diff_tr, var_l2_diff_bo, var_l2_diff_no,
		var_l2_adf_co, var_l2_adf_tr, var_l2_adf_bo, var_l2_adf_no)


names(var_all) <-c(
		"var_l1_co", "var_l1_tr", "var_l1_bo", "var_l1_no", 
		"var_l3_co", "var_l3_tr", "var_l3_bo", "var_l3_no",
		"var_l2_diff_co", "var_l2_diff_tr", "var_l2_diff_bo", "var_l2_diff_no",
		"var_l2_adf_co", "var_l2_adf_tr", "var_l2_adf_bo", "var_l2_adf_no")


## Check methods:
lapply(var_all, print)
lapply(var_all, summary)

lapply(var_all, function(x) head(residuals(x), 3))
lapply(var_all, function(x) head(fitted(x), 3))
sapply(var_all, deviance)


## logLik/AIC/BIC
sapply(var_all, logLik)
sapply(var_all, AIC)
sapply(var_all, AIC, fitMeasure="LL")
sapply(var_all, BIC)
sapply(var_all, BIC, fitMeasure="LL")




### fevd
var_all_level <- var_all[-grep("diff|adf", names(var_all))]
lapply(var_all_level , function(x) sapply(fevd(x, n.ahead=2), head))



## predict
var_all_pred <- var_all[-grep("bo|no|adf|diff", names(var_all))]
lapply(var_all_pred, predict, n.ahead=2)
lapply(var_all, function(x) try(predict(x, n.ahead=2), silent=TRUE))
lapply(var_all_pred, function(x) sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"]))
lapply(var_all, function(x) try(sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"]), silent=TRUE))

all.equal(lapply(var_all_pred, predict, n.ahead=2), lapply(var_all_pred, function(x) sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"])), check=FALSE)

lapply(var_all_level , predict_rolling, nroll=2)
lapply(var_all_level , predict_rolling, nroll=2, refit.every=1)

## boot
var_all_boot <- var_all[-grep("adf|diff", names(var_all))]
lapply(var_all_boot, function(x) tail(VAR.boot(x, seed=1234),2))

## sim 
comp_tvar_sim <- function(mod, serie){
  ns <- nrow(serie)
  sim_mod <- TVAR.sim(B=coef(mod), lag=mod$lag, include=mod$include,nthresh=0, n=ns-mod$lag, innov=residuals(mod), starting=serie[1:mod$lag,,drop=FALSE])
  all.equal(sim_mod, as.matrix(serie)[-c(1:mod$lag),], check.attr=FALSE)
}

lapply(var_all_level, comp_tvar_sim, serie=barry)



###################################
####### predict_rolling check
###################################
n_ca<- nrow(barry)

#### No refit lag=1
pred_roll<-predict_rolling(object=lineVar(barry, lag=1), nroll=10, n.ahead=1)
pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1, newdata=barry[n_ca-1,,drop=FALSE]), tail(pred_roll,1), check=FALSE)

#### No refit lag=3
pred_roll<-predict_rolling(object=lineVar(barry, lag=3), nroll=10, n.ahead=1)
pred1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-10,,drop=FALSE])
pred2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-9,,drop=FALSE])
all.equal(rbind(pred1, pred2), head(pred_roll,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1, newdata=barry[((n_ca-2):n_ca)-1,,drop=FALSE]), tail(pred_roll,1), check=FALSE)

### Refit lag=1
pred_ref<-predict_rolling(object=lineVar(barry, lag=1), nroll=10, n.ahead=1, refit.every=1)
pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=1), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=1), n.ahead=1), tail(pred_ref,1), check=FALSE)

### Refit lag=3
pred_ref<-predict_rolling(object=lineVar(barry, lag=3), nroll=10, n.ahead=1, refit.every=1)
pred_ref_1 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3), n.ahead=1)
pred_ref_2 <- predict(lineVar(tsDyn:::myHead(barry,n_ca-9), lag=3), n.ahead=1)
all.equal(rbind(pred_ref_1, pred_ref_2), head(pred_ref,2), check=FALSE)
all.equal(predict(lineVar(tsDyn:::myHead(barry,n_ca-1), lag=3), n.ahead=1), tail(pred_ref,1), check=FALSE)
