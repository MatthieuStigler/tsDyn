library(tsDyn)


data(zeroyld)
data(barry)




## Test a few VECM models
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
var_all_fevd <- var_all[-grep("diff|adf", names(var_all))]
lapply(var_all_fevd , function(x) sapply(fevd(x, n.ahead=2), head))



## predict
var_all_pred <- var_all[-grep("bo|no|adf|diff", names(var_all))]
lapply(var_all_pred, function(x) sapply(predict(x, n.ahead=2)$fcst, function(y) y[,"fcst"]))
lapply(var_all, function(x) try(sapply(predict(x, n.ahead=2)$fcst, function(y) y[,"fcst"]), silent=TRUE))

## boot
var_all_boot <- var_all[-grep("adf|diff", names(var_all))]
lapply(var_all_boot, function(x) head(VAR.boot(x, seed=1234),2))

